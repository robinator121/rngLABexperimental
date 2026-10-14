"""Local, bounded V0.3 jobs with independent stage identities and compact checkpoints."""
from __future__ import annotations
import copy,csv,gzip,hashlib,importlib.metadata,json,math,multiprocessing,os,platform,socket,time,traceback
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
from rnglab.v02.statistics import correlated_summary
from rnglab.v02.runner import canonical,digest,atomic_json,csv_write,input_files
from rnglab.v02.sources import research_provider
from rnglab.v03.physics import make_model,KERNEL_SOURCE_SHA256
from rnglab.v03.uniforms import UniformProcess
from rnglab.v03.monitor import Monitor,RULES


def engine_hash():
    h=hashlib.sha256();root=Path(__file__).parent
    for name in ['runner.py','physics.py','uniforms.py','monitor.py']:
        h.update(name.encode());h.update((root/name).read_bytes())
    for relative in ['v02/statistics.py','v02/sources.py','v02/physics.py','v02/pimc.py','v02/monitor.py','providers/__init__.py','degradation/__init__.py']:
        h.update(relative.encode());h.update((root.parent/relative).read_bytes())
    return h.hexdigest()


def settings(config,workload):
    return {key:workload.get(key,config.get(key,default)) for key,default in
            [('samples',65536),('burn_in',8192),('window',2048),('monitor_rows',4096),
             ('application_window',8192),('min_measurement',8192),('confirmation_windows',2)]}


def model_spec(workload):
    operational={'name','samples','burn_in','window','monitor_rows','application_window','min_measurement','confirmation_windows','calibration_key','initialization'}
    return {k:v for k,v in workload.items() if k not in operational}


def calibration_key(config,workload):
    return digest({'model':model_spec(workload),'settings':settings(config,workload),'initialization':workload.get('initialization','deterministic')})[:20]


def channels(kind,dimension):
    if kind=='ising':return ('site','acceptance')
    if dimension==3:return ('site','proposal','acceptance')
    if dimension==2:return ('proposal','acceptance')
    return ('proposal',)


class Source:
    def __init__(self,spec,seed,dimension,kind,steps):
        self.spec=spec;self.dimension=dimension;self.position=0;self.audits=[];self.base=None
        self.initial_key_hash=None
        if spec.get('kind')=='provider':
            provider_spec=copy.deepcopy(spec['provider_spec'])
            interval=provider_spec.get('reseed_interval',2048);entropy_bits=provider_spec.get('entropy_bits',256)
            def slices(node,role='direct'):
                if isinstance(node,dict):
                    if 'file' in node:
                        default=256 if role=='initial_source' else (steps//interval+1)*entropy_bits if role=='reseed_source' else math.ceil(steps*dimension*32/256)*256
                        stride=node.pop('trial_stride_bits',default)
                        node['offset_bits']=node.get('offset_bits',0)+node.pop('stage_offset_bits',0)+seed[-1]*stride
                        node['limit_bits']=node.get('limit_bits',stride)
                    for key,value in node.items():slices(value,key if key in {'initial_source','reseed_source'} else role)
                elif isinstance(node,list):
                    for value in node:slices(value,role)
            slices(provider_spec)
            self.source,self.base,self.audits=research_provider(provider_spec,seed,steps)
            if hasattr(self.base,'key'):self.initial_key_hash=hashlib.sha256(self.base.key).hexdigest()
            self.provenance={'classification':'existing bit-source/CSPRNG demonstration; see seed/reseed audit','spec':spec,'seed':seed}
        else:
            self.source=UniformProcess(spec,seed,dimension,channels(kind,dimension))
            self.provenance=self.source.provenance
    def sample(self,count):
        if self.base is not None:u,_=self.source.sample(count,self.dimension)
        else:u=self.source.sample(count)
        self.position+=count
        return u
    def details(self):
        if self.base is None:return {'latent_draws':self.source.draw_counters}
        return {'initial_key_sha256':self.initial_key_hash,'final_key_sha256':hashlib.sha256(self.base.key).hexdigest() if hasattr(self.base,'key') else None,
                'reseed_events':self.base.events,'source_audit':[e for a in self.audits for e in a.events],
                'seed_material_bits':sum(a.position for a in self.audits)}


def boundary_events(scores,times,threshold,k,onset=None,recovery=None):
    """Confirmation times are never backdated. Post-onset counters reset explicitly."""
    excursion=confirmation=post_excursion=post_confirmation=regain=None
    bad=post_bad=good=0
    for value,time_index in zip(scores,times):
        bad=bad+1 if value>threshold else 0
        if value>threshold and excursion is None:excursion=int(time_index)
        if bad>=k and confirmation is None:confirmation=int(time_index)
        if onset is not None and time_index>onset:
            post_bad=post_bad+1 if value>threshold else 0
            if value>threshold and post_excursion is None:post_excursion=int(time_index)
            if post_bad>=k and post_confirmation is None:post_confirmation=int(time_index)
        if recovery is not None and time_index>recovery:
            good=good+1 if value<=threshold else 0
            if good>=k and regain is None:regain=int(time_index)
    return {'first_excursion':excursion,'confirmation':confirmation,'post_onset_excursion':post_excursion,
            'post_onset_confirmation':post_confirmation,'post_recovery_within_band':regain,
            'pre_onset_confirmation':bool(onset is not None and confirmation is not None and confirmation<=onset)}


def stage_phase(start,end,onset,recovery,termination):
    if onset is None:return 'healthy'
    if end<=onset:return 'healthy_baseline'
    if start<onset:return 'onset_crossing'
    if termination is not None and start<termination<end:return 'termination_crossing'
    if recovery is not None and start<recovery<end:return 'recovery_crossing'
    if recovery is not None and start>=recovery:return 'post_recovery'
    if termination is not None and start>=termination:return 'source_recovery'
    return 'active_degradation'


def one_trial(config,workload,scenario,trial):
    tick=time.perf_counter();s=settings(config,workload);total_steps=s['burn_in']+s['samples']
    seed=[config['master_seed'],config['seed_namespace'],int.from_bytes(hashlib.sha256(workload['name'].encode()).digest()[:4],'big'),trial]
    model=make_model(model_spec(workload));source=Source(scenario['stream'],seed,model.dimension,workload['kind'],total_steps)
    initial_uniforms=0;initialization_hash=None
    if workload.get('initialization')=='equilibrium':
        initial_uniforms=model.initialization_uniform_count
        if initial_uniforms:
            initial=UniformProcess({'kind':'iid'},[*seed,439],1,('initial_state',)).sample(initial_uniforms).ravel()
            model.initialize_from_uniforms(initial)
            initialization_hash=hashlib.sha256(initial.astype('<f8').tobytes()).hexdigest()
    monitor=Monitor({'monitor_rows':s['monitor_rows'],'confirmation_windows':s['confirmation_windows']},model.dimension)
    onset=scenario.get('onset');recovery=scenario.get('source_recovery');termination=scenario.get('end')
    calibration=config.get('_calibration',{}).get(calibration_key(config,workload))
    retained=0;total=0.;recent=np.empty(0);chunks=[];aux_chunks=defaultdict(list);windows=[];monitors=[]
    uniform_hash=hashlib.sha256();u_count=pair_count=0;u_last=None
    u_mean=np.zeros(model.dimension);u_m2=np.zeros(model.dimension)
    pair_a=np.zeros(model.dimension);pair_b=np.zeros(model.dimension);pair_m2a=np.zeros(model.dimension);pair_m2b=np.zeros(model.dimension);pair_cov=np.zeros(model.dimension)
    u_min=np.full(model.dimension,np.inf);u_max=np.full(model.dimension,-np.inf)
    for start in range(0,total_steps,s['window']):
        count=min(s['window'],total_steps-start);end=start+count
        u=source.sample(count);uniform_hash.update(np.asarray(u,dtype='<f8').tobytes())
        block_mean=u.mean(axis=0);delta=block_mean-u_mean;combined=u_count+count
        u_m2+=np.sum((u-block_mean)**2,axis=0)+delta*delta*(u_count*count/combined)
        u_mean+=delta*(count/combined);u_count=combined
        u_min=np.minimum(u_min,u.min(axis=0));u_max=np.maximum(u_max,u.max(axis=0))
        left=np.vstack((u_last,u[:-1])) if u_last is not None else u[:-1]
        right=u if u_last is not None else u[1:]
        npairs=len(left)
        if npairs:
            ma,mb=left.mean(axis=0),right.mean(axis=0);da,db=ma-pair_a,mb-pair_b;combined_pairs=pair_count+npairs
            factor=pair_count*npairs/combined_pairs
            pair_m2a+=np.sum((left-ma)**2,axis=0)+da*da*factor
            pair_m2b+=np.sum((right-mb)**2,axis=0)+db*db*factor
            pair_cov+=np.sum((left-ma)*(right-mb),axis=0)+da*db*factor
            pair_a+=da*(npairs/combined_pairs);pair_b+=db*(npairs/combined_pairs);pair_count=combined_pairs
        u_last=u[-1].copy()
        values,aux=model.advance(u)
        mr=monitor.update(u,end)
        monitors.append({key:mr[key] for key in ('time',*RULES)})
        mean=recent_mean=error=recent_error=None
        if end>s['burn_in']:
            # Budgets enforce window-aligned burn-in, so no partial discard ambiguity.
            chunks.append(values);retained+=count;total+=float(values.sum());mean=total/retained
            recent=np.concatenate((recent,values))[-s['application_window']:];recent_mean=float(recent.mean())
            error=abs(mean-model.target);recent_error=abs(recent_mean-model.target)
            for name,array in aux.items():aux_chunks[name].append(array)
        window={'start':start,'time':end,'retained':retained,'phase':stage_phase(start,end,onset,recovery,termination),
                'window_mean':float(np.mean(values)),'cumulative_mean':mean,'recent_mean':recent_mean,
                'cumulative_error':error,'recent_error':recent_error,
                'acceptance':float(np.mean(aux['accepted'])) if 'accepted' in aux else None}
        if calibration and retained>=s['min_measurement']:
            for criterion in ['cumulative','recent']:
                scale=calibration['scale'][criterion][str(end)]
                window[criterion+'_score']=window[criterion+'_error']/scale
                window[criterion+'_tolerance']=scale*calibration['application_cutoffs'][criterion]
        windows.append(window)
    values=np.concatenate(chunks);stats=correlated_summary(values)
    diagnostics={}
    for name,arrays in aux_chunks.items():
        x=np.concatenate(arrays);summary={'mean':float(x.mean()),'sd':float(x.std(ddof=1))}
        if name in model.targets:summary.update({'target':model.targets[name],'bias':float(x.mean()-model.targets[name])})
        if name in {'magnetization','abs_magnetization','thermodynamic','centroid','x2'}:summary.update(correlated_summary(x))
        diagnostics[name]=summary
    outcome={'detectors':{},'applications':{}}
    times=[m['time'] for m in monitors]
    if calibration:
        for rule in RULES:
            outcome['detectors'][rule]=boundary_events([m[rule] for m in monitors],times,
                calibration['monitor_cutoffs'][rule],s['confirmation_windows'],onset,recovery)
        eligible=[w for w in windows if w['retained']>=s['min_measurement']]
        for criterion in ['cumulative','recent']:
            outcome['applications'][criterion]=boundary_events([w[criterion+'_score'] for w in eligible],
                [w['time'] for w in eligible],calibration['application_cutoffs'][criterion],s['confirmation_windows'],onset,recovery)
    u_var=u_m2/u_count;degenerate=u_min==u_max
    denominator=np.sqrt(pair_m2a*pair_m2b)
    rho=np.clip(np.divide(pair_cov,denominator,out=np.zeros_like(pair_cov),where=denominator>0),-1,1)
    correlations=[None if degenerate[j] or denominator[j]<=0 else float(rho[j]) for j in range(model.dimension)]
    u_var[degenerate]=0
    result={'workload':workload['name'],'scenario':scenario['name'],'trial':trial,'seed':seed,
            'seed_namespace':config['seed_namespace'],'purpose':config['purpose'],'calibration_key':calibration_key(config,workload),
            'status':'complete','kind':workload['kind'],'model_spec':model_spec(workload),'model_target':model.target,
            'sampler_backend':model.backend if hasattr(model,'backend') else 'numpy','dimension':model.dimension,
            'samples':s['samples'],'burn_in':s['burn_in'],'proposals':total_steps,'emitted_uniforms':total_steps*model.dimension,
            'initial_state':workload.get('initialization','deterministic'),'initialization_uniforms':initial_uniforms,
            'initialization_uniform_sha256':initialization_hash,'total_simulation_uniforms':total_steps*model.dimension+initial_uniforms,
            'reference_kind':'exact finite discretization' if workload['kind']=='pimc' else 'exact variational' if workload['kind']=='vmc' else 'exact finite lattice' if workload['kind']=='ising' else 'analytic integral',
            'mean':stats['mean'],'signed_error':stats['mean']-model.target,'absolute_error':abs(stats['mean']-model.target),
            'statistics':stats,'auxiliary':diagnostics,'stream_sha256':uniform_hash.hexdigest(),
            'uniform_means':u_mean.tolist(),'uniform_variances':u_var.tolist(),'uniform_lag1_correlations':correlations,
            'degenerate_uniform_channels':degenerate.tolist(),
            'source_provenance':source.provenance,'source_details':source.details(),'fault_onset':onset,'source_recovery':recovery,
            'nonstationary_input':bool(scenario.get('nonstationary',False)),'horizon':total_steps,
            'elapsed_seconds':time.perf_counter()-tick}
    payload={'result':result,'windows':windows,'monitor':monitors,'events':outcome}
    traces=values if trial<config.get('trace_trials',0) else None
    return payload,traces


def job_identity(workload,scenario,trial):return f"{workload['name']}__{scenario['name']}__{trial:06d}"


def applicable(workload,scenario):
    return (not scenario.get('kinds') or workload['kind'] in scenario['kinds']) and (not scenario.get('workloads') or workload['name'] in scenario['workloads'])


def execute(job):
    config,w,scenario,trial,output=job;out=Path(output);name=job_identity(w,scenario,trial)
    lock=out/'trials'/f'{name}.lock';path=out/'trials'/f'{name}.json.gz'
    fd=os.open(lock,os.O_WRONLY|os.O_CREAT|os.O_EXCL)
    with os.fdopen(fd,'w') as f:f.write(canonical({'pid':os.getpid(),'hostname':socket.gethostname(),'started':time.time()}))
    try:
        payload,trace=one_trial(config,w,scenario,trial)
        temporary=path.with_name(path.name+f'.{os.getpid()}.tmp')
        with gzip.open(temporary,'wt',compresslevel=3) as f:f.write(canonical({'checksum':digest(payload),'payload':payload}))
        os.replace(temporary,path)
        if trace is not None:
            (out/'traces').mkdir(exist_ok=True);np.savez_compressed(out/'traces'/f'{name}.npz',observable=trace)
        return name,payload['result']['elapsed_seconds']
    except BaseException as exc:
        atomic_json(out/'failures'/f'{name}.json',{'status':'failed','exception':repr(exc),'traceback':traceback.format_exc(),'time':time.time()})
        raise
    finally:lock.unlink(missing_ok=True)


def payloads(output):
    for path in sorted((Path(output)/'trials').glob('*.json.gz')):
        with gzip.open(path,'rt') as f:obj=json.load(f)
        if digest(obj['payload'])!=obj['checksum']:raise ValueError('Corrupt trial '+path.name)
        yield obj['payload']


def validate(config):
    for key in ['trials','master_seed','seed_namespace']:
        if not isinstance(config[key],int) or config[key]<0:raise ValueError('invalid '+key)
    if config['trials']<2:raise ValueError('at least two trials required')
    for group in ['workloads','scenarios']:
        names=[s['name'] for s in config[group]]
        if not names or len(set(names))!=len(names) or any(not n.replace('_','').replace('-','').isalnum() for n in names):raise ValueError('unique safe names required')
    for w in config['workloads']:
        s=settings(config,w);model=make_model(model_spec(w))
        if s['samples']<s['application_window'] or s['samples']%s['window'] or s['burn_in']%s['window']:raise ValueError('window-aligned positive budgets required')
        if s['window']<16 or s['monitor_rows']<256 or s['confirmation_windows']<1 or s['min_measurement']<s['window']:raise ValueError('invalid windows')
        if s['min_measurement']+(s['confirmation_windows']-1)*s['window']>s['samples']:raise ValueError('insufficient retained samples for application confirmation')
        if not np.isfinite(model.target):raise ValueError('finite validated target required')
        if w.get('initialization','deterministic') not in {'deterministic','equilibrium'}:raise ValueError('unknown chain initialization')
    if config.get('purpose') in {'scale_fitting','threshold_fitting','healthy_validation','confirmation'}:
        def check_file(node):
            if isinstance(node,dict):
                if 'file' in node and 'stage_offset_bits' not in node:raise ValueError('fresh file-source stages require explicit disjoint stage_offset_bits')
                for value in node.values():check_file(value)
            elif isinstance(node,list):
                for value in node:check_file(value)
        check_file(config['scenarios'])


def run(config,output,workers=4,resume=False,max_jobs=None,shard_index=0,shard_count=1):
    config=dict(config);validate(config);out=Path(output);begin=time.perf_counter()
    if workers<1 or workers>max(1,len(os.sched_getaffinity(0))//2):raise ValueError('worker count exceeds responsiveness policy')
    calibration=None
    if config.get('calibration_file'):
        calibration_path=Path(config['calibration_file'])
        if not calibration_path.is_absolute():calibration_path=Path(__file__).parents[3]/calibration_path
        calibration=json.loads(calibration_path.read_text())
        if calibration['engine_sha256']!=engine_hash():raise ValueError('calibration engine mismatch')
        used={calibration.get(key,{}).get('namespace') for key in ['scale_stage','threshold_stage']}
        if config['seed_namespace'] in used:raise ValueError('fitted calibration namespace cannot be reused as new evidence')
        config['_calibration']=calibration['strata']
        for w in config['workloads']:
            if calibration_key(config,w) not in config['_calibration']:raise ValueError('unmatched calibration stratum')
    public={k:v for k,v in config.items() if k!='_calibration'}
    identity={'config':public,'engine_sha256':engine_hash(),'kernel_sha256':KERNEL_SOURCE_SHA256,
              'python':platform.python_version(),'dependencies':{n:importlib.metadata.version(n) for n in ['numpy','scipy','cryptography']},
              'calibration_sha256':digest(calibration) if calibration else None,'input_files':input_files(public)}
    manifest=out/'manifest.json'
    if manifest.exists():
        if not resume:raise ValueError('existing experiment requires --resume')
        if json.loads(manifest.read_text())['identity']!=identity:raise ValueError('resume identity changed')
    else:
        if out.exists() and any(out.iterdir()):raise ValueError('output nonempty without manifest')
        out.mkdir(parents=True,exist_ok=True);atomic_json(manifest,{'identity':identity,'created_unix':time.time(),'host':platform.platform(),'seed_policy':'disjoint stage namespace; scenarios share workload/trial seeds for pairing','purpose':config['purpose']})
    for name in ['trials','failures']:(out/name).mkdir(exist_ok=True)
    for _ in payloads(out):pass
    jobs=[];index=0
    for w in config['workloads']:
        for scenario in config['scenarios']:
            if not applicable(w,scenario):continue
            for trial in range(config['trials']):
                name=job_identity(w,scenario,trial)
                if index%shard_count==shard_index and not (out/'trials'/f'{name}.json.gz').exists():jobs.append((config,w,scenario,trial,str(out)))
                index+=1
    if max_jobs is not None:jobs=jobs[:max_jobs]
    progress=out/f'progress-{shard_index}.jsonl'
    def record(item,n):
        name,seconds=item
        with progress.open('a') as f:f.write(canonical({'trial':name,'seconds':seconds,'finished_unix':time.time()})+'\n')
        if n%100==0 or n==len(jobs):print(f'{n}/{len(jobs)} trials, {time.perf_counter()-begin:.1f}s',flush=True)
    if workers==1:
        for n,job in enumerate(jobs,1):record(execute(job),n)
    else:
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
            futures=[pool.submit(execute,job) for job in jobs]
            for n,future in enumerate(as_completed(futures),1):record(future.result(),n)
    atomic_json(out/f'execution-{shard_index}.json',{'new_jobs':len(jobs),'completed_trials':len(list((out/'trials').glob('*.json.gz'))),
        'configured_trials':sum(applicable(w,sc) for w in config['workloads'] for sc in config['scenarios'])*config['trials'],'failed_attempt_records':len(list((out/'failures').glob('*.json'))),
        'workers':workers,'elapsed_seconds':time.perf_counter()-begin})
    return len(jobs)
