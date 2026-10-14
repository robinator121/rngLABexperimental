"""Atomic trial checkpoints, deterministic process jobs and monitored physical runs."""
import csv,hashlib,importlib.metadata,json,math,multiprocessing,os,platform,time,traceback
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor,as_completed
from copy import deepcopy
from pathlib import Path
import numpy as np
from rnglab.v02.physics import make_model,ising_exact
from rnglab.v02.sources import research_provider
from rnglab.v02.monitor import BitMonitor
from rnglab.v02.statistics import correlated_summary,mean_interval,proportion_interval,sustained_events


def canonical(data):return json.dumps(data,sort_keys=True,separators=(',',':'),allow_nan=False)
def digest(data):return hashlib.sha256(canonical(data).encode()).hexdigest()


def engine_hash():
    root=Path(__file__).parent
    h=hashlib.sha256()
    for name in ['physics.py','sources.py','monitor.py','statistics.py','runner.py']:
        h.update(name.encode());h.update((root/name).read_bytes())
    # Include reused V0.1 source/degradation implementation.
    for name in ['providers','degradation']:
        h.update((root.parent/name/'__init__.py').read_bytes())
    return h.hexdigest()


def input_files(config):
    files={}
    def visit(node):
        if isinstance(node,dict):
            if 'file' in node:
                path=Path(node['file'])
                if not path.is_absolute():path=Path(__file__).parents[3]/path
                files[node['file']]={'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'bytes':path.stat().st_size}
            for value in node.values():visit(value)
        elif isinstance(node,list):
            for value in node:visit(value)
    visit(config);return files


def atomic_json(path,data):
    path=Path(path);temp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    temp.write_text(canonical(data)+'\n');os.replace(temp,path)


def csv_write(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for row in rows for k in row))
    temp=Path(path).with_suffix('.csv.tmp')
    with temp.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
    os.replace(temp,path)


def validate(config):
    for name in ['trials','samples','window']:
        if not isinstance(config[name],int) or config[name]<1:raise ValueError(f'invalid {name}')
    if config['trials']<2 or config['samples']<16:raise ValueError('at least two trials and sixteen retained samples required')
    if config.get('burn_in',0)<0 or config.get('burn_in',0)%config['window']:raise ValueError('burn_in must be a nonnegative multiple of window')
    if not config.get('scenarios') or not config.get('workloads'):raise ValueError('empty experiment')
    for key in ['scenarios','workloads']:
        names=[s['name'] for s in config[key]]
        if len(set(names))!=len(names) or any(not n.replace('_','').replace('-','').isalnum() for n in names):raise ValueError('names must be unique safe identifiers')
    for w in config['workloads']:
        m=make_model(w)
        if m.target is None:raise ValueError('research detection requires an exact target; use L=3 or L=4')
        if not w.get('thresholds') or any(t<=0 for t in w['thresholds']):raise ValueError('positive predeclared thresholds required')
    if not 0<config.get('monitor_alpha',.005)<1:raise ValueError('invalid alpha')


def file_slices(spec,trial,steps):
    """File trials use disjoint slices; the same slices are paired across workloads."""
    spec=deepcopy(spec)
    def visit(node):
        if isinstance(node,dict):
            if 'file' in node:
                stride=node.pop('trial_stride_bits',steps*64)
                node['offset_bits']=node.get('offset_bits',0)+trial*stride
                node['limit_bits']=node.get('limit_bits',stride)
            for v in node.values():visit(v)
        elif isinstance(node,list):
            for v in node:visit(v)
    visit(spec);return spec


def one_trial(config,scenario,workload,trial):
    started=time.perf_counter();burn=config.get('burn_in',0);steps=burn+config['samples'];window=config['window']
    wid=int.from_bytes(hashlib.sha256(workload['name'].encode()).digest()[:4],'big')
    seed=[config['seed'],trial,wid]
    spec=file_slices(scenario,trial,steps);spec['monitor_stride']=window
    provider,base,audits=research_provider(spec,seed,steps)
    key_hash=hashlib.sha256(base.key).hexdigest() if hasattr(base,'key') else None
    model=make_model(workload)
    monitor=BitMonitor(config.get('monitor_alpha',.005),math.ceil(steps/window),config.get('monitor_bits',window*64))
    chunks=[];aux_chunks=defaultdict(list);windows=[];all_trace=[];consumed_hash=hashlib.sha256()
    total=0.;retained=0;accept_burn=0;burn_observations=0;start=0
    while start<steps:
        count=min(window,steps-start);end=start+count
        u,bits=provider.sample(count,2)
        consumed_hash.update(np.packbits(bits).tobytes())
        values,aux=model.advance(u)
        alarm=monitor.update(bits,end,u.ravel())
        if scenario.get('provider')=='triggered' and alarm and end<steps:
            base._reseed()
        phase='burn_in' if end<=burn else 'measurement'
        if phase=='measurement':
            chunks.append(values);total+=float(values.sum());retained+=len(values)
            for name,data in aux.items():aux_chunks[name].append(data)
            cumulative=total/retained
        else:
            cumulative=None
            if 'accepted' in aux:accept_burn+=int(aux['accepted'].sum());burn_observations+=count
        windows.append({'time':end,'retained':retained,'phase':phase,'window_mean':float(np.mean(values)),
                        'cumulative_mean':cumulative,'absolute_error':abs(cumulative-model.target) if cumulative is not None else None,
                        'acceptance':float(aux['accepted'].mean()) if 'accepted' in aux else None,
                        'alarm':alarm,'bit_bias':monitor.records[-1]['bit_bias']})
        if trial<config.get('trace_trials',0):all_trace.append(values)
        start=end
    values=np.concatenate(chunks);stats=correlated_summary(values)
    measured=[w for w in windows if w['phase']=='measurement' and w['retained']>=config.get('failure_min_samples',2048)]
    onset=scenario.get('fault_onset_step')
    after=[r['time'] for r in monitor.records if r['alarm'] and (onset is None or r['time']>onset)]
    detection=min(after) if after else None
    source_records=[];source_events=[]
    if audits:
        # The shared audit monitor records actual initial/reseed material, not an independent surrogate.
        source_records=audits[0].monitor.records
        source_events=[e for a in audits for e in a.events]
    source_alarm=min((r['time'] for r in source_records if r['alarm']),default=None)
    outcome=[]
    for threshold in workload['thresholds']:
        failure,recovery=sustained_events([w['absolute_error'] for w in measured],[w['time'] for w in measured],threshold,config.get('sustained_windows',3))
        status=('both' if detection is not None else 'failure_only') if failure is not None else ('detection_only' if detection is not None else 'neither')
        outcome.append({'threshold':threshold,'failure_time':failure,'recovery_time':recovery,'detection_time':detection,
                        'delta_t':failure-detection if failure is not None and detection is not None else None,
                        'status':status,'failure_before_detection':bool(failure is not None and (detection is None or failure<detection)),
                        'source_detection_time':source_alarm,
                        'delta_t_source':failure-source_alarm if failure is not None and source_alarm is not None else None,
                        'horizon':steps})
    aux_summary={}
    for name,arrays in aux_chunks.items():
        data=np.concatenate(arrays)
        aux_summary[name+'_mean']=float(data.mean())
        if name in {'magnetization','abs_magnetization','x2'}:
            aux_summary[name+'_tau']=correlated_summary(data)['tau']
    result={'scenario':scenario['name'],'workload':workload['name'],'kind':workload['kind'],'trial':trial,'seed':seed,
            'status':'complete','samples':config['samples'],'burn_in':burn,'steps':steps,'uniforms':steps*2,'output_bits':steps*64,
            'sweeps':steps/workload.get('length',4)**2 if workload['kind']=='ising' else None,
            'target':model.target,'ground_state_energy':.5 if workload['kind']=='vmc' else None,
            'variational_offset':model.target-.5 if workload['kind']=='vmc' else None,
            **stats,**aux_summary,'error':stats['mean']-model.target,
            'burn_acceptance':accept_burn/burn_observations if burn_observations else None,
            'first_alarm':monitor.first_alarm,'detection_after_onset':detection,'fault_onset_step':onset,
            'detection_delay':detection-onset if onset is not None and detection is not None else None,
            'pre_onset_alarm':bool(onset is not None and monitor.first_alarm is not None and monitor.first_alarm<=onset),
            'source_first_alarm':source_alarm,'initial_key_sha256':key_hash,'output_sha256':consumed_hash.hexdigest(),
            'entropy_bits_consumed':sum(a.position for a in audits),'reseed_count':len(base.events),
            'direct_input_bits':getattr(getattr(base,'source',None),'raw_bits_consumed',steps*64) if not audits else 0,
            'elapsed_seconds':time.perf_counter()-started,
            'stationary_uncertainty_applicable':not bool(scenario.get('nonstationary',False))}
    payload={'result':result,'windows':windows,'monitor':monitor.records,'source_monitor':source_records,
             'source_events':source_events,'reseed_events':base.events,'outcomes':outcome}
    return payload,np.concatenate(all_trace) if all_trace else None


def execute_job(job):
    config,scenario,workload,trial,output=job
    name=f"{scenario['name']}__{workload['name']}__{trial:05d}"
    output=Path(output);path=output/'trials'/f'{name}.json';lock=path.with_suffix('.lock')
    fd=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    os.close(fd)
    try:
        payload,trace=one_trial(config,scenario,workload,trial)
        if trace is not None:
            (output/'traces').mkdir(exist_ok=True)
            np.savez_compressed(output/'traces'/f'{name}.npz',observable=trace)
        atomic_json(path,{'checksum':digest(payload),'payload':payload})
        return name,payload['result']['elapsed_seconds']
    except Exception as e:
        atomic_json(output/'failures'/f'{name}.json',{'status':'failed','exception':repr(e),'traceback':traceback.format_exc()})
        raise
    finally:
        lock.unlink(missing_ok=True)


def read_checkpoints(output):
    data=[]
    for path in sorted((Path(output)/'trials').glob('*.json')):
        record=json.loads(path.read_text())
        if digest(record['payload'])!=record['checksum']:raise ValueError(f'corrupt trial checkpoint: {path.name}')
        data.append(record['payload'])
    return data


def summarize(output,config):
    data=read_checkpoints(output);rows=[p['result'] for p in data];groups=defaultdict(list)
    for row in rows:groups[(row['scenario'],row['workload'])].append(row)
    summary=[];detection_rows=[];paired=[]
    for (scenario,workload),g in groups.items():
        estimates=np.array([r['mean'] for r in g]);errors=np.array([r['error'] for r in g]);unique=len({r['output_sha256'] for r in g})
        mean,lo,hi=mean_interval(estimates)
        if unique<len(g):lo=hi=None
        summary.append({'scenario':scenario,'workload':workload,'trials':len(g),'unique_streams':unique,'target':g[0]['target'],
                        'mean':mean,'mean_ci95_low':lo,'mean_ci95_high':hi,'bias':float(errors.mean()),'variance':float(estimates.var(ddof=1)) if len(g)>1 else None,
                        'rmse':float(np.sqrt(np.mean(errors**2))),'mae':float(np.mean(np.abs(errors))),
                        'mean_ess':float(np.mean([r['ess'] for r in g if r['ess'] is not None])) if any(r['ess'] is not None for r in g) else None,
                        'mean_acceptance':float(np.mean([r['accepted_mean'] for r in g])) if 'accepted_mean' in g[0] else None,
                        'mean_seconds':float(np.mean([r['elapsed_seconds'] for r in g]))})
    outcomes=[{k:p['result'][k] for k in ['scenario','workload','trial']}|o for p in data for o in p['outcomes']]
    for scenario,workload in groups:
        thresholds=sorted({r['threshold'] for r in outcomes if r['scenario']==scenario and r['workload']==workload})
        g=groups[(scenario,workload)];unique=len({r['output_sha256'] for r in g})
        for threshold in thresholds:
            o=[r for r in outcomes if r['scenario']==scenario and r['workload']==workload and r['threshold']==threshold]
            detection=sum(r['detection_time'] is not None for r in o);failure=sum(r['failure_time'] is not None for r in o)
            dp,dl,dh=proportion_interval(detection,len(o));fp,fl,fh=proportion_interval(failure,len(o))
            if unique<len(g):dl=dh=fl=fh=None
            delta=[r['delta_t'] for r in o if r['delta_t'] is not None]
            detection_rows.append({'scenario':scenario,'workload':workload,'threshold':threshold,'trials':len(o),'unique_streams':unique,
                                   'detection_probability':dp,'detection_ci_low':dl,'detection_ci_high':dh,
                                   'failure_probability':fp,'failure_ci_low':fl,'failure_ci_high':fh,
                                   'failure_before_detection':sum(r['failure_before_detection'] for r in o),
                                   'both_observed':len(delta),'median_delta_t_conditional':float(np.median(delta)) if delta else None,
                                   **{f'n_{s}':sum(r['status']==s for r in o) for s in ['both','failure_only','detection_only','neither']}})
    control=config.get('paired_control','synthetic')
    for (scenario,workload),g in groups.items():
        baseline={r['trial']:r for r in groups.get((control,workload),[])}
        matches=[(r,baseline[r['trial']]) for r in g if r['trial'] in baseline]
        if not matches or scenario==control:continue
        if len({r['output_sha256'] for r,_ in matches})<len(matches):continue
        diff=[abs(r['error'])-abs(b['error']) for r,b in matches]
        mean,lo,hi=mean_interval(diff,confidence=1-.05/max(1,(len(config['scenarios'])-1)*len(config['workloads'])))
        paired.append({'scenario':scenario,'control':control,'workload':workload,'pairs':len(diff),'mean_absolute_error_difference':mean,
                       'familywise95_low':lo,'familywise95_high':hi,'method':'paired t; Bonferroni over configured scenario/workload contrasts'})
    csv_write(Path(output)/'trial_summary.csv',rows);csv_write(Path(output)/'summary.csv',summary)
    csv_write(Path(output)/'outcomes.csv',outcomes);csv_write(Path(output)/'detection_summary.csv',detection_rows);csv_write(Path(output)/'paired_differences.csv',paired)
    window_rows=[{k:p['result'][k] for k in ['scenario','workload','trial']}|w for p in data for w in p['windows']]
    csv_write(Path(output)/'windows.csv',window_rows)
    atomic_json(Path(output)/'execution.json',{'completed_trials':len(rows),'configured_trials':len(config['scenarios'])*len(config['workloads'])*config['trials'],
                                            'sum_trial_seconds':sum(r['elapsed_seconds'] for r in rows),'failed_records':len(list((Path(output)/'failures').glob('*.json')))})
    return summary


def run(config,output,workers=1,resume=False,shard_index=0,shard_count=1,max_jobs=None):
    validate(config);output=Path(output);started=time.perf_counter()
    if not 0<=shard_index<shard_count:raise ValueError('invalid shard')
    versions={n:importlib.metadata.version(n) for n in ['numpy','scipy','cryptography']}
    identity={'config':config,'engine_sha256':engine_hash(),'dependencies':versions,'python':platform.python_version(),'input_files':input_files(config)}
    manifest=output/'manifest.json'
    if manifest.exists():
        if not resume:raise ValueError('existing results require --resume')
        old=json.loads(manifest.read_text())
        if old['identity']!=identity:raise ValueError('resume rejected: configuration, code, Python or dependencies changed')
    else:
        if output.exists() and any(output.iterdir()):raise ValueError('output is not empty and has no valid manifest')
        output.mkdir(parents=True,exist_ok=True)
        atomic_json(manifest,{'identity':identity,'identity_sha256':digest(identity),'environment':platform.platform(),'created_unix':time.time(),
                             'seed_scheme':'[master,trial,sha256(workload_name)[:4] as big-endian integer]; scenarios paired',
                             'provenance':'synthetic models unless file source explicitly declared; all consumed output is hashed'})
    for name in ['trials','failures']:(output/name).mkdir(exist_ok=True)
    # Verify existing checkpoints before skipping any work.
    read_checkpoints(output)
    jobs=[];index=0
    for scenario in config['scenarios']:
        for workload in config['workloads']:
            for trial in range(config['trials']):
                name=f"{scenario['name']}__{workload['name']}__{trial:05d}.json"
                if index%shard_count==shard_index and not (output/'trials'/name).exists():jobs.append((config,scenario,workload,trial,str(output)))
                index+=1
    if max_jobs is not None:jobs=jobs[:max_jobs]
    log=output/f'progress-shard-{shard_index}.jsonl'
    def record(result,done):
        name,seconds=result
        with log.open('a') as f:f.write(canonical({'trial':name,'seconds':seconds,'finished_unix':time.time()})+'\n')
        if done%25==0 or done==len(jobs):print(f'{done}/{len(jobs)} jobs; {time.perf_counter()-started:.1f}s elapsed',flush=True)
    if workers==1:
        for i,job in enumerate(jobs,1):record(execute_job(job),i)
    else:
        with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
            futures=[pool.submit(execute_job,job) for job in jobs]
            for i,future in enumerate(as_completed(futures),1):record(future.result(),i)
    if shard_count==1:summarize(output,config)
    print(f'Completed this invocation: {len(jobs)} jobs in {time.perf_counter()-started:.1f}s',flush=True)
    return len(jobs)
