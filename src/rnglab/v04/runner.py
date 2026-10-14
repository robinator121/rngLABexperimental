"""V0.4 trials using the existing atomic checksum/checkpoint conventions.

All computation remains local. Historical modules are imported read-only.
Recorded input layers, supplied-draw identities and paired seed groups are
explicit; model names do not create artificial independence for matched paths.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ProcessPoolExecutor,as_completed
import ctypes
import gzip
import hashlib
import importlib.metadata
import json
import multiprocessing
import os
from pathlib import Path
import platform
import socket
import time
import traceback

import numpy as np
from scipy.special import ndtri,ndtr
from cryptography.hazmat.primitives.ciphers import Cipher,algorithms
from rnglab.v02.runner import canonical,digest,atomic_json
from rnglab.v02.statistics import correlated_summary
from rnglab.v03.runner import boundary_events,payloads
from rnglab.v03.monitor import Monitor
from rnglab.v04.models import Model,KERNEL_SHA256,COMPILER_FLAGS

ROOT=Path(__file__).resolve().parents[3]
DETECTORS=('marginal','normal_cross','normal_lag','composite_v03')
CORE_FILES=('v04/runner.py','v04/models.py','v04/pipelines.py','v03/monitor.py',
            'v03/runner.py','v02/statistics.py','v02/runner.py')

def engine_hash():
    h=hashlib.sha256()
    for name in CORE_FILES:
        path=Path(__file__).parents[1]/name
        if path.exists():h.update(name.encode());h.update(path.read_bytes())
    return h.hexdigest()

def settings(config):
    return {k:config.get(k,v) for k,v in dict(samples=65536,burn_in=8192,window=2048,
                monitor_rows=4096,application_window=8192,min_measurement=8192,confirmation_windows=2).items()}

def model_spec(workload):
    return {k:v for k,v in workload.items() if k not in {'name','seed_group'}}

def calibration_key(config,w):
    # Matched classical/quantum paths share the position and counterfactual
    # quantum-energy diagnostics. Observable coupling matters for energy.
    spec=model_spec(w)
    return digest({'alpha':spec['alpha'],'trial_quartic':spec.get('trial_quartic',0.),
                   'coupling':spec.get('coupling',0.),'input_family':config.get('input_family','pcg64'),
                   'proposal':spec.get('proposal',1.5),'settings':settings(config)})[:20]

def generate_inputs(scenario,seed,rows):
    if 'pipeline' in scenario:
        from rnglab.v04.pipelines import generate_pipeline
        u,meta=generate_pipeline(scenario['pipeline'],seed,rows)
    else:
        family=scenario.get('base_generator','pcg64')
        if family in {'pcg64','philox'}:
            rng=np.random.Generator((np.random.PCG64 if family=='pcg64' else np.random.Philox)(np.random.SeedSequence([*seed,902])))
            u=rng.random((rows,2));meta={'generator':family,'layer':'deterministic_generator',
                'base_uniforms':2*rows,'initial_seed':seed,'raw_physical_measurement':False}
        elif family=='chacha20':
            key=hashlib.sha256(np.asarray(seed,dtype='<u8').tobytes()+b'v04-chacha-key').digest()
            cipher=Cipher(algorithms.ChaCha20(key,bytes(16)),mode=None).encryptor()
            words=np.frombuffer(cipher.update(bytes(16*rows)),dtype='<u8').reshape(rows,2)
            u=(words>>11).astype(float)/2**53
            meta={'generator':family,'layer':'deterministic_generator','key_sha256':hashlib.sha256(key).hexdigest(),
                  'initial_seed':seed,'cipher_bytes':16*rows,'raw_physical_measurement':False,
                  'entropy_scope':'Public deterministic experimental key; not a measured entropy or security claim.'}
        else:raise ValueError('unknown base generator')
    u=np.array(u,dtype=np.float64,copy=True)
    copula=scenario.get('copula')
    if copula:
        rho=float(copula['rho']);onset=int(copula.get('onset',0));end=int(copula.get('end',rows))
        if abs(rho)>=1 or not 0<=onset<=end<=rows:raise ValueError('invalid copula schedule')
        part=u[onset:end];eps=np.finfo(float).eps
        z=ndtri(np.clip(part,eps,1-eps))
        if rho!=0:part[:,1]=np.minimum(np.nextafter(1.,0.),ndtr(rho*z[:,0]+np.sqrt(1-rho*rho)*z[:,1]))
        meta['application_copula']={**copula,'layer':'role_assignment_after_generator',
                                  'uniform_pearson':6/np.pi*np.arcsin(rho/2)}
    if u.shape!=(rows,2) or not np.isfinite(u).all() or np.any(u<0) or np.any(u>=1):
        raise ValueError('pipeline must supply [0,1) Nx2 uniforms')
    return u,meta

def one_trial(config,w,scenario,trial):
    tick=time.perf_counter();s=settings(config);rows=s['samples']+s['burn_in']
    group=w.get('seed_group',w['name'])
    seed=[config['master_seed'],config['seed_namespace'],int.from_bytes(hashlib.sha256(group.encode()).digest()[:4],'big'),trial]
    u,source=generate_inputs(scenario,seed,rows);model=Model(model_spec(w));values=model.advance(u)
    clean={k:np.asarray(v[s['burn_in']:]) for k,v in values.items()}
    metrics={}
    for key,v in clean.items():
        m=correlated_summary(v) if key in {'x','energy','x2','x3','x4'} else {'mean':float(v.mean()),'sd':float(v.std(ddof=1))}
        if key in model.targets:
            m.update(target=model.targets[key],error=m['mean']-model.targets[key])
            if m.get('mcse') is not None:m['intended_target_inside_mcse95']=bool(abs(m['error'])<=1.959963984540054*m['mcse'])
        metrics[key]=m
    monitor=Monitor({'monitor_rows':s['monitor_rows'],'confirmation_windows':s['confirmation_windows']},2)
    windows=[];monitors=[];total=energy_total=0.;retained=0
    calibration=config.get('_calibration',{}).get(calibration_key(config,w))
    for start in range(0,rows,s['window']):
        end=start+s['window'];part=u[start:end]
        old=monitor.update(part,end);latest=u[max(0,end-s['monitor_rows']):end]
        z=ndtri(np.clip(latest,np.finfo(float).eps,1-np.finfo(float).eps))
        cross=float(abs(np.mean(z[:,0]*z[:,1]))*np.sqrt(len(z)))
        lag=float(np.max(np.abs(np.mean(z[1:]*z[:-1],axis=0)))*np.sqrt(len(z)-1))
        monitors.append({'time':end,'marginal':max(old['monobit'],old['distribution']),
                         'normal_cross':cross,'normal_lag':lag,'composite_v03':old['composite']})
        row={'start':start,'time':end,'window_means':{k:float(v[start:end].mean()) for k,v in values.items()}}
        if end>s['burn_in']:
            total+=float(values['x'][start:end].sum());retained+=s['window']
            energy_total+=float(values['energy'][start:end].sum())
            recent=values['x'][max(s['burn_in'],end-s['application_window']):end]
            recent_energy=values['energy'][max(s['burn_in'],end-s['application_window']):end]
            row.update(retained=retained,cumulative_mean=total/retained,recent_mean=float(recent.mean()),
                       energy_cumulative_mean=energy_total/retained,energy_recent_mean=float(recent_energy.mean()))
            if retained>=s['min_measurement']:
                for criterion in ['cumulative','recent','energy_cumulative','energy_recent']:
                    target=model.targets['energy'] if criterion.startswith('energy_') else 0.
                    row[criterion+'_error']=abs(row[criterion+'_mean']-target)
                    if calibration:
                        scale=calibration['scale'][criterion][str(end)]
                        row[criterion+'_score']=row[criterion+'_error']/scale
                        row[criterion+'_tolerance']=scale*calibration['application_cutoffs'][criterion]
        else:row['retained']=0
        windows.append(row)
    onset=scenario.get('onset');recovery=scenario.get('source_recovery')
    events={'detectors':{},'applications':{}}
    if calibration:
        for key in DETECTORS:
            events['detectors'][key]=boundary_events([m[key] for m in monitors],[m['time'] for m in monitors],
                calibration['monitor_cutoffs'][key],s['confirmation_windows'],onset,recovery)
        eligible=[v for v in windows if v['retained']>=s['min_measurement']]
        for key in ['cumulative','recent','energy_cumulative','energy_recent']:
            events['applications'][key]=boundary_events([v[key+'_score'] for v in eligible],
                [v['time'] for v in eligible],calibration['application_cutoffs'][key],s['confirmation_windows'],onset,recovery)
    result={'workload':w['name'],'scenario':scenario['name'],'trial':trial,'seed':seed,'seed_group':group,
            'seed_namespace':config['seed_namespace'],'status':'complete','purpose':config['purpose'],
            'model_spec':model_spec(w),'targets':model.targets,'metrics':metrics,
            'mean':metrics['x']['mean'],'model_target':0.,'samples':s['samples'],'burn_in':s['burn_in'],
            'horizon':rows,'emitted_uniforms':2*rows,'sampler_backend':model.backend,
            'stream_sha256':hashlib.sha256(u.astype('<f8').tobytes()).hexdigest(),
            'position_sha256':hashlib.sha256(values['x'].astype('<f8').tobytes()).hexdigest(),
            'acceptance_sha256':hashlib.sha256(values['accepted'].astype('<f8').tobytes()).hexdigest(),
            'uniform_means':u.mean(axis=0).tolist(),'uniform_variances':u.var(axis=0).tolist(),
            'uniform_cross_correlation':float(np.corrcoef(u.T)[0,1]) if np.std(u,axis=0).min()>0 else None,
            'source_provenance':source,'fault_onset':onset,'source_recovery':recovery,
            'calibration_key':calibration_key(config,w),'elapsed_seconds':time.perf_counter()-tick}
    payload={'result':result,'windows':windows,'monitor':monitors,'events':events}
    return payload,clean['x'] if trial<config.get('trace_trials',0) else None

def applicable(w,scenario):return not scenario.get('workloads') or w['name'] in scenario['workloads']
def job_identity(w,s,t):return f"{w['name']}__{s['name']}__{t:06d}"

def execute(job):
    config,w,scenario,trial,output=job;out=Path(output);name=job_identity(w,scenario,trial)
    lock=out/'trials'/f'{name}.lock';path=out/'trials'/f'{name}.json.gz'
    fd=os.open(lock,os.O_WRONLY|os.O_CREAT|os.O_EXCL)
    with os.fdopen(fd,'w') as f:f.write(canonical({'pid':os.getpid(),'hostname':socket.gethostname(),'started':time.time()}))
    try:
        payload,trace=one_trial(config,w,scenario,trial)
        temp=path.with_name(path.name+f'.{os.getpid()}.tmp')
        with gzip.open(temp,'wt',compresslevel=3) as f:f.write(canonical({'checksum':digest(payload),'payload':payload}))
        os.replace(temp,path)
        if trace is not None:
            (out/'traces').mkdir(exist_ok=True);np.savez_compressed(out/'traces'/f'{name}.npz',x=trace)
        return name,payload['result']['elapsed_seconds']
    except BaseException as exc:
        atomic_json(out/'failures'/f'{name}.{time.time_ns()}.json',{'exception':repr(exc),'traceback':traceback.format_exc()})
        raise
    finally:lock.unlink(missing_ok=True)

def run(config,output,workers=4,resume=False):
    config=dict(config);s=settings(config);out=Path(output);begin=time.perf_counter()
    if config['trials']<2 or config['seed_namespace']<50000:raise ValueError('fresh V04 namespace and at least two trials required')
    if workers<1 or workers>min(4,max(1,len(os.sched_getaffinity(0))//2)):raise ValueError('local worker policy')
    if s['burn_in']%s['window'] or s['samples']%s['window']:raise ValueError('aligned budgets required')
    if any(not isinstance(v,int) or v<1 for k,v in s.items() if k!='burn_in') or s['burn_in']<0:
        raise ValueError('positive integer budgets required')
    if s['samples']<s['application_window'] or s['min_measurement']+(s['confirmation_windows']-1)*s['window']>s['samples']:
        raise ValueError('insufficient observation budget')
    for category in ['workloads','scenarios']:
        names=[x['name'] for x in config[category]]
        if len(set(names))!=len(names) or any(not n.replace('_','').replace('-','').isalnum() for n in names):raise ValueError('unique safe names required')
    for sc in config['scenarios']:
        if 'copula' in sc:
            if sc.get('onset')!=sc['copula'].get('onset',0):raise ValueError('copula and event onset clocks differ')
            if sc['copula'].get('end') is not None and sc.get('source_recovery')!=sc['copula']['end']:
                raise ValueError('copula termination and recovery clocks differ')
    calibration=None
    if config.get('calibration_file'):
        calibration=json.loads((ROOT/config['calibration_file']).read_text())
        if calibration['engine_sha256']!=engine_hash():raise ValueError('calibration engine mismatch')
        if config['seed_namespace'] in calibration['fitting_namespaces']:raise ValueError('cannot reuse fitted streams')
        config['_calibration']=calibration['strata']
        if any(calibration_key(config,w) not in config['_calibration'] for w in config['workloads']):raise ValueError('missing calibrated target')
    public={k:v for k,v in config.items() if k!='_calibration'}
    identity={'config':public,'engine_sha256':engine_hash(),'kernel_sha256':KERNEL_SHA256,'compiler_flags':COMPILER_FLAGS,
              'python':platform.python_version(),'dependencies':{n:importlib.metadata.version(n) for n in ['numpy','scipy','cryptography']},
              'calibration_sha256':digest(calibration) if calibration else None}
    manifest=out/'manifest.json'
    if manifest.exists():
        if not resume or json.loads(manifest.read_text())['identity']!=identity:raise ValueError('existing experiment requires identical resume')
    else:
        if out.exists() and any(out.iterdir()):raise ValueError('nonempty output without manifest')
        out.mkdir(parents=True,exist_ok=True);atomic_json(manifest,{'identity':identity,'created_unix':time.time(),
            'host':platform.platform(),'purpose':config['purpose'],'seed_policy':'Fresh stage namespaces; explicit seed_group pairs models/scenarios.'})
        snapshot=out/'source_snapshot';snapshot.mkdir()
        for name in CORE_FILES:
            path=Path(__file__).parents[1]/name
            if path.exists():shutil_path=snapshot/name;shutil_path.parent.mkdir(parents=True,exist_ok=True);shutil_path.write_bytes(path.read_bytes())
    for name in ['trials','failures']:(out/name).mkdir(exist_ok=True)
    for _ in payloads(out):pass
    jobs=[(config,w,sc,t,str(out)) for w in config['workloads'] for sc in config['scenarios'] if applicable(w,sc)
          for t in range(config['trials']) if not (out/'trials'/f'{job_identity(w,sc,t)}.json.gz').exists()]
    with ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:
        for n,f in enumerate(as_completed([pool.submit(execute,j) for j in jobs]),1):
            name,seconds=f.result()
            with (out/'progress.jsonl').open('a') as file:file.write(canonical({'trial':name,'seconds':seconds,'finished_unix':time.time()})+'\n')
            if n%100==0 or n==len(jobs):print(f'{n}/{len(jobs)} trials; {time.perf_counter()-begin:.1f}s',flush=True)
    execution={'new_jobs':len(jobs),'completed_trials':len(list((out/'trials').glob('*.json.gz'))),
        'configured_trials':sum(applicable(w,sc) for w in config['workloads'] for sc in config['scenarios'])*config['trials'],
        'failed_attempt_records':len(list((out/'failures').glob('*.json'))),'workers':workers,'elapsed_seconds':time.perf_counter()-begin}
    # Every execution has a separate record; completed metadata is not replaced.
    atomic_json(out/f'execution-{time.time_ns()}.json',execution)
    return execution

def main():
    p=argparse.ArgumentParser();p.add_argument('--config',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--resume',action='store_true')
    a=p.parse_args();run(json.loads(a.config.read_text()),a.output,a.workers,a.resume)

if __name__=='__main__':main()
