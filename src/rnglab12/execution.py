"""Finite, task-checkpointed research execution; no service or provider control."""
from pathlib import Path
from datetime import datetime, timezone
import csv, hashlib, json, math, os, time
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import numpy as np
from rnglab11.applications import make_application
from rnglab11.sources import SyntheticSource, UniformStream
from .seeds import derive_seed

ROOT = Path(__file__).resolve().parents[2]
APPLICATIONS = ('integration','pi','mcmc','ising','vmc','harmonic_pimc')
DESIGN = {'steps':8192,'burn_in':1024,'rolling':512,'stride':128,
          'initialization':'deterministic legacy constructor; fixed burn-in',
          'source':'task-specified source; default MT19937 noncryptographic synthetic control',
          'uniform_read_rows':128,'units':'attempted proposals; absolute zero-based fault rows',
          'application_parameters':'legacy defaults; Ising L3, harmonic finite M16 beta2'}
def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',',':'), allow_nan=False).encode()
def digest(value): return hashlib.sha256(canonical(value)).hexdigest()
CONFIG_ID = digest(DESIGN)
def now(): return datetime.now(timezone.utc).isoformat()
def write_new(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('xb') as f:
        f.write(json.dumps(value,indent=2,allow_nan=False).encode()+b'\n');f.flush();os.fsync(f.fileno())
def source_identity():
    files=sorted((ROOT/'src').rglob('*.py'))+sorted((ROOT/'scripts').glob('*.py'))+[ROOT/'pyproject.toml']
    mapping={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    return {'sha256':digest(mapping),'files':mapping}
def task_id(task):return digest(task)[:24]

def simulate(task):
    from .uncertainty import summarize_series
    started=time.monotonic()
    app=task['application']; design=task.get('design',DESIGN)
    steps=design['steps']; burn=design['burn_in']; window=design['rolling'];stride=design['stride']
    if not 1024<=steps<=16384 or not 0<=burn<steps-window or window%stride:
        raise ValueError('Bounded valid design required')
    seed=derive_seed(task['split'],task.get('seed_domain','healthy-v012')+'/'+app,task['group'])
    model=make_application({'kind':app}); dim=model.dimension
    source_spec=task.get('source_spec')
    if source_spec is None:
        source=SyntheticSource(seed,max_bytes=1048576)
        source_family='synthetic_fault'; source_label='MT19937 healthy noncryptographic control'
    else:
        from .sources import make_source,source_type
        seeded=source_type(source_spec) in {'deterministic_seeded_or_reseeded_csprng','synthetic_fault'}
        if task.get('reseed_material_fault'):
            from .sources import FixedRequestDRBG,ClassifiedSource
            from .faults import make_byte_fault
            if source_type(source_spec)!='deterministic_seeded_or_reseeded_csprng':raise ValueError('Research reseed control needs deterministic family')
            initial=SyntheticSource(seed,max_bytes=1048576)
            reseed_seed=derive_seed(task['split'],task.get('seed_domain','healthy-v012')+'/'+app,task['group'],stream='reseed')
            reseed=make_byte_fault(SyntheticSource(reseed_seed,max_bytes=1048576),task['reseed_material_fault'])
            source=ClassifiedSource(FixedRequestDRBG(initial,reseed_source=reseed,interval=2,max_bytes=1048576),
                'deterministic_seeded_or_reseeded_csprng',details={'intervention':'reseed material; initial seed intact'},max_bytes=1048576)
        else:
            extra={'seed':seed} if seeded else {}
            if source_spec.get('reseed_interval_generate_calls'):
                extra['reseed_seed']=derive_seed(task['split'],task.get('seed_domain','healthy-v012')+'/'+app,task['group'],stream='reseed')
            source=make_source({**source_spec,**extra},max_bytes=1048576)
        source_family=source_type(source_spec);source_label=source_family
    with source:
        stream=UniformStream(source)
        u=np.concatenate([stream.random(min(128,steps-i),dim) for i in range(0,steps,128)])
        provenance=source.provenance()
    fault=task.get('fault')
    if fault:
        from .faults import apply_uniform_fault
        u=apply_uniform_fault(u,0,fault,app)
    if task.get('linear_scale'):
        control=task['linear_scale']; rows=np.arange(steps)
        severity=np.clip((rows-control['start_row'])/control['duration_rows'],0,1)*control['final_severity']
        u=u*(1-severity[:,None])
    if task.get('positive_acceptance_control'):
        channel={'mcmc':1,'vmc':1,'harmonic_pimc':2}[app]
        first=fault['start_row']; last=first+fault['duration_rows']
        u[first:last,channel]=np.maximum(u[first:last,channel],2**-53)
    reset_at=task.get('reset_at')
    if reset_at is None:
        values,aux=model.advance(u)
    else:
        first,aux=model.advance(u[:reset_at]); model=make_application({'kind':app})
        second,_=model.advance(u[reset_at:]); values=np.concatenate([first,second])
    reference=float(model.target) if model.target is not None else None
    kept=np.asarray(values[burn:],dtype=float)
    if not np.isfinite(kept).all():raise ValueError('Nonfinite trajectory; record failed attempt explicitly')
    sums=np.concatenate([[0.],np.cumsum(values,dtype=float)])
    times=np.arange(burn+window,steps+1,stride,dtype=int)
    rolling=(sums[times]-sums[times-window])/window
    errors=np.abs(rolling-reference) if reference is not None else None
    # Delivered-MSB APT-style score on the same frozen look grid; not NIST certification.
    signs=(u>=.5).astype(float);cs=np.vstack([np.zeros(dim),np.cumsum(signs,axis=0)])
    proportions=(cs[times]-cs[times-window])/window
    source_scores=np.max(np.abs(proportions-.5),axis=1)*2*math.sqrt(window)
    raw_bits=(u[:,0]>=.5).astype(np.uint8)
    boundaries=np.r_[0,np.flatnonzero(raw_bits[1:]!=raw_bits[:-1])+1,len(raw_bits)]
    max_run=int(np.diff(boundaries).max())
    # Whole-horizon ideal fair-bit union bound; valid only under that reference law.
    rct_cutoff=math.ceil(math.log2(steps/.01))+1
    rct_event=None
    for a,b in zip(boundaries[:-1],boundaries[1:]):
        if b-a>=rct_cutoff:rct_event=int(a+rct_cutoff);break
    result={'task':task,'task_id':task_id(task),'seed':seed,
        'lineage_id':f"{task['split']}/{task.get('seed_domain','healthy-v012')}/{app}/{task['group']}",
        'application':app,'observable':{'integration':'x2','pi':'pi','mcmc':'x2','ising':'energy','vmc':'energy','harmonic_pimc':'centroid_virial'}[app],
        'reference':reference,'config_id':CONFIG_ID if design==DESIGN else digest(design),
        'times':times.tolist(),'errors':None if errors is None else errors.tolist(),
        'rolling_estimates':rolling.tolist(),'cumulative_estimates':((sums[times]-sums[burn])/(times-burn)).tolist(),
        'source_scores':source_scores.tolist(),'rct_max_run':max_run,'rct_reference_cutoff':rct_cutoff,'rct_event':rct_event,
        'rct_scope':'First delivered-uniform MSB, not physical raw bits; IID fair-bit union-bound reference only',
        'source_family':source_family,'source_label':source_label,'source_provenance':provenance,
        'uncertainty':summarize_series(kept,reference=reference),
        'elapsed_seconds':time.monotonic()-started,'completed_at_utc':now(),
        'status':'COMPLETE','reset_at':reset_at}
    if task.get('issued_at') is not None:
        issued=task['issued_at']
        if issued<2*window or issued>steps or issued not in times:raise ValueError('Issuance needs frozen past-only windows')
        past=u[:issued,0]
        result['issuance']={'at':issued,'features_through':issued,
            'features':[float(past[-window:].mean()),float(past[-window:].mean()-past[-2*window:-window].mean()),
                        float(errors[np.where(times==issued)[0][0]])],
            'feature_names':['past_uniform_mean','past_uniform_mean_change','current_application_absolute_error']}
    return result

def run_tasks(directory,tasks,*,workers=4,classification='PILOT',protocol=None):
    """An immutable result per task; retry after interruption recomputes only missing tasks."""
    import fcntl
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True)
    if not 1<=workers<=4:raise ValueError('This release execution is bounded to four workers')
    if len(tasks)>8000:raise ValueError('Large batches need a separate measured design')
    code=source_identity()
    manifest={'classification':classification,'tasks':tasks,'protocol':protocol,'source_identity':code,
              'workers':workers,'resumption':'task boundary; interrupted chains recomputed, not resumed in-memory'}
    lock=(directory/'.runner.lock').open('a')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    path=directory/'manifest.json'
    if path.exists():
        old=json.loads(path.read_text())
        if old!=manifest:raise ValueError('Frozen manifest/code mismatch; select a new run directory')
    else:write_new(path,manifest)
    pending=[];records=[]
    if (directory/'receipt.json').exists():
        previous=json.loads((directory/'receipt.json').read_text())
        for filename,expected in previous['result_hashes'].items():
            p=directory/'results'/filename
            if not p.exists() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
                raise ValueError('Preserved result hash mismatch; refusing to rebless modified evidence')
    for task in tasks:
        path=directory/'results'/(task_id(task)+'.json')
        if path.exists():
            r=json.loads(path.read_text())
            if r['task']!=task or r['status']!='COMPLETE':raise ValueError('Result/task mismatch or prior failure')
            records.append(r)
        else:pending.append(task)
    context=multiprocessing.get_context('spawn')
    start=time.monotonic()
    with ProcessPoolExecutor(max_workers=workers,mp_context=context) as pool:
        futures={pool.submit(simulate,t):t for t in pending}
        for future in as_completed(futures):
            task=futures[future]
            try:r=future.result()
            except Exception as exc:
                write_new(directory/'failures'/(task_id(task)+'.json'),{'task':task,'error':repr(exc),'at':now()})
                raise
            r['implementation_sha256']=code['sha256']
            write_new(directory/'results'/(task_id(task)+'.json'),r);records.append(r)
            if len(records)%128==0:print(f'{classification}: {len(records)}/{len(tasks)} complete',flush=True)
    if source_identity()!=code:raise RuntimeError('Source changed during execution; evidence retained but not confirmatory')
    receipt={'status':'COMPLETE','classification':classification,'records':len(records),'elapsed_seconds_this_invocation':time.monotonic()-start,
             'code_sha256':code['sha256'],'at':now(),'result_hashes':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((directory/'results').glob('*.json'))}}
    name='receipt.json' if not (directory/'receipt.json').exists() else f'resume-receipt-{time.time_ns()}.json'
    write_new(directory/name,receipt);lock.close()
    return records
