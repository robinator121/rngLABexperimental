from pathlib import Path
import hashlib,json,os,time
ROOT=Path(__file__).resolve().parents[1]

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb')as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def atomic(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_name(p.name+'.tmp')
    with t.open('w')as f:json.dump(value,f,sort_keys=True,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    os.replace(t,p)

def immutable(p,value):
    p=Path(p)
    if p.exists():
        if json.loads(p.read_text())!=value:raise ValueError('Refuse immutable overwrite:'+str(p))
    else:atomic(p,value)

def rng(namespace,index,channel=0):
    import numpy as np
    return np.random.Generator(np.random.PCG64DXSM(np.random.SeedSequence([2026101410,namespace,index,channel])))

def source_manifest():
    paths=[p for folder in ['rnglab10','scripts','tests','vendor']for p in(ROOT/folder).rglob('*')if p.is_file()and'__pycache__'not in p.parts]
    return {str(p.relative_to(ROOT)):sha(p)for p in sorted(paths)}
