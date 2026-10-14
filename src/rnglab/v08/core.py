"""One identity shares clean latent streams across strengths, never across cohorts."""
import hashlib
import json
import os
from pathlib import Path
import numpy as np
from scipy.special import ndtri
from rnglab.v02.statistics import correlated_summary
from rnglab.v03.physics import PIMC
from rnglab.v03.uniforms import UniformProcess
from rnglab.v07.calibration import event_time

COLUMNS = ('mean','variance','ess','mcse_low','mcse_high','mixing_warning',
           'endpoint_failure','nominal95_coverage','acceptance','cv_energy','th_energy',
           'cross_time','union_time','application_time','physical_time',
           'post_physical_time','cross_reset_time','union_reset_time','app_reset_time')
COL = {name:i for i,name in enumerate(COLUMNS)}
FROZEN_SHA = '4b0ec747573c1ff396c8e719bf4739d10359d1b78bec5635cf1cefd8383170b7'

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def atomic(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    temp.write_text(json.dumps(value,sort_keys=True,allow_nan=False)+'\n')
    os.replace(temp,path)

def frozen(root):
    path=Path(root)/'frozen/v07-cutoffs.json'
    if sha(path)!=FROZEN_SHA:raise ValueError('Archived v07 threshold checksum mismatch')
    return json.loads(path.read_text())['cutoffs']

def inputs(config,namespace,index,rho):
    seed=[config['master_seed'],namespace,index]
    spec={'kind':'cross_correlation','rho':rho,'pair':['proposal','acceptance'],
          'schedule':{'onset':config['onset'],'end':config['end']}}
    u=UniformProcess(spec,seed,3).sample(config['burn_in']+config['samples'])
    init=np.random.Generator(np.random.PCG64(np.random.SeedSequence([*seed,991]))).random(config['beads'])
    return u,init

def marginal_score(u):
    # Exactly the two v07 components used in its frozen union, no unused lag work.
    n=len(u);a=np.asfortranarray(u)
    ones=np.bitwise_count(np.floor(a*4294967296.).astype(np.uint32)).sum(axis=0,dtype=np.int64)
    mono=float(np.max(abs(2*ones-32*n))/np.sqrt(32*n))
    bins=np.minimum((a*16).astype(np.int64),15)
    counts=np.array([np.bincount(bins[:,j],minlength=16) for j in range(2)])
    dist=float(max(0.,np.max(((counts-n/16)**2/(n/16)).sum(axis=1)-15)/np.sqrt(30)))
    return max(mono,dist)

def events(looks,thresholds,onset):
    ins=[r for r in looks if 'cross'in r];apps=[r for r in looks if 'application'in r]
    def first(rows,key,threshold,reset=None):
        return event_time([r[key]for r in rows],[r['time']for r in rows],threshold,onset=reset)
    def union(reset=None):
        times=[first(ins,k,thresholds['source_union'],reset)for k in ['cross','marginal']]
        return min([v for v in times if v is not None],default=None)
    return {'cross_time':first(ins,'cross',thresholds['cross']),
        'union_time':union(),'application_time':first(apps,'application',thresholds['application']),
        'physical_time':first(apps,'physical_error',thresholds['physical']),
        'post_physical_time':first(apps,'physical_error',thresholds['physical'],onset),
        'cross_reset_time':first(ins,'cross',thresholds['cross'],onset),
        'union_reset_time':union(onset),'app_reset_time':first(apps,'application',thresholds['application'],onset)}

def simulate(config,namespace,index,rho,thresholds,backend='compiled',detail=False):
    u,init=inputs(config,namespace,index,rho)
    model=PIMC(beta=config['beta'],beads=config['beads'],sampler='normal_mode',
        mode_selection=config['mode_selection'],mode_proposal=config['mode_proposal'],backend=backend)
    model.initialize_from_uniforms(init);energy,aux=model.advance(u)
    x=aux['centroid'];csum=np.r_[0.,np.cumsum(x)]
    z=ndtri(np.clip(u[:,1:],2.**-54,np.nextafter(1.,0.)))
    psum=np.r_[0.,np.cumsum(z[:,0]*z[:,1])];h=len(u);looks=[]
    for end in range(config['look_rows'],h+1,config['look_rows']):
        r={'time':end}
        if end>=config['input_window_rows']:
            begin=end-config['input_window_rows']
            r.update(cross=float(abs(psum[end]-psum[begin])/np.sqrt(config['input_window_rows'])),
                marginal=marginal_score(u[begin:end,1:]))
        retained=end-config['burn_in']
        if retained>=config['application_min_retained']:
            mean=float((csum[end]-csum[config['burn_in']])/retained);error=abs(mean)
            r.update(retained=retained,centroid_mean=mean,physical_error=error,
                application=float(np.sqrt(retained)*error/config['healthy_centroid_sd'])if error>config['physical_tolerance']else 0.)
        looks.append(r)
    s=correlated_summary(x[config['burn_in']:])
    if s['ess'] is None:raise ArithmeticError('Degenerate chain; preserve failed identity')
    times=events(looks,{**thresholds,'physical':config['physical_tolerance']},config['onset'])
    b=config['burn_in']
    row=[s['mean'],s['variance'],s['ess'],s['mcse_ci_low'],s['mcse_ci_high'],int(s['mixing_warning']),
        int(abs(s['mean'])>config['physical_tolerance']),int(s['mcse_ci_low']<=0<=s['mcse_ci_high']),
        float(aux['accepted'][b:].mean()),float(energy[b:].mean()),float(aux['thermodynamic'][b:].mean())]
    row.extend(-1 if times[k]is None else times[k]for k in COLUMNS[11:])
    row=np.asarray(row,float)
    if not np.isfinite(row).all():raise ArithmeticError('Nonfinite simulation output')
    if detail:return row,{'looks':looks,'input_sha256':hashlib.sha256(u.astype('<f8').tobytes()).hexdigest(),
        'centroid_sha256':hashlib.sha256(x.astype('<f8').tobytes()).hexdigest(),'backend':model.backend}
    return row

def identity(config,namespace,index,thresholds):
    return np.stack([simulate(config,namespace,index,r,thresholds)for r in config['rhos']])
