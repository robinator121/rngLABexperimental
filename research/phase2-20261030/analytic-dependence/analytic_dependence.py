#!/usr/bin/env python3
"""Synthetic stationary harmonic-target dependence control, not RNG degradation.

Only the centroid follows AR(1); internal harmonic modes are independently
redrawn on every step. Uses independent PCG64 namespaces for each replica/role.
"""
from __future__ import annotations

import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse
import csv
import hashlib
import json
import math
import platform
import sys
import time
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import scipy
from scipy.stats import beta as beta_dist, t, norm, skew
import numba
from numba import njit

OUT=Path(__file__).resolve().parent
NAMESPACE='rnglab-phase2-analytic-dependence-v1'
M,BETA=16,2.
REPLICAS,N=1024,8192
PHIS=(0.,.997,.999)
BATCHES=(32,128,512)


def save(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(2**20),b''):h.update(chunk)
    return h.hexdigest()


def seed(condition,replica,role):
    return hashlib.sha256(f'{NAMESPACE}/{condition}/{replica:05d}/{role}'.encode()).hexdigest()


def generator(condition,replica,role):
    return np.random.default_rng(int(seed(condition,replica,role),16))


def parameters(m=M,beta=BETA):
    k=np.arange(m)
    q=beta/m+4*m/beta*np.sin(np.pi*k/m)**2
    weights=1/(m*q[1:])
    vc=.5/beta**2
    vs=float(2*np.sum(weights**2))
    mean=1/(2*beta)+.5/beta+float(weights.sum())
    return {'M':m,'beta':beta,'energy_mean':mean,'centroid_component_variance':vc,
            'shape_mean':float(weights.sum()),'shape_variance':vs,'CV_variance':vc+vs,
            'shape_third_cumulant':float(8*np.sum(weights**3)),
            'shape_weights':weights,'normal_mode_precision':q}


def exact_mc_variance(n,phi,p=None):
    """Finite-n stationary variance, including the (1-k/n) lag weights."""
    if p is None:p=parameters()
    k=np.arange(1,n,dtype=float)
    correlation=np.power(phi*phi,k)
    tr_r2=n+2*np.sum((n-k)*correlation)
    return float(p['shape_variance']/n+p['centroid_component_variance']*tr_r2/n**2)


def exact_mean_skewness(n,phi,p=None):
    if p is None:p=parameters()
    k=np.arange(1,n,dtype=float)
    tr_r3=n+6*np.sum((n-k)*k*np.power(phi*phi,k))
    third=tr_r3/(p['beta']**3*n**3)+p['shape_third_cumulant']/n**2
    return float(third/exact_mc_variance(n,phi,p)**1.5)


def exact_expected_bm_mcvariance(n,batch,phi,p=None):
    assert n%batch==0
    blocks=n//batch
    return (exact_mc_variance(batch,phi,p)-exact_mc_variance(n,phi,p))/(blocks-1)


@njit
def simulate(cg,sg,c0,n,phi,beta,weights):
    c=c0;scale=math.sqrt((1-phi*phi)/beta)
    result=np.empty(n)
    for j in range(n):
        c=phi*c+scale*cg.standard_normal()
        shape=0.
        for w in weights:
            z=sg.standard_normal();shape+=w*z*z
        result[j]=1/(2*beta)+.5*c*c+shape
    return result


def binomial_ci(k,n):
    return [0. if k==0 else float(beta_dist.ppf(.025,k,n-k+1)),
            1. if k==n else float(beta_dist.ppf(.975,k+1,n-k))]


def analyze(y,phi,condition,p):
    r,n=y.shape;means=y.mean(1);target=p['energy_mean'];emp=float(means.var(ddof=1))
    exact=exact_mc_variance(n,phi,p);bias=float(means.mean()-target)
    meancritical=float(t.ppf(.975,r-1));meanse=math.sqrt(emp/r)
    boot=generator(condition,0,'replica_bootstrap').integers(r,size=(2500,r))
    emp_boot=means[boot].var(1,ddof=1)
    base={'condition':condition,'phi':phi,'replicas':r,'retained':n,
          'CV_mean_target':target,'CV_marginal_variance_exact':p['CV_variance'],
          'empirical_variance_replica_means':emp,'exact_variance_replica_means':exact,
          'empirical_to_exact_variance_ratio':emp/exact,
          'empirical_to_exact_variance_ratio_bootstrap_ci95':(np.quantile(emp_boot,[.025,.975])/exact).tolist(),
          'bias':bias,'bias_ci95':[bias-meancritical*meanse,bias+meancritical*meanse],
          'exact_mean_skewness':exact_mean_skewness(n,phi,p),
          'empirical_mean_skewness':float(skew(means,bias=False))}
    rows=[];endpoints=[]
    for batch in BATCHES:
        assert n%batch==0
        blocks=n//batch;bm=y.reshape(r,blocks,batch).mean(2)
        estimated=bm.var(1,ddof=1)/blocks
        critical=float(t.ppf(.975,blocks-1))
        covered=np.abs(means-target)<=critical*np.sqrt(estimated)
        expected=exact_expected_bm_mcvariance(n,batch,phi,p)
        estimate_mean=float(estimated.mean());estimate_se=float(estimated.std(ddof=1)/math.sqrt(r))
        bootratio=estimated[boot].mean(1)/emp_boot
        row={**base,'method':f'BM{batch}','batch':batch,'batch_count':blocks,
             'critical':critical,'coverage_count':int(covered.sum()),'coverage':float(covered.mean()),
             'coverage_ci95':binomial_ci(int(covered.sum()),r),
             'mean_estimated_variance':estimate_mean,'exact_expected_estimated_variance':expected,
             'mean_estimated_variance_ci95':[estimate_mean-meancritical*estimate_se,estimate_mean+meancritical*estimate_se],
             'measured_to_exact_expected_BM_variance_ratio':estimate_mean/expected,
             'BM_expectation_validation_z_score':(estimate_mean-expected)/estimate_se,
             'estimated_to_empirical_variance_ratio':estimate_mean/emp,
             'estimated_to_empirical_variance_ratio_bootstrap_ci95':np.quantile(bootratio,[.025,.975]).tolist(),
             'exact_expected_BM_to_true_MCvariance_ratio':expected/exact}
        rows.append(row)
        for i in range(r):endpoints.append({'condition':condition,'phi':phi,'replica':i,'method':row['method'],
                  'energy_mean':float(means[i]),'mcvariance':float(estimated[i]),'covered':bool(covered[i])})
    critical=float(norm.ppf(.975));oracle=np.abs(means-target)<=critical*math.sqrt(exact)
    rows.append({**base,'method':'oracle_known_variance_normal','critical':critical,
                 'coverage_count':int(oracle.sum()),'coverage':float(oracle.mean()),
                 'coverage_ci95':binomial_ci(int(oracle.sum()),r),
                 'mean_estimated_variance':exact,'estimated_to_empirical_variance_ratio':exact/emp,
                 'scope':'exact finite-n variance; normal interval still has distributional approximation'})
    for i in range(r):endpoints.append({'condition':condition,'phi':phi,'replica':i,'method':'oracle_known_variance_normal',
                  'energy_mean':float(means[i]),'mcvariance':exact,'covered':bool(oracle[i])})
    return rows,endpoints


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',action='store_true');a=parser.parse_args()
    if not a.run:parser.error('Explicit --run required; results are preserved rather than overwritten')
    if (OUT/'results.json').exists():raise FileExistsError('Preserve existing results')
    started=time.monotonic();p=parameters();rows=[];endpoints=[];manifests=[]
    save(OUT/'protocol.json',{'id':NAMESPACE,'classification':'prespecified analytically solvable synthetic dependence control',
        'written_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':sha(__file__),
        'replicas':REPLICAS,'retained':N,'phis':PHIS,'batch_sizes':BATCHES,'M':M,'beta':BETA,
        'hypothesis':'BM estimated MC variance follows its exact stationary finite-n expectation; exact variance normal intervals distinguish variance estimation error from remaining finite-run distributional effects',
        'not_an_RNG_fault':'AR1 is an intentionally specified target-preserving transition; healthy PCG64 normal inputs',
        'initialization':'c0~N(0,1/beta), independent internal harmonic modes redrawn at every recorded step',
        'replicate_unit':'one complete independent trajectory; distinct namespaces across phi and replica',
        'resource_limit':'single process, BLAS one thread; at most one R1024 by n8192 CV array retained in memory',
        'primary_endpoints':['BM32/BM128/BM512 coverage','BM expected variance against analytic exact value','oracle known finite-n variance normal-interval coverage'],
        'sample_scope':'1024 replicas estimates 95% coverage with roughly ±1.34 percentage point precision; no guarantee at exactly 95%'})
    # Compile before timing per-condition generation.
    simulate(np.random.default_rng(1),np.random.default_rng(2),0.,2,.5,BETA,p['shape_weights'])
    for phi in PHIS:
        condition='phi'+str(phi).replace('.','p');s=time.monotonic();seeds=[];y=np.empty((REPLICAS,N))
        for i in range(REPLICAS):
            c0=generator(condition,i,'initial').standard_normal()/math.sqrt(BETA)
            y[i]=simulate(generator(condition,i,'centroid'),generator(condition,i,'shape'),c0,N,phi,BETA,p['shape_weights'])
            seeds.append({'replica':i,**{role:seed(condition,i,role) for role in ['initial','centroid','shape']}})
        path=OUT/f'{condition}-CV-traces.npz';np.savez_compressed(path,CV=y)
        result,individual=analyze(y,phi,condition,p);rows+=result;endpoints+=individual
        manifest={'condition':condition,'phi':phi,'timestamp_utc':datetime.now(timezone.utc).isoformat(),
                  'reference_git_revision':'dd96adba23b8100aba9147d71fffc63932fc6b0d',
                  'script_sha256':sha(__file__),'parameters':{'R':REPLICAS,'n':N,'M':M,'beta':BETA},
                  'seeds':seeds,'seed_rule':'SHA256(namespace/condition/replica:05d/role), full integer -> NumPy default_rng PCG64',
                  'software':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'numba':numba.__version__},
                  'platform':platform.platform(),'trace_file':path.name,'trace_sha256':sha(path),
                  'elapsed_seconds':time.monotonic()-s,'validation_status':'generated and analytically compared; inspect stated uncertainty'}
        save(OUT/f'{condition}-manifest.json',manifest);manifests.append(manifest)
        for row in result:print(condition,row['method'],row['coverage_count'],'/',REPLICAS,
                'variance ratio',row['estimated_to_empirical_variance_ratio'],flush=True)
        del y
    save(OUT/'results.json',{'protocol':NAMESPACE,'rows':rows,'elapsed_seconds':time.monotonic()-started,
         'reference':{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in p.items()}})
    with (OUT/'replica_endpoints.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(endpoints[0]));w.writeheader();w.writerows(endpoints)
    flat=[{k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in row.items()} for row in rows]
    with (OUT/'summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in flat for k in row)));w.writeheader();w.writerows(flat)
    save(OUT/'manifest.json',{'protocol':NAMESPACE,'written_utc':datetime.now(timezone.utc).isoformat(),
        'script_sha256':sha(__file__),'total_elapsed_seconds':time.monotonic()-started,
        'files':{p.name:sha(p) for p in sorted(OUT.iterdir()) if p.is_file() and p.name!='manifest.json'}})


if __name__=='__main__':main()
