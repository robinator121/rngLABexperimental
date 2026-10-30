"""RNGLab Phase II harmonic PIMC study. NumPy/SciPy/Numba; no fastmath.

Each statistical unit is a complete independently seeded trajectory. Conditions
within a stage share per-replica seeds intentionally for paired comparisons.
Stages have disjoint SHA-256 namespaces. Original Phase I artifacts are read-only.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import scipy
from scipy.stats import beta as beta_dist, t as student_t, norm
import numba
from numba import njit

ROOT = Path(__file__).resolve().parent
REVISION = 'dd96adba23b8100aba9147d71fffc63932fc6b0d'
M, BETA = 16, 2.0
OBSERVABLES = ('TD', 'CV', 'centroid', 'centroid_component', 'shape')


def save(path, value):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(2**20), b''):
            h.update(chunk)
    return h.hexdigest()


def seed(stage, role, replica):
    return hashlib.sha256(f'rnglab-phase2-v1/{stage}/{role}/{replica:05d}'.encode()).hexdigest()


def generator(stage, role, replica):
    return np.random.default_rng(int(seed(stage, role, replica), 16))


def reference(m=M, beta=BETA):
    eye = np.eye(m)
    difference = np.roll(eye, -1, axis=0) - eye
    laplacian = difference.T @ difference
    q = beta/m*eye + m/beta*laplacian
    covariance = np.linalg.inv(q)
    a_td = eye/(2*m) - m/(2*beta**2)*laplacian
    a_cv = eye/m - np.ones((m,m))/(2*m*m)
    a_shape = eye/m - np.ones((m,m))/(m*m)
    a_centroid = np.ones((m,m))/(2*m*m)
    z = beta/(2*m)
    finite = .5/(math.tanh(m*math.asinh(z))*math.sqrt(1+z*z))
    matrices = {'TD':a_td,'CV':a_cv,'centroid_component':a_centroid,'shape':a_shape}
    means = {'TD':m/(2*beta)+float(np.trace(a_td@covariance)),
             'CV':1/(2*beta)+float(np.trace(a_cv@covariance)),
             'centroid':0.,'centroid_component':.5/beta,
             'shape':float(np.trace(a_shape@covariance))}
    variances = {name:float(2*np.trace(a@covariance@a@covariance)) for name,a in matrices.items()}
    variances['centroid'] = 1/beta
    k = np.arange(m)
    eigenvalues = beta/m + 4*m/beta*np.sin(np.pi*k/m)**2
    return {'M':m,'beta':beta,'finite_energy':finite,'continuum_energy':.5/math.tanh(beta/2),
            'means':means,'variances':variances,'centroid_marginal_CV_variance_fraction':variances['centroid_component']/variances['CV'],
            'precision':q,'covariance':covariance,'cholesky':np.linalg.cholesky(covariance),
            'normal_mode_precision':eigenvalues}


@njit
def local_kernel(g, centroid_g, scan_g, initial, steps, beta, delta, refresh=False,
                 reverse_mix=False, fault=0, onset=0, severity=.08, source_record=False):
    """One proposal and acceptance input per bead; explicit optional centroid heat bath.

    Fault 1 shifts acceptance uniform, fault 2 couples acceptance to proposal,
    fault 3 shifts proposal uniform. All interventions affect an offline copy.
    """
    m = len(initial)
    x = initial.copy()
    dt = beta/m
    raw = np.empty((steps,3))  # TD, centroid, radius of gyration squared
    sources = np.empty((steps,3)) if source_record else np.empty((0,3))
    accepted = 0
    for sweep in range(steps):
        backwards = reverse_mix and scan_g.random() < .5
        source_p, source_a, source_cross = 0., 0., 0.
        for position in range(m):
            j = m-1-position if backwards else position
            up, ua = g.random(), g.random()
            if sweep >= onset:
                if fault == 1:
                    ua = min(ua+severity, np.nextafter(1.,0.))
                elif fault == 2:
                    ua = up
                elif fault == 3:
                    up = min(up+severity, np.nextafter(1.,0.))
            old = x[j]
            new = old + delta*(2*up-1)
            left, right = x[(j-1)%m],x[(j+1)%m]
            ds = ((new-left)**2+(right-new)**2-(old-left)**2-(right-old)**2)/(2*dt)+dt*(new*new-old*old)/2
            if ds <= 0 or ua < math.exp(-ds):
                x[j] = new
                accepted += 1
            source_p += up
            source_a += ua
            source_cross += (up-.5)*(ua-.5)
        if refresh:
            old_centroid = np.mean(x)
            new_centroid = centroid_g.standard_normal()/math.sqrt(beta)
            x += new_centroid-old_centroid
        c = np.mean(x)
        shape = np.mean((x-c)**2)
        springs = 0.
        for j in range(m):
            springs += (x[(j+1)%m]-x[j])**2
        raw[sweep,0] = m/(2*beta)-m*springs/(2*beta*beta)+.5*(shape+c*c)
        raw[sweep,1] = c
        raw[sweep,2] = shape
        if source_record:
            sources[sweep,0] = source_p/m
            sources[sweep,1] = source_a/m
            sources[sweep,2] = source_cross/m
    return raw, sources, accepted/(steps*m)


def numpy_kernel(g, cg, sg, initial, steps, beta, delta, refresh=False, reverse_mix=False):
    """Independent slow kernel for exact finite-seed regression comparison."""
    x = initial.copy(); m=len(x); dt=beta/m; rows=[]; accepted=0
    for _ in range(steps):
        order = range(m-1,-1,-1) if reverse_mix and sg.random()<.5 else range(m)
        for j in order:
            up,ua=g.random(),g.random();old=x[j];new=old+delta*(2*up-1)
            left,right=x[(j-1)%m],x[(j+1)%m]
            ds=((new-left)**2+(right-new)**2-(old-left)**2-(right-old)**2)/(2*dt)+dt*(new*new-old*old)/2
            if ds<=0 or ua<math.exp(-ds):x[j]=new;accepted+=1
        if refresh:x+=cg.standard_normal()/math.sqrt(beta)-x.mean()
        c=x.mean();shape=np.mean((x-c)**2);springs=np.sum((np.roll(x,-1)-x)**2)
        rows.append([m/(2*beta)-m*springs/(2*beta*beta)+.5*(shape+c*c),c,shape])
    return np.asarray(rows),accepted/(steps*m)


def observables(raw):
    td,c,shape = raw[...,0],raw[...,1],raw[...,2]
    centroid_component=.5*c*c
    return {'TD':td,'CV':1/(2*BETA)+centroid_component+shape,
            'centroid':c,'centroid_component':centroid_component,'shape':shape}


def run_condition(directory, stage, config):
    directory=Path(directory); directory.mkdir(parents=True,exist_ok=True)
    path=directory/f"{config['id']}.npz"
    if path.exists():
        raise FileExistsError(f'Preserve prior outputs: {path} already exists')
    started=time.monotonic(); r=config['replicas'];steps=config['steps'];ref=reference()
    data=np.empty((r,steps,3)); sources=np.empty((r,steps,3)) if config.get('source_record') else np.empty((0,0,0));accept=[];seeds=[]
    for i in range(r):
        gi=generator(stage,'initial',i);gu=generator(stage,'update',i);gc=generator(stage,'centroid',i);gs=generator(stage,'scan',i)
        if config['initial']=='stationary':initial=gi.standard_normal(M)@ref['cholesky'].T
        elif config['initial']=='uniform':initial=gi.random(M)-.5
        elif config['initial']=='zero':initial=np.zeros(M)
        elif config['initial']=='displaced':initial=np.full(M,3.)
        else:raise ValueError(config['initial'])
        if config.get('iid'):
            paths=generator(stage,'independent_paths',i).standard_normal((steps,M))@ref['cholesky'].T
            c=paths.mean(1);shape=((paths-c[:,None])**2).mean(1)
            td=M/(2*BETA)-M*np.sum((np.roll(paths,-1,axis=1)-paths)**2,axis=1)/(2*BETA**2)+.5*(shape+c*c)
            data[i]=np.column_stack([td,c,shape]);acc=1.
        else:
            raw,src,acc=local_kernel(gu,gc,gs,initial,steps,BETA,config.get('delta',.4),
                config.get('refresh',False),config.get('reverse_mix',False),config.get('fault',0),
                config.get('onset',0),config.get('severity',.08),config.get('source_record',False))
            data[i]=raw
            if config.get('source_record'):sources[i]=src
        accept.append(acc)
        seeds.append({'replica':i,**{role:seed(stage,role,i) for role in ['initial','update','centroid','scan','independent_paths']}})
    np.savez_compressed(path, raw=data, sources=sources, acceptance=np.asarray(accept))
    manifest={'experiment_id':config['id'],'stage':stage,'timestamp_utc':datetime.now(timezone.utc).isoformat(),
              'git_revision_reference':REVISION,'script_sha256':digest(__file__),'parameters':config,
              'seed_rule':'SHA256(rnglab-phase2-v1/stage/role/replica:05d), full 256bit integer seed for NumPy PCG64',
              'pairing':'same replica/role seeds deliberately reused across conditions within stage; stages disjoint',
              'independent_replicate_unit':'whole trajectory, not window or sweep','seeds':seeds,
              'versions':{'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,'numba':numba.__version__},
              'platform':platform.platform(),'output':path.name,'sha256':digest(path),'elapsed_seconds':time.monotonic()-started,
              'validation_status':'generated; analyze separately; no stationarity inference from ESS'}
    save(directory/f"{config['id']}.manifest.json",manifest)
    print(f"generated {stage}/{config['id']} R={r} n={steps} in {manifest['elapsed_seconds']:.2f}s",flush=True)
    return path


def binomial_interval(k,n,alpha=.05):
    return [0. if k==0 else float(beta_dist.ppf(alpha/2,k,n-k+1)),
            1. if k==n else float(beta_dist.ppf(1-alpha/2,k+1,n-k))]


def fft_autocovariance(y):
    centered=y-y.mean(1,keepdims=True);n=y.shape[1];fft_n=1<<(2*n-1).bit_length()
    f=np.fft.rfft(centered,n=fft_n,axis=1)
    return np.fft.irfft(f*f.conj(),n=fft_n,axis=1)[:,:n]/n


def geyer_variance(gamma,n):
    pairs=gamma[:,:2*(n//2)].reshape(len(gamma),-1,2).sum(2)
    estimates=[];stops=[]
    for g,p in zip(gamma,pairs):
        nonpositive=np.flatnonzero(p<=0);stop=int(nonpositive[0]) if len(nonpositive) else len(p)
        if stop==0:estimates.append(np.nan);stops.append(0);continue
        lrv=-g[0]+2*np.minimum.accumulate(p[:stop]).sum()
        estimates.append(lrv/n if lrv>0 else np.nan);stops.append(stop)
    return np.asarray(estimates),np.asarray(stops)


def variance_methods(y,batches=(32,128,512,1024,2048)):
    r,n=y.shape;means=y.mean(1);v=y.var(1,ddof=1);result={'naive_IID':(v/n,float(student_t.ppf(.975,n-1)),{})}
    for b in batches:
        blocks=n//b
        if blocks>=8:
            bm=y[:,:blocks*b].reshape(r,blocks,b).mean(2)
            result[f'BM{b}']=(bm.var(1,ddof=1)/blocks,float(student_t.ppf(.975,blocks-1)),{'batches':blocks,'batch':b})
    gamma=fft_autocovariance(y)
    gv,stops=geyer_variance(gamma,n)
    result['Geyer_IMS_sensitivity']=(gv,float(norm.ppf(.975)),{'reversibility_required_for_standard_theory':True,'median_positive_pairs':float(np.median(stops))})
    b=min(512,n//8)
    if b>=2:
        cs=np.concatenate([np.zeros((r,1)),np.cumsum(y,axis=1)],axis=1)
        block=(cs[:,b:]-cs[:,:-b])/b
        obm=b/((n-b)*(n-b+1))*np.sum((block-means[:,None])**2,axis=1)
        result[f'OBM{b}']=(obm,float(norm.ppf(.975)),{'batch':b,'interval':'asymptotic normal'})
        bart=(gamma[:,0]+2*np.sum(gamma[:,1:b]*(1-np.arange(1,b)/b),axis=1))/n
        result[f'Bartlett{b}']=(bart,float(norm.ppf(.975)),{'bandwidth':b,'interval':'asymptotic normal'})
    return result


def analyze_window(path,burn,n,label,all_observables=False):
    loaded=np.load(path);raw=loaded['raw'][:,burn:burn+n];obs=observables(raw);r=len(raw);ref=reference();rows=[];replicate_rows=[]
    if raw.shape[1]!=n:raise ValueError('window extends beyond trace')
    for name in (OBSERVABLES if all_observables else ('TD','CV')):
        y=obs[name];target=ref['means'][name];means=y.mean(1);emp=float(means.var(ddof=1));bias=float(means.mean()-target)
        bias_se=math.sqrt(emp/r);bias_ci=[bias-float(student_t.ppf(.975,r-1))*bias_se,bias+float(student_t.ppf(.975,r-1))*bias_se]
        boot=generator('analysis',f'{path.stem}/{label}/{name}/bootstrap',0).integers(r,size=(1500,r))
        empirical_boot=means[boot].var(1,ddof=1)
        for method,(mcvar,critical,extra) in variance_methods(y).items():
            valid=np.isfinite(mcvar)&(mcvar>0);covered=valid&(np.abs(means-target)<=critical*np.sqrt(np.maximum(mcvar,0)))
            k=int(covered.sum());ratio=float(np.mean(mcvar[valid])/emp) if emp>0 else None
            ratio_boot=np.nanmean(mcvar[boot],axis=1)/empirical_boot
            row={'condition':path.stem,'window':label,'burn':burn,'retained':n,'observable':name,'method':method,
                'replicas':r,'covered':k,'coverage':k/r,'coverage_ci95':binomial_interval(k,r),'abstentions':int((~valid).sum()),
                'mean':float(means.mean()),'target':target,'bias':bias,'bias_ci95':bias_ci,'empirical_variance_chain_means':emp,
                'mean_estimated_variance':float(np.mean(mcvar[valid])),'estimated_to_empirical_variance_ratio':ratio,
                'variance_ratio_bootstrap_ci95':np.quantile(ratio_boot,[.025,.975]).tolist(),
                'median_estimated_tau':float(np.median(n*mcvar[valid]/y.var(1,ddof=1)[valid])),
                'median_estimated_ESS':float(np.median(y.var(1,ddof=1)[valid]/mcvar[valid])),**extra}
            rows.append(row)
            for i in range(r):
                replicate_rows.append({'condition':path.stem,'window':label,'observable':name,'method':method,'replica':i,
                    'estimate':float(means[i]),'estimated_variance':None if not valid[i] else float(mcvar[i]),'covered':bool(covered[i])})
    cm=obs['centroid_component'].mean(1);sm=obs['shape'].mean(1);vc=float(cm.var(ddof=1));vs=float(sm.var(ddof=1));cov=float(np.cov(cm,sm,ddof=1)[0,1])
    decomposition={'condition':path.stem,'window':label,'centroid_component_variance':vc,'shape_variance':vs,
                   'twice_covariance':2*cov,'sum':vc+vs+2*cov,'CV_empirical_variance':float(obs['CV'].mean(1).var(ddof=1)),
                   'scope':'between-replica variance of time averages, not same-time marginal decomposition'}
    acf={}
    for name in ['CV','centroid','centroid_component','shape']:
        a=obs[name]-ref['means'][name];nfft=1<<(2*n-1).bit_length();f=np.fft.rfft(a,n=nfft,axis=1)
        sums=np.fft.irfft(f*f.conj(),n=nfft,axis=1)[:,:min(n//2,2048)+1].sum(0)
        gamma=sums/(r*(n-np.arange(len(sums))));rho=gamma/gamma[0]
        acf[name]={'lag':list(range(len(rho))),'pooled_known_mean_acf':rho.tolist(),
                   'marginal_variance_empirical':float(gamma[0]),'marginal_variance_exact':ref['variances'][name],
                   'cumulative_tau':{str(k):float(1+2*rho[1:k+1].sum()) for k in [32,128,512,2048] if k<len(rho)},
                   'scope':'known-target-mean pooled stationary diagnostic; nonstationary starts are descriptive only'}
    return rows,replicate_rows,decomposition,acf


def write_analysis(directory,configs,windows):
    rows=[];rep_rows=[];decomp=[];acfs={}
    for config in configs:
        path=Path(directory)/f"{config['id']}.npz"
        for label,burn,n in windows[config['id']]:
            a,b,c,d=analyze_window(path,burn,n,label)
            rows+=a;rep_rows+=b;decomp.append(c);acfs[f"{config['id']}/{label}"]=d
    save(Path(directory)/'analysis.json',rows);save(Path(directory)/'variance_decomposition.json',decomp);save(Path(directory)/'acf.json',acfs)
    for name,data in [('summary.csv',rows),('replicate_endpoints.csv',rep_rows)]:
        flat=[{k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in row.items()} for row in data]
        keys=list(dict.fromkeys(k for row in flat for k in row))
        with (Path(directory)/name).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(flat)
    return rows


def exploration():
    out=ROOT/'exploratory-v1';out.mkdir(exist_ok=True)
    configs=[{'id':'local_stationary_d04','initial':'stationary','steps':32768,'replicas':128,'delta':.4},
       *[{'id':f'local_{init}_d04','initial':init,'steps':12288,'replicas':128,'delta':.4} for init in ['uniform','zero','displaced']],
       *[{'id':f'local_stationary_d{tag}','initial':'stationary','steps':8192,'replicas':128,'delta':d} for tag,d in [('015',.15),('10',1.)]],
       {'id':'centroid_refresh','initial':'stationary','steps':8192,'replicas':128,'delta':.4,'refresh':True},
       {'id':'reversible_scan','initial':'stationary','steps':8192,'replicas':128,'delta':.4,'reverse_mix':True},
       {'id':'iid_gaussian','initial':'stationary','steps':8192,'replicas':128,'iid':True}]
    windows={c['id']:[('n8192',0,8192)] for c in configs}
    windows['local_stationary_d04']=[('n2048',0,2048),('n8192',0,8192),('n32768',0,32768)]
    for init in ['uniform','zero','displaced']:
        windows[f'local_{init}_d04']=[(f'burn{b}_n8192',b,8192) for b in [0,1024,4096]]
    save(out/'protocol.json',{'classification':'exploratory; no final holdout used','timestamp_utc':datetime.now(timezone.utc).isoformat(),'configs':configs,'windows':windows,'budget':'128 replicas; <1GiB traces expected; single-process CPU; BLAS threads1'})
    for c in configs:run_condition(out,'exploratory-v1',c)
    rows=write_analysis(out,configs,windows)
    for row in rows:
        if row['observable']=='CV' and row['method']=='BM512':print(row['condition'],row['window'],row['covered'],row['replicas'],row['estimated_to_empirical_variance_ratio'],flush=True)


def freeze_confirmation():
    path=ROOT/'confirmation-v1'/'protocol.json'
    if path.exists():raise FileExistsError('Protocol already frozen')
    protocol={'id':'confirmation-v1','timestamp_utc':datetime.now(timezone.utc).isoformat(),'script_sha256':digest(__file__),
        'classification':'locked healthy validation after exploration; no confirmation seeds generated yet',
        'replicas':1024,'M':M,'beta':BETA,'delta':.4,'retained':8192,'burn':0,'initial':'exact stationary finiteM Gaussian',
        'primary_observable':'CV','primary_method':'BM512 (16 batch means; t15)',
        'configs':[{'id':'local_stationary','initial':'stationary','steps':8192,'replicas':1024,'delta':.4},
                   {'id':'centroid_refresh','initial':'stationary','steps':8192,'replicas':1024,'delta':.4,'refresh':True},
                   {'id':'iid_gaussian','initial':'stationary','steps':8192,'replicas':1024,'iid':True}],
        'healthy_gate':{'coverage_CP95_inside':[.92,.98],'variance_ratio_bootstrap_CI95_inside':[.75,1.25],'bias_CI95_inside':[-.01,.01],
                        'required_conditions':['centroid_refresh','iid_gaussian']},
        'mechanism_comparison':{'paired_coverage_improvement_lower95_exceeds':.10,'stationary_initialization_excludes_transient_bias_as_necessary_cause':True},
        'multiplicity':'one primary CV interval method; other methods/TD secondary descriptive; gate is conjunction, no post-test retuning',
        'fault_study_condition':'run only if both required healthy conditions pass frozen gate; fault cal/test use new namespaces',
        'sample_precision':'1024 independent replicas; nominal95% binomial half-width about1.34 percentage points',
        'not_a_general_guarantee':'finiteM harmonic workload, exact stationary initialization and fixed run length only'}
    save(path,protocol);print('Frozen protocol',path,digest(path))


def confirmation():
    out=ROOT/'confirmation-v1';p=json.loads((out/'protocol.json').read_text())
    if p['script_sha256']!=digest(__file__):raise RuntimeError('script changed after protocol freeze')
    for c in p['configs']:run_condition(out,'confirmation-v1',c)
    rows=write_analysis(out,p['configs'],{c['id']:[('primary',0,8192)] for c in p['configs']})
    primary=[r for r in rows if r['observable']=='CV' and r['method']=='BM512'];gates={}
    for r in primary:
        ci=r['coverage_ci95'];vr=r['variance_ratio_bootstrap_ci95'];bias=r['bias_ci95']
        gates[r['condition']]={'coverage':ci[0]>=.92 and ci[1]<=.98,'variance_ratio':vr[0]>=.75 and vr[1]<=1.25,'bias':bias[0]>=-.01 and bias[1]<=.01,'row':r}
        gates[r['condition']]['pass']=all(gates[r['condition']][k] for k in ['coverage','variance_ratio','bias'])
    endpoints=list(csv.DictReader((out/'replicate_endpoints.csv').open()))
    def covers(condition):return np.array([x['covered']=='True' for x in endpoints if x['condition']==condition and x['observable']=='CV' and x['method']=='BM512'])
    d=covers('centroid_refresh').astype(int)-covers('local_stationary').astype(int)
    boot=generator('analysis','paired-confirm',0).integers(len(d),size=(10000,len(d)))
    diff={'paired_coverage_difference':float(d.mean()),'bootstrap_ci95':np.quantile(d[boot].mean(1),[.025,.975]).tolist(),
          'improved_only':int((d==1).sum()),'local_only':int((d==-1).sum()),'passes_margin':bool(np.quantile(d[boot].mean(1),.025)>.10)}
    result={'gates':gates,'comparison':diff,'fault_study_allowed':all(gates[k]['pass'] for k in ['centroid_refresh','iid_gaussian'])}
    save(out/'primary_decision.json',result);print(json.dumps(result,indent=2))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('command',choices=['explore','freeze','confirm']);a=parser.parse_args()
    if a.command=='explore':exploration()
    elif a.command=='freeze':freeze_confirmation()
    else:confirmation()


if __name__=='__main__':main()
