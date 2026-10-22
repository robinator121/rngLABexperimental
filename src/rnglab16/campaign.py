"""Execute frozen, finite TRIDENT studies; atomic cell checkpoints and full traces."""
import argparse,hashlib,json,os,time,platform,resource
from pathlib import Path
import numpy as np
from scipy.stats import t as student_t,beta as beta_dist
from .core import *

CONFIG={'version':'0.16.0.dev3','status':'bounded_confirmation_with_exploratory_extensions',
 'A':{'seed':180001,'replicas':256,'n':2048,'parameters':[[.1,.2],[.01,.02],[.7,.8]],'epsilons':[.0001,.001,.01,.05]},
 'B':{'seed':181000,'replicas':32,'steps':1200,'burn':300,'lattices':[2,3,4,8],'temperatures':[1.5,float(TC),3.5],'faults':['healthy','bias','coupled'], 'acceptance_bias':.15},
 'quantum':{'seed':182000,'replicas':32,'steps':3000,'burn':500,'beads':[4,16],'faults':['healthy','bias','coupled']},
 'C':{'seed':183000,'L':3,'T':float(TC),'steps':500,'cutoff':150,'horizon':200,'window':64,'stride':20,'tolerance':.15,'persistent_checks':2,'severity':.3,'development_per_family':64,'validation_per_family':64,'test_per_family':128,'families':['healthy','abrupt','drift'],'OOD':['periodic','coupled','correlation']},
 'limits':{'threads':1,'wall_seconds':900,'retries':0,'no_cloud':True},
 'hypotheses':{'H1':'exact sensitivity remainder is quadratic away from singular a+b=0','H2':'bias changes Ising energy; observable response differs by model','H3':'theory features improve heldout Brier over source+application by at least .02'},
 'theory_domain':{'delta_max':.025,'C_amendment':'run002 had zero future events at tolerance .5; new tolerance .15 and cutoff150 frozen with fresh seeds before run004', 'reason':'attempt001 lost stationary uniqueness at larger acceptance bias; subsequent replay is amended development, not independent confirmation'},
 'inference':{'unit':'independent SeedSequence child chain; paired arms within same base namespace','alpha':.05,'multiple_testing':'Holm across Ising paired energy tests; Brier primary contrast paired trajectory bootstrap','bootstrap_repetitions':2000,'coverage':'pointwise exact binomial; no universal guarantee','confirmation':'no threshold or features tuned on final test','power':'32 paired B units detect about 0.51 paired SD at unadjusted 80% power; Holm reduces power. H3 n384 IID rows gives rough Brier difference SE<=0.051 worst-case; small gains may remain inconclusive.'}}

def encode(x):
    if isinstance(x,np.ndarray):return x.tolist()
    if isinstance(x,np.generic):return x.item()
    raise TypeError(type(x))
def save(path,obj):
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(obj,default=encode,indent=2,allow_nan=False));tmp.replace(path)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def interval(values):
    v=np.asarray(values);mean=v.mean();se=v.std(ddof=1)/np.sqrt(len(v));q=student_t.ppf(.975,len(v)-1)
    return {'n':len(v),'mean':float(mean),'se':float(se),'ci95':[float(mean-q*se),float(mean+q*se)]}
def chain_metrics(trace,burn,L,T):
    x=trace[:,burn:];E=x[:,:,0];m=x[:,:,1];N=L*L
    return np.column_stack([E.mean(1),m.mean(1),x[:,:,2].mean(1),N*(x[:,:,3].mean(1)-m.mean(1)**2)/T,N*(x[:,:,4].mean(1)-E.mean(1)**2)/T**2])
def tau(x):
    x=x-x.mean();n=len(x);f=np.fft.rfft(x,n=2*n);a=np.fft.irfft(f*f.conj())[:n];a/=max(a[0],1e-30)
    pairs=a[1:-1:2]+a[2::2];stop=np.where(pairs<=0)[0];k=stop[0] if len(stop) else len(pairs)
    return float(max(1,1+2*pairs[:k].sum()))
def uq(x,reference):
    ts=np.array([tau(y) for y in x]);n=x.shape[1];means=x.mean(1);se=np.sqrt(x.var(1,ddof=1)*ts/n)
    covered=abs(means-reference)<=1.96*se;k=int(covered.sum());R=len(x)
    ci=[0 if k==0 else beta_dist.ppf(.025,k,R-k+1),1 if k==R else beta_dist.ppf(.975,k+1,R-k)]
    return {'chain_means':means,'tau':None if np.all(x.var(1)==0) else ts,'ess':None if np.all(x.var(1)==0) else n/ts,'coverage_count':k,'chains':R,'coverage_ci95':ci,'method':'truncated initial-positive adjacent ACF pairs; exploratory finite-run MCSE; no theorem for systematic scan','bias':interval(means-reference),'mse':float(np.mean((means-reference)**2))}

def exact_ising_transition(L,T,delta):
    n=L*L;states=((np.arange(2**n)[:,None]>>np.arange(n))&1)*2-1;grid=states.reshape(-1,L,L);E=energy(grid);P=np.zeros((len(E),len(E)))
    for k in range(n):
        j=np.arange(len(E))^(1<<k);de=E[j]-E;acc=np.where(de<=0,1,np.clip(np.exp(-np.maximum(0,de)/T)-delta,0,1))/n
        P[np.arange(len(E)),j]+=acc;P[np.arange(len(E)),np.arange(len(E))]+=1/n-acc
    return P,E/n

def run(root):
    start=time.monotonic();root.mkdir(parents=True,exist_ok=True);(root/'traces').mkdir(exist_ok=True)
    cfg=root/'config.json'
    if cfg.exists():assert json.loads(cfg.read_text())==CONFIG
    else:save(cfg,CONFIG)
    source={str(p.relative_to(Path(__file__).parents[2])):digest(p) for p in Path(__file__).parent.glob('*.py')}
    frozen=root/'source_hashes.json'
    if frozen.exists():assert json.loads(frozen.read_text())==source,'Source changed: new run required'
    else:save(frozen,source)
    def checkpoint(name,fn):
        if time.monotonic()-start>CONFIG['limits']['wall_seconds']:raise TimeoutError('campaign budget exhausted')
        p=root/(name+'.json')
        if p.exists():return json.loads(p.read_text())
        a=time.monotonic();obj=fn();obj['wall_seconds']=time.monotonic()-a;save(p,obj);print(name,round(obj['wall_seconds'],2),flush=True);return obj
    def A():
        rows=[];g=np.random.default_rng(CONFIG['A']['seed']);R=256;n=2048
        for a,b in CONFIG['A']['parameters']:
            P=two_state(a,b);p=a/(a+b);x=(g.random(R)<p).astype(int);sums=np.zeros(R)
            for _ in range(n):u=g.random(R);x=np.where(x==0,u<a,~(u<b)).astype(int);sums+=x
            rows.append({'a':a,'b':b,'pi':stationary(P),'finite_n_theory_var':variance_mean(a,b,n),'empirical_variance':float((sums/n).var(ddof=1)),'means':sums/n,'eigenvalue':1-a-b,
             'perturbations':[{'eps':e,'delta_exact':(a+e)/(a+b+e)-p,'linear':b*e/(a+b)**2,'remainder':-b*e*e/((a+b)**2*(a+b+e))} for e in CONFIG['A']['epsilons']]})
        augmented_rows=[]
        for rho in [0,.5,.9]:
            P=augmented(rho=rho);pi=stationary(P);power=np.ones(len(P))/len(P)
            for _ in range(10000):power=power@P
            augmented_rows.append({'rho':rho,'pi_x1':float(pi.reshape(2,8)[1].sum()),'stationary_residual':float(abs(pi@P-pi).max()),'power_difference':float(abs(power-pi).max()),'naive_original_x_prediction':1/3})
        return {'rows':rows,'augmented':augmented_rows,'three_state':stationary([[.7,.2,.1],[.1,.6,.3],[.2,.3,.5]])}
    checkpoint('A',A)
    B=CONFIG['B'];names=['energy','magnetization','abs_magnetization','susceptibility','heat_capacity']
    counter=0
    for L in B['lattices']:
        for T in B['temperatures']:
            base=B['seed']+counter;counter+=1
            for fault in B['faults']:
                name=f'B_L{L}_T{T:.5f}_{fault}'
                def cell(L=L,T=T,base=base,fault=fault,name=name):
                    tr,u,acc=ising(L,T,B['replicas'],B['steps'],base,fault,B['acceptance_bias'])
                    np.savez_compressed(root/'traces'/f'{name}.npz',observables=tr,source_mean=u,acceptance=acc)
                    metrics=chain_metrics(tr,B['burn'],L,T);ref=enumeration(L,T) if L<=4 else None
                    return {'L':L,'T':T,'fault':fault,'base_seed':base,'metrics':{k:interval(metrics[:,i]) for i,k in enumerate(names)},'chain_metrics':metrics,'exact':ref,'energy_uq':uq(tr[:,B['burn']:,0],ref['energy']) if ref else None,'consumed_uniform_mean':float(u.mean()),'trace':f'traces/{name}.npz'}
                checkpoint(name,cell)
    for L in [4,8]:
        for T in B['temperatures']:
            def clusters(L=L,T=T):
                traces=[];sizes=[]
                for i in range(8):
                    tr,sz=wolff(L,T,1000,184000+L*100+i+int(T*10));traces.append(tr);sizes.append(sz)
                traces=np.array(traces);np.savez_compressed(root/'traces'/f'wolff{L}_{T:.5f}.npz',observables=traces,cluster_sizes=sizes)
                return {'L':L,'T':T,'unit':'cluster flips, not sweeps; unequal computational work','metrics':{k:interval(chain_metrics(traces,200,L,T)[:,i]) for i,k in enumerate(names)},'mean_cluster_fraction':float(np.mean(sizes)/(L*L)),'energy_tau':[tau(x) for x in traces[:,200:,0]]}
            checkpoint(f'W_L{L}_T{T:.5f}',clusters)
    Q=CONFIG['quantum']
    for model in ['vmc','pimc4','pimc16']:
        for fault in Q['faults']:
            def quantum(model=model,fault=fault):
                if model=='vmc':tr=vmc(Q['replicas'],Q['steps'],Q['seed'],fault);refs=[.5,.5,.75];obs=['energy','x2','x4']
                else:
                    M=int(model[4:]);tr=pimc(Q['replicas'],Q['steps'],M,Q['seed']+M,fault);e=exact_finite_bead_energy(M,2.);refs=[e,e,e,.5];obs=['TD','CV','x2','centroid2']
                np.savez_compressed(root/'traces'/f'{model}_{fault}.npz',observables=tr)
                return {'model':model,'fault':fault,'metrics':{k:uq(tr[:,Q['burn']:,i],refs[i]) for i,k in enumerate(obs)},'reset_vs_cumulative_means':{'last500':tr[:,-500:].mean(1),'retained':tr[:,Q['burn']:].mean(1)}}
            checkpoint(f'Q_{model}_{fault}',quantum)
    # Theory table is deterministic, independent of all forecast training/test trajectories.
    def theory():
        grid=[]
        for delta in [0,.005,.01,.015,.02,.025]:
            P,E=exact_ising_transition(3,TC,delta);grid.append([delta,float(stationary(P)@E)])
        return {'grid':grid,'scope':'stationary IID acceptance-bias response only; not transient or coupled-source law'}
    theory_table=checkpoint('C_theory',theory)['grid'];c=CONFIG['C'];ref=enumeration(3,TC)['energy']
    for split,R,offset in [('development',64,0),('validation',64,100),('test',128,200),('OOD',128,300)]:
        families=c['OOD'] if split=='OOD' else c['families']
        for family in families:
            def forecasts(family=family,R=R,split=split,offset=offset):
                base=c['seed']+offset+c['families'].index(family) if family in c['families'] else c['seed']+offset+c['OOD'].index(family)
                tr,u,acc=ising(3,TC,R,c['steps'],base,family,.8 if family=='correlation' else c['severity'])
                E=tr[:,:,0];cut=c['cutoff'];window=c['window'];times=np.arange(window,c['steps']+1,c['stride'])
                rolling=np.column_stack([E[:,t-window:t].mean(1) for t in times]);exceed=abs(rolling-ref)>c['tolerance'];event=np.full(R,-1)
                for j in range(1,len(times)):
                    hit=exceed[:,j-1]&exceed[:,j]&(event<0);event[hit]=times[j]
                eligible=(event<0)|(event>cut);label=(event>cut)&(event<=cut+c['horizon'])
                src=u[:,cut-window:cut].mean(1)-.5;drift=u[:,cut-32:cut].mean(1)-u[:,cut-64:cut-32].mean(1)
                app=E[:,cut-window:cut].mean(1)-ref;appdrift=E[:,cut-32:cut].mean(1)-E[:,cut-64:cut-32].mean(1)
                inferred=np.maximum(0,2*src);table=np.array(theory_table);theory_error=np.interp(inferred,table[:,0],table[:,1])-ref
                X=np.column_stack([src,drift,app,appdrift,abs(theory_error),abs(theory_error)*np.maximum(0,drift)*c['horizon']])
                np.savez_compressed(root/'traces'/f'C_{split}_{family}.npz',observables=tr,source_mean=u,features=X,events=event,eligible=eligible,labels=label)
                return {'split':split,'family':family,'base_seed':base,'independent_units':R,'eligible':eligible,'events':event,'labels':label.astype(int),'features':X,'theory_domain_exceeded':inferred>.025,'future_oracle_used_in_features':False,'right_censored_at':c['steps']}
            checkpoint(f'C_{split}_{family}',forecasts)
    save(root/'execution.json',{'elapsed_seconds':time.monotonic()-start,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':platform.python_version(),'numpy':np.__version__,'platform':platform.platform(),'config_sha256':digest(cfg),'completed':True})

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('output',type=Path);args=p.parse_args();run(args.output)
