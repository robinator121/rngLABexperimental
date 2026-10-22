import argparse,json,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import ttest_1samp,beta
from rnglab16.core import fit_logistic,predict
from rnglab16.campaign import save,interval

def run(root):
    def read(n):return json.loads((root/(n+'.json')).read_text())
    cfg=read('config');output={};tests=[]
    for p in sorted(root.glob('B_*.json')):
        r=json.loads(p.read_text())
        if r['fault']=='healthy':continue
        h=json.loads(p.with_name(p.name.replace('_'+r['fault']+'.json','_healthy.json')).read_text());diff=np.array(r['chain_metrics'])[:,0]-np.array(h['chain_metrics'])[:,0]
        tests.append({'cell':p.stem,'paired_energy_effect':interval(diff),'p':float(ttest_1samp(diff,0).pvalue) if np.std(diff)>0 else (0. if diff.mean()!=0 else 1.)})
    order=np.argsort([r['p'] for r in tests]);running=0
    for rank,i in enumerate(order):running=max(running,min(1.,tests[i]['p']*(len(tests)-rank)));tests[i]['holm_p']=running
    output['ising_paired_energy']=tests
    def data(split):
        X=[];y=[];families=[];excluded={}
        for p in sorted(root.glob(f'C_{split}_*.json')):
            r=json.loads(p.read_text());good=np.array(r['eligible'],bool);X.extend(np.array(r['features'])[good]);y.extend(np.array(r['labels'])[good]);families.extend([r['family']]*sum(good));excluded[r['family']]=int(sum(~good))
        return np.array(X),np.array(y),np.array(families),excluded
    X,y,_,_=data('development');models={};features={'source':[0,1],'application':[2,3],'combined':[0,1,2,3],'theory':[0,1,2,3,4,5]}
    for key,cols in features.items():models[key]=fit_logistic(X[:,cols],y)
    output['models']=models;output['development_eligible']=len(y);output['forecast']={}
    predictions={}
    for split in ['validation','test','OOD']:
        X,y,families,excluded=data(split);predictions[split]={};summary={'n':len(y),'events':int(sum(y)),'excluded_before_cutoff':excluded,'models':{}}
        for key,cols in features.items():
            p=predict(models[key],X[:,cols]);predictions[split][key]=p
        predictions[split]['constant']=np.full(len(y),np.mean(data('development')[1]));predictions[split]['never_fail']=np.zeros(len(y))
        for key,p in predictions[split].items():
            bins=[]
            for a,b in zip(np.arange(0,1,.2),np.arange(.2,1.01,.2)):
                mask=(p>=a)&(p<b if b<1 else p<=b)
                if mask.any():bins.append({'n':int(mask.sum()),'predicted':float(p[mask].mean()),'observed':float(y[mask].mean())})
            # Descriptive operating decision at fixed risk threshold .1; abstain outside training range.
            out=np.any((X < np.min(data('development')[0],axis=0))|(X>np.max(data('development')[0],axis=0)),axis=1)
            theory_out=(2*np.maximum(X[:,0],0)>.025) if key=='theory' else np.zeros(len(y),bool)
            out=out|theory_out
            permit=(p<.1)&~out
            summary['models'][key]={'brier':float(np.mean((p-y)**2)),'calibration':bins,'permitted':int(permit.sum()),'failures_permitted':int(y[permit].sum()),'abstained':int(out.sum()),'by_family':{family:{'n':int(sum(families==family)),'brier':float(np.mean((p[families==family]-y[families==family])**2))} for family in np.unique(families)}}
        diff=(predictions[split]['theory']-y)**2-(predictions[split]['combined']-y)**2
        rng=np.random.default_rng(166000+len(y));boots=np.array([diff[rng.integers(len(y),size=len(y))].mean() for _ in range(2000)])
        summary['primary_theory_minus_combined_brier']={'difference':float(diff.mean()),'paired_bootstrap_ci95':np.quantile(boots,[.025,.975]),'H3_supported':bool(np.quantile(boots,.975)<-.02),'unit':'one eligible trajectory; IID seed groups, conditional on fixed fitted model'}
        output['forecast'][split]=summary
        np.savez_compressed(root/f'predictions_{split}.npz',labels=y,families=families,**predictions[split])
    save(root/'analysis.json',output)
    print(json.dumps({'test':output['forecast']['test'],'OOD':output['forecast']['OOD']},default=lambda x:x.tolist(),indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);run(p.parse_args().root)
