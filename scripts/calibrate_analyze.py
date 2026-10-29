from pathlib import Path
import json
import numpy as np
from scipy.stats import binomtest
from rnglab16.campaign import save
from rnglab17.campaign import CONFIG,cell_specs,seed
from rnglab17.endpoints import rank_threshold,rate,causal_loss

def holm(rows,key='p'):
 order=np.argsort([r[key] for r in rows]);value=0
 for j,i in enumerate(order):value=max(value,min(1,rows[i][key]*(len(rows)-j)));rows[i]['holm_p']=value

def run(root):
 root=Path(root);summary={'cells':{},'H1':{},'H2':{},'coverage_tests':[]};primary=[]
 for spec in cell_specs():
  name=spec['id'];cal=json.loads((root/f'{name}__calibration__healthy.json').read_text());h=json.loads((root/f'{name}__confirmation__healthy.json').read_text());cs=np.array(cal['scores']);hs=np.array(h['scores']);tols=np.array(h['tolerance']);threshold=np.maximum([rank_threshold(cs[:,i]) for i in range(cs.shape[1])],tols)
  records=[]
  for i,obs in enumerate(h['observables']):
   k=int(sum(hs[:,i]>threshold[i]));rec={'observable':obs,'threshold':threshold[i],'practical_tolerance':tols[i],'healthy_calibrated':rate(k,len(hs)),'healthy_accuracy_persistent':rate(int(np.sum(np.array(h['accuracy_persistent'])[:,i])),len(hs)),'mean_error':h['diagnostics'][obs]['mean_error'],'mse':h['diagnostics'][obs]['mse']}
   if i==spec['primary']:
    rec['simultaneous_one_sided_upper']=rate(k,len(hs),alpha=.05/12,one_sided=True)['interval'][1];primary.append({'cell':name,**rec});rec['calibration_gate']=rec['simultaneous_one_sided_upper']<=.10
   for batch in ['batch32','batch128']:
    if batch in h['diagnostics'][obs]:
     covered=h['diagnostics'][obs][batch]['covered'];kk=sum(covered);rec[batch]=rate(kk,len(covered))
     if i==spec['primary']:summary['coverage_tests'].append({'cell':name,'method':batch,'covered':int(kk),'n':len(covered),'p':float(binomtest(kk,len(covered),.95,alternative='less').pvalue)})
   records.append(rec)
  faults=[]
  for p in root.glob(f'{name}__confirmation__*.json'):
   if p.name.endswith('__healthy.json'):continue
   f=json.loads(p.read_text());R=f['replicas'];fs=np.array(f['scores']);means=np.array(f['chain_means']);hh=np.array(h['chain_means'])[:R];dc=(fs>threshold);hc=(hs[:R]>threshold);obsrecords=[]
   for i,obs in enumerate(h['observables']):
    b=int(sum(dc[:,i]&~hc[:,i]));c=int(sum(~dc[:,i]&hc[:,i]));pval=1. if b+c==0 else float(binomtest(b,b+c,.5).pvalue)
    row={'observable':obs,'detected':rate(int(sum(dc[:,i])),R),'paired_detection_improvement':float(np.mean(dc[:,i].astype(float)-hc[:,i])),'discordant_fault_only':b,'discordant_healthy_only':c,'p':pval,'causal_loss':causal_loss(hh[:,i],means[:,i],h['reference'][i],tols[i],seed('bootstrap',name,f['fault'],i))}
    if i==spec['primary'] and f['fault']=='mild_bias':summary['H2'].setdefault('primary_tests',[]).append({'cell':name,**row})
    obsrecords.append(row)
   faults.append({'fault':f['fault'],'structural_invalid_analytic':f['structural_invalid_analytic'],'source_diagnostic_is_consumed':f['source_diagnostic_is_consumed'],'observables':obsrecords})
  summary['cells'][name]={'spec':spec,'primary':spec['primary'],'thresholds':threshold,'observables':records,'faults':faults}
 holm(summary['H2']['primary_tests']);holm(summary['coverage_tests'])
 summary['H1']={'primary_cells':primary,'all_pass':all(r['simultaneous_one_sided_upper']<=.10 for r in primary),'nominal_rank_target':.05,'validation_upper_ceiling':.10}
 summary['H2']['supported_cells']=[r['cell'] for r in summary['H2']['primary_tests'] if r['holm_p']<.05 and r['paired_detection_improvement']>=.20]
 summary['coverage_tests_failed']=[r for r in summary['coverage_tests'] if r['holm_p']<.05]
 summary['forecast_gate']=summary['cells']['ising_L3_T1']['observables'][0]['calibration_gate']
 # Prespecified near-critical interaction of paired excess energy MSE, distinct from raw severity.
 interactions=[]
 for L in [2,3,4]:
  samples=[]
  for ti in range(3):
   h=json.loads((root/f'ising_L{L}_T{ti}__confirmation__healthy.json').read_text());f=json.loads((root/f'ising_L{L}_T{ti}__confirmation__mild_bias.json').read_text());ref=h['reference'][0];a=np.array(h['chain_means'])[:128,0]-ref;b=np.array(f['chain_means'])[:,0]-ref;samples.append(b*b-a*a)
  g=np.random.default_rng(seed('criticality',L));delta=samples[1].mean()-max(samples[0].mean(),samples[2].mean());boot=[]
  for _ in range(2000):
   draws=[x[g.integers(len(x),size=len(x))].mean() for x in samples];boot.append(draws[1]-max(draws[0],draws[2]))
  interactions.append({'L':L,'near_minus_max_other_excess_MSE':float(delta),'approx_family98_33pct_ci':np.quantile(boot,[.05/6,1-.05/6]),'caveat':'Lattice-specific severity is min(.02,.5exp(-8/T)), so this is vulnerability under a common positive-transition safety rule, not equal physical severity.'})
 summary['near_criticality']=interactions
 save(root/'calibration_analysis.json',summary)
 print(json.dumps({'H1_all_pass':summary['H1']['all_pass'],'H1_cells':[(r['cell'],r['healthy_calibrated']['events'],r['simultaneous_one_sided_upper']) for r in primary],'H2_supported_cells':summary['H2']['supported_cells'],'coverage_failed':summary['coverage_tests_failed'],'forecast_gate':summary['forecast_gate']},indent=2))
if __name__=='__main__':
 import argparse;p=argparse.ArgumentParser();p.add_argument('root');run(p.parse_args().root)
