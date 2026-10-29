"""CALIBRATE: frozen whole-trajectory endpoints, independent confirmation, paired faults."""
from pathlib import Path
import argparse,json,hashlib,time,resource,platform,os,multiprocessing
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from scipy.stats import t as student_t
from rnglab16.core import ising,vmc,pimc,enumeration,TC
from rnglab16.campaign import save,interval,tau
from rnglab14.models.pimc_harmonic import exact_finite_bead_energy
from .endpoints import *

CONFIG={'version':'0.17.0.dev1','root_namespace':'CALIBRATE-20261022-v1','workers':2,'threads_per_worker':1,
 'Ising':{'lattices':[2,3,4],'temperatures':[1.5,float(TC),3.5],'steps':2048,'burn':512,'window':512,'stride':128,'accuracy_tolerance':[.1,.1,.1,.2,.2]},
 'quantum':{'models':['vmc','pimc4','pimc16'],'steps':4096,'burn':1024,'window':512,'stride':128},
 'counts':{'pilot':32,'calibration':399,'healthy_confirmation':512,'fault_confirmation':128},
 'faults':['mild_bias','coupled','correlation','periodic','structural_bias','unused_source'],
 'hypotheses':{'H1':'Whole-path healthy false declaration target .05; validation simultaneous upper95<=.10 for twelve primary cells',
 'H2':'Paired harmful mild-bias declarations exceed healthy at Holm alpha.05 AND detection improvement>=.20; excess loss independently reported',
 'H3':'Physical references independently agree; interval coverage validated separately, never assumed from formula validity',
 'H4':'Frozen combined predictor beats train prevalence constant in test Brier by .01 with upper adjusted CI below -.01, and useful continuation>=.20 with risk gate and lower cost than always stop',
 'H5':'Same combined candidate retains H4 score/utility advantage in excluded OOD families'},
 'statistics':{'alpha':.05,'primary_cells':12,'calibration_unit':'complete independent trajectory, one max-persistent score per observable','primary_H1_bound':'one-sided exact Bonferroni .05/12 upper','H2_family':'twelve mild-bias primary endpoints, Holm exact paired McNemar','bootstrap':2000,'rank_target':.05,'tolerance_floor':True,'within_path_points_are_not_replicates':True},
 'forecast':{'L':3,'T':float(TC),'steps':2048,'burn':512,'window':512,'stride':128,'cutoff':1408,'horizon':512,'families':['healthy','abrupt','drift'],'OOD':['periodic','correlation'],'counts':[128,128,256,256],'gate':'L3/Tc whole-path H1 upper<=.10; otherwise do not run forecasting','primary_candidate':'combined','primary_baseline':'development prevalence constant','H4_H5_CI':'paired trajectory percentile bootstrap98.75% per endpoint; approximate family adjustment','continued_risk_validation_ceiling':.1,'min_validation_continued':32,'continue_fraction_target':.2,'costs':{'continue_failure':10,'stop':1,'review_and_stop':1.3,'intervention':.5,'intervention_failure_loss_no_mitigation':10,'intervention_failure_loss_hypothetical50percent':5},'intervention_benefit_measured':False},
 'limits':{'workers':2,'expected_peak_memory_GiB':3,'wall_budget_seconds':1200,'cloud':False,'GPU':False,'retry':0},
 'references':['https://proceedings.mlr.press/v25/vovk12/vovk12.pdf','https://arxiv.org/abs/0811.1729']}

def seed(*parts):return int.from_bytes(hashlib.sha256(('/'.join([CONFIG['root_namespace'],*map(str,parts)])).encode()).digest()[:8],'big')
def cell_specs():
 out=[]
 for L in CONFIG['Ising']['lattices']:
  for ti,T in enumerate(CONFIG['Ising']['temperatures']):out.append({'id':f'ising_L{L}_T{ti}','model':'ising','L':L,'T':T,'primary':0})
 for m in CONFIG['quantum']['models']:out.append({'id':m,'model':m,'primary':1})
 return out

def simulate(spec,R,base,fault='healthy'):
 model=spec['model'];s=CONFIG['Ising'] if model=='ising' else CONFIG['quantum'];start=time.monotonic();source=None;consumed=None
 if model=='ising':
  T=spec['T'];mild=min(.02,.5*np.exp(-8/T));f={'mild_bias':'bias','structural_bias':'bias','unused_source':'healthy'}.get(fault,fault)
  severity=.8 if fault=='correlation' else .15 if fault=='structural_bias' else mild
  trace,source,accept=ising(spec['L'],T,R,s['steps'],base,f,severity)
  consumed=source.copy()
  if fault=='unused_source':source=np.clip(source+.1,0,1)
  exact=enumeration(spec['L'],T);reference=np.array([exact[k] for k in ['energy','magnetization','abs_magnetization','susceptibility','heat_capacity']]);tols=np.array(s['accuracy_tolerance']);observables=['energy','magnetization','abs_magnetization','susceptibility','heat_capacity']
  # Retain primitive sufficient statistics; nonlinear observables are computed per window.
  N=spec['L']**2;times=np.arange(s['burn']+s['window'],s['steps']+1,s['stride']);means=[]
  for t in times:
   a=trace[:,t-s['window']:t];e=a[:,:,0];m=a[:,:,1]
   means.append(np.column_stack([e.mean(1),m.mean(1),a[:,:,2].mean(1),N*(a[:,:,3].mean(1)-m.mean(1)**2)/T,N*(a[:,:,4].mean(1)-e.mean(1)**2)/T**2]))
  means=np.stack(means,axis=1);error=abs(means-reference);persistent=np.minimum(error[:,:-1],error[:,1:]);endpoint={'times':times[1:],'means':means,'error':error,'persistent_error':persistent,'score':persistent.max(1),'accuracy_any':(error>tols).any(1),'accuracy_persistent':(persistent>tols).any(1)}
  retained=trace[:,s['burn']:];e=retained[:,:,0];m=retained[:,:,1];chainmeans=np.column_stack([e.mean(1),m.mean(1),retained[:,:,2].mean(1),N*(retained[:,:,3].mean(1)-m.mean(1)**2)/T,N*(retained[:,:,4].mean(1)-e.mean(1)**2)/T**2]);structural=(fault=='structural_bias' and .15>=np.exp(-8/T))
 else:
  fault2='bias' if fault=='mild_bias' else 'healthy' if fault in ['unused_source','structural_bias','periodic','correlation'] else fault
  if model=='vmc':trace=vmc(R,s['steps'],base,fault2);reference=np.array([.5,.5,.75]);tols=np.array([.02,.1,.25]);observables=['energy','x2','x4']
  else:
   M=int(model[4:]);trace=pimc(R,s['steps'],M,base,fault2);e=exact_finite_bead_energy(M,2.);reference=np.array([e,e,e,.5]);tols=np.array([.1,.1,.1,.1]);observables=['TD','CV','x2','centroid2']
  endpoint=window_metrics(trace,s['burn'],s['window'],s['stride'],reference,tols);retained=trace[:,s['burn']:];chainmeans=retained.mean(1);structural=False
 result={'spec':spec,'base_seed':base,'replicas':R,'fault':fault,'reference':reference,'tolerance':tols,'observables':observables,'chain_means':chainmeans,'scores':endpoint['score'],'accuracy_persistent':endpoint['accuracy_persistent'],'structural_invalid_analytic':bool(structural),'source_diagnostic_is_consumed':fault!='unused_source','wall_seconds':time.monotonic()-start,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
 summaries={}
 for i,name in enumerate(observables):
  q={'mean_error':interval(chainmeans[:,i]-reference[i]),'mse':float(np.mean((chainmeans[:,i]-reference[i])**2))}
  if model!='ising' or i<3:
   q['batch32']=batch_intervals(retained[:,:,i],reference[i],32);q['batch128']=batch_intervals(retained[:,:,i],reference[i],128)
   q['tau']=None if np.all(retained[:,:,i].var(1)==0) else [tau(x) for x in retained[:,:,i]]
   q['ess']=None if q['tau'] is None else (retained.shape[1]/np.array(q['tau']))
  summaries[name]=q
 result['diagnostics']=summaries
 return result,{'observables':trace,'source_diagnostic_mean':np.empty((0,)) if source is None else source,'consumed_source_mean':np.empty((0,)) if consumed is None else consumed,'persistent':endpoint['persistent_error'],'times':endpoint['times'],'window_means':endpoint['means']}

def execute_task(task):
 root=Path(task['root']);spec=task['spec'];phase=task['phase'];fault=task['fault'];name=f"{spec['id']}__{phase}__{fault}";path=root/(name+'.json')
 if path.exists():return name
 R=CONFIG['counts']['calibration' if phase=='calibration' else 'healthy_confirmation' if fault=='healthy' else 'fault_confirmation'];base=seed(spec['id'],phase)
 result,traces=simulate(spec,R,base,fault)
 np.savez_compressed(root/'traces'/(name+'.npz'),**traces);save(path,result)
 return name

def main(root,stage):
 root=Path(root);root.mkdir(parents=True,exist_ok=True);(root/'traces').mkdir(exist_ok=True)
 if stage=='pilot':
  out=[]
  for spec in [cell_specs()[4],cell_specs()[-1]]:
   r,z=simulate(spec,32,seed('pilot',spec['id']));out.append(r);np.savez_compressed(root/'traces'/('pilot_'+spec['id']+'.npz'),**z)
  save(root/'pilot.json',{'classification':'feasibility only, excluded from calibration and confirmation','runs':out});return
 cfg=root/'config.json'
 if cfg.exists():assert json.loads(cfg.read_text())==CONFIG
 else:save(cfg,CONFIG)
 sourcefiles=list(Path(__file__).parent.glob('*.py'))+[Path(__file__).parents[1]/'rnglab16/core.py']
 hashes={str(p.relative_to(Path(__file__).parents[1])):hashlib.sha256(p.read_bytes()).hexdigest() for p in sourcefiles}
 locked=root/'source_hashes.json'
 if locked.exists():assert json.loads(locked.read_text())==hashes,'source mismatch; new run required'
 else:save(locked,hashes)
 start=time.monotonic()
 for phase in ['calibration','confirmation']:
  tasks=[]
  for spec in cell_specs():
   faults=['healthy'] if phase=='calibration' else ['healthy','mild_bias','coupled','unused_source']+(['correlation','periodic','structural_bias'] if spec['model']=='ising' else [])
   tasks +=[{'root':str(root),'spec':spec,'phase':phase,'fault':fault} for fault in faults]
  with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context("fork")) as pool:
   for name in pool.map(execute_task,tasks):
    print(name,flush=True)
    if time.monotonic()-start>CONFIG['limits']['wall_budget_seconds']:raise TimeoutError('stage wall budget')
 save(root/'execution.json',{'elapsed_seconds':time.monotonic()-start,'python':platform.python_version(),'numpy':np.__version__,'workers':2,'threads':1,'completed':True})

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root');p.add_argument('--stage',choices=['pilot','confirmation'],default='confirmation');a=p.parse_args();main(a.root,a.stage)
