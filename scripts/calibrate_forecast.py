"""Incident risk only after the independently tested healthy-calibration gate."""
from pathlib import Path
import argparse,json,time,hashlib
import numpy as np
from rnglab16.core import ising,enumeration,TC,fit_logistic,predict
from rnglab16.campaign import save
from rnglab17.campaign import CONFIG,seed
from rnglab17.endpoints import window_metrics,first_events,rank_threshold,rate

FEATURES={'source':[0,1,6],'application':[2,3],'combined':[0,1,2,3,6],'theory':[0,1,2,3,4,5,6]}
def generate(root,split,R,family,threshold,table):
 p=root/f'forecast_{split}_{family}.json'
 if p.exists():return json.loads(p.read_text())
 c=CONFIG['forecast'];base=seed('forecast',split,family);T=c['T'];mild=.02
 tr,u,acc=ising(3,T,R,c['steps'],base,family,.8 if family=='correlation' else mild)
 ref=enumeration(3,T)['energy'];ep=window_metrics(tr[:,:,:1],c['burn'],c['window'],c['stride'],[ref],[.1]);event=first_events(ep['persistent_error'],ep['times'],[threshold])[:,0];cut=c['cutoff'];h=c['horizon'];eligible=(event<0)|(event>cut);y=(event>cut)&(event<=cut+h)
 source=u[:,cut-512:cut];E=tr[:,:,0];s=source.mean(1)-.5;drift=u[:,cut-256:cut].mean(1)-u[:,cut-512:cut-256].mean(1);app=E[:,cut-512:cut].mean(1)-ref;trend=E[:,cut-256:cut].mean(1)-E[:,cut-512:cut-256].mean(1)
 inferred=1-np.sqrt(np.maximum(0,1-2*np.maximum(s,0)));delta=np.clip(inferred,0,.025);theory=np.interp(delta,table[:,0],table[:,1])-ref;future=np.clip(inferred+2*np.maximum(drift,0),0,.025);projected=np.interp(future,table[:,0],table[:,1])-ref
 centered=source-source.mean(1)[:,None];ac=np.sum(centered[:,:-1]*centered[:,1:],1)/np.maximum(1e-20,np.sum(centered*centered,1));X=np.column_stack([s,drift,app,trend,abs(theory),abs(projected)-abs(app),ac])
 result={'split':split,'family':family,'n':R,'base_seed':base,'eligible':eligible,'events':event,'labels':y.astype(int),'features':X,'theory_domain_exceeded':inferred>.025,'complete_horizon_followup':True,'cutoff':cut,'horizon':h,'threshold':threshold}
 np.savez_compressed(root/'traces'/f'forecast_{split}_{family}.npz',observables=tr,source_mean=u,acceptance=acc,features=X,events=event)
 save(p,result);print('forecast',split,family,'eligible',int(sum(eligible)),'events',int(sum(y&eligible)),flush=True);return result

def run(root):
 root=Path(root);analysis=json.loads((root/'calibration_analysis.json').read_text());c=CONFIG['forecast']
 if not analysis['forecast_gate']:
  save(root/'forecast_analysis.json',{'status':'NOT_EXECUTED_CALIBRATION_GATE_FAILED','reason':'Ising L3/Tc adjusted healthy upper exceeds .10'});return
 threshold=analysis['cells']['ising_L3_T1']['thresholds'][0]
 table=np.array(json.loads((Path(__file__).resolve().parents[2]/'reference/v016_theory_grid.json').read_text())['grid'])
 splitdata={};registry=[]
 for split,R in zip(['development','validation','test','OOD'],c['counts']):
  arr=[]
  for family in (c['OOD'] if split=='OOD' else c['families']):
   row=generate(root,split,R,family,threshold,table);arr.append(row);registry.append({'split':split,'family':family,'n':R,'seed':row['base_seed'],'eligible':int(sum(row['eligible'])),'precutoff':int(sum(~np.array(row['eligible'],bool))),'incident_failures':int(sum(np.array(row['labels'])[np.array(row['eligible'],bool)]))})
  X=[];y=[];family=[]
  for r in arr:
   good=np.array(r['eligible'],bool);X.extend(np.array(r['features'])[good]);y.extend(np.array(r['labels'])[good]);family.extend([r['family']]*sum(good))
  splitdata[split]=(np.array(X),np.array(y),np.array(family))
 trainX,trainy,_=splitdata['development']
 if len(np.unique(trainy))<2:
  save(root/'forecast_analysis.json',{'status':'NO_TWO_CLASS_DEVELOPMENT_TARGET','registry':registry});return
 models={k:fit_logistic(trainX[:,cols],trainy) for k,cols in FEATURES.items()}
 if not all(m['success'] for m in models.values()):
  save(root/'forecast_analysis.json',{'status':'BLOCKED_OPTIMIZER_NONCONVERGENCE','models':models,'registry':registry});return
 # Fixed static warning boundaries derived exclusively from healthy calibration prefix data.
 z=np.load(root/'traces/ising_L3_T1__calibration__healthy.npz');cut=c['cutoff'];ref=enumeration(3,TC)['energy'];source_boundary=rank_threshold(abs(z['source_diagnostic_mean'][:,cut-512:cut].mean(1)-.5),.10);app_boundary=rank_threshold(abs(z['observables'][:,cut-512:cut,0].mean(1)-ref),.10)
 flags={'static_source':abs(trainX[:,0])>source_boundary,'static_application':abs(trainX[:,2])>app_boundary};static={}
 for k,flag in flags.items():static[k]=[(float(trainy[flag==b].sum())+1)/(int(sum(flag==b))+2) for b in [False,True]]
 summaries={};allpred={};validation_safe={}
 for split in ['validation','test','OOD']:
  X,y,family=splitdata[split];preds={k:predict(models[k],X[:,cols]) for k,cols in FEATURES.items()};preds['constant']=np.full(len(X),trainy.mean())
  # Boolean array lookup must be integer indexing, not boolean selection.
  preds['static_source']=np.take(static['static_source'],(abs(X[:,0])>source_boundary).astype(int));preds['static_application']=np.take(static['static_application'],(abs(X[:,2])>app_boundary).astype(int))
  summary={'n':len(y),'events':int(sum(y)),'prevalence':float(y.mean()),'models':{}};out_of_support=np.any((X<trainX.min(0))|(X>trainX.max(0)),axis=1)
  for key,p in preds.items():
   unsupported=out_of_support.copy()
   if key=='theory':unsupported|=(1-np.sqrt(np.maximum(0,1-2*np.maximum(X[:,0],0))))>.025
   proposed=(p<=.1)&~unsupported
   if split=='validation':validation_safe[key]=bool(sum(proposed)>=32 and rate(int(sum(y[proposed])),int(sum(proposed)),one_sided=True)['interval'][1]<=.1)
   policy=np.where(unsupported,'insufficient',np.where(p>=.3,'stop',np.where(p>.1,'intervene','continue')))
   if not validation_safe.get(key,False):policy=np.where(policy=='continue','insufficient',policy)
   cont=policy=='continue';intervene=policy=='intervene';stop=policy=='stop';insuff=policy=='insufficient'
   costs=10*y*cont+stop+1.3*insuff+(.5+10*y)*intervene;mitigated=10*y*cont+stop+1.3*insuff+(.5+5*y)*intervene
   bins=[]
   for low,high in zip(np.arange(0,1,.2),np.arange(.2,1.01,.2)):
    mask=(p>=low)&(p<high if high<1 else p<=high)
    if mask.any():bins.append({'n':int(sum(mask)),'predicted':float(p[mask].mean()),'observed':rate(int(sum(y[mask])),int(sum(mask)))})
   summary['models'][key]={'brier':float(np.mean((p-y)**2)),'calibration':bins,'ECE':float(sum(b['n']/len(y)*abs(b['predicted']-b['observed']['rate']) for b in bins)),'policy_counts':{k:int(sum(policy==k)) for k in ['continue','intervene','stop','insufficient']},'continued_risk':rate(int(sum(y[cont])),int(sum(cont))),'validation_continuation_gate':validation_safe[key],'no_mitigation_cost':float(costs.mean()),'hypothetical50percent_mitigation_cost':float(mitigated.mean()),'cost_minus_stop_ci95':np.quantile([np.mean(costs[np.random.default_rng(seed('cost',split,key,j)).integers(len(costs),size=len(costs))]-1) for j in range(1000)],[.025,.975]),'baseline_always_stop_cost':1.,'baseline_always_continue_cost':float(10*y.mean()),'by_family':{f:{'n':int(sum(family==f)),'events':int(sum(y[family==f])),'brier':float(np.mean((p[family==f]-y[family==f])**2))} for f in np.unique(family)}}
  primary=(preds['combined']-y)**2-(preds['constant']-y)**2;ablate=(preds['theory']-y)**2-(preds['combined']-y)**2;g=np.random.default_rng(seed('forecast-bootstrap',split));b=np.array([primary[g.integers(len(y),size=len(y))].mean() for _ in range(2000)]);ab=np.array([ablate[g.integers(len(y),size=len(y))].mean() for _ in range(2000)])
  model=summary['models']['combined'];summary['primary_combined_minus_constant']={'difference':float(primary.mean()),'approx_adjusted_ci':np.quantile(b,[.00625,.99375]),'score_gate':bool(np.quantile(b,.99375)<-.01),'usefulness_gate':bool(model['policy_counts']['continue']/len(y)>=.2 and model['no_mitigation_cost']<.9 and model['validation_continuation_gate'])}
  summary['theory_minus_combined_secondary']={'difference':float(ablate.mean()),'pointwise_ci95':np.quantile(ab,[.025,.975])};summaries[split]=summary
  np.savez_compressed(root/f'forecast_predictions_{split}.npz',labels=y,family=family,**preds);allpred[split]=preds
 save(root/'forecast_analysis.json',{'status':'EXECUTED','registry':registry,'features':FEATURES,'models':models,'static_thresholds':{'source':source_boundary,'application':app_boundary},'splits':summaries,'selection':'combined candidate fixed before test; other comparisons secondary','intervention_mitigation_empirically_tested':False})
 print(json.dumps({s:{'n':v['n'],'events':v['events'],'primary':v['primary_combined_minus_constant']} for s,v in summaries.items()},default=lambda x:x.tolist(),indent=2))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('root');run(p.parse_args().root)
