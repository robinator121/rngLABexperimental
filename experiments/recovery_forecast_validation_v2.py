"""Finite exploratory frozen-design validation, one issuance per trajectory.

No survival guarantee, repeated-look significance or novel predictor claim.
Run --freeze first, then --run with the same new output directory.
"""
from pathlib import Path
import argparse, dataclasses, hashlib, json, resource, sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from rnglab14.forecast import (NoChangeBaselinePredictor,LinearExtrapolationPredictor,
    CUSUMChangePointPredictor,FixedHorizonPredictor,evaluate_forecast_on_trajectory,summarize_forecasting_performance)
from rnglab14.forecast.survival import KaplanMeierBaseline,CensoringKaplanMeier,ipcw_brier

SPEC={'version':2,'classification':'prospectively_frozen_exploratory_validation_not_confirmation',
      'question':'Can truthful censoring labels and a development-fit survival baseline produce interpretable heldout scores?',
      'hypothesis':'Original unresolved/all-abstention/post-event controls fail old evaluation; corrected controls match independent logic',
      'independent_unit':'one seeded trajectory, one landmark at100',
      'development_seeds':[851000,851255],'heldout_seeds':[969000,969511],
      'latent_steps':500,'landmark':100,'horizons':[100,200],
      'censoring':'independent administrative last-observed step uniform200..500; conditional on eligibility',
      'families':['linear_random_rate','step_random_onset','constant_healthy'],
      'primary_endpoints':['exact defect controls','IPCW Brier versus fixed0 and0.5','resolved/unsafe/unresolved/abstained counts'],
      'acceptance':['unresolved labels never counted safe','no post-event forecasts','undefined MAE for no predictions','KM/IPCW matches independently checked arithmetic'],
      'negative_controls':['all abstaining flat signals','insufficient censored followup','zero-horizon utility baseline'],
      'stop_conditions':['one run only','CPU60seconds','memory2GB','no retuning after heldout'],
      'responsible_verifier':'principal independent discrete-time logic and separate test arithmetic',
      'references':['https://scikit-survival.readthedocs.io/en/stable/api/generated/sksurv.metrics.brier_score.html'],
      'limits':'Scores assume this synthetic population and independent censoring. Path IDs do not by themselves establish independence.'}

def canonical(v):return json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
def sample(seed):
 r=np.random.default_rng(seed);t=np.arange(1,501);family=int(r.integers(3))
 if family==0:v=t*r.uniform(.0005,.002)+r.normal(0,.004,500)
 elif family==1:v=np.where(t>=int(r.integers(140,491)),.5,0.)
 else:v=np.full(500,.12)
 censor=int(r.integers(200,501));hits=np.flatnonzero(v>=.4);tau=int(hits[0]+1) if len(hits) else None
 if tau is not None and tau<=100:return None
 duration=(min(tau,censor) if tau is not None else censor)-100
 return v[:censor].tolist(),duration,bool(tau is not None and tau<=censor),str(seed),family

def main():
 p=argparse.ArgumentParser();p.add_argument('phase',choices=['freeze','run']);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 if a.phase=='freeze':
  a.output.mkdir(parents=True,exist_ok=False);(a.output/'spec.json').write_bytes(canonical(SPEC));print(hashlib.sha256(canonical(SPEC)).hexdigest());return
 assert json.loads((a.output/'spec.json').read_text())==SPEC
 if (a.output/'results.json').exists():raise ValueError('immutable results already exist')
 resource.setrlimit(resource.RLIMIT_CPU,(60,65));resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
 dev=[sample(s) for s in range(851000,851256)];test=[sample(s) for s in range(969000,969512)]
 dev=[x for x in dev if x];test=[x for x in test if x]
 km=KaplanMeierBaseline.fit([x[1] for x in dev],[x[2] for x in dev])
 g=CensoringKaplanMeier.fit([x[1] for x in dev],[x[2] for x in dev])
 scores={}
 for h in SPEC['horizons']:
  scores[str(h)]={name:ipcw_brier([x[1] for x in test],[x[2] for x in test],[prob]*len(test),h,g,[x[3] for x in test]) for name,prob in [('KM_development_fit',km.failure_probability(h)),('never_fail',0.),('constant_half',.5)]}
 summaries={}
 for name,model in [('linear',LinearExtrapolationPredictor(.4)),('cusum',CUSUMChangePointPredictor(.4)),('fixed50',FixedHorizonPredictor(50,.4)),('zero',FixedHorizonPredictor(0,.4)),('no_change',NoChangeBaselinePredictor(.4))]:
  fs=[]
  for values,d,e,key,family in test:fs.extend(evaluate_forecast_on_trajectory(model,values,.4,[100],name,key))
  summaries[name]=dataclasses.asdict(summarize_forecasting_performance(fs))
 def controls():
  f=evaluate_forecast_on_trajectory(NoChangeBaselinePredictor(.4),[.1]*50,.4,[10])[0]
  flat=summarize_forecasting_performance(evaluate_forecast_on_trajectory(LinearExtrapolationPredictor(.4),[0.]*50,.4,[10,20,30]))
  post=evaluate_forecast_on_trajectory(NoChangeBaselinePredictor(.4),[0.]*49+[.5]*51,.4,[20,50,60])
  return {'censor_unresolved':f.covered_by_lower_bound is None,'empty_MAE_undefined':flat.event_mae is None,'empty_calibration_undefined':flat.is_well_calibrated is None,'postevent_excluded':[x.current_t for x in post]==[20]}
 result={'spec_sha256':hashlib.sha256(canonical(SPEC)).hexdigest(),'development_trajectories':len(dev),'heldout_trajectories':len(test),'one_landmark_per_trajectory':True,'controls':controls(),'point_forecast_summaries':summaries,'survival_scores':scores,'no_predictive_superiority_claim':True}
 assert all(result['controls'].values())
 (a.output/'results.json').write_text(json.dumps(result,indent=2,allow_nan=False))
 src=Path(__file__).resolve().parents[1];paths=[Path(__file__),*sorted((src/'src/rnglab14/forecast').glob('*.py'))]
 (a.output/'manifest.json').write_text(json.dumps({'spec_sha256':result['spec_sha256'],'results_sha256':hashlib.sha256((a.output/'results.json').read_bytes()).hexdigest(),'source_sha256':{str(f.relative_to(src)):hashlib.sha256(f.read_bytes()).hexdigest() for f in paths},'classification':SPEC['classification']},indent=2))
 print(json.dumps({'controls':result['controls'],'heldout':len(test),'Brier':{h:{k:v['brier_ipcw'] for k,v in sc.items()} for h,sc in scores.items()}},indent=2))
if __name__=='__main__':main()
