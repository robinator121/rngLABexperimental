"""Interpretable prospective RUS baseline; source exposure is a separate quantity.

No historical confirmations are loaded. Each fit/calibration/evaluation split is
explicit, immutable and disjoint by trajectory identity. Prediction defaults to
unavailable until a held-out, pre-registered validation has passed.
"""
from __future__ import annotations
from dataclasses import asdict,dataclass
import hashlib
import json
import math
from typing import Any
import numpy as np
from scipy.stats import beta

UNITS={'raw_bits','conditioned_bits','generated_uniforms','drbg_bytes','generator_calls','simulation_steps','effective_samples'}

def stable_hash(value: Any) -> str:
 return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()

def integer(value,name,minimum=0):
 if type(value) is not int or value<minimum:raise ValueError(f'{name} must be an integer >= {minimum}')
 return value

def finite(value,name):
 if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError(f'{name} must be finite')
 return float(value)

def _safe_operations_before_failure(lower):
 """Completed integer operations strictly before a first-failure time bound."""
 return max(0,math.ceil(lower)-1)

@dataclass(frozen=True)
class ApplicationContract:
 application: str
 observable: str
 absolute_error_threshold: float
 operations_unit: str='simulation_steps'
 feature_definition: str='source severity and severity change per generator call'
 max_overprediction_risk: float=.05
 version: int=1
 def __post_init__(self):
  if not self.application or not self.observable or not self.feature_definition:raise ValueError('Specify application, observable and feature definition')
  if self.operations_unit not in UNITS:raise ValueError('Unknown operation unit; no implicit conversion')
  if finite(self.absolute_error_threshold,'error threshold')<=0:raise ValueError('positive error threshold required')
  if not 0<finite(self.max_overprediction_risk,'risk')<1:raise ValueError('risk must lie in (0,1)')
  integer(self.version,'version',1)
 @property
 def digest(self):return stable_hash(asdict(self))

@dataclass(frozen=True)
class Trajectory:
 trajectory_id: str
 severity: float
 degradation_rate: float
 failure_operations: float|None
 censor_operations: float|None=None
 def __post_init__(self):
  if not self.trajectory_id:raise ValueError('trajectory_id required')
  finite(self.severity,'severity');finite(self.degradation_rate,'degradation_rate')
  if self.failure_operations is None:
   if self.censor_operations is None or finite(self.censor_operations,'censor')<0:raise ValueError('Censored trajectory requires follow-up length')
  else:
   if finite(self.failure_operations,'failure operations')<0:raise ValueError('failure operations must be nonnegative')
   if self.censor_operations is not None:raise ValueError('Trajectory cannot simultaneously declare failure and right censoring')

@dataclass(frozen=True)
class Dataset:
 dataset_id: str
 split: str
 contract_hash: str
 trajectories: tuple[Trajectory,...]
 provenance: str
 synthetic: bool=False
 def __post_init__(self):
  if self.split not in {'development','calibration','confirmation'}:raise ValueError('Unknown split')
  if not self.dataset_id or not self.provenance or not self.contract_hash:raise ValueError('Dataset identity, provenance and contract hash required')
  if not isinstance(self.trajectories,tuple):raise ValueError('Freeze trajectories as a tuple')
  ids=[t.trajectory_id for t in self.trajectories]
  if not ids or len(ids)!=len(set(ids)):raise ValueError('Nonempty, unique trajectory identities required')
  if type(self.synthetic) is not bool:raise ValueError('synthetic must be boolean')
 @property
 def digest(self):return stable_hash(asdict(self))

@dataclass(frozen=True)
class ValidationRule:
 max_unsafe_overprediction: float
 confidence: float=.95
 minimum_trajectories: int=30
 def __post_init__(self):
  if not 0<finite(self.max_unsafe_overprediction,'max unsafe risk')<1:raise ValueError('risk in (0,1) required')
  if not 0<finite(self.confidence,'confidence')<1:raise ValueError('confidence in (0,1) required')
  integer(self.minimum_trajectories,'minimum_trajectories',1)

class ReliabilityPredictor:
 """Linear failure-delay baseline + split-conformal one-sided time bound.

 Failure_operations is the first threshold-crossing time in operations *after
 the feature observation*. A crossing at operation N excludes N from safe work. The model
 never maps physical entropy or reseed exposure directly to application failure.
 Right-censored data are retained but unsupported by this simple time-regression
 calibration: fitting is refused rather than silently discarding them.
 """
 def __init__(self,contract:ApplicationContract):
  self.contract=contract;self.coefficients=None;self.development=None;self.calibration=None
  self.residuals=None;self.correction=None;self.rule=None;self.validation=None;self.confirmation=None;self._ids=set();self._dataset_ids=set();self.domain=None
 def _accept(self,data,split):
  if data.split!=split or data.contract_hash!=self.contract.digest:raise ValueError('Dataset split/application contract mismatch')
  if self.development is not None and data.synthetic!=self.development.synthetic:raise ValueError('Synthetic and nonsynthetic validation populations cannot be silently interchanged')
  ids={t.trajectory_id for t in data.trajectories}
  if ids&self._ids or data.dataset_id in self._dataset_ids:raise ValueError('Train/calibration/confirmation leakage detected')
  if any(t.failure_operations is None for t in data.trajectories):raise ValueError('INSUFFICIENT_EVIDENCE: censored trajectories require a separately validated survival model; no silent omission')
  return ids
 def fit(self,data:Dataset):
  if self.development is not None:raise ValueError('Model already frozen; fit a new model version')
  ids=self._accept(data,'development')
  if len(ids)<4:raise ValueError('At least four independent development trajectories required')
  features=np.array([[t.severity,t.degradation_rate] for t in data.trajectories])
  x=np.column_stack([np.ones(len(features)),features]);y=np.array([t.failure_operations for t in data.trajectories])
  coefficients,_,rank,_=np.linalg.lstsq(x,y,rcond=None)
  if rank!=3 or not np.isfinite(coefficients).all():raise ValueError('Unidentified severity/rate response; more development variation required')
  self.coefficients=coefficients;self.domain=np.array([features.min(axis=0),features.max(axis=0)])
  self.development=data;self._ids.update(ids);self._dataset_ids.add(data.dataset_id)
  return self
 def _point(self,severity,rate):
  values=np.array([finite(severity,'severity'),finite(rate,'degradation_rate')])
  if self.coefficients is None:raise ValueError('Unfitted model')
  if np.any(values<self.domain[0]) or np.any(values>self.domain[1]):raise ValueError('OUT_OF_DOMAIN: no extrapolation beyond development bounds')
  return max(0.,float(np.r_[1.,values]@self.coefficients))
 def calibrate(self,data:Dataset):
  if self.development is None:raise ValueError('Fit development data first')
  if self.calibration is not None or self.validation is not None:raise ValueError('Calibration is frozen; create a new predictor version')
  ids=self._accept(data,'calibration');alpha=self.contract.max_overprediction_risk
  n=len(data.trajectories);rank=math.ceil((n+1)*(1-alpha))
  if rank>n:raise ValueError('INSUFFICIENT_CALIBRATION: finite-sample conformal rank exceeds calibration size')
  residuals=np.array([self._point(t.severity,t.degradation_rate)-t.failure_operations for t in data.trajectories])
  self.residuals=np.sort(residuals);self.correction=float(self.residuals[rank-1]);self.calibration=data
  self._ids.update(ids);self._dataset_ids.add(data.dataset_id)
  return self
 def freeze_validation(self,rule:ValidationRule):
  if self.calibration is None or self.rule is not None:raise ValueError('Calibrate once, then freeze one held-out acceptance rule')
  self.rule=rule
 def evaluate_heldout(self,data:Dataset):
  if self.rule is None or self.validation is not None:raise ValueError('Freeze validation rule before the first held-out evaluation; no post-outcome tuning')
  ids=self._accept(data,'confirmation')
  n=len(ids);over=0;errors=[];point_errors=[]
  for trajectory in data.trajectories:
   point=self._point(trajectory.severity,trajectory.degradation_rate);lower=max(0.,point-self.correction)
   excess=lower-trajectory.failure_operations;over+=excess>1e-10;errors.append(max(0.,excess));point_errors.append(point-trajectory.failure_operations)
  tail=1-self.rule.confidence
  upper=1. if over==n else float(beta.ppf(1-tail,over+1,n-over))
  self.validation={'dataset_id':data.dataset_id,'dataset_hash':data.digest,'n':n,'unsafe_overpredictions':over,
   'observed_unsafe_rate':over/n,'one_sided_binomial_upper':upper,'confidence':self.rule.confidence,
   'accepted':n>=self.rule.minimum_trajectories and upper<=self.rule.max_unsafe_overprediction,
   'rule':asdict(self.rule),'max_unsafe_excess_operations':max(errors),'point_rmse':float(np.sqrt(np.mean(np.square(point_errors)))),
   'synthetic':data.synthetic,'censoring':'No censored trajectories accepted by this baseline',
   'limits':'Pointwise binomial validation on independent held-out trajectories; no guarantee against future population shift.'}
  self.confirmation=data;self._ids.update(ids);self._dataset_ids.add(data.dataset_id)
  return dict(self.validation)
 def predict(self,severity,rate,*,contract:ApplicationContract|None=None,health_state='unknown',risk_horizon=None,require_confirmation=True,residual_invariance_assumed=False):
  base={'source_health_state':health_state,'application':self.contract.application,'observable':self.contract.observable,
   'unit':self.contract.operations_unit,'failure_threshold':self.contract.absolute_error_threshold,
   'remaining_safe_operations':None,'predicted_failure_risk':None,'risk_upper':None,'confidence':None,
   'not_interchangeable_with':'source exposure, entropy or cryptographic security',
   'future_degradation_assumption':'Feature severity/rate remains represented by the development/calibration trajectory population.'}
  if contract is not None and contract.digest!=self.contract.digest:return {**base,'status':'NO_PREDICTION','reason':'application contract mismatch'}
  if self.calibration is None:return {**base,'status':'INSUFFICIENT_EVIDENCE','reason':'Calibration unavailable; no supported application reliability prediction'}
  if require_confirmation and (self.validation is None or not self.validation['accepted']):return {**base,'status':'INSUFFICIENT_EVIDENCE','reason':'Held-out prospective validation unavailable or failed'}
  try:point=self._point(severity,rate)
  except ValueError as exc:return {**base,'status':'NO_PREDICTION','reason':str(exc)}
  lower=max(0.,point-self.correction)
  result={**base,'status':'VALIDATED_PROTOTYPE' if require_confirmation else 'EXPLORATORY_CALIBRATED',
   'point_failure_delay_operations':point,'conservative_failure_delay_lower':lower,'remaining_safe_operations':_safe_operations_before_failure(lower),
   'marginal_overprediction_target':self.contract.max_overprediction_risk,'calibration_n':len(self.residuals),
   'bound_assumptions':['Independent trajectories','Exchangeable calibration/new residuals','Identical application, threshold, unit and feature definition'],
   'model_hash':self.artifact()['model_hash'],'independent_review_required':True}
  result['validation_is_synthetic']=None if self.validation is None else self.validation['synthetic']
  result['bound_scope']='Marginal across exchangeable new trajectories; not a conditional guarantee for each feature value'
  if risk_horizon is not None:
   horizon=finite(risk_horizon,'risk horizon')
   if horizon<0:raise ValueError('risk horizon must be nonnegative')
   empirical=float(np.mean(self.residuals>=point-horizon));delta=.05
   band=math.sqrt(math.log(2/delta)/(2*len(self.residuals)))
   result.update(residual_tail_proxy=empirical,residual_tail_upper=min(1.,empirical+band),residual_band_confidence=1-delta,
    risk_horizon_operations=horizon,risk_method='Empirical residual tail with DKW uniform CDF band under IID representative residuals; conditional risk additionally requires feature-invariant residual law',
    conditional_risk_status='UNAVAILABLE: feature-invariant residual model not asserted')
   if residual_invariance_assumed is True:
    result.update(predicted_failure_risk=empirical,risk_upper=min(1.,empirical+band),confidence=1-delta,
     conditional_risk_status='MODEL_ASSUMED: feature-invariant residual law asserted by caller; not established by DKW')
   elif residual_invariance_assumed is not False:raise ValueError('Residual invariance assertion must be an explicit boolean')
  return result
 def artifact(self):
  value={'schema_version':1,'algorithm':'linear-severity-rate/split-conformal-lower/v1','contract':asdict(self.contract),
   'coefficients':self.coefficients.tolist() if self.coefficients is not None else None,
   'domain':self.domain.tolist() if self.domain is not None else None,'correction':self.correction,
   'residuals':self.residuals.tolist() if self.residuals is not None else None,
   'development_hash':self.development.digest if self.development else None,'calibration_hash':self.calibration.digest if self.calibration else None,
   'development_dataset':asdict(self.development) if self.development else None,
   'calibration_dataset':asdict(self.calibration) if self.calibration else None,
   'validation_rule':asdict(self.rule) if self.rule else None,'validation':self.validation,
   'confirmation_dataset':asdict(self.confirmation) if self.confirmation else None}
  return {**value,'model_hash':stable_hash(value)}
 @classmethod
 def restore(cls,artifact):
  value=dict(artifact);expected=value.pop('model_hash')
  if stable_hash(value)!=expected:raise ValueError('Model artifact checksum mismatch')
  if value.get('schema_version')!=1 or value.get('algorithm')!='linear-severity-rate/split-conformal-lower/v1':raise ValueError('Unknown predictor artifact schema/algorithm')
  model=cls(ApplicationContract(**value['contract']))
  def dataset(raw):
   return Dataset(**{**raw,'trajectories':tuple(Trajectory(**t) for t in raw['trajectories'])})
  if value['development_dataset']:model.fit(dataset(value['development_dataset']))
  if value['calibration_dataset']:model.calibrate(dataset(value['calibration_dataset']))
  if value['validation_rule']:model.freeze_validation(ValidationRule(**value['validation_rule']))
  if value['validation']:
   if not value.get('confirmation_dataset'):raise ValueError('Held-out evidence required to restore validated predictor')
   model.evaluate_heldout(dataset(value['confirmation_dataset']))
  elif value.get('confirmation_dataset'):raise ValueError('Confirmation dataset without validation outcome')
  # Integrity hashes do not authenticate scientific claims. Recompute every
  # fitted/calibrated/validated quantity from the embedded immutable splits.
  if stable_hash(model.artifact())!=stable_hash(artifact):raise ValueError('Restored model differs from recomputed frozen evidence')
  return model


def source_exposure_horizon(*,current_cycle:int,interval:int|None,source_offset:int,reseed_bytes:int,defect_segments:list[dict]):
 """Exact fixed-consumption future-reseed exposure; never an application prediction.

 A cycle is one configured generator call. The initial state serves N calls; reseeding occurs before call N+1.
 current_cycle is completed calls.
 source_offset is the next unread source byte, after initialization/previous
 reseeds. Segments are half-open byte intervals [start,end).
 """
 integer(current_cycle,'current_cycle');integer(source_offset,'source_offset');integer(reseed_bytes,'reseed_bytes',1)
 for segment in defect_segments:
  integer(segment['start'],'start');integer(segment['end'],'end')
  if segment['end']<=segment['start']:raise ValueError('Nonempty ordered half-open defect segment required')
 common={'quantity':'source_exposure_horizon','cycle_unit':'generator_calls','source_unit':'bytes',
  'application_reliability_horizon':None,'warning':'Incorporating degraded reseed input does not itself establish defective DRBG output or numerical failure.'}
 if interval is None:return {**common,'status':'NO_FUTURE_RESEED','next_reseed_cycle':None,'remaining_calls_to_exposure':None,'scope':'Future reseeding only; initial-state defects are separate.'}
 integer(interval,'interval',1);next_cycle=(((max(1,current_cycle)-1)//interval)+1)*interval
 future=[s for s in defect_segments if s['end']>source_offset]
 if not future:return {**common,'status':'NO_IDENTIFIED_FUTURE_EXPOSURE','next_reseed_cycle':next_cycle,'next_reseed_before_generate_call':next_cycle+1,'calls_to_next_reseed':next_cycle-current_cycle,'remaining_calls_to_exposure':None}
 start=min(max(source_offset,s['start']) for s in future);index=(start-source_offset)//reseed_bytes
 cycle=next_cycle+index*interval;segment_start=source_offset+index*reseed_bytes
 return {**common,'status':'SCHEDULED_EXPOSURE','next_reseed_cycle':next_cycle,'next_reseed_before_generate_call':next_cycle+1,'calls_to_next_reseed':next_cycle-current_cycle,
  'exposure_cycle':cycle,'exposure_before_generate_call':cycle+1,'remaining_calls_to_exposure':cycle-current_cycle,'reseed_source_segment':[segment_start,segment_start+reseed_bytes]}
