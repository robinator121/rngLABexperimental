"""Source-interchangeable consumers reusing the preserved RNGLab physics code."""
from __future__ import annotations
from dataclasses import dataclass
import math
import time
from typing import Any
import numpy as np
from scipy.stats import t as student_t
from scipy.stats import beta
from rnglab.v03.physics import Independent,VMC,Ising,PIMC
from rnglab10.polynomial import PolynomialPIMC,double_well,NAMES
from rnglab10.faults import inject
from .sources import ByteSource,UniformStream,SourceExhausted,SourceBudgetExceeded,BytesSource,FileSource
from .drbg import ReseedUnavailable,HMACDRBG,ReseedingDRBG

_OPTIONS={
 'integration':{'integrand'},'pi':set(),'mcmc':{'proposal','x0','observable'},
 'vmc':{'alpha','proposal','x0','observable'},'ising':{'length','temperature','initial','observable'},
 'harmonic_pimc':{'beta','beads','centroid_probability','bead_proposal','centroid_proposal','estimator','path0','mode_proposal'},
 'double_well_pimc':{'beta','beads','width','path','observable'},
 'tilted_double_well_pimc':{'beta','beads','width','path','observable'}}

def _preflight_application(spec):
 if not isinstance(spec,dict) or spec.get('kind') not in _OPTIONS:raise ValueError('Known application object required')
 if set(spec)-{'kind'}-_OPTIONS[spec['kind']]:raise ValueError('Unknown or unsupported application options')
 beads=spec.get('beads',16)
 if type(beads) is not int or not 3<=beads<=256:raise PermissionError('APPROVAL_REQUIRED: development PIMC state is bounded to3..256 beads before allocation')
 for name in ('beta','alpha','proposal','width','temperature','bead_proposal','centroid_proposal','mode_proposal'):
  if name in spec and spec[name] is not None:
   value=spec[name]
   if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<=0:raise ValueError('Finite positive '+name+' required')
 for name in ('x0','centroid_probability'):
  if name in spec:
   value=spec[name]
   if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError('Finite numeric '+name+' required')
 if 'centroid_probability' in spec and not 0<=spec['centroid_probability']<=1:raise ValueError('centroid_probability in[0,1] required')
 if 'length' in spec and (type(spec['length']) is not int or spec['length'] not in (3,4)):raise ValueError('Ising length3 or4 required')
 for name in ('path','path0'):
  if name not in spec or spec[name] is None:continue
  value=spec[name]
  if np.isscalar(value):
   if name=='path' or isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value):raise ValueError('Invalid initial path')
  else:
   if not isinstance(value,(list,tuple,np.ndarray)) or len(value)!=beads:raise ValueError('Initial path must match bounded bead count')
   array=np.asarray(value,dtype=float)
   if array.shape!=(beads,) or not np.isfinite(array).all():raise ValueError('Finite one-dimensional initial path required')
 return 1 if spec['kind']=='integration' else 3 if 'pimc' in spec['kind'] else 2

def _source_nodes(source):
 nodes=[];seen=set()
 def visit(node):
  if id(node) in seen:return
  if not isinstance(node,ByteSource) or len(nodes)>=32:raise ValueError('Bounded ByteSource graph required')
  seen.add(id(node));nodes.append(node)
  for name in ('source','seed_source','reseed_source'):
   child=getattr(node,name,None)
   if child is not None:visit(child)
 visit(source);return nodes

def _read_path(source):
 path=[];seen=set();node=source
 while node is not None:
  if id(node) in seen or len(path)>=32:raise ValueError('Cyclic or excessive source graph')
  seen.add(id(node));path.append(node)
  node=getattr(node,'source',None) if not isinstance(node,ReseedingDRBG) else (node.reseed_source if node.interval is not None else None)
 return path

def _preflight_source(source,limits,request_bytes):
 for node in _read_path(source):
  bound=node.max_bytes
  if bound is None:
   bound=len(node._data) if isinstance(node,BytesSource) else node._identity[0] if isinstance(node,FileSource) else None
  if type(bound) is not int or not 0<=bound<=limits.max_source_bytes:raise PermissionError('APPROVAL_REQUIRED: every active development source needs a finite modest byte budget')
  raw_bound=getattr(node,'max_input_bytes',None)
  if raw_bound is not None:
   identity=getattr(node,'identity',getattr(node,'_identity',None))
   if identity is not None:raw_bound=min(raw_bound,identity[0])
   if type(raw_bound) is not int or not 0<=raw_bound<=limits.max_source_bytes:raise PermissionError('APPROVAL_REQUIRED: raw-input processing exceeds development budget')
  cap=getattr(node,'max_request_bytes',None)
  if isinstance(node,ReseedingDRBG):cap=HMACDRBG.max_request_bytes
  # Byte-preserving wrappers pass the same request through. A conditioner has
  # its own frozen upstream chunk size; do not silently split Generate calls.
  if cap is not None and request_bytes>cap:raise ValueError('Frozen chunk exceeds generator per-request limit; select a smaller chunk_rows')
  if hasattr(node,'chunk_bytes'):request_bytes=node.chunk_bytes
  elif isinstance(node,ReseedingDRBG):request_bytes=node.reseed_bytes

def _raw_count(node):
 if isinstance(node,HMACDRBG):return None
 if any(getattr(node,name,None) is not None for name in ('source','seed_source','reseed_source')):return None
 if hasattr(node,'input_bytes'):return node.input_bytes
 if hasattr(node,'input_bytes_consumed'):return node.input_bytes_consumed
 return node.bytes_consumed

def _snapshot(nodes):
 return {id(node):{'bytes':node.bytes_consumed,'raw':_raw_count(node),
  'calls':node.cycles if isinstance(node,ReseedingDRBG) else node.generate_calls if isinstance(node,HMACDRBG) else None} for node in nodes}

@dataclass(frozen=True)
class DevelopmentLimits:
 max_steps:int=16384
 max_source_bytes:int=1_048_576
 max_wall_seconds:float=30
 def __post_init__(self):
  for name in ('max_steps','max_source_bytes'):
   value=getattr(self,name)
   if type(value) is not int or value<1:raise ValueError('Positive integer development limits required')
  if not math.isfinite(self.max_wall_seconds) or self.max_wall_seconds<=0:raise ValueError('Finite wall limit required')

def development_gate(steps:int,channels:int,limits:DevelopmentLimits):
 if type(steps) is not int or steps<1:raise ValueError('steps must be a positive integer')
 if steps>limits.max_steps or steps*channels*8>limits.max_source_bytes:
  raise PermissionError('APPROVAL_REQUIRED: substantial scientific execution exceeds development-only limits; freeze a protocol and obtain explicit human authorization')

class _Pi:
 dimension=2;initialization_uniform_count=0;target=math.pi
 mixing_note='Independent unit-square indicator under IID uniforms; paired source dependence can alter the integral.'
 def advance(self,u):return 4*((u[:,0]**2+u[:,1]**2)<1).astype(float),{}

class _Polynomial:
 dimension=3;initialization_uniform_count=0;target=None
 mixing_note='Finite-bead polynomial PIMC. Burn-in, inter-well mixing and independent reference convergence remain application-specific checks.'
 def __init__(self,tilt,observable='centroid',**options):
  if tilt not in (0.,.1):raise ValueError('Only historical tested tilt=0 or0.1 available through this adapter')
  if observable not in NAMES:raise ValueError('Unknown polynomial observable')
  self.model=PolynomialPIMC(double_well(tilt=tilt),backend='python',**options);self.observable=observable
 def advance(self,u):
  out=self.model.advance(u);return out[:,NAMES.index(self.observable)],{name:out[:,i] for i,name in enumerate(NAMES)}

def make_application(spec:dict):
 _preflight_application(spec)
 spec=dict(spec);kind=spec.pop('kind')
 if kind=='integration':return Independent(backend='numpy',**spec)
 if kind=='pi':
  if spec:raise ValueError('Pi workload has no adjustable options')
  return _Pi()
 if kind=='mcmc':return VMC(alpha=.5,observable=spec.pop('observable','x2'),backend='python',**spec)
 if kind=='vmc':return VMC(backend='python',**spec)
 if kind=='ising':return Ising(length=spec.pop('length',3),backend='python',**spec)
 if kind=='harmonic_pimc':return PIMC(sampler='normal_mode',mode_selection='cv_optimal',backend='python',**spec)
 if kind in {'double_well_pimc','tilted_double_well_pimc'}:return _Polynomial(0. if kind=='double_well_pimc' else .1,**spec)
 raise ValueError('Unknown application: '+kind)

class _Moments:
 def __init__(self):self.n=0;self.mean=0.;self.m2=0.
 def add(self,x):
  self.n+=1;delta=x-self.mean;self.mean+=delta/self.n;self.m2+=delta*(x-self.mean)
 @property
 def variance(self):return self.m2/(self.n-1) if self.n>1 else None

class StreamingSummary:
 """Welford moments + fixed nonoverlapping batch means; no unbounded series."""
 def __init__(self,batch_size=32):
  if type(batch_size) is not int or batch_size<2:raise ValueError('batch_size>=2 required')
  self.batch_size=batch_size;self.values=_Moments();self.batches=_Moments();self.pending=0.;self.pending_n=0;self.prefix=_Moments();self.pending_values=[]
 def add(self,values):
  for value in np.asarray(values).ravel():
   if not np.isfinite(value):raise ValueError('Nonfinite application output')
   value=float(value);self.values.add(value);self.pending+=value;self.pending_n+=1;self.pending_values.append(value)
   if self.pending_n==self.batch_size:
    self.batches.add(self.pending/self.batch_size)
    for completed_value in self.pending_values:self.prefix.add(completed_value)
    self.pending=self.pending_n=0;self.pending_values=[]
 def report(self,reference=None):
  n=self.values.n;variance=self.values.variance;mcse=interval=ess=None
  reason='insufficient complete nonoverlapping batches'
  if self.batches.n>=8 and self.batches.variance is not None and self.batches.variance>0:
   mcse=math.sqrt(self.batches.variance/self.batches.n)
   half=float(student_t.ppf(.975,self.batches.n-1))*mcse
   interval=[self.batches.mean-half,self.batches.mean+half]
   ess=min(float(self.prefix.n),max(1.,self.prefix.variance/mcse**2)) if self.prefix.variance is not None and self.prefix.variance>0 else None
   reason='approximate batch-means uncertainty; stationarity/mixing and batch-size adequacy unverified'
  return {'n':n,'estimate':self.values.mean if n else None,'sample_variance':variance,
   'mcse':mcse,'interval_95':interval,'interval_estimate':self.batches.mean if interval else None,'interval_sample_count':self.prefix.n,'interval_uses_complete_batch_prefix':True,'ess_estimate':ess,'effective_samples_are_estimated':True,
   'batch_size':self.batch_size,'complete_batches':self.batches.n,'incomplete_tail_samples':self.pending_n,
   'uncertainty_status':reason,'reference':reference,
   'signed_error':self.values.mean-reference if n and reference is not None else None,
   'interval_contains_reference':interval[0]<=reference<=interval[1] if interval and reference is not None else None}

def run_application(source:ByteSource,spec:dict,*,steps:int=1024,burn_in:int=0,chunk_rows:int=128,
                    limits:DevelopmentLimits|None=None,monitor=None,post_uniform_fault:dict|None=None,
                    run_id:str|None=None,independent_group_id:str|None=None,observation_callback=None):
 """Small development execution only. A future confirmation runner needs approval.

 DRBG calls occur once per random chunk; chunk_rows is therefore frozen and
 reported as part of a reseeded-DRBG design. Input byte consumption and sampler
 iterations remain separate. No raw physical data or private keys are exported.
 """
 limits=limits or DevelopmentLimits();channels=_preflight_application(spec)
 development_gate(steps,channels,limits)
 if type(burn_in) is not int or not 0<=burn_in<steps:raise ValueError('burn_in must be integer in [0,steps)')
 if type(chunk_rows) is not int or not 1<=chunk_rows<=4096:raise ValueError('chunk_rows must be in [1,4096]')
 for name,value in [('run_id',run_id),('independent_group_id',independent_group_id)]:
  if value is not None and (not isinstance(value,str) or not value.strip() or len(value)>512):raise ValueError('Bounded nonempty '+name+' required')
 if (run_id is None)!=(independent_group_id is None):raise ValueError('Supply both run and independent-group IDs')
 if observation_callback is not None and not callable(observation_callback):raise ValueError('observation_callback must be callable')
 nodes=_source_nodes(source);_preflight_source(source,limits,min(steps,chunk_rows)*channels*8)
 if monitor is not None and (getattr(monitor,'channels',None)!=channels or getattr(monitor,'horizon',0)-getattr(monitor,'rows',0)<steps):raise ValueError('Monitor channels/horizon must cover the run before input consumption')
 started=time.monotonic();model=make_application(spec)
 if post_uniform_fault is not None and post_uniform_fault.get('kind')=='temporal_acceptance':
  raise ValueError('Legacy temporal filter resets per call; use a chunk-invariant source-fault adapter or freeze a separately reviewed stateful filter')
 if post_uniform_fault is not None and model.dimension!=3:raise ValueError('Legacy numerical uniform faults explicitly require three channels')
 stream=UniformStream(source,max_draws=steps*model.dimension);summary=StreamingSummary();events=[]
 before=_snapshot(nodes);initial_bytes=source.bytes_consumed;done=0;status='DEVELOPMENT_COMPLETE';error=None;event_count=0
 while done<steps:
  if time.monotonic()-started>limits.max_wall_seconds:status='TIME_LIMIT';break
  count=min(chunk_rows,steps-done)
  if source.bytes_consumed-initial_bytes+count*model.dimension*8>limits.max_source_bytes:status='SOURCE_BUDGET';break
  try:u=stream.random(count,model.dimension)
  except SourceExhausted as exc:status='SOURCE_EXHAUSTED';error=str(exc);break
  except ReseedUnavailable as exc:status='RESEED_SOURCE_EXHAUSTED';error=str(exc);break
  except SourceBudgetExceeded as exc:status='SOURCE_BUDGET';error=str(exc);break
  if post_uniform_fault:u=inject(u,post_uniform_fault['kind'],post_uniform_fault.get('value',0.))
  if monitor:
   new=monitor.update(u);event_count+=len(new);events=(events+new)[-64:]
  values,_=model.advance(u);discard=max(0,min(count,burn_in-done));summary.add(values[discard:]);done+=count
  if observation_callback is not None:observation_callback(done,values[discard:],getattr(model,'target',None))
  for node in _read_path(source):
   if hasattr(node,'record_application_iterations'):
    node.record_application_iterations(count);break
 report=summary.report(getattr(model,'target',None))
 after=_snapshot(nodes)
 accounting=[{'layer':i,'source_id':node.source_id,'type':type(node).__name__,
  'output_bytes_before':before[id(node)]['bytes'],'output_bytes_after':after[id(node)]['bytes'],
  'output_bytes_during_run':after[id(node)]['bytes']-before[id(node)]['bytes'],
  'raw_bytes_during_run':None if after[id(node)]['raw'] is None else after[id(node)]['raw']-before[id(node)]['raw'],
  'generator_calls_during_run':None if after[id(node)]['calls'] is None else after[id(node)]['calls']-before[id(node)]['calls']} for i,node in enumerate(nodes)]
 calls=[item['generator_calls_during_run'] for item in accounting if item['generator_calls_during_run'] is not None]
 return {'schema_version':1,'classification':'development_smoke_only','scientific_confirmation':False,
  'run_id':run_id,'independent_group_id':independent_group_id,
  'application':spec,'status':status,'error':error,'requested_steps':steps,'simulation_steps':done,'burn_in':burn_in,
  'generated_uniforms':stream.uniforms_generated,'uniform_channels':model.dimension,'uniform_bytes_consumed':source.bytes_consumed-initial_bytes,
  'raw_source_bytes_consumed':sum(item['raw_bytes_during_run'] or 0 for item in accounting),
  'raw_source_accounting_scope':'Leaf input bytes during this run; initial seeds/preexisting reads excluded; direct HMAC seed material is outside ByteSource accounting',
  'layered_source_accounting':accounting,
  'generator_calls':calls[0] if calls else None,'all_generator_calls_during_run':sum(calls) if calls else None,'chunk_rows_frozen':chunk_rows,
  'source_provenance':source.provenance(),'elapsed_seconds':time.monotonic()-started,'summary':report,
  'monitor_event_count':event_count,'recent_monitor_events':events,'mixing_note':model.mixing_note,
  'monitor_stage':'post_uniform_fault' if post_uniform_fault else 'consumer_uniforms',
  'initialization':{'method':'configured deterministic constructor state','random_initialization_uniforms_consumed':0,'burn_in_steps':burn_in},
  'reference_status':'finite discretization/analytic legacy reference' if getattr(model,'target',None) is not None else 'not supplied; bias/MSE/coverage cannot be inferred',
  'limitations':['This short run does not establish source health, security or a prospective reliability horizon.',
    'Batch-means uncertainty and ESS require independent stationarity/mixing checks.',
    'Changes in target, sampling efficiency and uncertainty failure are distinct mechanisms.']}


def aggregate_replicates(records:list[dict],*,reference:float,tolerance:float):
 if not records or not math.isfinite(reference) or not math.isfinite(tolerance) or tolerance<=0:raise ValueError('Finite reference and positive tolerance required')
 if len({id(r) for r in records})!=len(records):raise ValueError('Duplicate replicate object is not an independent run')
 identities=[(r.get('run_id'),r.get('independent_group_id')) for r in records]
 independent=all(isinstance(a,str) and a.strip() and isinstance(b,str) and b.strip() for a,b in identities)
 if independent and (len({a for a,b in identities})!=len(records) or len({b for a,b in identities})!=len(records)):raise ValueError('Duplicate run/group IDs imply pseudoreplication')
 estimates=[];cover=[]
 for record in records:
  if record['summary']['estimate'] is None or not math.isfinite(record['summary']['estimate']):raise ValueError('Incomplete/nonfinite replicate estimate')
  estimates.append(record['summary']['estimate'])
  if record['summary']['interval_95'] is not None:
   lo,hi=record['summary']['interval_95'];cover.append(lo<=reference<=hi)
 errors=np.asarray(estimates)-reference;n=len(errors);fail=int(np.sum(abs(errors)>tolerance))
 return {'replicates':n,'reference':reference,'tolerance':tolerance,'estimated_bias':float(np.mean(errors)),
  'variance_of_estimates':float(np.var(estimates,ddof=1)) if n>1 else None,'mse':float(np.mean(errors**2)),
  'coverage_rate':sum(cover)/len(cover) if cover else None,'coverage_denominator':len(cover),
  'failure_rate':fail/n,'failure_count':fail,
  'failure_rate_pointwise_95_interval':[0. if fail==0 else float(beta.ppf(.025,fail,n-fail+1)),1. if fail==n else float(beta.ppf(.975,fail+1,n-fail))] if independent else None,
  'inference_status':'DECLARED_INDEPENDENT_GROUPS' if independent else 'UNAVAILABLE_INDEPENDENCE_METADATA; descriptive summaries only',
  'assumptions':'Replicates must be independent; reference/tolerance are frozen externally. Descriptive pointwise intervals do not establish calibrated prediction.'}
