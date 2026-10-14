"""Bounded development CLI. No command here authorizes confirmation or provisioning."""
from __future__ import annotations
import argparse
from dataclasses import asdict
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import sys
import time
import uuid
import numpy as np
from . import __version__
from .configuration import load_json,source_from_spec
from .applications import run_application,DevelopmentLimits
from .sources import SyntheticSource,BytesSource,iter_chunks
from .statistics import evaluate_bits,HealthTests
from .monitoring import StreamingMonitor
from .conditioning import compare_treatments
from .crypto_workloads import nonce_collision_experiment,stream_repetition_experiment,ed25519_key_experiment
from .reliability import ApplicationContract,ReliabilityPredictor,source_exposure_horizon
from .events import ApplicationObservation,RiskForecast,ResearchRun,evaluate_run,aggregate_runs


def environment():
 versions={}
 for name in ['numpy','scipy','cryptography','pytest']:
  try:versions[name]=importlib.metadata.version(name)
  except importlib.metadata.PackageNotFoundError:versions[name]=None
 return {'version':__version__,'python':sys.version,'platform':platform.platform(),'packages':versions}


def write_output(path,value):
 encoded=json.dumps(value,indent=2,allow_nan=False)+'\n'
 path=Path(path).resolve()
 historical=Path(__file__).resolve().parents[2]/'reports/HISTORICAL_BASELINE_MANIFEST.json'
 if historical.exists():
  root=Path(load_json(historical)['root']).resolve()
  if path.is_relative_to(root):raise PermissionError('Historical research tree is read-only')
 if path.exists():raise FileExistsError('Evidence already exists; select a new output path')
 path.parent.mkdir(parents=True,exist_ok=True)
 # O_EXCL reserves the destination; no existing research record is replaced.
 with path.open('x') as stream:stream.write(encoded)
 return str(path)


def analyze(source,*,total_bytes=8192,window_bytes=1024,min_entropy=None):
 if type(total_bytes) is not int or not 1<=total_bytes<=1_048_576:raise PermissionError('APPROVAL_REQUIRED: analysis limited to1MiB during development')
 if type(window_bytes) is not int or not 16<=window_bytes<=8192:raise ValueError('window_bytes in [16,8192] required')
 health=HealthTests(min_entropy=min_entropy) if min_entropy is not None else None
 windows=[];health_events=[];health_event_count=0;start=source.bytes_consumed
 for data in iter_chunks(source,chunk_bytes=window_bytes,total_bytes=total_bytes):
  bits=np.unpackbits(np.frombuffer(data,dtype=np.uint8));tests=evaluate_bits(bits)
  if health:
   new_events=health.update(bits);health_event_count+=len(new_events);health_events=(health_events+new_events)[-100:]
  windows.append({'source_end_byte':source.offset,'bytes':len(data),'statistical_tests':tests})
 return {'classification':'development_diagnostics','source':source.provenance(),'bytes_analyzed':source.bytes_consumed-start,
  'windows':windows,'health_test_events':health_events,'health_event_count':health_event_count,'health_events_omitted':health_event_count-len(health_events),'health_status':'90B_INFORMED_ASSUMED_ENTROPY' if health else 'NOT_RUN_NO_ENTROPY_ASSUMPTION',
  'health_assumed_min_entropy_per_bit':min_entropy,'cryptographic_security':'NOT_ESTABLISHED','physical_entropy':'NOT_ESTABLISHED',
  'limits':'Three SP800-22 tests only. Window p-values are not a calibrated whole-run alarm or full entropy assessment.'}


def evaluate_event_spec(path):
 if Path(path).stat().st_size>8_388_608:raise ValueError('Event import limited to8MiB')
 config=load_json(path)
 if config.get('schema_version')!=1:raise ValueError('Event schema_version1 required')
 runs=[]
 if not isinstance(config.get('runs'),list) or not 1<=len(config['runs'])<=10000:raise ValueError('Bounded event runs list required')
 for row in config['runs']:
  fields=dict(row)
  fields['observations']=tuple(ApplicationObservation(**o) for o in fields['observations'])
  fields['forecasts']=tuple(RiskForecast(**f) for f in fields.get('forecasts',[]))
  runs.append(ResearchRun(**fields))
 result={'classification':'imported_observed_event_analysis','scientific_confirmation':False,
  'spec_sha256':hashlib.sha256(Path(path).read_bytes()).hexdigest(),
  'runs':[evaluate_run(run,config.get('warning_probability',.5)) for run in runs],
  'provenance_warning':'Imported calibration and independence declarations require independent audit; software does not authenticate them.'}
 if config.get('aggregate',False):
  result['aggregate']=aggregate_runs(runs,config.get('warning_probability',.5),
   forecast_issued_at=config.get('forecast_issued_at'),forecast_horizon=config.get('forecast_horizon'))
 return result


def smoke(output_dir):
 directory=Path(output_dir);directory.mkdir(parents=True,exist_ok=False)
 source_spec={'kind':'synthetic','seed':20261016,'max_bytes':131072}
 applications=[];started=time.monotonic()
 for kind in ['integration','pi','mcmc','ising','vmc','harmonic_pimc','double_well_pimc','tilted_double_well_pimc']:
  with source_from_spec(source_spec) as source:
   model_channels=1 if kind=='integration' else 2 if kind in {'pi','mcmc','ising','vmc'} else 3
   monitor=StreamingMonitor(window=64,step=64,channels=model_channels,horizon=256)
   applications.append(run_application(source,{'kind':kind},steps=256,chunk_rows=64,monitor=monitor))
 with source_from_spec(source_spec) as source:statistics=analyze(source,total_bytes=2048,window_bytes=1024,min_entropy=1.)
 with source_from_spec(source_spec) as source:
  treatment=compare_treatments(source,input_bytes=2048)
  comparisons={name:{'origin':item['origin'],'transformation':item['report'],
   'diagnostics':analyze(item['source'],total_bytes=max(1,item['report']['output_bytes']),window_bytes=256)} for name,item in treatment.items() if item['report']['output_bytes']}
 with source_from_spec(source_spec) as source:nonces=nonce_collision_experiment(source,count=32)
 with source_from_spec(source_spec) as source:keys=ed25519_key_experiment(source,count=4)
 c=ApplicationContract('harmonic_pimc','centroid',.05)
 prediction=ReliabilityPredictor(c).predict(.1,.01,health_state='uncalibrated diagnostics')
 report={'status':'DEVELOPMENT_SMOKE_COMPLETE','environment':environment(),'seed':20261016,'source_spec':source_spec,
  'scientific_confirmation':False,'protocol_status':'DRAFT_NOT_FROZEN','elapsed_seconds':time.monotonic()-started,
  'applications':applications,'diagnostics':statistics,'matched_raw_debiased_conditioned':comparisons,
  'crypto_test_only':{'nonces':nonces,'ed25519':keys},'reliability':prediction,'api_calls':0,'gpu_used':False}
 write_output(directory/'report.json',report)
 return {'report':str((directory/'report.json').resolve()),'status':report['status'],'applications':len(applications),'api_calls':0}


def main(argv=None):
 p=argparse.ArgumentParser(prog='rnglab11',description='RNGLab v0.11 development framework; confirmation requires separate human authorization')
 sub=p.add_subparsers(dest='command',required=True)
 sub.add_parser('status')
 sm=sub.add_parser('smoke');sm.add_argument('--output-dir',type=Path,default=None)
 run=sub.add_parser('run');run.add_argument('spec',type=Path);run.add_argument('--output',type=Path,required=True)
 an=sub.add_parser('analyze');an.add_argument('source_spec',type=Path);an.add_argument('--bytes',type=int,default=8192);an.add_argument('--window-bytes',type=int,default=1024);an.add_argument('--min-entropy',type=float);an.add_argument('--output',type=Path,required=True)
 pred=sub.add_parser('predict');pred.add_argument('--artifact',type=Path);pred.add_argument('--severity',type=float,required=True);pred.add_argument('--rate',type=float,default=0.);pred.add_argument('--risk-horizon',type=float)
 re=sub.add_parser('reseed-status');re.add_argument('--interval',type=int);re.add_argument('--cycle',type=int,default=0);re.add_argument('--source-offset',type=int,default=48);re.add_argument('--reseed-bytes',type=int,default=32);re.add_argument('--defect-start',type=int);re.add_argument('--defect-end',type=int)
 infra=sub.add_parser('infrastructure');infra.add_argument('--output',type=Path)
 ev=sub.add_parser('events');ev.add_argument('spec',type=Path);ev.add_argument('--output',type=Path,required=True)
 dash=sub.add_parser('dashboard');dash.add_argument('--config',type=Path,required=True)
 args=p.parse_args(argv)
 try:
  if args.command=='status':result={'environment':environment(),'implementation':'DEVELOPMENT_RELEASE','protocol_status':'DRAFT_NOT_FROZEN','scientific_confirmation_authorized':False,'review_status':'STANDALONE_LOCAL','substantial_experiments':'APPROVAL_REQUIRED'}
  elif args.command=='smoke':result=smoke(args.output_dir or Path('runs')/('smoke-'+time.strftime('%Y%m%d-%H%M%S')+'-'+uuid.uuid4().hex[:6]))
  elif args.command=='run':
   config=load_json(args.spec)
   if config.get('purpose','development')!='development':raise PermissionError('APPROVAL_REQUIRED: no frozen/confirmatory campaign is authorized')
   if args.output.exists():raise FileExistsError('Existing output is immutable')
   with source_from_spec(config['source'],base_dir=args.spec.parent) as source:
    result=run_application(source,config['application'],steps=config.get('steps',1024),burn_in=config.get('burn_in',0),chunk_rows=config.get('chunk_rows',128),post_uniform_fault=config.get('post_uniform_fault'),run_id=config.get('run_id'),independent_group_id=config.get('independent_group_id'))
   result.update(environment=environment(),config_sha256=hashlib.sha256(args.spec.read_bytes()).hexdigest(),protocol_status='DRAFT_NOT_FROZEN')
   result={'report':write_output(args.output,result),'status':result['status']}
  elif args.command=='analyze':
   with source_from_spec(load_json(args.source_spec),base_dir=args.source_spec.parent) as source:report=analyze(source,total_bytes=args.bytes,window_bytes=args.window_bytes,min_entropy=args.min_entropy)
   result={'report':write_output(args.output,report),'bytes_analyzed':report['bytes_analyzed']}
  elif args.command=='predict':
   model=ReliabilityPredictor.restore(load_json(args.artifact)) if args.artifact else ReliabilityPredictor(ApplicationContract('harmonic_pimc','centroid',.05))
   result=model.predict(args.severity,args.rate,risk_horizon=args.risk_horizon)
  elif args.command=='reseed-status':
   segments=[] if args.defect_start is None else [{'start':args.defect_start,'end':args.defect_end if args.defect_end is not None else args.defect_start+1}]
   result=source_exposure_horizon(current_cycle=args.cycle,interval=args.interval,source_offset=args.source_offset,reseed_bytes=args.reseed_bytes,defect_segments=segments)
  elif args.command=='infrastructure':
   from .infrastructure import inventory_resources,backend_capabilities
   result={'resources':inventory_resources(),'backends':backend_capabilities(),'execution_policy':'LOCAL_ONLY_USER_REQUEST','vm_provisioning':'DOCUMENTED_NOT_AUTHORIZED'}
   if args.output:write_output(args.output,result)
  elif args.command=='events':
   result={'report':write_output(args.output,evaluate_event_spec(args.spec))}
  elif args.command=='dashboard':
   raise NotImplementedError('Standalone dashboard available via web monitor')
  print(json.dumps(result,indent=2,allow_nan=False));return 0
 except (ValueError,TypeError,KeyError,OSError,PermissionError,RuntimeError) as exc:
  print('rnglab11: '+str(exc),file=sys.stderr);return 2

if __name__=='__main__':raise SystemExit(main())
