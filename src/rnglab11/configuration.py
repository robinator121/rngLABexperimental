"""Strict local source configuration; no acquisition or remote transmission."""
from pathlib import Path
import json
from .sources import FileSource,SyntheticSource,SystemCSPRNGSource,TimestampCSVSource
from .measurements import MeasurementCSVSource
from .conditioning import ConditionedSource,VonNeumannConditioner,HashConditioner
from .drbg import ReseedingDRBG
from .degradation import DegradationSource


def load_json(path):
 def unique(pairs):
  result={}
  for key,value in pairs:
   if key in result:raise ValueError('Duplicate JSON key: '+key)
   result[key]=value
  return result
 return json.loads(Path(path).read_text(),object_pairs_hook=unique,parse_constant=lambda x:(_ for _ in ()).throw(ValueError('Nonfinite JSON')))


def source_from_spec(spec,*,base_dir=None,depth=0):
 if not isinstance(spec,dict) or depth>8:raise ValueError('Bounded source object required')
 cfg=dict(spec);kind=cfg.pop('kind');base=Path(base_dir or '.').resolve()
 def nested(child):return source_from_spec(child,base_dir=base,depth=depth+1)
 def path(value):return (base/Path(value).expanduser()).resolve()
 if kind=='synthetic':return SyntheticSource(**cfg)
 if kind=='os_csprng':
  if cfg.get('capture_path'):cfg['capture_path']=path(cfg['capture_path'])
  return SystemCSPRNGSource(**cfg)
 if kind=='file':
  cfg['path']=path(cfg['path']);cfg['kind']=cfg.pop('source_type','binary');return FileSource(**cfg)
 if kind=='timestamps':cfg['path']=path(cfg['path']);return TimestampCSVSource(**cfg)
 if kind=='measurements':cfg['path']=path(cfg['path']);return MeasurementCSVSource(**cfg)
 if kind=='conditioned':
  child=nested(cfg.pop('source'));method=cfg.pop('method');block=cfg.pop('hash_block_bytes',64)
  try:
   conditioner=VonNeumannConditioner() if method=='von_neumann' else HashConditioner(block) if method=='sha256' else None
   if conditioner is None:raise ValueError('Unknown conditioning method')
   return ConditionedSource(child,conditioner,**cfg)
  except BaseException:child.close();raise
 if kind=='drbg':
  seed=nested(cfg.pop('seed_source'));reseed=None
  try:
   reseed=nested(cfg.pop('reseed_source')) if 'reseed_source' in cfg else None
   return ReseedingDRBG(seed,reseed_source=reseed,**cfg)
  except BaseException:
   seed.close()
   if reseed:reseed.close()
   raise
 if kind=='degraded':
  child=nested(cfg.pop('source'))
  try:
   if 'substitution_hex' in cfg:cfg['substitution']=bytes.fromhex(cfg.pop('substitution_hex'))
   return DegradationSource(child,**cfg)
  except BaseException:child.close();raise
 raise ValueError('Unknown source kind: '+str(kind))
