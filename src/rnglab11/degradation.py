"""Controlled source interventions; stochastic faults and attacks stay distinct."""
from __future__ import annotations
import hashlib
import math
import random
from .sources import ByteSource,SourceExhausted,SourceFaultSource

class DegradationSource(ByteSource):
 """Byte-position interventions with bounded bursts and linear severity drift.

 Counter/substitution/periodic streams are intentionally predictable controls,
 not physical source models. Underlying bytes are consumed and retained in their
 original file, even when an adversary replaces delivered output. No forecast
 law for an arbitrary attacker follows from a programmed intervention schedule.
 """
 MODES={'bias','autocorrelation','periodicity','drift','predictability','adversarial_substitution','bit_flip','stuck_zero'}
 def __init__(self,source:ByteSource,*,mode:str,start:int=0,end:int|None=None,severity:float=1.,
              final_severity:float|None=None,seed:int=0,period:int=32,substitution:bytes=b'RNGLab controlled substitution',max_bytes:int=1_048_576):
  if mode not in self.MODES:raise ValueError('Unknown controlled degradation')
  if type(start) is not int or start<0 or (end is not None and (type(end) is not int or end<=start)):raise ValueError('Invalid half-open fault interval')
  final=severity if final_severity is None else final_severity
  if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) or not 0<=v<=1 for v in [severity,final]):raise ValueError('Severity must be finite in [0,1]')
  if final!=severity and end is None:raise ValueError('Drift needs a finite interval')
  if type(period) is not int or not 1<=period<=4096 or type(seed) is not int:raise ValueError('Bounded period and integer seed required')
  if not isinstance(substitution,bytes) or not 1<=len(substitution)<=4096:raise ValueError('Bounded substitution bytes required')
  self.source=source;self.mode=mode;self.start=start;self.end=end;self.severity=severity;self.final=final
  self.period=period;self.substitution=substitution;self.rng=random.Random(seed);self.previous=None;self.cycle=[]
  self.exposed=0;self.modified=0
  super().__init__({'adapter':'controlled_degradation_v1','mode':mode,'start_byte':start,'end_byte_exclusive':end,
   'severity':severity,'final_severity':final,'control_seed':seed,'period_bytes':period,
   'substitution_sha256':hashlib.sha256(substitution).hexdigest(),'origin_source_id':source.source_id,
   'degradation_kind':'adversarial' if mode in {'predictability','adversarial_substitution'} else 'programmed_stochastic',
   'physical_entropy_assessed':False},max_bytes=max_bytes)
 def _read(self,n):
  try:raw=self.source.read(n)
  except SourceExhausted as exc:raw=exc.partial
  result=bytearray()
  for i,value in enumerate(raw):
   index=self.offset+i;changed=value
   if index>=self.start and (self.end is None or index<self.end):
    self.exposed+=1;fraction=0 if self.end is None else (index-self.start)/max(1,self.end-self.start-1)
    p=self.severity+fraction*(self.final-self.severity)
    if self.mode in {'bias','drift'}:
     for bit in range(8):
      if self.rng.random()<p:changed&=~(1<<bit)
    else:
     if self.mode=='periodicity' and len(self.cycle)<self.period:self.cycle.append(value)
     if self.rng.random()<p:
      if self.mode=='autocorrelation' and self.previous is not None:changed=self.previous
      elif self.mode=='periodicity':changed=self.cycle[(index-self.start)%len(self.cycle)]
      elif self.mode=='predictability':changed=(index-self.start)%256
      elif self.mode=='adversarial_substitution':changed=self.substitution[(index-self.start)%len(self.substitution)]
      elif self.mode=='bit_flip':changed^=1
      elif self.mode=='stuck_zero':changed=0
   self.modified+=changed!=value;self.previous=changed;result.append(changed)
  return bytes(result)
 def provenance(self):return {**super().provenance(),'exposed_bytes':self.exposed,'modified_bytes':self.modified,'origin':self.source.provenance()}
 def close(self):self.source.close()
