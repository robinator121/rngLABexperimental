"""Explicit measurement-to-bit adapter; acquisition and entropy claims are absent."""
import csv
import hashlib
import math
import os
from pathlib import Path
from .sources import ByteSource,SourceBudgetExceeded

class MeasurementCSVSource(ByteSource):
 """Finite numeric CSV samples thresholded and packed MSB-first, eight per byte.

 This deterministic lossy transformation is not a validated entropy extractor.
 Raw lines are untouched and accounted separately. A final incomplete output
 byte is discarded and reported; output exhaustion never invents padding bits.
 """
 def __init__(self,path,*,column:str,threshold:float,measurement_unit:str,max_bytes:int=1_048_576,max_input_bytes:int=16_777_216):
  if not isinstance(column,str) or not column or not isinstance(measurement_unit,str) or not measurement_unit or isinstance(threshold,bool) or not math.isfinite(threshold):raise ValueError('Document measurement column, unit and finite threshold')
  if type(max_input_bytes) is not int or max_input_bytes<1 or type(max_bytes) is not int or max_bytes<0:raise ValueError('Finite positive raw-input and nonnegative output budgets required')
  self.path=Path(path).resolve();self.file=self.path.open('rb');self.column=column;self.threshold=threshold;self.unit=measurement_unit
  self.input_bytes=0;self.input_samples=0;self.discarded_tail_bits=0;self.input_hash=hashlib.sha256();self.columns=None
  stat=os.fstat(self.file.fileno());self.identity=(stat.st_size,stat.st_mtime_ns,stat.st_ino)
  self.max_input_bytes=max_input_bytes;self.pending=0;self.pending_n=0;self.eof=False
  super().__init__({'adapter':'numeric_csv_threshold_msb','path':str(self.path),'column':column,'threshold':threshold,
   'measurement_unit':measurement_unit,'file_identity':list(self.identity),'entropy_claim':'none; threshold mapping requires source-specific validation'},max_bytes=max_bytes)
 def _line(self):
  remaining=self.max_input_bytes-self.input_bytes
  if remaining<=0:raise SourceBudgetExceeded('raw measurement input limit exhausted')
  data=self.file.readline(min(4097,remaining));self.input_bytes+=len(data);self.input_hash.update(data)
  if len(data)>4096:raise ValueError('Measurement CSV line too long')
  if data and not data.endswith(b'\n') and self.input_bytes>=self.max_input_bytes:raise SourceBudgetExceeded('input limit ended inside measurement line')
  return data
 def _read(self,n):
  stat=os.fstat(self.file.fileno())
  if (stat.st_size,stat.st_mtime_ns,stat.st_ino)!=self.identity:raise RuntimeError('Raw measurement file changed')
  out=bytearray()
  while len(out)<n and not self.eof:
   line=self._line()
   if not line:
    self.eof=True;self.discarded_tail_bits=self.pending_n;break
   fields=next(csv.reader([line.decode('utf-8')]))
   if self.columns is None:
    if len(set(fields))!=len(fields) or self.column not in fields:raise ValueError('Unique header with measurement column required')
    self.columns=fields;continue
   if len(fields)!=len(self.columns):raise ValueError('Measurement row width mismatch')
   value=float(fields[self.columns.index(self.column)])
   if not math.isfinite(value):raise ValueError('Nonfinite physical measurement')
   bit=int(value>=self.threshold);self.input_samples+=1;self.pending=(self.pending<<1)|bit;self.pending_n+=1
   if self.pending_n==8:out.append(self.pending);self.pending=self.pending_n=0
  return bytes(out)
 def provenance(self):return {**super().provenance(),'raw_input_bytes':self.input_bytes,'raw_input_prefix_sha256':self.input_hash.hexdigest(),'measurement_samples':self.input_samples,'discarded_final_bits':self.discarded_tail_bits}
 def clone(self):
  stat=self.path.stat()
  if (stat.st_size,stat.st_mtime_ns,stat.st_ino)!=self.identity:raise RuntimeError('Raw measurement file changed before replay')
  replay=type(self)(self.path,column=self.column,threshold=self.threshold,measurement_unit=self.unit,max_bytes=self.max_bytes,max_input_bytes=self.max_input_bytes)
  if replay.identity!=self.identity:
   replay.close();raise RuntimeError('Raw measurement file changed while opening replay')
  return replay
 def close(self):self.file.close()
