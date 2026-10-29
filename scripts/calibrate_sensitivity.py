"""Post-confirmation stratified bootstrap sensitivity, preserving family counts."""
from pathlib import Path
import json,numpy as np
from rnglab16.campaign import save
from rnglab17.campaign import seed
ROOT=Path(__file__).resolve().parents[2];RUN=ROOT/'confirmation-v2';rows={}
for split in ['validation','test','OOD']:
 z=np.load(RUN/f'forecast_predictions_{split}.npz');y=z['labels'];families=z['family'];d=(z['combined']-y)**2-(z['constant']-y)**2;g=np.random.default_rng(seed('stratified-sensitivity',split));groups=[np.flatnonzero(families==f) for f in np.unique(families)];boot=[]
 for _ in range(2000):
  ix=np.concatenate([g.choice(group,size=len(group),replace=True) for group in groups]);boot.append(float(d[ix].mean()))
 rows[split]={'difference':float(d.mean()),'approx_adjusted_stratified_ci':np.quantile(boot,[.00625,.99375]),'eligible_family_counts':{str(f):int(sum(families==f)) for f in np.unique(families)}}
save(ROOT/'STRATIFIED_SENSITIVITY.json',{'classification':'post-confirmation sensitivity; does not replace frozen primary analysis','results':rows})
