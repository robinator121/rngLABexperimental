"""Post-hoc descriptive endpoint/monitor audit; not confirmation evidence."""
from pathlib import Path
import json,numpy as np
from rnglab16.core import enumeration,TC
from rnglab16.campaign import save
root=Path(__file__).parents[2];run=root/'confirmation-run004';ref=enumeration(3,TC)['energy'];rows=[]
for split in ['development','validation','test','OOD']:
 for p in run.glob(f'C_{split}_*.json'):
  r=json.loads(p.read_text());ok=np.array(r['eligible'],bool);y=np.array(r['labels']);x=np.array(r['features']);family=r['family']
  z=np.load(run/'traces'/f'C_{split}_{family}.npz');E=z['observables'][:,:,0];u=z['source_mean'];events=np.array(r['events']);alarm=[]
  for typ in ['source','application']:
   times=np.arange(64,501,20);values=u if typ=='source' else E
   score=np.column_stack([abs(values[:,t-64:t].mean(1)-(.5 if typ=='source' else ref)) for t in times]);hit=score>(.05 if typ=='source' else .15)
   first=np.array([times[np.flatnonzero(h)[0]] if np.any(h) else -1 for h in hit]);det=(first>=0)&(events>=0)
   alarm.append({'monitor':typ,'any_alarm_count':int(sum(first>=0)),'median_event_minus_first_alarm_conditional':float(np.median(events[det]-first[det])) if det.any() else None,'joint_events_and_alarms':int(sum(det))})
  rows.append({'split':split,'family':family,'total':len(y),'eligible':int(sum(ok)),'precutoff_failures':int(sum(~ok)),'future_failures':int(sum(y[ok])),'posthoc_source_threshold_brier':float(np.mean(((abs(x[ok,0])>.05)-y[ok])**2)) if ok.any() else None,'posthoc_application_threshold_brier':float(np.mean(((abs(x[ok,2])>.15)-y[ok])**2)) if ok.any() else None,'monitor':alarm})
save(root/'FORECAST_UNIT_TABLE.json',{'classification':'post-hoc descriptive diagnostic; thresholds .05/.15 not independent confirmation','rows':rows})
print(json.dumps([r for r in rows if r['split']=='test'],indent=2))
