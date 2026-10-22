from pathlib import Path
import json,hashlib
import numpy as np
r=Path(__file__).resolve().parents[2];run=r/'confirmation-run004';cfg=json.loads((run/'config.json').read_text());C=cfg['C'];checks=[]
ref=json.loads((run/'B_L3_T2.26919_healthy.json').read_text())['exact']['energy']
for p in run.glob('C_*.json'):
 if p.name=='C_theory.json':continue
 q=json.loads(p.read_text());z=np.load(run/'traces'/f"C_{q['split']}_{q['family']}.npz");E=z['observables'][:,:,0];u=z['source_mean'];cut=C['cutoff'];X=np.array(q['features'])
 expected=np.column_stack([u[:,cut-64:cut].mean(1)-.5,u[:,cut-32:cut].mean(1)-u[:,cut-64:cut-32].mean(1),E[:,cut-64:cut].mean(1)-ref,E[:,cut-32:cut].mean(1)-E[:,cut-64:cut-32].mean(1)])
 assert np.allclose(expected,X[:,:4]);checks.append(p.stem+':prefix_features')
 # Recompute labels from retained future observations using a separate scan.
 times=list(range(64,C['steps']+1,20));event=[]
 for seq in E:
  flags=[abs(float(np.mean(seq[t-64:t]))-ref)>.15 for t in times]
  hits=[times[j] for j in range(1,len(times)) if flags[j] and flags[j-1]];event.append(hits[0] if hits else -1)
 event=np.array(event);assert np.array_equal(event,np.array(q['events']));checks.append(p.stem+':event_reconstruction')
 assert np.array_equal((event>cut)&(event<=cut+C['horizon']),q['labels']);checks.append(p.stem+':horizon')
a=json.loads((run/'A.json').read_text());ratios=[]
for q in a['rows']:
 means=np.array(q['means']);g=np.random.default_rng(188000);v=np.array([means[g.integers(len(means),size=len(means))].var(ddof=1) for _ in range(2000)])
 ratios.append({'a':q['a'],'b':q['b'],'empirical_to_theoretical_variance_ratio':q['empirical_variance']/q['finite_n_theory_var'],'empirical_variance_bootstrap_ci95':np.quantile(v,[.025,.975]).tolist(),'theoretical_variance':q['finite_n_theory_var']})
result={'status':'PASS','deterministic_trace_checks':len(checks),'checks':checks,'A_variance_controls':ratios,'note':'Post-hoc alternative reduction of saved data. No additional simulations; bootstrap uncertainty approximate.'}
(r/'FINAL_TRACE_CHECKS.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'checks':len(checks),'A_variance_controls':ratios}))
