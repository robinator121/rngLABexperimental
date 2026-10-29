from pathlib import Path
import json,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OLD=ROOT.parent/'rnglab-v016-trident-20261022';run=OLD/'confirmation-run004'
checks=[]
for line in (OLD/'SHA256SUMS').read_text().splitlines():
 h,rel=line.split(maxsplit=1);p=OLD/rel
 assert hashlib.sha256(p.read_bytes()).hexdigest()==h,rel
checks.append('all'+str(len((OLD/'SHA256SUMS').read_text().splitlines()))+'_existing_sealed_entries_verified')
r=json.loads((run/'C_test_healthy.json').read_text());z=np.load(run/'traces/C_test_healthy.npz');E=z['observables'][:,:,0]
ref=json.loads((run/'B_L3_T2.26919_healthy.json').read_text())['exact']['energy'];times=np.arange(64,501,20)
rolling=np.column_stack([E[:,t-64:t].mean(1) for t in times]);flags=abs(rolling-ref)>.15;event=np.full(len(E),-1)
for j in range(1,len(times)):event=np.where((event<0)&flags[:,j-1]&flags[:,j],times[j],event)
assert np.array_equal(event,r['events']);eligible=(event<0)|(event>150);y=(event>150)&(event<=350)
assert sum(~eligible)==49 and sum(y[eligible])==57 and sum(eligible)==79
npz=np.load(run/'predictions_test.npz');diff=(npz['theory']-npz['labels'])**2-(npz['combined']-npz['labels'])**2
analysis=json.loads((run/'analysis.json').read_text());assert np.isclose(diff.mean(),analysis['forecast']['test']['primary_theory_minus_combined_brier']['difference'])
result={'checks':checks,'healthy_pre_cutoff_failures':int(sum(~eligible)),'healthy_eligible':int(sum(eligible)),'healthy_future_failures':int(sum(y[eligible])),'total_healthy_trajectories':len(E),'theory_minus_combined_brier':float(diff.mean()),'source_checksums_unchanged':True,'causal_interpretation':'v016 label records repeated finite-window deviation, not a causal defect finding'}
(ROOT/'FORENSIC_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')
(ROOT/'V017_FORENSIC_BASELINE.md').write_text('''# v0.16 forensic reproduction

All entries in the existing v0.16 seal verified without changing any old file. Reconstructed the original rolling64/two-crossing endpoint from the saved healthy Ising traces. The reported49/128 pre-cutoff events and57/79 subsequent eligible events are exact. They arise from an oversensitive operational definition applied to correlated finite windows, not evidence that healthy PRNGs broke.

The saved test Brier contrast is recomputed from raw probabilities; theory adds no improvement. OOD deterioration and zero permitted continuations are directly preserved in analysis.json. These are failed predictive-usefulness/utility tests, not a software arithmetic defect.

Strong acceptance bias can make two magnetized states absorbing: loss of ergodicity is genuine structural damage caused by that intervention. Interpreting it through a unique stationary equilibrium is a model-assumption failure. Corrected harmonic PIMC formulas remain valid, but26/32 healthy CV coverage in the last v0.16 sample is below the nominal95% target; this is an uncertainty/mixing concern, not an RNG defect. VMC constant local energy at the exact state is known mathematics and cannot certify position sampling.

The two-state perturbation identity and augmented-driver controls passed prior checks within their stated domain. v0.15 remains a proposal in the inspected archive, not an accepted experimental release. New CALIBRATE evidence uses separate sources and seeds; no old result is reclassified as independent confirmation.

'''+json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
