"""Exploratory estimator/mixing attribution using preserved healthy traces."""
from pathlib import Path
import numpy as np,json
from rnglab16.campaign import save,interval
ROOT=Path(__file__).resolve().parents[2];RUN=ROOT/'confirmation-v2';rows=[]
for L in [2,3,4]:
 for ti in [0,1,2]:
  name=f'ising_L{L}_T{ti}';h=json.loads((RUN/f'{name}__confirmation__healthy.json').read_text());z=np.load(RUN/'traces'/f'{name}__confirmation__healthy.npz');x=z['observables'][:,512:];T=h['spec']['T'];N=L*L;ref=h['reference'][3];within=N*(x[:,:,3].mean(1)-x[:,:,1].mean(1)**2)/T;symmetric=N*x[:,:,3].mean(1)/T
  signs=np.sign(x[:,:,1]);switches=np.sum((signs[:,:-1]*signs[:,1:])<0,axis=1)
  rows.append({'cell':name,'reference':ref,'within_chain_plugin':interval(within),'known_zero_mean_symmetry_estimator':interval(symmetric),'within_chain_error':interval(within-ref),'symmetry_error':interval(symmetric-ref),'paired_symmetry_minus_plugin':interval(symmetric-within),'chains_without_observed_direct_sign_switch':int(sum(switches==0)),'qualification':'Strict direct-sign switches omit transitions via zero. Known E[M]=0 uses exact zero-field finite-lattice spin-flip symmetry, not general-purpose repair. Correct stationary second-moment sampling still needs validation.'})
save(ROOT/'SUSCEPTIBILITY_PROBE.json',{'classification':'post-confirmation exploratory estimator attribution; original endpoints unchanged','rows':rows})
print(json.dumps(rows[-3:],indent=2))
