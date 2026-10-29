"""Adversarial alternate reductions and dual-pipeline verification."""
from pathlib import Path
import json,hashlib
import numpy as np
from scipy.stats import beta,binomtest
from scipy.special import expit
from rnglab17.campaign import CONFIG,cell_specs,seed
from rnglab16.core import enumeration,TC
from rnglab16.campaign import save
from rnglab14.models.pimc_harmonic import exact_finite_bead_energy
ROOT=Path(__file__).resolve().parents[2];RUN=ROOT/'confirmation-v2'
def audit():
 checks=[]
 def verify(name,value):
  assert bool(value),name
  checks.append(name)
 a=json.loads((RUN/'calibration_analysis.json').read_text());primary=[];controlled=[];hashes=json.loads((RUN/'source_hashes.json').read_text())
 for rel,h in hashes.items():verify('source_'+rel,hashlib.sha256((ROOT/'repo/src'/rel).read_bytes()).hexdigest()==h)
 namespaces=[seed(s['id'],p) for s in cell_specs() for p in ['calibration','confirmation']];verify('unique_calibration_confirmation_namespaces',len(set(namespaces))==len(namespaces))
 for spec in cell_specs():
  name=spec['id'];i=spec['primary'];cal=json.loads((RUN/f'{name}__calibration__healthy.json').read_text());conf=json.loads((RUN/f'{name}__confirmation__healthy.json').read_text());obs=conf['observables'][i];scores=[]
  for phase,q in [('calibration',cal),('confirmation',conf)]:
   z=np.load(RUN/'traces'/f'{name}__{phase}__healthy.npz');x=z['observables'][:,:,i];s=CONFIG['Ising'] if spec['model']=='ising' else CONFIG['quantum'];ref=q['reference'][i];ts=list(range(s['burn']+s['window'],s['steps']+1,s['stride']));windows=np.stack([np.mean(x[:,t-s['window']:t],axis=1) for t in ts],axis=1);err=abs(windows-ref);score=np.max(np.minimum(err[:,:-1],err[:,1:]),axis=1);scores.append(score)
   verify(name+phase+'_score_reconstruction',np.allclose(score,np.array(q['scores'])[:,i],rtol=0,atol=1e-12))
  threshold=max(float(np.partition(scores[0],379)[379]),conf['tolerance'][i]);verify(name+'_threshold',abs(threshold-a['cells'][name]['thresholds'][i])<1e-12)
  k=int(sum(scores[1]>threshold));upper=1 if k==512 else float(beta.ppf(1-.05/12,k+1,512-k));verify(name+'_bound',abs(upper-a['cells'][name]['observables'][i]['simultaneous_one_sided_upper'])<1e-12);primary.append({'cell':name,'false_events':k,'n':512,'upper':upper,'passed':upper<=.10})
  # Identical trajectories under an unused upstream source demonstrate a causal null, not entropy quality.
  healthy=np.load(RUN/'traces'/f'{name}__confirmation__healthy.npz')['observables'];null=np.load(RUN/'traces'/f'{name}__confirmation__unused_source.npz')['observables'];verify(name+'_negative_control_identity',np.array_equal(healthy[:128],null))
  f=json.loads((RUN/f'{name}__confirmation__mild_bias.json').read_text());p=np.array(f['scores'])[:,i]>threshold;h=scores[1][:128]>threshold;b=int(sum(p&~h));c=int(sum(~p&h));pval=1 if b+c==0 else binomtest(b,b+c,.5).pvalue
  row=next(x for x in a['H2']['primary_tests'] if x['cell']==name);verify(name+'_paired_test',abs(pval-row['p'])<1e-14)
  if upper<=.10 and row['holm_p']<.05 and row['paired_detection_improvement']>=.20:controlled.append(name)
 # Independent Ising references: explicit bond sums and thermodynamic derivative of partition function.
 for L in [2,3,4]:
  n=L*L;x=((np.arange(2**n)[:,None]>>np.arange(n))&1)*2-1;E=np.zeros(len(x),float)
  for row in range(L):
   for col in range(L):
    j=row*L+col;E-=x[:,j]*x[:,row*L+(col+1)%L]+x[:,j]*x[:,((row+1)%L)*L+col]
  for T in [1.5,TC,3.5]:
   w=np.exp(-(E-E.min())/T);w/=w.sum();M=x.sum(1);ref=enumeration(L,T);expected={'energy':w@E/n,'magnetization':w@M/n,'abs_magnetization':w@abs(M)/n,'susceptibility':(w@(M*M)-(w@M)**2)/(n*T),'heat_capacity':(w@(E*E)-(w@E)**2)/(n*T*T)}
   for key,value in expected.items():verify(f'Ising_exact_{L}_{T}_{key}',abs(value-ref[key])<1e-11)
 for M in [4,8,16,32]:
  beta_=2.;modes=beta_/M+2*M/beta_*(1-np.cos(2*np.pi*np.arange(M)/M));reference=np.sum(1/modes)/M;verify('PIMC_spectral_reference_'+str(M),abs(reference-exact_finite_bead_energy(M,2.))<1e-11)
 # Incident forecast labels and all seven features are reconstructed from prefix-only arrays.
 f=json.loads((RUN/'forecast_analysis.json').read_text());reg=f.get('registry',[]);theory=np.array(json.loads((ROOT/'reference/v016_theory_grid.json').read_text())['grid']);ref=enumeration(3,TC)['energy'];threshold=a['cells']['ising_L3_T1']['thresholds'][0];allseeds=namespaces.copy()
 for row in reg:
  q=json.loads((RUN/f"forecast_{row['split']}_{row['family']}.json").read_text());z=np.load(RUN/'traces'/f"forecast_{row['split']}_{row['family']}.npz");E=z['observables'][:,:,0];u=z['source_mean'];X=np.array(q['features']);cut=1408
  prefix=u[:,cut-512:cut];s=prefix.mean(1)-.5;drift=u[:,cut-256:cut].mean(1)-u[:,cut-512:cut-256].mean(1);app=E[:,cut-512:cut].mean(1)-ref;trend=E[:,cut-256:cut].mean(1)-E[:,cut-512:cut-256].mean(1);d=1-np.sqrt(np.maximum(0,1-2*np.maximum(s,0)));current=np.interp(np.clip(d,0,.025),theory[:,0],theory[:,1])-ref;future=np.interp(np.clip(d+2*np.maximum(drift,0),0,.025),theory[:,0],theory[:,1])-ref;v=prefix-prefix.mean(1)[:,None];ac=np.sum(v[:,:-1]*v[:,1:],1)/np.maximum(1e-20,np.sum(v*v,1));expected=np.column_stack([s,drift,app,trend,abs(current),abs(future)-abs(app),ac]);verify('prefix_'+row['split']+row['family'],np.allclose(X,expected))
  times=list(range(1024,2049,128));means=np.column_stack([E[:,t-512:t].mean(1) for t in times]);flags=abs(means-ref)>threshold;event=np.full(len(E),-1)
  for j in range(1,len(times)):event=np.where((event<0)&flags[:,j-1]&flags[:,j],times[j],event)
  verify('events_'+row['split']+row['family'],np.array_equal(event,q['events']));verify('horizon_'+row['split']+row['family'],np.array_equal((event>1408)&(event<=1920),q['labels']));verify('eligible_'+row['split']+row['family'],np.array_equal((event<0)|(event>1408),q['eligible']));allseeds.append(q['base_seed'])
 verify('all_dataset_namespaces_unique',len(set(allseeds))==len(allseeds))
 if f['status']=='EXECUTED':
  feature_lists=f['features'];train=[];y=[]
  for row in reg:
   if row['split']!='development':continue
   q=json.loads((RUN/f"forecast_development_{row['family']}.json").read_text());good=np.array(q['eligible'],bool);train.extend(np.array(q['features'])[good]);y.extend(np.array(q['labels'])[good])
  train=np.array(train)
  for name,model in f['models'].items():verify('normalization_'+name,np.allclose(np.mean(train[:,feature_lists[name]],axis=0),model['center']))
  for split in ['validation','test','OOD']:
   z=np.load(RUN/f'forecast_predictions_{split}.npz');lab=z['labels'];verify('Brier_'+split,abs(np.mean((z['combined']-lab)**2)-f['splits'][split]['models']['combined']['brier'])<1e-14)
 save(ROOT/'AUDIT_RESULTS.json',{'passed':True,'checks':checks,'count':len(checks),'H1_primary':primary,'H2_effect_gates':a['H2']['supported_cells'],'H2_with_validated_false_alarm_control':controlled,'independence':'Separate alternative calculations within the validation suite, not external personnel or blinded scientific review','interpretation_correction':'Nine H2 effect tests may pass, but only those also passing H1 healthy validation qualify as controlled discrimination.'})
 print('Adversarial deterministic checks passed:',len(checks));print('H2 controlled:',controlled)
if __name__=='__main__':audit()
