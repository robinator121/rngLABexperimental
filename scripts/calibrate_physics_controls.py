"""Additional explicitly exploratory mixing/work controls, not pooled confirmation."""
from pathlib import Path
import json,time,resource
import numpy as np
from rnglab16.core import energy,wolff,pimc,Draws,TC,enumeration
from rnglab16.campaign import interval,save
from rnglab17.campaign import seed
from rnglab17.endpoints import batch_intervals,rate
from rnglab14.models.pimc_harmonic import harmonic_energy_estimators,exact_finite_bead_energy
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'physics-controls';OUT.mkdir(exist_ok=True)

def instrumented(L,T,steps,base):
 g=np.random.default_rng(base);s=g.choice([-1,1],size=(L,L));out=[];sizes=[];bonds=[];p=1-np.exp(-2/T)
 for _ in range(steps):
  start=(int(g.integers(L)),int(g.integers(L)));cluster={start};todo=[start];sign=s[start];work=0
  while todo:
   x,y=todo.pop()
   for z in [((x+1)%L,y),((x-1)%L,y),(x,(y+1)%L),(x,(y-1)%L)]:
    work+=1
    if z not in cluster and s[z]==sign and g.random()<p:cluster.add(z);todo.append(z)
  for z in cluster:s[z]*=-1
  e=energy(s)/(L*L);m=s.mean();out.append([e,m,abs(m),m*m,e*e]);sizes.append(len(cluster));bonds.append(work)
 return np.array(out),np.array(sizes),np.array(bonds)

def run():
 start=time.monotonic();rows=[];L=4
 for ti,T in enumerate([1.5,float(TC),3.5]):
  pilot=[]
  for i in range(8):
   tr,size,bonds=instrumented(L,T,256,seed('cluster-pilot',ti,i));pilot.extend(size[64:])
  count=int(np.ceil(2048*L*L/np.mean(pilot)));burn=int(np.ceil(512*L*L/np.mean(pilot)));values=[];work=[];traces=[]
  # Freeze operation counts from pilot before independent cluster confirmations.
  save(OUT/f'cluster_protocol_T{ti}.json',{'pilot_groups':8,'pilot_mean_cluster_size':float(np.mean(pilot)),'fixed_operations':count,'burn_operations':burn,'replicas':32,'attempted_spin_budget_reference':2048*L*L,'seed_rule':'SHA256(cluster-confirmation,temperature-index,replica-index)','no_random_stopping_time':True})
  t=time.monotonic()
  for i in range(32):
   trace,size,bonds=instrumented(L,T,count,seed('cluster-confirmation',ti,i));a=trace[burn:];values.append([a[:,0].mean(),a[:,1].mean(),a[:,2].mean(),L*L*(a[:,3].mean()-a[:,1].mean()**2)/T,L*L*(a[:,4].mean()-a[:,0].mean()**2)/T**2]);work.append([int(sum(size)),int(sum(bonds))]);traces.append(trace)
  np.savez_compressed(OUT/f'cluster_T{ti}.npz',observables=traces,work=work)
  rows.append({'T':T,'L':L,'operations':count,'burn':burn,'seed_namespace':'cluster-confirmation','wall_seconds':time.monotonic()-t,'spin_flips_mean':float(np.mean(work,0)[0]),'bond_visits_mean':float(np.mean(work,0)[1]),'metrics':{k:interval(np.array(values)[:,j]) for j,k in enumerate(['energy','magnetization','abs_magnetization','susceptibility','heat_capacity'])},'reference':enumeration(L,T),'limits':'Fixed op counts approximately match flipped-spin work; not identical to attempted-spin or wall cost. No failure-threshold transfer from Metropolis.'})
 # Independent Gaussian reference and longer local chains separate estimator/formula from dynamics.
 M=16;beta=2.;R=64;N=2048;D=np.roll(np.eye(M),-1,axis=0)-np.eye(M);K=beta/M*np.eye(M)+M/beta*(D.T@D);cov=np.linalg.inv(K);chol=np.linalg.cholesky(cov);g=np.random.default_rng(seed('PIMC-IID-reference'));paths=g.normal(size=(R,N,M))@chol.T;td,cv=harmonic_energy_estimators(paths,beta);ref=exact_finite_bead_energy(M,beta)
 exact=[]
 for q,k in [(td,'TD'),(cv,'CV')]:
  means=q.mean(1);se=q.std(1,ddof=1)/np.sqrt(N);from scipy.stats import t as student_t
  covered=abs(means-ref)<=student_t.ppf(.975,N-1)*se;exact.append({'observable':k,'covered':rate(int(sum(covered)),R),'between_replica_bias':interval(means-ref)})
 np.savez_compressed(OUT/'PIMC_IID.npz',TD=td,CV=cv)
 longer=pimc(R,12288,M,seed('PIMC-longer'),'healthy');np.savez_compressed(OUT/'PIMC_longer.npz',observables=longer)
 controls=[]
 for label,trace,burn in [('short_same_prefix',longer[:,:4096],1024),('longer',longer,4096)]:
  for j,obs in [(0,'TD'),(1,'CV')]:
   for batch in [32,128,512]:
    u=batch_intervals(trace[:,burn:,j],ref,batch);controls.append({'label':label,'observable':obs,'batch':batch,'coverage':rate(int(sum(u['covered'])),R),'mean_error':interval(u['mean']-ref),'note':'same64 chains across lengths; not independent replications or confirmation. Longer burn and length change together.'})
 # Initialization sensitivity: different halves of each chain; does not prove stationarity.
 memory={'early_retained':longer[:,1024:4096,1].mean(1),'late':longer[:,8192:,1].mean(1),'paired_late_minus_early':interval(longer[:,8192:,1].mean(1)-longer[:,1024:4096,1].mean(1))}
 testA=instrumented(3,TC,100,1700000001)[0];testB=wolff(3,TC,100,1700000001)[0];assert np.array_equal(testA,testB)
 save(OUT/'results.json',{'classification':'Exploratory controls after main confirmation; no primary hypotheses retuned','cluster':rows,'PIMC_IID':exact,'PIMC_lengths':controls,'memory':memory,'instrumentation_replay_exact':True,'elapsed_seconds':time.monotonic()-start,'peak_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'PIMC_reference':ref})
 print('Exploratory physics controls completed',time.monotonic()-start,flush=True)
if __name__=='__main__':run()
