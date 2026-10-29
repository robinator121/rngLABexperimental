import math
import numpy as np
from scipy.stats import beta,t as student_t

def rank_threshold(scores,alpha=.05):
    s=np.asarray(scores,float)
    if s.ndim!=1 or not np.isfinite(s).all() or not 0<alpha<1:raise ValueError('finite scores and alpha in(0,1) required')
    rank=math.ceil((len(s)+1)*(1-alpha))
    if rank>len(s):raise ValueError('calibration sample too small')
    return float(np.sort(s)[rank-1])

def rate(k,n,alpha=.05,one_sided=False):
    if n<1:return {'n':n,'events':k,'rate':None,'interval':None}
    return {'n':n,'events':int(k),'rate':float(k/n),'interval':[0 if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),1 if k==n else float(beta.ppf(1-(alpha if one_sided else alpha/2),k+1,n-k))]}

def window_metrics(trace,burn,window,stride,reference,tolerances):
    x=np.asarray(trace,float);times=np.arange(burn+window,x.shape[1]+1,stride)
    means=np.stack([x[:,t-window:t].mean(1) for t in times],axis=1)
    error=abs(means-np.asarray(reference));tol=np.asarray(tolerances)
    # Persistence score: both successive errors must exceed the threshold.
    persistent=np.minimum(error[:,:-1],error[:,1:])
    return {'times':times[1:],'means':means,'error':error,'persistent_error':persistent,'score':persistent.max(1),'accuracy_any':(error>tol).any(1),'accuracy_persistent':(persistent>tol).any(1)}

def first_events(persistent,times,thresholds):
    hit=persistent>np.asarray(thresholds);result=np.full((len(hit),hit.shape[2]),-1,dtype=int)
    for j,t in enumerate(times):result=np.where((result<0)&hit[:,j],int(t),result)
    return result

def batch_intervals(series,reference,batch):
    x=np.asarray(series,float);blocks=x.shape[1]//batch
    if blocks<4:raise ValueError('at least four batches')
    means=x[:,:blocks*batch].reshape(len(x),blocks,batch).mean(2);center=means.mean(1);se=means.std(1,ddof=1)/np.sqrt(blocks);q=student_t.ppf(.975,blocks-1)
    return {'mean':center,'se':se,'covered':abs(center-reference)<=q*se,'lo':center-q*se,'hi':center+q*se,'batches':blocks,'batch':batch}

def causal_loss(healthy,faulty,reference,tolerance,seed=170900):
    a=np.asarray(healthy)-reference;b=np.asarray(faulty)-reference;delta=b*b-a*a
    rng=np.random.default_rng(seed);boots=np.array([np.mean(delta[rng.integers(len(delta),size=len(delta))]) for _ in range(2000)])
    return {'paired_excess_squared_loss':float(delta.mean()),'ci95':np.quantile(boots,[.025,.975]).tolist(),'healthy_mse':float(np.mean(a*a)),'fault_mse':float(np.mean(b*b)),'practical_squared_loss_margin':tolerance*tolerance,'practically_supported':bool(np.quantile(boots,.025)>tolerance*tolerance)}
