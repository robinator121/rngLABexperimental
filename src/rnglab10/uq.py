"""Fixed-sample scalar MCSE methods; none certifies the physical target."""
import math
import numpy as np
from scipy.stats import norm,t,beta

METHODS=('naive','batch_means','overlapping_batch_means','bartlett_spectral')
LEVELS=np.array([.68,.90,.95,.99])

def estimates(x):
    x=np.asarray(x,dtype=float)
    if x.ndim!=1 or len(x)<32 or not np.isfinite(x).all():raise ValueError('Finite scalar trace n>=32 required')
    n=len(x);mean=float(x.mean());var=float(x.var(ddof=1));b=int(np.sqrt(n));a=n//b
    if var<=0:raise ValueError('Degenerate trace; a zero MCSE is not a reliability certificate')
    batch=x[:a*b].reshape(a,b).mean(axis=1)
    bm=float(b*batch.var(ddof=1))
    pref=np.r_[0.,np.cumsum(x)]
    overlapping=(pref[b:]-pref[:-b])/b
    obm=float(n*b/((n-b+1)*(n-b))*np.sum((overlapping-mean)**2))
    y=x-mean;size=1<<(2*n-1).bit_length();ft=np.fft.rfft(y,n=size)
    acov=np.fft.irfft(ft*ft.conj(),n=size)[:b]/n
    spec=float(acov[0]+2*np.dot(1-np.arange(1,b)/b,acov[1:]))
    lrv=np.array([var,bm,obm,max(spec,0.)])
    return {'mean':mean,'variance':var,'n':n,'batch_size':b,'batch_count':a,
        'mcse':np.sqrt(lrv/n),'lrv':lrv,'ess':np.divide(n*var,lrv,out=np.full(4,np.nan),where=lrv>0)}

def coverage(mean,se,target=0.,reference_error=0.):
    radius=np.asarray(se)[:,None]*norm.ppf((1+LEVELS)/2)
    distance=abs(mean-target)
    return {'definite':distance+reference_error<=radius,
            'possible':distance<=radius+reference_error}

def exact_rate(k,n,alpha=.05):
    if n<=0:raise ValueError('No independent identities')
    return [k/n,0. if k==0 else float(beta.ppf(alpha/2,k,n-k+1)),
            1. if k==n else float(beta.ppf(1-alpha/2,k+1,n-k))]

def mean_ci(x,alpha=.05):
    x=np.asarray(x,dtype=float);mean=float(x.mean())
    h=float(t.ppf(1-alpha/2,len(x)-1)*x.std(ddof=1)/math.sqrt(len(x)))
    return [mean,mean-h,mean+h]

def normal_target_coverage(bias,var_mean,level=.95):
    se=np.sqrt(var_mean);z=norm.ppf((1+level)/2)
    return float(norm.cdf(z-bias/se)-norm.cdf(-z-bias/se))

def ar1_mean_variance(n,phi,stationary_variance=1.):
    lag=np.arange(1,n,dtype=float)
    return float(stationary_variance/n*(1+2*np.sum((1-lag/n)*phi**lag)))
