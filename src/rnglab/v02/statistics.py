"""Correlated-sample uncertainty; diagnostics are not a stationarity guarantee."""
import numpy as np
from scipy.stats import t, binomtest


def correlated_summary(values):
    x = np.asarray(values,dtype=float)
    n = len(x)
    mean = float(x.mean())
    centered = x-mean
    variance = float(np.var(x,ddof=1)) if n > 1 else 0.
    if n < 16 or np.ptp(x) == 0:
        return {"mean":mean,"variance":variance,"tau":None,"ess":None,"mcse":None,
                "batch_mcse":None,"mcse_ci_low":None,"mcse_ci_high":None,"mixing_warning":True}
    size = 1 << (2*n-1).bit_length()
    ft = np.fft.rfft(centered,n=size)
    acov = np.fft.irfft(ft*ft.conj(),n=size)[:n]/n
    rho = acov/acov[0]
    pairs = rho[:n//2*2].reshape(-1,2).sum(axis=1)
    stop = np.flatnonzero(pairs <= 0)
    pairs = pairs[:stop[0]] if len(stop) else pairs
    pairs = np.minimum.accumulate(pairs) if len(pairs) else pairs
    # tau = 1 + 2 sum rho; conservatively do not report ESS>N for anticorrelation.
    tau = float(max(1.,-1+2*pairs.sum()))
    ess = n/tau
    mcse = float(np.sqrt(variance*tau/n))
    batch_size = max(1,int(np.sqrt(n)))
    batches = n//batch_size
    batch_means = x[:batches*batch_size].reshape(batches,batch_size).mean(axis=1)
    bse = float(batch_means.std(ddof=1)/np.sqrt(batches))
    half = float(t.ppf(.975,batches-1)*max(mcse,bse))
    return {"mean":mean,"variance":variance,"tau":tau,"ess":float(ess),"mcse":mcse,
            "batch_mcse":bse,"mcse_ci_low":mean-half,"mcse_ci_high":mean+half,
            "mixing_warning":bool(ess<50 or batch_size<5*tau)}


def mean_interval(values, confidence=.95):
    values = np.asarray(values,dtype=float)
    mean = float(values.mean())
    half = float(t.ppf((1+confidence)/2,len(values)-1)*values.std(ddof=1)/np.sqrt(len(values))) if len(values)>1 else None
    return mean, None if half is None else mean-half, None if half is None else mean+half


def proportion_interval(successes,total):
    if total == 0:
        return None,None,None
    ci = binomtest(successes,total).proportion_ci(confidence_level=.95,method="exact")
    return successes/total,float(ci.low),float(ci.high)


def sustained_events(errors, times, threshold, consecutive=3):
    """Times are confirmation times; never backdate to first threshold excursion."""
    failure,recovery=None,None
    bad,good=0,0
    for error,time in zip(errors,times):
        bad = bad+1 if error>threshold else 0
        good = good+1 if error<=threshold else 0
        if failure is None and bad>=consecutive:
            failure=int(time)
        elif failure is not None and recovery is None and good>=consecutive:
            recovery=int(time)
    return failure,recovery
