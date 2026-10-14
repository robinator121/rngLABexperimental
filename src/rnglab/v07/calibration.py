"""Strict sustained events and nonparametric whole-history tolerance limits."""
import numpy as np
from scipy.stats import beta,binomtest,t

METRICS=('cross','source_union','application')


def sustained_score(scores,k=2):
    a=np.asarray(scores,dtype=float)
    if k<1 or not np.isfinite(a).all():raise ValueError('Finite scores and k>=1 required')
    if len(a)<k:return 0.
    return float(np.max(np.lib.stride_tricks.sliding_window_view(a,k).min(axis=1)))


def event_time(scores,times,threshold,k=2,onset=None):
    count=0
    for score,time in zip(scores,times):
        if onset is not None and time<=onset:continue
        count=count+1 if score>threshold else 0
        if count>=k:return int(time)
    return None


def order_certificate(n,rank,target=.01,family_size=3):
    if not 1<=rank<=n or not 0<target<1 or family_size<1:raise ValueError('Invalid order certificate')
    failure=float(beta.sf(target,n+1-rank,rank))
    return {'training_histories':n,'one_based_rank':rank,'realized_probability_target_upper':target,
        'per_metric_certificate_failure_bound':failure,'simultaneous_confidence_lower':max(0.,1-family_size*failure),
        'marginal_predictive_false_alarm_bound':(n+1-rank)/(n+1),
        'assumptions':'Independent IID complete healthy-history scores under frozen settings; strict > handles ties conservatively. Certificate is conditional population tolerance coverage, not exact equality of realized FPR.'}


def cutoff(values,rank):
    values=np.asarray(values,float)
    if not 1<=rank<=len(values)or not np.isfinite(values).all():raise ValueError('Invalid calibration data')
    return float(np.sort(values)[rank-1])


def proportion(successes,total,confidence=.95):
    if not 0<=successes<=total or total<1:raise ValueError('Invalid binomial counts')
    interval=binomtest(successes,total).proportion_ci(confidence_level=confidence,method='exact')
    return {'count':successes,'n':total,'rate':successes/total,'low':float(interval.low),'high':float(interval.high),'confidence':confidence}


def tost(values,target,margin,alpha=.05):
    a=np.asarray(values,float)
    if len(a)<2 or margin<=0 or not np.isfinite(a).all():raise ValueError('Invalid equivalence data')
    se=float(a.std(ddof=1)/np.sqrt(len(a)))
    if se==0:raise ValueError('Degenerate variance requires separate investigation')
    difference=float(a.mean()-target)
    lower=float(t.sf((difference+margin)/se,len(a)-1));upper=float(t.cdf((difference-margin)/se,len(a)-1))
    return {'difference':difference,'p':max(lower,upper),'equivalent':max(lower,upper)<alpha,
        'interval':[difference-float(t.ppf(1-alpha,len(a)-1))*se,difference+float(t.ppf(1-alpha,len(a)-1))*se]}
