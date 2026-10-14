"""Independent-identity inference, with explicit paired discordances and censoring."""
import numpy as np
from scipy.stats import binomtest,t

def rate(values,alpha=.05):
    a=np.asarray(values,bool);n=len(a)
    if not n:return {'n':0,'count':0,'rate':None,'low':None,'high':None}
    k=int(a.sum());ci=binomtest(k,n).proportion_ci(confidence_level=1-alpha,method='exact')
    return {'n':n,'count':k,'rate':k/n,'low':float(ci.low),'high':float(ci.high)}

def paired_risk(fault,healthy,alpha=.05):
    f=np.asarray(fault,bool);h=np.asarray(healthy,bool)
    if f.shape!=h.shape or not len(f):raise ValueError('Nonempty paired vectors required')
    plus=rate(f&~h,alpha/2);minus=rate(h&~f,alpha/2)
    return {'n':len(f),'difference':float(f.mean()-h.mean()),'low':max(-1.,plus['low']-minus['high']),
        'high':min(1.,plus['high']-minus['low']),'discordant_increase':plus['count'],
        'discordant_decrease':minus['count'],'method':'difference of two discordant probabilities; Bonferroni exact binomial bounds'}

def mean_ci(values,alpha=.05):
    a=np.asarray(values,float)
    if len(a)<2 or not np.isfinite(a).all():raise ValueError('Two finite identities required')
    m=float(a.mean());sd=float(a.std(ddof=1));se=sd/np.sqrt(len(a));hw=float(t.ppf(1-alpha/2,len(a)-1)*se)
    return {'n':len(a),'mean':m,'sd':sd,'se':float(se),'low':m-hw,'high':m+hw,
        'p_zero':float(2*t.sf(abs(m/se),len(a)-1))if se else (1. if m==0 else 0.)}

def latency(times,onset,horizon,alpha=.05):
    a=np.asarray(times,int);pre=(a>=0)&(a<=onset);risk=~pre;post=(a>onset)&risk
    durations=np.where(post,a-onset,horizon-onset)[risk]
    # Administrative censoring only, after excluding monitor-specific prior alarms.
    cdf=[{'delay':int(v),'count_detected':int(np.sum(post&(a-onset<=v))),
          'risk_n':int(risk.sum()),'cdf':float(np.sum(post&(a-onset<=v))/risk.sum())if risk.sum()else None}
         for v in range(0,horizon-onset+1,2048)]
    return {'pre_onset_count':int(pre.sum()),'risk_n':int(risk.sum()),'post_onset':rate(post[risk],alpha),
        'right_censored':int(np.sum(risk&~post)),'restricted_mean_delay':float(durations.mean())if len(durations)else None,
        'observed_only_median_delay':float(np.median(a[post]-onset))if post.any()else None,
        'cdf':cdf,'censor_time':horizon-onset,'semantics':'First whole-history confirmation; retain streak at onset. Pre-onset alarms excluded from detector-specific latency risk set.'}
