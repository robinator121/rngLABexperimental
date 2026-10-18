"""Simple development-fit survival baseline and explicit IPCW horizon scoring.

These functions require independent trajectories and independent censoring,
with the censoring distribution estimated on development data. They do not turn
point-horizon heuristics into probabilities or validate unknown future faults.
Discrete convention: censor time is the last confirmed nonfailure time.
"""
from dataclasses import dataclass
import math
import numpy as np


def _outcomes(durations, observed):
    d=np.asarray(durations,dtype=float);e=np.asarray(observed)
    if d.ndim!=1 or len(d)==0 or e.shape!=d.shape or not np.all(np.isfinite(d)) or np.any(d<0):
        raise ValueError("finite nonnegative one-dimensional durations required")
    if e.dtype.kind!='b':raise ValueError("observed events must be booleans")
    return d,e


@dataclass(frozen=True)
class KaplanMeierBaseline:
    durations: tuple[float,...]
    observed: tuple[bool,...]

    @classmethod
    def fit(cls,durations,observed):
        d,e=_outcomes(durations,observed)
        return cls(tuple(float(v) for v in d),tuple(bool(v) for v in e))

    def survival(self,time,before=False):
        if not math.isfinite(time) or time<0:raise ValueError("finite nonnegative time required")
        d=np.asarray(self.durations);e=np.asarray(self.observed);s=1.
        for at in np.unique(d):
            if at>time or (before and at==time):break
            risk=int(np.sum(d>=at));events=int(np.sum((d==at)&e))
            s*=1-events/risk
        if time>max(self.durations) and s>0:
            raise ValueError("horizon beyond observed training support; do not extrapolate censoring")
        return float(s)

    def failure_probability(self,horizon):
        return 1-self.survival(horizon)


class CensoringKaplanMeier(KaplanMeierBaseline):
    """Reverse KM with event-first ties, fitted from ORIGINAL event indicators.

    If T=C the failure is observed. Original failures at t leave the censoring
    risk set before its hazard at t; merely replacing event by ~event gives a
    biased censoring estimator for discrete/tied times.
    """
    def survival(self,time,before=False):
        if not math.isfinite(time) or time<0:raise ValueError("finite nonnegative time required")
        d=np.asarray(self.durations);e=np.asarray(self.observed);s=1.
        for at in np.unique(d):
            if at>time or (before and at==time):break
            risk=int(np.sum(d>=at))-int(np.sum((d==at)&e))
            censored=int(np.sum((d==at)&~e))
            if censored:
                if risk<=0:raise ValueError("invalid censoring risk set")
                s*=1-censored/risk
        if time>max(self.durations) and s>0:
            raise ValueError("horizon beyond observed training support")
        return float(s)


def ipcw_brier(durations,observed,failure_probabilities,horizon,censoring_model,trajectory_ids=None):
    """One prediction per independent trajectory; repeated landmarks forbidden.

    Failure by h contributes (1-p)^2/G(T-); survival through h contributes
    p^2/G(h-); censoring before h contributes zero to the IPCW sum. G is fitted
    with CensoringKaplanMeier on ORIGINAL event indicators in development data. The denominator includes ALL
    evaluation trajectories. Positivity and support are required, not assumed
    away by dropping unresolved outcomes.
    """
    if not isinstance(censoring_model,CensoringKaplanMeier):
        raise ValueError("fit an event-first CensoringKaplanMeier on original outcomes")
    d,e=_outcomes(durations,observed);p=np.asarray(failure_probabilities,dtype=float)
    if p.shape!=d.shape or not np.all(np.isfinite(p)) or np.any((p<0)|(p>1)):
        raise ValueError("one finite failure probability in [0,1] per trajectory required")
    if not math.isfinite(horizon) or horizon<0:raise ValueError("invalid horizon")
    if trajectory_ids is not None and (len(trajectory_ids)!=len(d) or len(set(trajectory_ids))!=len(d)):
        raise ValueError("use one landmark per independent trajectory for this score")
    g_h=censoring_model.survival(horizon,before=True)
    if g_h<=0:raise ValueError("censoring positivity/support failure")
    scores=[];unresolved=0
    for time,event,prob in zip(d,e,p):
        if event and time<=horizon:
            g=censoring_model.survival(float(time),before=True)
            if g<=0:raise ValueError("censoring positivity failure at event")
            scores.append((1-prob)**2/g)
        elif time>=horizon:
            scores.append(prob**2/g_h)
        else:
            scores.append(0.);unresolved+=1
    return {'brier_ipcw':float(np.mean(scores)),'trajectory_count':len(d),
            'censored_before_horizon':unresolved,'horizon':horizon,
            'assumptions':'independent_trajectories_and_censoring; development-fit_G; positive_support',
            'calibration_status':'SCOPED_SCORE_NOT_A_GUARANTEE','per_trajectory_scores':scores}
