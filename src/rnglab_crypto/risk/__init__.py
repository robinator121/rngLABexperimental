"""Observation windows, administratively censored survival and model risk."""
import math

import numpy as np
from rnglab14.forecast.survival import KaplanMeierBaseline
from rnglab17.endpoints import rate


def event_rate(events, alpha=.05):
    events = list(events)
    if not events or any(type(x) is not bool for x in events):
        raise ValueError("nonempty boolean experimental units required")
    return {**rate(sum(events), len(events), alpha=alpha), "method": "Clopper-Pearson exact two-sided"}


def birthday_probability(draws: int, bits: int):
    """Exact uniform finite-support birthday risk (stable log product).

    This is conditional on IID draws into exactly 2**bits equally likely
    values. It is not a predictor for arbitrary correlated or drifting sources.
    """
    if type(draws) is not int or not 0 <= draws <= 100_000 or type(bits) is not int or not 0 <= bits <= 256:
        raise ValueError("bounded draw count and bit width required")
    space = 1 << bits
    if draws > space:
        return 1.0
    return -math.expm1(sum(math.log1p(-i/space) for i in range(draws)))


def survival_summary(first_events, horizon):
    """Operation indices are zero-based; horizon is number of observed operations."""
    if type(horizon) is not int or horizon < 1 or not first_events:
        raise ValueError("positive horizon and nonempty observations required")
    if any(x is not None and (type(x) is not int or not 0 <= x < horizon) for x in first_events):
        raise ValueError("event outside observation horizon")
    durations = [horizon if x is None else x+1 for x in first_events]
    observed = [x is not None for x in first_events]
    km = KaplanMeierBaseline.fit(durations, observed)
    times = sorted(set([0, horizon]+durations))
    curve = [{"operations": t, "survival": km.survival(t)} for t in times]
    # Administrative common horizon means event indicators need no censoring weights.
    return {"censored": observed.count(False), "event_rate": event_rate(observed),
            "curve": curve, "median_operations": next((v["operations"] for v in curve if v["survival"] <= .5), None),
            "uncertainty": "pointwise exact binomial bands at the common observation horizon; no tail extrapolation"}


def paired_difference(a, b, seed=919, resamples=2000):
    """Difference a-b, bootstrap whole independent paired trajectories."""
    x = np.asarray(a, dtype=float)-np.asarray(b, dtype=float)
    if x.ndim != 1 or not len(x) or not np.isfinite(x).all():
        raise ValueError("finite paired observations required")
    rng = np.random.default_rng(seed)
    values = x[rng.integers(len(x), size=(resamples, len(x)))].mean(1)
    return {"mean_difference": float(x.mean()), "interval95": np.quantile(values, [.025, .975]).tolist(),
            "method": "descriptive paired trajectory percentile bootstrap; boundary degeneracy possible"}
