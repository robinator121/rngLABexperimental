"""Scores for probability forecasts of validity and for failure-time distributions.

Binary outcomes: the Brier score and the logarithmic score are strictly proper
(Gneiting and Raftery 2007). Reliability is reported in bins with
Clopper-Pearson intervals; a forecaster is not called calibrated on the basis
of a score alone.

Failure times are right-censored at a fixed, known sweep (administrative
censoring). A predictive distribution therefore carries an explicit last cell
'not failed inside the horizon', and
- the ranked probability score is taken over the thresholds inside the
  horizon, where the event 'failed by the threshold' is observed for every run;
- a central prediction interval may extend into the censored cell, and whether
  it covers the outcome is then still decidable for every run;
- a point error is reported only for runs where both the forecast median and
  the outcome are inside the horizon, with the other three cases counted.
Nothing here imputes a failure time beyond the horizon.
"""
from __future__ import annotations

from typing import Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy import stats

FloatArray = NDArray[np.float64]


def brier(q: Sequence[float], y: Sequence[bool]) -> float:
    q, y = np.asarray(q, dtype=np.float64), np.asarray(y, dtype=np.float64)
    return float(np.mean((q - y)**2))


def log_score(q: Sequence[float], y: Sequence[bool], floor: float = 1e-6) -> tuple[float, int]:
    """Mean negative log-likelihood with probabilities clipped to [floor, 1-floor]; also the number clipped."""
    q, y = np.asarray(q, dtype=np.float64), np.asarray(y, dtype=bool)
    clipped = np.clip(q, floor, 1 - floor)
    return float(-np.mean(np.where(y, np.log(clipped), np.log1p(-clipped)))), int(np.count_nonzero(clipped != q))


def skill(score: float, baseline: float) -> float | None:
    """1 - score/baseline for negatively oriented scores; None when the baseline is zero."""
    return None if baseline == 0 else 1 - score / baseline


def clopper_pearson(k: int, n: int, level: float = 0.95) -> tuple[float, float] | None:
    if n == 0:
        return None
    a = (1 - level) / 2
    return (0.0 if k == 0 else float(stats.beta.ppf(a, k, n - k + 1)),
            1.0 if k == n else float(stats.beta.ppf(1 - a, k + 1, n - k)))


def reliability(q: Sequence[float], y: Sequence[bool], bins: int = 10) -> list[dict]:
    """Equal-width bins of the forecast probability; the last bin includes 1."""
    q, y = np.asarray(q, dtype=np.float64), np.asarray(y, dtype=bool)
    index = np.minimum((q * bins).astype(np.int64), bins - 1)
    rows = []
    for b in range(bins):
        inside = index == b
        count, hits = int(inside.sum()), int(y[inside].sum())
        rows.append({"bin": b, "low": b / bins, "high": (b + 1) / bins, "count": count,
                     "mean_forecast": float(q[inside].mean()) if count else None,
                     "observed": hits / count if count else None, "cp95": clopper_pearson(hits, count)})
    return rows


def calibration_gap(rows: Sequence[Mapping], minimum_count: int = 100) -> dict:
    """Largest |observed - mean forecast| over bins with at least `minimum_count` runs, and how many
    of those bins have a Clopper-Pearson interval that excludes their mean forecast."""
    used = [r for r in rows if r["count"] >= minimum_count]
    gaps = [abs(r["observed"] - r["mean_forecast"]) for r in used]
    outside = sum(1 for r in used if not r["cp95"][0] <= r["mean_forecast"] <= r["cp95"][1])
    return {"bins_used": len(used), "max_gap": max(gaps) if gaps else None, "bins_excluding_forecast": outside,
            "weighted_gap": (float(np.average(gaps, weights=[r["count"] for r in used])) if gaps else None)}


def spiegelhalter_z(q: Sequence[float], y: Sequence[bool]) -> float | None:
    """Spiegelhalter's calibration statistic: standard normal when every forecast is the true probability.

    Z = sum (y - q)(1 - 2q) / sqrt(sum (1 - 2q)^2 q (1 - q)). None when the variance is zero
    (all forecasts 0, 1/2 or 1).
    """
    q, y = np.asarray(q, dtype=np.float64), np.asarray(y, dtype=np.float64)
    variance = float(np.sum((1 - 2 * q)**2 * q * (1 - q)))
    return None if variance <= 0 else float(np.sum((y - q) * (1 - 2 * q)) / np.sqrt(variance))


def paired_bootstrap_difference(y: Sequence[bool], first: Sequence[float], second: Sequence[float], seed: int,
                                resamples: int = 10000, level: float = 0.95) -> dict:
    """Percentile interval of mean Brier(first) - mean Brier(second), resampling runs. Negative favours `first`."""
    y = np.asarray(y, dtype=np.float64)
    difference = (np.asarray(first) - y)**2 - (np.asarray(second) - y)**2
    rng = np.random.default_rng(seed)
    values = np.concatenate([difference[rng.integers(0, len(y), size=(min(500, resamples - s), len(y)))].mean(axis=1)
                             for s in range(0, resamples, 500)])
    return {"difference": float(difference.mean()), "level": level, "n": len(y),
            "interval": np.quantile(values, [(1 - level) / 2, 1 - (1 - level) / 2]).tolist()}


def paired_bootstrap_skill(y: Sequence[bool], model: Sequence[float], baseline: Sequence[float], seed: int,
                           resamples: int = 10000, level: float = 0.95) -> dict:
    """Percentile interval of the Brier skill of `model` over `baseline`, resampling runs."""
    y = np.asarray(y, dtype=np.float64)
    loss_m, loss_b = (np.asarray(model) - y)**2, (np.asarray(baseline) - y)**2
    rng = np.random.default_rng(seed)
    values, undefined = np.empty(resamples), 0
    for start in range(0, resamples, 500):
        index = rng.integers(0, len(y), size=(min(500, resamples - start), len(y)))
        base = loss_b[index].mean(axis=1)
        undefined += int(np.count_nonzero(base == 0))
        values[start:start + len(index)] = 1 - loss_m[index].mean(axis=1) / np.where(base == 0, np.nan, base)
    point = skill(float(loss_m.mean()), float(loss_b.mean()))
    interval = None if undefined else np.quantile(values, [(1 - level) / 2, 1 - (1 - level) / 2]).tolist()
    return {"skill": point, "interval": interval, "level": level, "resamples": resamples, "undefined_resamples": undefined,
            "n": len(y)}


# ---------------------------------------------------------------------------
# Failure-time distributions with administrative right-censoring
# ---------------------------------------------------------------------------
def outcome_column(time: NDArray[np.int64], observed: NDArray[np.bool_], ends: NDArray[np.int64]) -> NDArray[np.int64]:
    """Column of the outcome in a pmf over (checkpoints..., censored)."""
    return np.where(observed, np.searchsorted(ends, time), len(ends)).astype(np.int64)


def ranked_probability_score(pmf: FloatArray, column: NDArray[np.int64]) -> FloatArray:
    """Per-run mean over the in-horizon thresholds of (P(failed by threshold) - 1{failed by threshold})^2."""
    cdf = np.cumsum(pmf[:, :-1], axis=1)
    happened = np.arange(pmf.shape[1] - 1)[None, :] >= column[:, None]
    return np.mean((cdf - happened)**2, axis=1)


def central_interval(pmf: FloatArray, level: float = 0.8) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Columns of the lower and upper (1-level)/2 quantiles; the last column means 'beyond the horizon'."""
    cdf = np.cumsum(pmf, axis=1)
    tail = (1 - level) / 2
    lower = np.argmax(cdf >= tail - 1e-12, axis=1)
    upper = np.argmax(cdf >= 1 - tail - 1e-12, axis=1)
    return lower, upper


def interval_coverage(pmf: FloatArray, column: NDArray[np.int64], level: float = 0.8) -> dict:
    """Coverage of the central interval. Decidable for censored runs too: a censored outcome is covered
    exactly when the interval's upper end is the censored cell."""
    lower, upper = central_interval(pmf, level)
    covered = (lower <= column) & (column <= upper)
    k, n = int(covered.sum()), len(column)
    return {"level": level, "covered": k, "n": n, "coverage": k / n if n else None, "cp95": clopper_pearson(k, n),
            "upper_end_censored": int(np.count_nonzero(upper == pmf.shape[1] - 1))}


def median_error(pmf: FloatArray, column: NDArray[np.int64], ends: NDArray[np.int64]) -> dict:
    """Absolute error of the predictive median, in sweeps, where both median and outcome lie inside the horizon."""
    median = np.argmax(np.cumsum(pmf, axis=1) >= 0.5 - 1e-12, axis=1)
    last = pmf.shape[1] - 1
    both = (median < last) & (column < last)
    error = np.abs(ends[median[both]] - ends[column[both]]) if both.any() else np.zeros(0)
    return {"both_inside": int(both.sum()), "median_absolute_error_sweeps": float(np.median(error)) if len(error) else None,
            "mean_absolute_error_sweeps": float(np.mean(error)) if len(error) else None,
            "forecast_inside_outcome_censored": int(np.count_nonzero((median < last) & (column == last))),
            "forecast_censored_outcome_inside": int(np.count_nonzero((median == last) & (column < last))),
            "both_censored": int(np.count_nonzero((median == last) & (column == last)))}


def kaplan_meier(time: Sequence[float], observed: Sequence[bool]) -> tuple[FloatArray, FloatArray]:
    """Product-limit estimate of the survival function with right-censoring.

    Returns the distinct event times and S(t) just after each. Censored runs leave the risk
    set at their censoring time without counting as events; ties put events before censorings.
    """
    time, observed = np.asarray(time, dtype=np.float64), np.asarray(observed, dtype=bool)
    events = np.unique(time[observed])
    survival, current = np.empty(len(events)), 1.0
    for i, t in enumerate(events):
        at_risk = np.count_nonzero(time >= t)
        current *= 1 - np.count_nonzero((time == t) & observed) / at_risk
        survival[i] = current
    return events, survival


def km_median(time: Sequence[float], observed: Sequence[bool]) -> float | None:
    """Smallest event time at which the Kaplan-Meier survival is at most 1/2; None if never reached."""
    events, survival = kaplan_meier(time, observed)
    below = np.flatnonzero(survival <= 0.5)
    return float(events[below[0]]) if len(below) else None


def minimax_excess_brier(p_first: float, p_second: float) -> float:
    """Smallest worst-case excess Brier score of any forecast that must be the same in two worlds
    whose histories have the same law but whose outcome probabilities are p_first and p_second.

    A common forecast q has excess (q - p_i)^2 in world i; the maximum over the two worlds
    is least at the midpoint, where it equals (p_first - p_second)^2 / 4.
    """
    return (p_first - p_second)**2 / 4
