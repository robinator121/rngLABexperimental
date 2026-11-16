"""RNGLab v0.20 nuclear module: statistical checks against exact laws.

Each check returns a plain dict with the observed and expected quantities, the
statistic and the p-value, so that a result file can be recounted
independently. The checks test a simulator against an analytic law; with a
healthy generator a p-value below the declared level is either a chance
rejection (rate alpha) or a defect. A non-finite input, statistic or p-value
raises ValueError: an invalid calculation must never read as a passing check.
"""
from __future__ import annotations

import math
from typing import Any, Callable, Sequence

import numpy as np
from scipy import stats

from rnglab.stats import holm_bonferroni_correction

from .decay import MODE_SPECS, DecayRun, Mode


def _finite(*values: float) -> None:
    if not all(math.isfinite(float(v)) for v in values):
        raise ValueError("non-finite statistic or p-value")


def binomial_check(observed: int, trials: int, probability: float) -> dict[str, Any]:
    """Two-sided exact binomial test (scipy minlike convention); p = 1 when there are no trials."""
    observed, trials, probability = int(observed), int(trials), float(probability)
    if not 0 <= observed <= trials or not 0.0 <= probability <= 1.0:
        raise ValueError("need 0 <= observed <= trials and a finite probability in [0, 1]")
    expected = trials * probability
    sd = math.sqrt(trials * probability * (1.0 - probability))
    if trials == 0:
        p_value = 1.0
    elif probability in (0.0, 1.0):
        p_value = 1.0 if observed == expected else 0.0
    else:
        p_value = float(stats.binomtest(observed, trials, probability).pvalue)
    _finite(p_value)
    return {"test": "binomial_exact_two_sided", "observed": observed, "trials": trials,
            "probability": probability, "expected": expected,
            "z": (observed - expected) / sd if sd > 0 else 0.0, "p_value": p_value}


def multinomial_check(observed: Sequence[int], probabilities: Sequence[float], *,
                      min_expected: float = 5.0) -> dict[str, Any]:
    """Pearson chi-square goodness of fit to Multinomial(N, p), with a fixed pooling rule.

    Categories with probability exactly 0 must be empty (otherwise p = 0 and
    ``impossible`` is set) and are left out. Categories with expected count
    below ``min_expected`` are pooled into one cell, in their original order;
    if that cell is still below ``min_expected`` it is merged into the
    remaining category with the smallest expected count (lowest index on ties).
    Degrees of freedom: cells - 1.
    """
    counts = np.asarray(observed, dtype=np.int64)
    probs = np.asarray(probabilities, dtype=np.float64)
    if counts.shape != probs.shape or counts.ndim != 1 or np.any(counts < 0) or not np.all(np.isfinite(probs)) \
            or np.any(probs < 0):
        raise ValueError("observed and probabilities must be matching finite nonnegative vectors")
    if abs(math.fsum(probs) - 1.0) > 1e-9:
        raise ValueError("probabilities must sum to 1")
    total = int(counts.sum())
    impossible = bool(np.any(counts[probs == 0.0] > 0))
    keep = np.flatnonzero(probs > 0.0)
    expected = total * probs[keep]
    small = expected < min_expected
    cells = [[int(i)] for i, s in zip(keep, small) if not s]
    pooled = [int(i) for i, s in zip(keep, small) if s]
    if pooled:
        pooled_expected = total * probs[pooled].sum()
        if pooled_expected < min_expected and cells:
            target = min(range(len(cells)), key=lambda c: (total * probs[cells[c]].sum(), cells[c][0]))
            cells[target] = sorted(cells[target] + pooled)
        else:
            cells.append(pooled)
    if len(cells) < 2:
        raise ValueError("fewer than two cells after pooling; the check is undefined")
    f_obs = np.array([counts[c].sum() for c in cells], dtype=np.float64)
    f_exp = np.array([total * probs[c].sum() for c in cells])
    f_exp *= f_obs.sum() / f_exp.sum()
    statistic, p_value = stats.chisquare(f_obs, f_exp)
    _finite(statistic, p_value)
    return {"test": "multinomial_chisquare", "observed": counts.tolist(), "probabilities": probs.tolist(),
            "total": total, "cells": cells, "cell_observed": f_obs.tolist(), "cell_expected": f_exp.tolist(),
            "statistic": float(statistic), "dof": len(cells) - 1,
            "p_value": 0.0 if impossible else float(p_value), "impossible": impossible}


def ks_check(sample: np.ndarray, cdf: Callable[[np.ndarray], np.ndarray]) -> dict[str, Any]:
    """One-sample Kolmogorov-Smirnov test against a continuous CDF."""
    sample = np.asarray(sample, dtype=np.float64)
    if sample.size == 0 or not np.all(np.isfinite(sample)):
        raise ValueError("need a nonempty, finite sample")
    result = stats.kstest(sample, cdf)
    _finite(result.statistic, result.pvalue)
    return {"test": "kolmogorov_smirnov", "n": int(sample.size), "statistic": float(result.statistic),
            "p_value": float(result.pvalue)}


def variance_ratio_check(values: Sequence[float], variance: float) -> dict[str, Any]:
    """Two-sided test of a sample variance against a known variance.

    (R - 1) s^2 / variance is referred to chi-square(R - 1); exact for normal
    data, an approximation otherwise (stated wherever it is used).
    """
    values = np.asarray(values, dtype=np.float64)
    if values.size < 2 or not np.all(np.isfinite(values)) or not (math.isfinite(variance) and variance > 0):
        raise ValueError("need at least two finite values and a finite positive variance")
    dof = values.size - 1
    statistic = float(((values - values.mean()) ** 2).sum() / variance)
    lower = stats.chi2.cdf(statistic, dof)
    _finite(statistic, lower)
    return {"test": "variance_ratio_chisquare", "n": int(values.size), "sample_variance": statistic * variance / dof,
            "variance": float(variance), "ratio": statistic / dof, "statistic": statistic, "dof": dof,
            "p_value": float(min(1.0, 2.0 * min(lower, stats.chi2.sf(statistic, dof))))}


def decay_invariants(run: DecayRun) -> dict[str, bool]:
    """Deterministic consistency checks of a Layer A run; every value must be True."""
    scheme, ev = run.scheme, run.events
    z = np.array([n.z for n in scheme.nuclides])
    a = np.array([n.a for n in scheme.nuclides])
    level = np.array([n.level for n in scheme.nuclides])
    depth = max(scheme.depth, 1)
    ids = ev["event_id"]
    out: dict[str, bool] = {}
    out["times_within_horizon"] = bool(np.all((ev["time"] >= 0.0) & (ev["time"] <= run.horizon)))
    out["event_ids_follow_layout"] = bool(np.array_equal(ids, ev["nucleus"] * scheme.depth + ev["step"])
                                          and np.all(np.diff(ids) > 0)
                                          and np.all((ev["step"] >= 0) & (ev["step"] < depth)))
    later = ev["step"] > 0
    position = np.searchsorted(ids, ev["previous_event"][later])
    found = bool(np.all(position < ids.size)) and bool(np.array_equal(ids[np.minimum(position, ids.size - 1)],
                                                                       ev["previous_event"][later]))
    out["lineage_links"] = found and bool(
        np.array_equal(ev["daughter"][position], ev["parent"][later])
        and np.all(ev["time"][position] <= ev["time"][later])
        and np.array_equal(ev["nucleus"][position], ev["nucleus"][later]))
    first = ~later
    out["first_events_leave_initial_state"] = bool(np.array_equal(ev["previous_event"][first], np.full(first.sum(), -1))
                                                   and np.array_equal(ev["parent"][first],
                                                                      run.initial_state[ev["nucleus"][first]]))
    out["channel_defines_daughter_and_mode"] = bool(
        np.array_equal(ev["daughter"], scheme.daughter_index[ev["parent"], ev["channel"]])
        and np.array_equal(ev["mode"], scheme.mode_code[ev["parent"], ev["channel"]]))
    delta = {spec.code: (spec.delta_z, spec.delta_a) for spec in MODE_SPECS.values()}
    want_z = np.array([delta[int(m)][0] for m in ev["mode"]], dtype=np.int64)
    want_a = np.array([delta[int(m)][1] for m in ev["mode"]], dtype=np.int64)
    out["mode_bookkeeping"] = bool(np.array_equal(z[ev["daughter"]] - z[ev["parent"]], want_z)
                                   and np.array_equal(a[ev["daughter"]] - a[ev["parent"]], want_a))
    isomeric = ev["mode"] == MODE_SPECS[Mode.ISOMERIC_TRANSITION].code
    out["isomeric_transitions_lower_level"] = bool(np.all(level[ev["daughter"][isomeric]] < level[ev["parent"][isomeric]]))
    out["neutrons_only_from_beta_delayed_neutron"] = bool(
        np.array_equal(ev["neutron"], ev["mode"] == MODE_SPECS[Mode.BETA_DELAYED_NEUTRON].code))
    photon_yield = scheme.photon_yield[ev["parent"], ev["channel"]]
    out["photon_yield_extremes"] = bool(not np.any(ev["photon"][photon_yield == 0.0])
                                        and np.all(ev["photon"][photon_yield == 1.0]))
    last = np.flatnonzero(np.r_[ev["nucleus"][1:] != ev["nucleus"][:-1], True]) if ids.size else np.empty(0, int)
    final = run.initial_state.copy()
    final[ev["nucleus"][last]] = ev["daughter"][last]
    out["final_states_match_last_events"] = bool(np.array_equal(final, run.final_state))
    out["populations_at_horizon_match_final_states"] = bool(np.array_equal(
        run.populations([run.horizon])[0], np.bincount(run.final_state, minlength=len(scheme))))
    return out


def holm(p_values: Sequence[float], alpha: float) -> dict[str, Any]:
    """Holm step-down family-wise control (rnglab.stats.holm_bonferroni_correction)."""
    p_values = [float(p) for p in p_values]
    if not all(math.isfinite(p) and 0.0 <= p <= 1.0 for p in p_values):
        raise ValueError("every p-value must be finite and in [0, 1]")
    if isinstance(alpha, bool) or not (isinstance(alpha, (int, float)) and 0.0 < alpha < 1.0):
        raise ValueError("alpha must be a number in (0, 1)")
    result = holm_bonferroni_correction(p_values, alpha=alpha)
    return {"method": "holm", "alpha": alpha, "tests": result.total_tests, "rejected": result.rejected_count,
            "adjusted_p_values": result.adjusted_p_values, "reject": result.reject_decisions}
