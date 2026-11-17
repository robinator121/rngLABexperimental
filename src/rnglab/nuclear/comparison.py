"""RNGLab v0.20 nuclear module, milestone B2: moment tests and two-sample comparisons.

Like ``checks``, every function returns a plain dict with its inputs, statistic
and p-value, and raises ValueError on a non-finite input or result.
"""
from __future__ import annotations

import math
from typing import Any, Sequence

import numpy as np
from scipy import stats


def _finite(*values: float) -> None:
    if not all(math.isfinite(float(v)) for v in values):
        raise ValueError("non-finite input, statistic or p-value")


def _counts(values: Any, what: str) -> list[int]:
    """A vector of nonnegative integer counts as Python integers: nothing is truncated, no sum can overflow."""
    items = values.tolist() if isinstance(values, np.ndarray) else list(values)
    out = []
    for value in items:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 0:
            raise ValueError(f"{what} must be nonnegative integers, got {value!r}")
        out.append(int(value))
    return out


def _satterthwaite(variance: float, variance_of_sample_variance: float) -> float:
    """nu = 2 variance**2 / Var(s^2), refused unless positive and finite (an input underflowed or overflowed)."""
    try:
        nu = 2.0 * float(variance) ** 2 / float(variance_of_sample_variance)
    except (OverflowError, ZeroDivisionError):
        nu = math.inf
    if not (math.isfinite(nu) and nu > 0):
        raise ValueError(f"degrees of freedom {nu!r} are not positive and finite (an input underflowed or overflowed)")
    return nu


def _two_sided(lower: float, upper: float) -> float:
    """2 min(lower, upper), capped at 1, after checking that both tails are finite (a NaN is never read as 1)."""
    _finite(lower, upper)
    return float(min(1.0, 2.0 * min(lower, upper)))


def z_check(estimate: float, expected: float, variance: float, label: str = "z") -> dict[str, Any]:
    """Two-sided normal test of an estimate against its expected value and exact variance."""
    _finite(estimate, expected, variance)
    if not variance > 0:
        raise ValueError("variance must be positive")
    z = (float(estimate) - float(expected)) / math.sqrt(variance)
    p = float(2.0 * stats.norm.sf(abs(z)))
    _finite(z, p)
    return {"test": label, "estimate": float(estimate), "expected": float(expected), "variance": float(variance),
            "z": z, "p_value": p}


def scaled_chisquare_variance_check(sample_variance: float, variance: float, variance_of_sample_variance: float,
                                    label: str = "scaled_chisquare_variance") -> dict[str, Any]:
    """Two-sided test of a sample variance against its exact variance, keeping the right skew of s^2.

    s^2 is approximated by (variance / nu) chi-square(nu) with nu = 2 variance^2 / Var(s^2), which
    matches the exact mean and variance of s^2 (Satterthwaite). A normal approximation of s^2 would
    understate its right tail, by about a factor of 2 at z = 4.
    """
    _finite(sample_variance, variance, variance_of_sample_variance)
    if not (variance > 0 and variance_of_sample_variance > 0):
        raise ValueError("variances must be positive")
    nu = _satterthwaite(variance, variance_of_sample_variance)
    statistic = nu * sample_variance / variance
    _finite(statistic)
    p = _two_sided(stats.chi2.cdf(statistic, nu), stats.chi2.sf(statistic, nu))
    return {"test": label, "sample_variance": float(sample_variance), "variance": float(variance),
            "variance_of_sample_variance": float(variance_of_sample_variance), "dof": nu, "statistic": statistic,
            "ratio": float(sample_variance / variance), "p_value": p}


def variance_ratio_f_check(sample_variance_a: float, sample_variance_b: float, variance: float,
                           variance_of_sample_variance: float, label: str = "variance_ratio_f") -> dict[str, Any]:
    """Two-sided test that two independent sample variances estimate the same variance.

    Each s^2 is approximated by (variance / nu) chi-square(nu) with nu = 2 variance^2 / Var(s^2)
    (as in scaled_chisquare_variance_check), so their ratio is F(nu, nu) and the common scale
    cancels: only nu comes from the analytic reference. A z-test of the difference s_a^2 - s_b^2
    with the exact null variance was tried first and rejected: its null rejection rate at the Holm
    level was 1.66 times nominal (research/v020-nuclear-b2-verify/calibration.json).
    """
    _finite(sample_variance_a, sample_variance_b, variance, variance_of_sample_variance)
    if not (sample_variance_a > 0 and sample_variance_b > 0 and variance > 0 and variance_of_sample_variance > 0):
        raise ValueError("variances must be positive")
    nu = _satterthwaite(variance, variance_of_sample_variance)
    ratio = float(sample_variance_a) / float(sample_variance_b)
    _finite(ratio)
    p = _two_sided(stats.f.cdf(ratio, nu, nu), stats.f.sf(ratio, nu, nu))
    return {"test": label, "sample_a": float(sample_variance_a), "sample_b": float(sample_variance_b),
            "dof": nu, "statistic": ratio, "p_value": p}


def null_difference_check(a: float, b: float, variance_each: float, label: str = "null_difference") -> dict[str, Any]:
    """z-test of a - b when both are independent draws of one statistic with known null variance.

    Under the null the difference of two independent, identically distributed statistics is symmetric,
    so a normal reference is accurate even when each statistic is skewed.
    """
    result = z_check(float(a) - float(b), 0.0, 2.0 * float(variance_each), label)
    result.update({"sample_a": float(a), "sample_b": float(b)})
    return result


def sample_covariance(x: Sequence[float], y: Sequence[float]) -> float:
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if x.size != y.size or x.size < 2:
        raise ValueError("need two samples of equal size >= 2")
    return float(((x - x.mean()) * (y - y.mean())).sum() / (x.size - 1))


def plugin_covariance_variance(x: Sequence[float], y: Sequence[float]) -> float:
    """Reference-free estimate of Var(s_xy) from the sample's own second and fourth moments."""
    x, y = np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)
    r = x.size
    dx, dy = x - x.mean(), y - y.mean()
    m22 = float(np.mean(dx * dx * dy * dy))
    sxx, syy, sxy = float(dx @ dx / (r - 1)), float(dy @ dy / (r - 1)), float(dx @ dy / (r - 1))
    return (m22 - sxy ** 2) / r + (sxx * syy + sxy ** 2) / (r * (r - 1))


def covariance_difference_check(xa, ya, xb, yb, label: str = "covariance_difference") -> dict[str, Any]:
    """Two-sample test of equal covariance (or variance when x is y), with plug-in variances."""
    sa, sb = sample_covariance(xa, ya), sample_covariance(xb, yb)
    va, vb = plugin_covariance_variance(xa, ya), plugin_covariance_variance(xb, yb)
    result = z_check(sa - sb, 0.0, va + vb, label)
    result.update({"sample_a": sa, "sample_b": sb})
    return result


def homogeneity_check(a: Sequence[int], b: Sequence[int], *, min_expected: float = 5.0,
                      cells: Sequence[Sequence[int]] | None = None) -> dict[str, Any]:
    """Chi-square test that two count vectors come from one multinomial law (2 x K table).

    With ``cells`` (the form the B2 verification uses), the columns are pooled
    into cells declared before any data is seen, so the degrees of freedom
    (cells - 1) and the applicability never depend on the data; a cell that is
    empty in both rows contributes 0 and keeps its degree of freedom
    (conservative), and a count outside the cells is an error.

    Without ``cells``, columns with an expected count below ``min_expected`` in
    either row are pooled into one column, in their original order; a pooled
    column still below the limit joins the remaining column with the smallest
    total. Columns that are empty in both rows are dropped. This rule depends
    on the observed margins.
    """
    a, b = _counts(a, "counts a"), _counts(b, "counts b")
    if len(a) != len(b) or not a:
        raise ValueError("need two nonnegative count vectors of equal length")
    table = np.array([a, b], dtype=np.float64)
    if cells is not None:
        return _homogeneity_fixed(table, cells, a, b)
    keep = np.flatnonzero(table.sum(axis=0) > 0)
    rows = table.sum(axis=1)
    if np.any(rows == 0):
        raise ValueError("each row needs at least one count")
    expected_min = lambda cols: float(table[:, cols].sum() * rows.min() / rows.sum())  # noqa: E731
    small = [int(c) for c in keep if expected_min([c]) < min_expected]
    cells = [[int(c)] for c in keep if int(c) not in small]
    if small:
        if expected_min(small) < min_expected and cells:
            target = min(range(len(cells)), key=lambda i: (table[:, cells[i]].sum(), cells[i][0]))
            cells[target] = sorted(cells[target] + small)
        else:
            cells.append(small)
    if len(cells) < 2:
        raise ValueError("fewer than two columns after pooling; the check is undefined")
    pooled = np.array([[table[r, c].sum() for c in cells] for r in range(2)])
    statistic, p, dof, _ = stats.chi2_contingency(pooled, correction=False)
    _finite(statistic, p)
    return {"test": "homogeneity_chisquare", "a": [int(v) for v in a], "b": [int(v) for v in b], "cells": cells,
            "statistic": float(statistic), "dof": int(dof), "p_value": float(p)}


def _homogeneity_fixed(table: np.ndarray, cells: Sequence[Sequence[int]], a, b) -> dict[str, Any]:
    columns = table.shape[1]
    members = [int(i) for cell in cells for i in cell]
    if len(cells) < 2 or len(set(members)) != len(members) or any(not 0 <= i < columns for i in members):
        raise ValueError("need two or more disjoint cells of valid column indices")
    outside = np.ones(columns, dtype=bool)
    outside[members] = False
    if np.any(table[:, outside] > 0):
        raise ValueError("a count lies outside the declared cells")
    pooled = np.array([[table[r, list(cell)].sum() for cell in cells] for r in range(2)])
    rows = pooled.sum(axis=1)
    if np.any(rows == 0):
        raise ValueError("each row needs at least one count")
    totals = pooled.sum(axis=0)
    expected = np.outer(rows, totals) / rows.sum()
    used = totals > 0
    statistic = float((((pooled - expected) ** 2)[:, used] / expected[:, used]).sum())
    dof = len(cells) - 1
    p = float(stats.chi2.sf(statistic, dof))
    _finite(statistic, p)
    return {"test": "homogeneity_chisquare_fixed_cells", "a": [int(v) for v in a], "b": [int(v) for v in b],
            "cells": [list(map(int, cell)) for cell in cells], "empty_cells": int((~used).sum()),
            "statistic": statistic, "dof": dof, "p_value": p}


def consecutive_cells(expected: Sequence[float], min_expected: float) -> list[list[int]]:
    """Cells of consecutive indices (grid intervals), decided before any data is seen.

    Walks from the left, closing a cell once its expected count reaches
    ``min_expected``; a remainder below the limit joins the last closed cell.
    Returns a single cell when no cell reaches the limit (then no per-interval
    test is defined).
    """
    expected = np.asarray(expected, dtype=np.float64)
    if expected.ndim != 1 or expected.size == 0 or not np.all(np.isfinite(expected)) or np.any(expected < 0):
        raise ValueError("expected counts must be a finite nonnegative nonempty vector")
    cells: list[list[int]] = []
    current: list[int] = []
    for i, value in enumerate(expected.tolist()):
        current.append(i)
        if float(expected[current].sum()) >= min_expected:
            cells.append(current)
            current = []
    if current:
        if cells:
            cells[-1] += current
        else:
            cells.append(current)
    return cells


def allocation_check(successes: Sequence[int], trials: Sequence[int], probability: float,
                     cells: Sequence[Sequence[int]]) -> dict[str, Any]:
    """Allocation of successes over intervals, given the trials of each interval.

    Under the null, successes of interval k are Binomial(trials_k, p), independent
    across intervals given the trials (photons given decays). The Pearson
    statistic sum_c (S_c - p T_c)^2 / (p (1 - p) T_c) over pre-declared cells of
    intervals (S_c, T_c cell sums) is referred to chi-square(number of cells):
    given the trials each term has mean 1 and is asymptotically chi-square(1).
    A cell without trials contributes 0 and keeps its degree of freedom
    (conservative). It sees how successes are spread over time, which a test of
    the totals cannot.
    """
    s, n = _counts(successes, "successes"), _counts(trials, "trials")
    if len(s) != len(n) or any(x > m for x, m in zip(s, n)):
        raise ValueError("successes and trials must be matching vectors with 0 <= successes <= trials")
    _finite(probability)
    if not 0.0 < probability < 1.0:
        raise ValueError("probability must lie strictly between 0 and 1")
    members = sorted(int(i) for cell in cells for i in cell)
    if len(cells) < 2 or members != list(range(len(s))):
        raise ValueError("cells must partition the intervals into two or more parts")
    statistic = 0.0
    cell_successes, cell_trials = [], []
    for cell in cells:
        sc, nc = sum(s[int(i)] for i in cell), sum(n[int(i)] for i in cell)  # Python integers: exact
        cell_successes.append(sc)
        cell_trials.append(nc)
        if nc > 0:
            statistic += (sc - probability * nc) ** 2 / (probability * (1.0 - probability) * nc)
    dof = len(cells)
    p = float(stats.chi2.sf(statistic, dof))
    _finite(statistic, p)
    return {"test": "binomial_allocation_chisquare", "probability": float(probability),
            "cells": [list(map(int, cell)) for cell in cells], "cell_successes": cell_successes,
            "cell_trials": cell_trials, "statistic": statistic, "dof": dof, "p_value": p}


def pooled_cells(expected: Sequence[float], min_expected: float) -> list[list[int]]:
    """Pooling of categories by expected count, decided before any data is seen.

    Same rule as checks.multinomial_check, applied to expected counts: a
    category with expected count 0 is left out; categories below
    ``min_expected`` are pooled into one cell, in their original order; a pooled
    cell still below the limit joins the remaining cell with the smallest
    expected count (lowest index on ties). The result may have fewer than two
    cells, in which case no chi-square test is defined.
    """
    expected = np.asarray(expected, dtype=np.float64)
    if expected.ndim != 1 or not np.all(np.isfinite(expected)) or np.any(expected < 0):
        raise ValueError("expected counts must be a finite nonnegative vector")
    keep = [int(i) for i in np.flatnonzero(expected > 0)]
    small = [i for i in keep if expected[i] < min_expected]
    cells = [[i] for i in keep if expected[i] >= min_expected]
    if small:
        if expected[small].sum() < min_expected and cells:
            target = min(range(len(cells)), key=lambda c: (expected[cells[c]].sum(), cells[c][0]))
            cells[target] = sorted(cells[target] + small)
        else:
            cells.append(small)
    return cells


def fixed_cells_chisquare(observed: Sequence[int], probabilities: Sequence[float],
                          cells: Sequence[Sequence[int]]) -> dict[str, Any]:
    """Pearson statistic of counts against Multinomial(n, probabilities) on pre-declared cells.

    n is the observed total. Every category with positive probability must lie
    in exactly one cell; a category outside the cells must have probability 0
    and count 0. Degrees of freedom: cells - 1. Returns the statistic and its
    degrees of freedom (no p-value: statistics of several intervals are summed).
    """
    exact = _counts(observed, "observed counts")
    counts = np.array(exact, dtype=object)
    probs = np.asarray(probabilities, dtype=np.float64)
    if counts.shape != probs.shape or counts.ndim != 1 or not np.all(np.isfinite(probs)) or np.any(probs < 0):
        raise ValueError("observed and probabilities must be matching finite nonnegative vectors")
    members = [int(i) for cell in cells for i in cell]
    if len(cells) < 2 or len(set(members)) != len(members) or any(not 0 <= i < probs.size for i in members):
        raise ValueError("need two or more disjoint cells of valid category indices")
    outside = np.ones(probs.size, dtype=bool)
    outside[members] = False
    if np.any(probs[outside] > 0) or np.any(counts[outside] > 0):
        raise ValueError("a category outside the cells has positive probability or a nonzero count")
    total = sum(exact)
    if total <= 0:
        raise ValueError("no observations: the statistic is undefined")
    cell_probability = np.array([probs[list(cell)].sum() for cell in cells])
    if abs(float(cell_probability.sum()) - 1.0) > 1e-9 or np.any(cell_probability <= 0):
        raise ValueError("cell probabilities must be positive and sum to 1")
    cell_observed = np.array([float(sum(exact[int(i)] for i in cell)) for cell in cells])
    cell_expected = total * cell_probability / cell_probability.sum()
    statistic = float(((cell_observed - cell_expected) ** 2 / cell_expected).sum())
    _finite(statistic)
    return {"test": "fixed_cells_pearson", "total": total, "cells": [list(map(int, cell)) for cell in cells],
            "cell_observed": cell_observed.tolist(), "cell_expected": cell_expected.tolist(),
            "statistic": statistic, "dof": len(cells) - 1}


def planned_path_cells(plan, expected_start_counts: Sequence[Sequence[float]], *,
                       min_expected: float = 10.0) -> dict[int, dict[int, list[list[int]]]]:
    """Pre-declared cells of the per-interval path-category tests of an AggregatePlan.

    expected_start_counts[k][i] is the expected number of nuclei in nuclide i
    at the start of interval k (from the exact law, before any data). Returns
    {source nuclide: {interval: cells}}, keeping the intervals whose expected
    category counts leave two or more cells. Because the cells are fixed in
    advance, the degrees of freedom of each interval's Pearson statistic are
    constants, and the statistics of successive intervals, each conditionally
    chi-square given the past, sum to a chi-square with the summed degrees of
    freedom (asymptotically).
    """
    out: dict[int, dict[int, list[list[int]]]] = {}
    for i in plan.sources:
        per_interval = {}
        for k, delta in enumerate(plan.deltas.tolist()):
            q = plan.interval_plans[delta][i].law.as_float()
            cells = pooled_cells(float(expected_start_counts[k][i]) * q, min_expected)
            if len(cells) >= 2:
                per_interval[k] = cells
        if per_interval:
            out[int(i)] = per_interval
    return out


def summed_chisquare(parts: Sequence[dict[str, Any]], label: str = "summed_conditional_chisquare") -> dict[str, Any]:
    """Sum of Pearson statistics with fixed degrees of freedom, referred to chi-square(summed dof)."""
    if not parts:
        raise ValueError("nothing to sum")
    statistic = float(sum(part["statistic"] for part in parts))
    dof = int(sum(part["dof"] for part in parts))
    p = float(stats.chi2.sf(statistic, dof))
    _finite(statistic, p)
    return {"test": label, "statistic": statistic, "dof": dof, "parts": len(parts), "p_value": p}
