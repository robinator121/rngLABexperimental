"""Numerical limits of the B2 aggregated decay engine: a measurement, not a hypothesis test.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-numerics/numerics.py

Writes results/numerics.json. Everything is deterministic: references are
exact (rational or 60-digit decimal) or asymptotic with a stated error order;
the only random draws are B2 runs in the namespace rnglab-v020-nuclear-b2-numerics.
Sections:
  cdf_exact       Boost incomplete beta (used by B2) and cephes bdtr against exact sums, n <= 1e7
  cdf_small_mean  large n with small means: Boost F and S, the recurrence table and the sampler's CDF
  cdf_large_n     n = 1e9 ... PART_MAX - 1: exact centre anchor (odd n, p = 1/2) and an Edgeworth reference, error O(1/n)
  cdf_huge_n_exact exact F at n = 1e12 + 1 ... PART_MAX - 1 (p = 1/2) from the centre value and exact pmf sums
  cdf_finiteness  Boost's F and S at every point the bisection can reach, n up to PART_MAX (and the NaN region above)
  table_truncation the mass beyond the recurrence table, exact (decimal), means 1e-3 ... 256
  inversion       B2 sampler against exact rational inversion, small n, every u of a grid and +-1e-12 around steps
  table_boundaries the uniforms adjacent to every distribution-function step: floating-point CDF against exact inversion
  inversion_tails B2 sampler at extreme u (0 ... 1 - 2**-53) in both regimes, against exact decimal quantiles
  inversion_anchor exact symmetry anchors up to PART_MAX - 1 and invariance of draws to the initial bracket
  midpoint_convention  the discretisation bound with tails (exact), the exact variance of the midpoint law, and
                  the variance cost of inverting at u itself
  small_p         P(X >= 1) = 1 - (1 - p)**n for tiny p, relative error against 80-digit decimal
  paths_close     close decay constants: digits used, decimal against matrix exponential, float closed form
  paths_extreme   decay constants spanning 24 orders of magnitude, intervals 1e-15 ... 1e15
  paths_tiny      tiny intervals: relative error of the matrix exponential for tiny path probabilities
  expm_stiff      float64 matrix exponential against 150-digit decimal laws on the stiff, close-rate chain S5
  populations     B2 runs from 1e6 to 2**62 nuclei: parts, time, conservation, z against exact means
"""
from __future__ import annotations

import datetime
import json
import math
import sys
import time
from decimal import Context, Decimal, localcontext
from fractions import Fraction
from pathlib import Path

import warnings

import numpy as np
from scipy import special

from rnglab.nuclear import (
    Channel, DecayScheme, Mode, Nuclide, WordStream, bateman_chain, labelled_source, state_probabilities,
)
from rnglab.nuclear.aggregate import INT64_MAX, simulate_aggregate
from rnglab.nuclear.binomial import (
    PART_MAX, TABLE_MEAN_MAX, binomial_inverse, sampler_tails, table_tails, uses_upper_tail,
)
from rnglab.nuclear.manifest import code_identity, environment
from rnglab.nuclear.paths import bounds_hold, interval_law, tree_probabilities

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
NAMESPACE = "rnglab-v020-nuclear-b2-numerics"


def exact_cdfs(ks: list[int], n: int, p: float, digits: int = 60) -> dict[int, Decimal]:
    """P(X <= k) for every k in ks, Binomial(n, p), from one 60-digit decimal sum over mean +- 60 sd."""
    with localcontext(Context(prec=digits, Emin=-10 ** 9, Emax=10 ** 9)):
        P = Decimal(p)
        Q = 1 - P
        sd = math.sqrt(n * p * (1 - p))
        lo = max(0, int(n * p - 60 * sd) - 5)
        hi = min(n, int(n * p + 60 * sd) + 5)
        wanted = sorted(ks)
        term, total, running, below = Decimal(1), Decimal(0), Decimal(0), {}
        for x in range(lo, hi + 1):
            if x > lo:
                term = term * (n - x + 1) / x * P / Q
            total += term
        position = 0
        for x in range(lo, hi + 1):
            if x > lo:
                term = term * (n - x + 1) / x * P / Q
            else:
                term = Decimal(1)
            running += term
            while position < len(wanted) and wanted[position] == x:
                below[x] = running / total
                position += 1
        for k in wanted[position:]:
            below[k] = Decimal(1) if k >= hi else Decimal(0)
        return below


def cdf_exact_section() -> dict:
    rows = []
    cases = [(n, p) for n in (10, 100, 1000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 7) for p in (0.5, 0.1, 1e-3, 1e-6)]
    cases += [(10 ** 8, 0.5), (10 ** 8, 0.1), (10 ** 9, 0.5), (10 ** 9, 0.1)]
    for n, p in cases:
            sd = math.sqrt(n * p * (1 - p))
            ks = sorted({min(n - 1, max(0, int(n * p + z * sd))) for z in (-6, -4, -2, -1, 0, 1, 2, 4, 6)})
            worst_f = worst_s = worst_cephes = worst_used = 0.0
            references = exact_cdfs(ks, n, p)
            for k in ks:
                ref = references[k]
                f = Decimal(float(special.betaincc(k + 1, n - k, p)))
                s = Decimal(float(special.betainc(k + 1, n - k, p)))
                c = float(special.bdtr(k, n, p))
                used_f, used_s = sampler_tails(k, n, p)
                consulted = (abs(float(Decimal(used_s) - (1 - ref))) if uses_upper_tail(n, p, float(ref))
                             else abs(float(Decimal(used_f) - ref)))
                worst_used = max(worst_used, consulted)
                worst_f = max(worst_f, abs(float(f - ref)))
                worst_s = max(worst_s, abs(float(s - (1 - ref))))
                worst_cephes = max(worst_cephes, abs(c - float(ref)) if math.isfinite(c) else math.inf)
            rows.append({"n": n, "p": p, "points": len(ks), "max_abs_error_F_boost": worst_f,
                         "max_abs_error_S_boost": worst_s, "max_abs_error_cephes_bdtr": worst_cephes,
                         "max_abs_error_sampler_tails": worst_used})
    return {"reference": "60-digit decimal sum over mean +- 60 sd", "rows": rows,
            "note": "sampler errors use the tail the sampler actually consults (binomial.uses_upper_tail)"}


def exact_table(n: int, p: float, kmax: int, digits: int = 50) -> list[Decimal]:
    """F(0..kmax) from the pmf recurrence started at (1 - p)**n, in decimal (small means)."""
    with localcontext(Context(prec=digits, Emin=-10 ** 9, Emax=10 ** 9)):
        P = Decimal(p)
        Q = 1 - P
        pmf, total, out = (Decimal(n) * Q.ln()).exp(), Decimal(0), []
        for x in range(kmax + 1):
            total += pmf
            out.append(total)
            if x < n:
                pmf = pmf * (n - x) / (x + 1) * P / Q
        return out


def cdf_small_mean_section() -> dict:
    rows = []
    for n in (10 ** 6, 10 ** 7, 10 ** 8, 10 ** 9, 10 ** 10, 10 ** 12, 10 ** 15, PART_MAX):
        for mean in (1, 3, 10, 30, 100, 256, 300, 1000):
            p = mean / n
            sd = math.sqrt(mean)
            kmax = int(mean + 12 * sd + 10)
            ref = exact_table(n, p, kmax)
            ks = sorted({int(k) for k in np.linspace(0, kmax - 1, 60)})
            errors = {"F_boost": 0.0, "S_boost": 0.0, "table": 0.0, "sampler_tails": 0.0}
            table = table_tails(n, p) if mean <= TABLE_MEAN_MAX else None
            for k in ks:
                errors["F_boost"] = max(errors["F_boost"],
                                        abs(float(Decimal(float(special.betaincc(k + 1, n - k, p))) - ref[k])))
                errors["S_boost"] = max(errors["S_boost"],
                                        abs(float(Decimal(float(special.betainc(k + 1, n - k, p))) - (1 - ref[k]))))
                if table is not None:
                    f_t, s_t = (float(table[0][k]), float(table[1][k])) if k < table[0].size else (1.0, 0.0)
                    errors["table"] = max(errors["table"], abs(float(Decimal(f_t) - ref[k])),
                                          abs(float(Decimal(s_t) - (1 - ref[k]))))
                f_u, s_u = sampler_tails(k, n, p)
                consulted = (abs(float(Decimal(s_u) - (1 - ref[k]))) if uses_upper_tail(n, p, float(ref[k]))
                             else abs(float(Decimal(f_u) - ref[k])))
                errors["sampler_tails"] = max(errors["sampler_tails"], consulted)
            rows.append({"n": n, "mean": mean, "regime": "table" if mean <= TABLE_MEAN_MAX else "boost",
                         **{f"max_abs_error_{k}": v for k, v in errors.items()}})
    return {"reference": "50-digit decimal pmf recurrence from (1 - p)**n", "table_mean_max": TABLE_MEAN_MAX,
            "rows": rows}


def edgeworth_cdf(k: int, n: int, p: float) -> float:
    """Normal approximation with continuity and skewness corrections; error O(1/n).

    The standardised point is computed with an exact numerator (decimal):
    k + 1/2 - n p would lose about log10(n) digits in floating point.
    """
    with localcontext(Context(prec=60)):
        x = float((Decimal(k) + Decimal("0.5") - Decimal(n) * Decimal(p))
                  / (Decimal(n) * Decimal(p) * (1 - Decimal(p))).sqrt())
    sd = math.sqrt(n * p * (1 - p))
    gamma = (1 - 2 * p) / sd
    return float(special.ndtr(x) - math.exp(-x * x / 2) / math.sqrt(2 * math.pi) * gamma * (x * x - 1) / 6)


def cdf_large_n_section() -> dict:
    rows = []
    for n in (10 ** 9 + 1, 10 ** 12 + 1, 10 ** 15 + 1, PART_MAX - 1):
        centre = (n - 1) // 2
        f = float(special.betaincc(centre + 1, n - centre, 0.5))
        s = float(special.betainc(centre + 1, n - centre, 0.5))
        t0 = time.perf_counter()
        worst = {}
        for p in (0.5, 0.1):
            sd = math.sqrt(n * p * (1 - p))
            ks = [int(n * p + z * sd) for z in (-4, -2, -1, 0, 1, 2, 4)]
            worst[str(p)] = max(abs(float(special.betaincc(k + 1, n - k, p)) - edgeworth_cdf(k, n, p)) for k in ks)
        seconds = (time.perf_counter() - t0) / 14
        rows.append({"n": n, "F_at_centre_p_half": f, "S_at_centre_p_half": s,
                     "exact_centre_value": 0.5, "centre_error_F": abs(f - 0.5), "centre_error_S": abs(s - 0.5),
                     "max_abs_boost_minus_edgeworth": worst, "edgeworth_error_order_1_over_n": 1.0 / n,
                     "seconds_per_evaluation": seconds})
    return {"reference": "exact centre value 1/2 for odd n and p = 1/2; Edgeworth (skewness, continuity) with error O(1/n)",
            "rows": rows}


_BERNOULLI = [(1, 6), (-1, 30), (1, 42), (-1, 30), (5, 66), (-691, 2730), (7, 6)]
_PI = Decimal("3.14159265358979323846264338327950288419716939937510582097494459230781640628620899")


def _log_gamma(z: int) -> Decimal:
    """ln Gamma(z) by Stirling's series in decimal; for z >= 1e6 seven terms give over 60 digits."""
    z = Decimal(z)
    total = (z - Decimal("0.5")) * z.ln() - z + (2 * _PI).ln() / 2
    for k, (num, den) in enumerate(_BERNOULLI, start=1):
        total += Decimal(num) / Decimal(den) / (2 * k * (2 * k - 1) * z ** (2 * k - 1))
    return total


def cdf_huge_n_exact_section() -> dict:
    """Exact F at n = 1e12 + 1 ... 2**53 - 1, p = 1/2: F(m) = 1/2 for odd n, plus exact pmf sums.

    pmf(m) comes from Stirling's series (60 digits), later terms from the exact recurrence.
    """
    rows = []
    for n, z in ((10 ** 12 + 1, 1.0), (10 ** 12 + 1, 3.0), (10 ** 15 + 1, 0.1), (10 ** 15 + 1, 0.25),
                 (PART_MAX - 1, 0.05)):
        m = (n - 1) // 2
        d = int(z * math.sqrt(n / 4.0))
        k = m + d
        with localcontext(Context(prec=60, Emin=-10 ** 9, Emax=10 ** 9)):
            pmf = (_log_gamma(n + 1) - _log_gamma(m + 1) - _log_gamma(n - m + 1) - n * Decimal(2).ln()).exp()
            exact, x = Decimal("0.5"), m
            for _ in range(d):
                exact += pmf
                pmf = pmf * (n - x - 1) / (x + 2)
                x += 1
        boost_f, boost_s = sampler_tails(k, n, 0.5)
        consulted_error = (float(Decimal(boost_s) - (1 - exact)) if uses_upper_tail(n, 0.5, float(exact))
                           else float(Decimal(boost_f) - exact))
        rows.append({"n": n, "z": z, "k": k, "exact_F": float(exact),
                     "boost_F_error": float(Decimal(boost_f) - exact),
                     "boost_S_error": float(Decimal(boost_s) - (1 - exact)),
                     "sampler_consulted_error": consulted_error,
                     "edgeworth_error": float(Decimal(edgeworth_cdf(k, n, 0.5)) - exact)})
    return {"reference": "exact: centre value 1/2 for odd n, p = 1/2, plus 60-digit pmf sums (Stirling start)",
            "rows": rows}


def exact_quantiles(n: int, p: float, us) -> np.ndarray:
    """min{x : F(x) > u + 2**-54}: exact rational inversion at the midpoint of u's cell (the sampler's convention)."""
    P = Fraction(p)
    half_cell = Fraction(1, 2 ** 54)
    cdf, total = [], Fraction(0)
    for x in range(n + 1):
        total += math.comb(n, x) * P ** x * (1 - P) ** (n - x)
        cdf.append(total)
    return np.array([next(x for x in range(n + 1) if cdf[x] > Fraction(u) + half_cell) for u in us])


def inversion_section() -> dict:
    grid = np.linspace(0.0, 1.0 - 2.0 ** -53, 2001)
    draws = mismatches = near_draws = near_mismatches = 0
    for n in (1, 2, 3, 5, 8, 13, 21, 34, 55):
        for p in (1e-9, 1e-3, 0.05, 0.25, 0.5):
            got = binomial_inverse(np.full(grid.size, n), np.full(grid.size, p), grid)
            mismatches += int(np.sum(got != exact_quantiles(n, p, grid)))
            draws += grid.size
            P = Fraction(p)
            total, steps = Fraction(0), []
            for x in range(n):
                total += math.comb(n, x) * P ** x * (1 - P) ** (n - x)
                steps.append(float(total))
            near = np.array([c + d for c in steps for d in (-1e-12, 1e-12) if 0.0 <= c + d < 1.0])
            if near.size:
                got = binomial_inverse(np.full(near.size, n), np.full(near.size, p), near)
                near_mismatches += int(np.sum(got != exact_quantiles(n, p, near)))
                near_draws += near.size
    return {"reference": "exact rational inversion", "grid_draws": draws, "grid_mismatches": mismatches,
            "near_step_draws": near_draws, "near_step_mismatches": near_mismatches}


def exact_tail_tables(n: int, p: float, digits: int = 60):
    """Exact (F, S) dictionaries over the window mean +- 60 sd (decimal), for exact quantiles at extreme u."""
    with localcontext(Context(prec=digits, Emin=-10 ** 9, Emax=10 ** 9)):
        P = Decimal(p)
        Q = 1 - P
        sd = math.sqrt(n * p * (1 - p))
        a = max(0, int(n * p - 60 * sd) - 5)
        b = min(n, int(n * p + 60 * sd) + 5)
        term, terms = Decimal(1), []
        for x in range(a, b + 1):
            if x > a:
                term = term * (n - x + 1) / x * P / Q
            terms.append(term)
        total = sum(terms)
        F, S, run = {}, {}, Decimal(0)
        for i, x in enumerate(range(a, b + 1)):
            run += terms[i]
            F[x] = run / total
        run = Decimal(0)
        for i in range(len(terms) - 1, -1, -1):
            S[a + i] = run / total
            run += terms[i]
        return a, b, F, S


def inversion_tails_section() -> dict:
    us = [0.0, 2.0 ** -53, 1e-15, 1e-10, 1e-5, 0.5, 1 - 1e-5, 1 - 1e-10, 1 - 1e-15, 1 - 2.0 ** -53]
    rows = []
    for n, p in ((50, 0.3), (600, 0.4), (10 ** 6, 1e-5), (10 ** 9, 3e-8), (600, 0.45), (10 ** 4, 0.5),
                 (10 ** 6, 3e-4), (10 ** 8, 0.1)):
        a, b, F, S = exact_tail_tables(n, p)
        exact = []
        with localcontext(Context(prec=60)):
            half_cell = Decimal(2) ** -54
            for u in us:
                if u < 0.5:
                    exact.append(next(x for x in range(a, b + 1) if F[x] > Decimal(u) + half_cell))
                else:
                    exact.append(next(x for x in range(a, b + 1) if S[x] < 1 - Decimal(u) - half_cell))
        got = binomial_inverse(np.full(len(us), n), np.full(len(us), p), np.array(us)).tolist()
        rows.append({"n": n, "p": p, "regime": "table" if n * p <= TABLE_MEAN_MAX else "boost", "u": us,
                     "exact": exact, "sampler": got, "mismatches": sum(g != e for g, e in zip(got, exact))})
    return {"reference": "exact quantiles at the cell midpoints u + 2**-54, from 60-digit decimal tails",
            "rows": rows, "total_mismatches": sum(r["mismatches"] for r in rows)}


def _chernoff_tail(n: int, p: float, k: int, lower: bool) -> float:
    """Chernoff bound on P(X <= k) (lower) or P(X >= k) (upper) for Binomial(n, p), k beyond the mean."""
    x = k / n
    if lower and x <= 0.0:
        return (1.0 - p) ** n
    if not lower and x >= 1.0:
        return p ** n
    kl = x * math.log(x / p) + (1.0 - x) * math.log((1.0 - x) / (1.0 - p))
    return math.exp(-n * kl)


def midpoint_convention_section() -> dict:
    """The discretisation of the midpoint rule with exact F, evaluated exactly.

    G(x) = clip(ceil(2**53 F(x) - 1/2), 0, 2**53) / 2**53 is the distribution function of the midpoint rule
    with exact F. Bound: |Var_G - Var_F| <= B2 + B1**2 with B1 = sum_x e(x), B2 = sum_x |2 (x - mean) + 1| e(x),
    e(x) = min(2**-54, F(x), 1 - F(x)) (summation by parts, tails included). F is exact (60-digit decimal) on
    mean +- 60 sd; outside, e(x) <= F(a - 1) or 1 - F(b), bounded by Chernoff, contributes at most
    a (2 mean + 1) F(a - 1) and (n - b) (2 (n - mean) + 1) (1 - F(b)) to B2. Var_G is computed exactly beside it.
    'left_endpoint_share' is the variance that inverting at u = 0 itself (probability 2**-53 of x = 0) would add.
    The old bound 2**-54 (z_low**2 + z_high**2) left out the tails and fails for tiny means (review finding).
    """
    rows = []
    cases = ((1, 2.0 ** -55), (1, 2.0 ** -50), (10, 1e-17), (10 ** 6, 1e-20), (10 ** 6, 1e-14), (10 ** 6, 1e-9),
             (10 ** 3, 0.01), (10 ** 3, 0.5), (10 ** 6, 0.5), (10 ** 9, 1e-6), (10 ** 8, 0.5), (10 ** 9, 0.5))
    for n, p in cases:
        a, b, F, _ = exact_tail_tables(n, p)
        with localcontext(Context(prec=60, Emin=-10 ** 9, Emax=10 ** 9)):
            mean, var = Decimal(n) * Decimal(p), Decimal(n) * Decimal(p) * (1 - Decimal(p))
            half, cells = Decimal(2) ** -54, Decimal(2) ** 53
            b1 = b2 = Decimal(0)
            g_prev, g_mean, g_second = Decimal(0), Decimal(0), Decimal(0)
            for x in range(a, b + 1):
                fx = F[x] if x < n else Decimal(1)
                e = min(half, fx, 1 - fx)
                if x < n:
                    b1 += e
                    b2 += abs(2 * (x - mean) + 1) * e
                g = min(cells, max(Decimal(0), (cells * fx - Decimal("0.5")).to_integral_value(rounding="ROUND_CEILING"))) / cells
                mass = g - g_prev
                g_mean += x * mass
                g_second += Decimal(x) * x * mass
                g_prev = g
            low_tail = Decimal(_chernoff_tail(n, p, a - 1, True)) if a > 0 else Decimal(0)
            high_tail = Decimal(_chernoff_tail(n, p, b + 1, False)) if b < n else Decimal(0)
            b1 += a * low_tail + (n - b) * high_tail
            b2 += a * (2 * mean + 1) * low_tail + (n - b) * (2 * (n - mean) + 1) * high_tail
            var_g = g_second - g_mean ** 2
            bound = (b2 + b1 ** 2) / var
            actual = abs(var_g - var) / var
        low = int(binomial_inverse(n, p, 0.0))
        high = int(binomial_inverse(n, p, 1.0 - 2.0 ** -53))
        sd = math.sqrt(n * p * (1 - p))
        rows.append({"n": n, "p": p, "mean": n * p, "x_at_u_0": low, "x_at_u_max": high,
                     "z_at_u_0": (low - n * p) / sd, "z_at_u_max": (high - n * p) / sd,
                     "B1": float(b1), "B2": float(b2), "relative_variance_bound": float(bound),
                     "relative_variance_error_exact": float(actual), "bound_holds": bool(actual <= bound),
                     "old_bound_without_tails": 2.0 ** -54 * (((low - n * p) / sd) ** 2 + ((high - n * p) / sd) ** 2),
                     "left_endpoint_share": 2.0 ** -53 * (n * p) ** 2 / (n * p * (1 - p))})
    return {"rows": rows, "bound_holds_everywhere": all(r["bound_holds"] for r in rows)}


def cdf_finiteness_section() -> dict:
    """Boost's F and S at every point the bisection can reach, n up to PART_MAX, and the NaN region above it.

    Points: mean + z sd for 73 values of z in [-9, 9] (draws stay within the quantiles of 2**-54 and
    1 - 2**-54, about +-8.3 sd) and every k within 50 of the mean; only n p > TABLE_MEAN_MAX (Boost regime).
    Rows with n > PART_MAX document why PART_MAX is 2**50: the sampler refuses those n.
    """
    rows = []
    for n in (2 ** 30, 2 ** 36, 2 ** 40, 2 ** 44, 2 ** 46, 2 ** 48, 2 ** 49, PART_MAX - 1, PART_MAX,
              2 ** 52, 2 ** 53 - 1, 2 ** 53):
        for p in (0.5, 0.45, 0.4, 0.3, 0.2, 0.1, 1e-2, 1e-4, 1e-6):
            if n * p <= TABLE_MEAN_MAX:
                continue
            mean, sd = n * p, math.sqrt(n * p * (1 - p))
            ks = {int(math.floor(mean + sd * z)) for z in np.linspace(-9.0, 9.0, 73)}
            centre = int(math.floor(mean))
            ks.update(range(centre - 50, centre + 51))
            ks = np.array(sorted(k for k in ks if 0 <= k < n), dtype=np.float64)
            with warnings.catch_warnings(), np.errstate(all="ignore"):
                warnings.simplefilter("ignore")
                f = special.betaincc(ks + 1.0, n - ks, p)
                s = special.betainc(ks + 1.0, n - ks, p)
            rows.append({"n": n, "p": p, "within_domain": n <= PART_MAX, "points": int(ks.size),
                         "nonfinite_F": int(np.sum(~np.isfinite(f))), "nonfinite_S": int(np.sum(~np.isfinite(s)))})
    inside = [r for r in rows if r["within_domain"]]
    outside = [r for r in rows if not r["within_domain"]]
    return {"part_max": PART_MAX, "rows": rows,
            "points_within_domain": sum(r["points"] for r in inside),
            "nonfinite_within_domain": sum(r["nonfinite_F"] + r["nonfinite_S"] for r in inside),
            "nonfinite_beyond_domain": sum(r["nonfinite_F"] + r["nonfinite_S"] for r in outside)}


def table_truncation_section() -> dict:
    """Mass beyond the recurrence table, P(X > L) with L = min(n, ceil(mean + 40 sqrt(mean) + 60)): exact decimal."""
    rows = []
    for n in (10 ** 3, 10 ** 6, 10 ** 9, PART_MAX):
        for mean in (1e-3, 1e-2, 0.1, 0.3, 0.5, 0.8, 1.0, 1.3, 1.6, 2.0, 3.0, 5.0, 10.0, 30.0, 100.0, 256.0):
            p = mean / n
            if p > 0.5:
                continue
            length = int(min(n, math.ceil(mean + 40.0 * math.sqrt(mean) + 60.0)))
            with localcontext(Context(prec=50, Emin=-10 ** 9, Emax=10 ** 9)):
                P = Decimal(p)
                Q = 1 - P
                pmf = (Decimal(n) * Q.ln()).exp()
                for x in range(min(length, n) + 1):
                    if x < n:
                        pmf = pmf * (n - x) / (x + 1) * P / Q
                tail, x = Decimal(0), length + 1
                while x <= n and pmf > 0:
                    tail += pmf
                    if pmf < tail * Decimal("1e-45"):
                        break
                    pmf = pmf * (n - x) / (x + 1) * P / Q
                    x += 1
                log10_tail = float(tail.log10()) if tail > 0 else None
            rows.append({"n": n, "mean": mean, "table_length": length, "log10_omitted_mass": log10_tail})
    values = [r["log10_omitted_mass"] for r in rows if r["log10_omitted_mass"] is not None]
    return {"rows": rows, "max_log10_omitted_mass": max(values) if values else None}


def table_boundaries_section() -> dict:
    """The grid uniforms adjacent to every step F(x): the sampler (floating-point F) against exact inversion.

    For every x with 2**-54 < F(x) < 1 - 2**-54, u- is the largest grid point with u + 2**-54 < F(x) and u+
    the next one (exact inversion gives X <= x at u- and X > x at u+); 'far' points are 8 cells further out.
    A mismatch at u-/u+ means the evaluated F put the step on the other side of the midpoint; the review
    found one at n = 11, p = 1/2, u = 1/2.
    """
    rows = []
    cases = [(n, p) for n in range(1, 101) for p in (0.5, 0.25, 0.1, 0.03)] + [(600, 0.5), (1000, 0.45), (3000, 0.3)]
    for n, p in cases:
        with localcontext(Context(prec=80, Emin=-10 ** 9, Emax=10 ** 9)):
            P = Decimal(p)
            Q = 1 - P
            pmf, run, exact = Q ** n, Decimal(0), []
            for x in range(n + 1):
                run += pmf
                exact.append(run)
                if x < n:
                    pmf = pmf * (n - x) / (x + 1) * P / Q
            cells = Decimal(2) ** 53
            near, far = [], []
            for x in range(n):
                fx = exact[x]
                if not (Decimal(2) ** -54 < fx < 1 - Decimal(2) ** -54):
                    continue
                k = int((cells * fx - Decimal("0.5")).to_integral_value(rounding="ROUND_CEILING")) - 1
                near += [k, k + 1]
                far += [k - 8, k + 9]
        def exact_draw(k: int) -> int:
            v = (Decimal(k) + Decimal("0.5")) / Decimal(2) ** 53
            return next(x for x in range(n + 1) if exact[x] > v)
        out = {}
        for label, ks in (("adjacent", near), ("far", far)):
            ks = [k for k in ks if 0 <= k < 2 ** 53]
            if not ks:
                out[label] = (0, 0)
                continue
            us = np.array([k * 2.0 ** -53 for k in ks])
            got = binomial_inverse(np.full(us.size, n), np.full(us.size, p), us)
            with localcontext(Context(prec=80)):
                want = [exact_draw(k) for k in ks]
            out[label] = (len(ks), int(sum(int(g) != w for g, w in zip(got.tolist(), want))))
        rows.append({"n": n, "p": p, "regime": "table" if n * p <= TABLE_MEAN_MAX else "boost",
                     "adjacent_draws": out["adjacent"][0], "adjacent_mismatches": out["adjacent"][1],
                     "far_draws": out["far"][0], "far_mismatches": out["far"][1]})
    return {"rows": [r for r in rows if r["adjacent_mismatches"] or r["far_mismatches"] or r["n"] >= 100],
            "cases": len(rows),
            "adjacent_draws": sum(r["adjacent_draws"] for r in rows),
            "adjacent_mismatches": sum(r["adjacent_mismatches"] for r in rows),
            "far_draws": sum(r["far_draws"] for r in rows),
            "far_mismatches": sum(r["far_mismatches"] for r in rows)}


def inversion_anchor_section() -> dict:
    """Exact anchors at huge n and invariance of draws to the bisection's initial bracket.

    For odd n and p = 1/2, F(m) = 1/2 exactly (m = (n - 1) / 2) and pmf(m + 1) > 2**-54, so exact midpoint
    inversion gives m + 1 at u = 1/2 and m at u = 1/2 - 2**-53. Draws must not depend on the initial
    half-width (1, 4, 64): before the NaN guard, a NaN moved a draw by 3,842 at n = 2**53 - 1.
    """
    anchors = []
    for n in (10 ** 9 + 1, 10 ** 12 + 1, 10 ** 15 + 1, PART_MAX - 1):
        m = (n - 1) // 2
        got = binomial_inverse(np.full(2, n), np.full(2, 0.5), np.array([0.5, 0.5 - 2.0 ** -53])).tolist()
        anchors.append({"n": n, "draws": got, "exact": [m + 1, m], "ok": got == [m + 1, m]})
    us = np.concatenate([[0.0, 2.0 ** -53, 1e-12], np.linspace(1e-6, 1 - 1e-6, 21), [1 - 1e-12, 1 - 2.0 ** -53]])
    rows = []
    for n in (10 ** 12, 10 ** 15, PART_MAX - 1, PART_MAX):
        for p in (0.5, 0.3, 0.1, 1e-3, 1e-9):
            draws = [binomial_inverse(np.full(us.size, n), np.full(us.size, p), us, half_width=w) for w in (1, 4, 64)]
            rows.append({"n": n, "p": p, "uniforms": int(us.size),
                         "differences": int(sum(int(np.sum(d != draws[1])) for d in draws))})
    return {"anchors": anchors, "anchors_ok": all(a["ok"] for a in anchors), "bracket_rows": rows,
            "bracket_differences": sum(r["differences"] for r in rows)}


def _log_one_minus(p: Decimal) -> Decimal:
    """ln(1 - p), by its series when p is tiny (1 - p would round to 1)."""
    if p > Decimal("1e-8"):
        return (1 - p).ln()
    total, term, k = Decimal(0), p, 1
    while term > Decimal(10) ** -100 * p:
        total -= term / k
        k += 1
        term *= p
    return total


def _one_minus_exp(x: Decimal) -> Decimal:
    """1 - exp(x) for x <= 0, by its series when |x| is small."""
    if abs(x) > Decimal("1e-6"):
        return 1 - x.exp()
    total, term, k = Decimal(0), x, 1
    while abs(term) > Decimal(10) ** -100 * abs(x):
        total -= term
        k += 1
        term = term * x / k
    return total


def small_p_section() -> dict:
    rows = []
    for n in (1, 10 ** 6, 10 ** 12, PART_MAX):
        for p in (1e-300, 1e-200, 1e-100, 1e-30, 1e-16, 1e-10, 1e-5):
            with localcontext(Context(prec=80, Emin=-10 ** 9, Emax=10 ** 9)):
                exact = _one_minus_exp(Decimal(n) * _log_one_minus(Decimal(p)))
            boost = float(special.betainc(1.0, float(n), p))
            relative = abs(float((Decimal(boost) - exact) / exact)) if exact != 0 else None
            rows.append({"n": n, "p": p, "P_at_least_one_exact": float(exact), "boost": boost,
                         "relative_error": relative})
    return {"quantity": "P(X >= 1) = S(0) = betainc(1, n, p)", "reference": "80-digit decimal 1 - exp(n ln(1 - p))",
            "rows": rows}


def chain(half_lives) -> DecayScheme:
    nuclides = [Nuclide(f"C{i}", 150 + i, 400, h, (Channel(f"C{i + 1}", Mode.BETA_MINUS, 1.0),))
                for i, h in enumerate(half_lives)]
    nuclides.append(Nuclide(f"C{len(half_lives)}", 150 + len(half_lives), 400, math.inf))
    return DecayScheme("chain", nuclides)


def float_closed_form(scheme: DecayScheme, delta: float):
    rates = [0.0 if n.stable else math.log(2) / n.half_life for n in scheme.nuclides]
    try:
        return float(bateman_chain(rates, [delta])[0])
    except ValueError as error:
        return f"refused: {error}"


def paths_close_section() -> dict:
    rows = []
    for spacing in (1e-1, 1e-3, 1e-6, 1e-9, 1e-12, 1e-14):
        scheme = chain([1.0 * (1 + i * spacing) for i in range(4)])
        law = interval_law(scheme, 0, 1.0)
        tree = tree_probabilities(scheme, 0, 1.0)
        rows.append({"relative_spacing": spacing, "digits_used": law.digits_used,
                     "max_abs_decimal_minus_expm": float(np.max(np.abs(law.as_float() - tree))),
                     "deepest_path_probability": float(law.probabilities[-1]),
                     "float_closed_form_deepest": float_closed_form(scheme, 1.0),
                     "bounds_hold": bounds_hold(scheme, law)})
    return {"scheme": "four-member chain, half-lives 1 + i * spacing, then stable; interval 1", "rows": rows}


def paths_extreme_section() -> dict:
    rows = []
    for half_lives in ([1e-12, 1e12], [1e12, 1e-12], [1e-12, 1.0, 1e12], [1e6, 1e-6, 1.0], [1e-9, 1e-3, 1e3, 1e9]):
        scheme = chain(half_lives)
        for delta in (1e-15, 1e-6, 1.0, 1e6, 1e15):
            law = interval_law(scheme, 0, delta)
            with localcontext(Context(prec=300)):
                total = sum(law.probabilities, Decimal(0))
            rows.append({"half_lives": half_lives, "delta": delta, "digits_used": law.digits_used,
                         "sum_minus_one": float(total - 1),
                         "max_abs_decimal_minus_expm": float(np.max(np.abs(law.as_float()
                                                                          - tree_probabilities(scheme, 0, delta)))),
                         "bounds_hold": bounds_hold(scheme, law)})
    return {"rows": rows}


def expm_stiff_section() -> dict:
    """Float64 matrix exponential against the decimal law on the stiff, close-rate chain of the B2 study (S5).

    The decimal law is confirmed by recomputing it with 150 minimum and 100 extra digits.
    """
    protocol = json.loads((REPO / "research/v020-nuclear-b2-verify/protocol.json").read_text())
    s5 = DecayScheme.from_dict(protocol["schemes"]["S5"])
    rows = []
    for delta in sorted(set(np.diff(protocol["cases"]["S5"]["times"]).tolist())):
        for i, nuclide in enumerate(s5.nuclides):
            if nuclide.stable:
                continue
            law = interval_law(s5, i, delta)
            high = interval_law(s5, i, delta, min_digits=150, extra_digits=100)
            with localcontext(Context(prec=300)):
                precision_check = max(abs(a - b) for a, b in zip(law.probabilities, high.probabilities))
            reference = np.array([float(x) for x in high.probabilities])
            rows.append({"start": nuclide.name, "delta": delta,
                         "decimal_default_vs_150_digits": float(precision_check),
                         "max_abs_expm_minus_decimal": float(np.max(np.abs(tree_probabilities(s5, i, delta) - reference)))})
    return {"scheme": "S5 of research/v020-nuclear-b2-verify (rate ratio 1e6, two rates 1e-3 apart)", "rows": rows}


def paths_tiny_section(s2: DecayScheme) -> dict:
    rows = []
    for delta in (1e-3, 1e-6, 1e-9, 1e-12, 1e-15):
        law = interval_law(s2, 0, delta)
        tree = tree_probabilities(s2, 0, delta)
        deepest = max(range(len(law.categories)), key=lambda c: len(law.categories[c].steps))
        exact = float(law.probabilities[deepest])
        rows.append({"delta": delta, "deepest_path_steps": len(law.categories[deepest].steps),
                     "deepest_path_probability": exact,
                     "expm_relative_error_on_it": abs(tree[deepest] - exact) / exact if exact > 0 else None,
                     "bounds_hold": bounds_hold(s2, law)})
    return {"scheme": "S2 from Xa-400 (Layer A protocol)", "rows": rows}


def populations_section(s2: DecayScheme) -> dict:
    times = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]
    rows = []
    for n0 in (10 ** 6, 10 ** 9, 10 ** 12, PART_MAX, 2 ** 56, 2 ** 60, 2 ** 62):
        stream = WordStream(labelled_source(NAMESPACE, "populations", n0))
        t0 = time.perf_counter()
        run = simulate_aggregate(s2, {"Xa-400": n0}, times, stream)
        seconds = time.perf_counter() - t0
        z = []
        for t, counts in zip(times[1:], run.populations[1:]):
            p = state_probabilities(s2, "Xa-400", t)[0]
            keep = (p > 0) & (p < 1)
            z.extend(((counts[keep] - n0 * p[keep]) / np.sqrt(n0 * p[keep] * (1 - p[keep]))).tolist())
        rows.append({"n0": n0, "parts": run.plan.parts, "seconds": seconds, "words_read": run.words_read,
                     "conserved": bool(set(run.populations.sum(axis=1).tolist()) == {n0}),
                     "max_abs_z_against_exact_means": float(np.max(np.abs(z))), "z_values": len(z)})
    try:
        simulate_aggregate(s2, {"Xa-400": INT64_MAX, "Xb-396": 1}, times,
                           WordStream(labelled_source(NAMESPACE, "populations", "overflow")))
        refused = False
    except ValueError:
        refused = True
    return {"scheme": "S2 from Xa-400, grid 0, 0.5, 1, 2, 4, 8, 12", "rows": rows,
            "population_above_int64_refused": refused,
            "note": "one replicate per size: the z values are descriptive, about 54 per size"}


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=HERE / "results" / "numerics.json",
                        help="output file (development runs write elsewhere)")
    out = parser.parse_args().out
    if out.exists():
        raise SystemExit("results/numerics.json exists; results are never overwritten")
    s2 = DecayScheme.from_dict(json.loads((REPO / "research/v020-nuclear-a-verify/protocol.json").read_text())
                               ["schemes"]["S2"])
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    clock = time.perf_counter()
    result = {"generated_by": "research/v020-nuclear-b2-numerics/numerics.py", "started_utc": started,
              "code": code_identity(REPO), "environment": environment(), "namespace": NAMESPACE}
    for name, section in (("cdf_exact", cdf_exact_section), ("cdf_small_mean", cdf_small_mean_section),
                          ("cdf_large_n", cdf_large_n_section), ("cdf_huge_n_exact", cdf_huge_n_exact_section),
                          ("cdf_finiteness", cdf_finiteness_section), ("table_truncation", table_truncation_section),
                          ("inversion", inversion_section), ("table_boundaries", table_boundaries_section),
                          ("inversion_tails", inversion_tails_section), ("inversion_anchor", inversion_anchor_section),
                          ("midpoint_convention", midpoint_convention_section),
                          ("small_p", small_p_section),
                          ("paths_close", paths_close_section), ("paths_extreme", paths_extreme_section),
                          ("paths_tiny", lambda: paths_tiny_section(s2)), ("expm_stiff", expm_stiff_section),
                          ("populations", lambda: populations_section(s2))):
        t0 = time.perf_counter()
        result[name] = section()
        result[name]["section_seconds"] = time.perf_counter() - t0
        print(f"{name}: {result[name]['section_seconds']:.1f} s", flush=True)
    result["elapsed_seconds"] = time.perf_counter() - clock
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1, default=str) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
