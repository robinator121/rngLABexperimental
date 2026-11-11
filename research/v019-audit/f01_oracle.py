"""F-01: numerical analysis of the Stage D oracle probability. Writes results/f01_oracle.json.

    PYTHONPATH=src python research/v019-audit/f01_oracle.py

The frozen quantity (research/v019-d-pilot/protocol.json, oracle.probability) is

    q(delta) = P(|T| <= c),  T ~ noncentral t(15, delta),  c = scipy.stats.t.ppf(0.975, df=15),

evaluated there as scipy.stats.nct.cdf(c, 15, delta) - scipy.stats.nct.cdf(-c, 15, delta).

This script (1) maps where that call is not finite and which operation fails,
(2) computes q with a 70-digit reference that uses only Python's decimal
module, (3) evaluates candidate double-precision formulations against it.
It changes nothing in the pilot.

Reference. With W = sqrt(V/15), V ~ chi-square(15), T = (Z + delta)/W, so for any real t

    F(t; delta) = P(T <= t) = integral_0^inf Phi(t w - delta) f(w) dw,
    f(w) = 2 (15/2)^(15/2) / Gamma(15/2) * w^14 * exp(-15 w^2 / 2).

Both tails are integrated directly (no subtraction from 1), with the substitution
w = exp(u) and the trapezoidal rule in u, which converges exponentially for this
integrand. Phi comes from erfc: a positive-term series below 2, a continued
fraction above. c and delta enter as the exact values of their doubles.
"""
from __future__ import annotations

import json
import math
import platform
import sys
import warnings
from decimal import Decimal, getcontext, localcontext
from pathlib import Path

import numpy as np
import scipy
import scipy.special as sc
from scipy import stats

_HERE = Path(__file__).resolve().parent
DF = 15
C = float(stats.t.ppf(0.975, df=DF))
PRECISION = 70


# --------------------------------------------------------------------------- 70-digit reference
def _pi() -> Decimal:
    """Gauss-Legendre (Brent-Salamin) iteration."""
    a, b, t, p = Decimal(1), 1 / Decimal(2).sqrt(), Decimal(1) / 4, Decimal(1)
    for _ in range(9):
        a_next = (a + b) / 2
        b = (a * b).sqrt()
        t -= p * (a - a_next)**2
        a, p = a_next, 2 * p
    return (a + b)**2 / (4 * t)


def _erfc(z: Decimal, root_pi: Decimal) -> Decimal:
    if z < 0:
        return 2 - _erfc(-z, root_pi)
    if z <= 2:
        # erf(z) = 2/sqrt(pi) * exp(-z^2) * sum_n 2^n z^(2n+1) / (1*3*...*(2n+1)): all terms positive
        term, total, n = z, z, 0
        threshold = Decimal(10)**(-(getcontext().prec + 5))
        while term > threshold * total:
            n += 1
            term = term * 2 * z * z / (2 * n + 1)
            total += term
        return 1 - 2 / root_pi * (-z * z).exp() * total
    # erfc(z) = exp(-z^2)/sqrt(pi) * 1/(z + (1/2)/(z + (2/2)/(z + (3/2)/(z + ...)))), evaluated bottom-up
    depth = int(6000 / float(z * z)) + 40
    value = z
    for k in range(depth, 0, -1):
        value = z + Decimal(k) / 2 / value
    return (-z * z).exp() / (root_pi * value)


def _trapezoid(delta: float, step_denominator: int, pi: Decimal) -> tuple[Decimal, Decimal, int]:
    """Trapezoidal sums of the upper and lower tail integrands in u = ln w with step 1/step_denominator."""
    root_pi, root_two = pi.sqrt(), Decimal(2).sqrt()
    c, d = Decimal(C), Decimal(delta)  # exact values of the doubles
    nu = Decimal(DF)
    gamma = Decimal(135135) / 128 * root_pi  # Gamma(15/2) = 13!! / 2^7 * sqrt(pi)
    constant = 2 * (nu / 2)**7 * (nu / 2).sqrt() / gamma
    h = Decimal(1) / step_denominator
    low = -14 * step_denominator
    high = int(math.log(abs(delta) / C + 6) * step_denominator) + 1
    upper = lower = Decimal(0)
    for k in range(low, high + 1):
        w = (h * k).exp()
        weight = constant * w**DF * (-nu * w * w / 2).exp()  # f(w) * dw/du
        upper += weight * _erfc(-(c * w - d) / root_two, root_pi) / 2
        lower += weight * _erfc((c * w + d) / root_two, root_pi) / 2
    return upper * h, lower * h, high - low + 1


def reference(delta: float, target: str = "1e-45", largest_denominator: int = 2048) -> dict:
    """q, the upper term F(c; |delta|) and the lower term F(-c; |delta|), with the step halved until
    two successive values of q agree to the relative target (the integrand narrows as delta grows).

    q is even in delta (the law of -T under delta is the law of T under -delta), so the reference is
    taken at |delta|: there both terms are small or moderate and their difference loses no digits.
    At negative delta both terms are close to 1 and the same subtraction would cancel.
    """
    delta = abs(delta)
    with localcontext() as context:
        context.prec = PRECISION
        pi = _pi()
        denominator = 16
        upper, lower, nodes = _trapezoid(delta, denominator, pi)
        while True:
            denominator *= 2
            new_upper, new_lower, nodes = _trapezoid(delta, denominator, pi)
            q, previous = new_upper - new_lower, upper - lower
            change = abs(q - previous) / q if q > 0 else Decimal(0)
            upper, lower = new_upper, new_lower
            if change < Decimal(target) or denominator >= largest_denominator:
                return {"q": q, "upper": upper, "lower": lower, "nodes": nodes, "step": f"1/{denominator}",
                        "last_relative_change": float(change), "converged": bool(change < Decimal(target))}


# --------------------------------------------------------------------------- double-precision formulations
def frozen(delta: float) -> float:
    """protocol.json oracle.probability, verbatim."""
    c = stats.t.ppf(0.975, df=15)
    return float(stats.nct.cdf(c, df=15, nc=delta) - stats.nct.cdf(-c, df=15, nc=delta))


def frozen_terms(delta: float) -> tuple[float, float]:
    return float(stats.nct.cdf(C, df=DF, nc=delta)), float(stats.nct.cdf(-C, df=DF, nc=delta))


def reflected(delta: float) -> float:
    """F(-c; delta) = 1 - F(c; -delta): the lower term taken from the library's survival function."""
    return float(stats.nct.cdf(C, df=DF, nc=delta) - stats.nct.sf(C, df=DF, nc=-delta))


def bounded(delta: float) -> float:
    """The frozen expression at |delta| (q is even in delta), with a non-finite lower term set to 0.

    Justification: T <= -c < 0 needs Z + |delta| < 0, so 0 <= F(-c; |delta|) <= Phi(-|delta|).
    Where the frozen expression is finite at |delta| this returns it unchanged.
    """
    upper, lower = frozen_terms(abs(delta))
    return upper - (lower if math.isfinite(lower) else 0.0)


_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(200)


def mixture(delta: float) -> float:
    """q as the chi mixture of normal probabilities, by Gauss-Legendre quadrature on fixed panels.

    The integrand is Phi(c w - |delta|) - Phi(-c w - |delta|) times the density of W. Panels cover
    [0, |delta|/c + 6]; beyond it the density is below exp(-270) of its value at |delta|/c.
    The difference of the two normal probabilities is formed as ndtr(a) - ndtr(b) with a > b,
    which has no cancellation because the second term is the smaller one.
    """
    d = abs(float(delta))
    top = d / C + 6.0
    edges = np.linspace(0.0, top, max(8, int(math.ceil(top / 0.75))) + 1)
    total = 0.0
    log_constant = math.log(2.0) + 7.5 * math.log(7.5) - math.lgamma(7.5)
    for a, b in zip(edges[:-1], edges[1:]):
        w = (a + b) / 2 + (b - a) / 2 * _NODES
        with np.errstate(divide="ignore"):
            density = np.exp(log_constant + 14 * np.log(w) - 7.5 * w * w)
        inside = sc.ndtr(C * w - d) - sc.ndtr(-C * w - d)
        total += float(np.sum(_WEIGHTS * inside * density)) * (b - a) / 2
    return min(max(total, 0.0), 1.0)


CANDIDATES = {"frozen": frozen, "reflected": reflected, "bounded": bounded, "mixture": mixture}


# --------------------------------------------------------------------------- analysis
def failure_map() -> dict:
    """Where each term of the frozen expression is not finite, and what the library reports."""
    grid = np.round(np.arange(-60.0, 60.0001, 0.01), 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        upper = stats.nct.cdf(C, df=DF, nc=grid)
        lower = stats.nct.cdf(-C, df=DF, nc=grid)
        direct = sc.nctdtr(DF, grid, -C)
        survival = stats.nct.sf(C, df=DF, nc=-grid)
    bad = ~np.isfinite(lower)
    out = {"grid": [-60.0, 60.0, 0.01], "points": int(len(grid)),
           "upper_term_not_finite": int(np.sum(~np.isfinite(upper))),
           "lower_term_not_finite": int(bad.sum()),
           "lower_term_first_last": [float(grid[bad].min()), float(grid[bad].max())] if bad.any() else None,
           "lower_term_not_finite_for_negative_delta": int(np.sum(bad & (grid < 0))),
           "share_not_finite_inside_that_range": (float(bad[(grid >= grid[bad].min()) & (grid <= grid[bad].max())].mean())
                                                  if bad.any() else None),
           "nct_cdf_is_special_nctdtr": bool(np.array_equal(np.isnan(direct), np.isnan(lower))
                                             and np.array_equal(direct[~bad], lower[~bad])),
           "reflected_survival_not_finite": int(np.sum(~np.isfinite(survival))),
           "frozen_q_not_finite": int(np.sum(~np.isfinite(upper - lower)))}
    # the library's own error category at failing and at neighbouring points
    reports = {}
    for delta in (8.54, 8.55, 9.5, 9.75, 10.0, 20.0, 43.6, 43.61):
        try:
            with sc.errstate(all="raise"):
                reports[str(delta)] = repr(float(sc.nctdtr(DF, delta, -C)))
        except Exception as error:  # noqa: BLE001 - the message is the evidence
            reports[str(delta)] = f"{type(error).__name__}: {error}"
    out["library_report_for_lower_term"] = reports
    # other arguments: the failure is a property of the negative-argument branch, not of df = 15 or of c
    region = {}
    deltas = np.round(np.arange(0.0, 60.0001, 0.05), 2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for df in (5, 15, 30):
            for x in (-0.5, -C, -5.0):
                values = sc.nctdtr(df, deltas, x)
                nan = ~np.isfinite(values)
                region[f"df={df},x={x:.4g}"] = {"not_finite": int(nan.sum()),
                                               "first": float(deltas[nan].min()) if nan.any() else None,
                                               "last": float(deltas[nan].max()) if nan.any() else None}
    out["other_arguments"] = region
    return out


def against_reference() -> dict:
    points = [0.0, 1e-300, 1e-8, 0.25, 1.0, 2.0, 3.0, 5.0, 7.0, 8.0, 8.5, 8.54, 8.55, 8.6, 9.0, 9.5, 9.75, 10.0, 12.0,
              15.0, 20.0, 25.0, 30.0, 37.0, 40.0, 43.6, 45.0, 60.0, 100.0, -1.0, -5.0, -10.0, -20.0]
    rows = []
    worst = {name: {"abs": 0.0, "rel_where_q_above_1e-12": 0.0, "not_finite": 0} for name in CANDIDATES}
    step_agreement = 0.0
    for delta in points:
        fine = reference(delta)
        q = fine["q"]
        step_agreement = max(step_agreement, fine["last_relative_change"])
        row = {"delta": delta, "reference_q": f"{q:.40E}", "reference_upper": f"{fine['upper']:.20E}",
               "reference_lower": f"{fine['lower']:.20E}", "normal_bound_on_lower": float(stats.norm.cdf(-abs(delta))),
               "nodes": fine["nodes"], "step": fine["step"], "reference_converged": fine["converged"],
               "reference_last_relative_change": fine["last_relative_change"], "candidates": {}}
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            terms = frozen_terms(delta)
            row["frozen_upper"], row["frozen_lower"] = (None if not math.isfinite(v) else v for v in terms)
            for name, function in CANDIDATES.items():
                value = function(delta)
                if not math.isfinite(value):
                    row["candidates"][name] = None
                    worst[name]["not_finite"] += 1
                    continue
                error = abs(Decimal(value) - q)
                relative = float(error / q) if q > 0 else None
                row["candidates"][name] = {"value": value, "abs_error": float(error), "rel_error": relative}
                worst[name]["abs"] = max(worst[name]["abs"], float(error))
                if q > Decimal("1e-12"):
                    worst[name]["rel_where_q_above_1e-12"] = max(worst[name]["rel_where_q_above_1e-12"], relative)
        rows.append(row)
    return {"precision_digits": PRECISION, "reference_last_relative_change_max": step_agreement,
            "reference_converged_everywhere": all(r["reference_converged"] for r in rows),
            "points": rows, "worst": worst}


def properties() -> dict:
    """Bounds, evenness, monotonicity and continuity of each formulation on a dense grid."""
    grid = np.round(np.arange(0.0, 60.0001, 0.01), 2)
    out = {}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for name, function in CANDIDATES.items():
            values = np.array([function(float(d)) for d in grid])
            mirror = np.array([function(-float(d)) for d in grid[::25]])
            finite = np.isfinite(values)
            v = values[finite]
            steps = np.diff(v)
            both = np.isfinite(mirror) & finite[::25]
            out[name] = {"not_finite": int((~finite).sum()), "below_0": int(np.sum(v < 0)), "above_1": int(np.sum(v > 1)),
                         "minimum": float(v.min()), "maximum": float(v.max()),
                         "increases_between_adjacent_finite_points": int(np.sum(steps > 0)),
                         "largest_increase": float(steps.max()) if len(steps) else None,
                         "largest_step_down": float(-steps.min()) if len(steps) else None,
                         "evenness_max_abs_difference": float(np.max(np.abs(values[::25][both] - mirror[both]))) if both.any() else None,
                         "mirror_not_finite": int(np.sum(~np.isfinite(mirror)))}
        dense = np.array([mixture(float(d)) for d in grid])
        other = np.array([bounded(float(d)) for d in grid])
        positive = np.array([frozen(float(d)) for d in grid])
        patched = np.array([bounded(float(d)) for d in grid])
        swapped = np.array([reflected(float(d)) for d in grid])
        negative = np.array([frozen(-float(d)) for d in grid])
    available = np.isfinite(positive)
    out["bounded_equals_frozen_bitwise_where_frozen_is_finite"] = bool(np.array_equal(positive[available], patched[available]))
    out["reflected_against_frozen_where_finite"] = {
        "bitwise_equal_points": int(np.sum(positive[available] == swapped[available])), "points": int(available.sum()),
        "max_abs_difference": float(np.max(np.abs(positive[available] - swapped[available])))}
    out["frozen_at_negative_delta"] = {"not_finite": int(np.sum(~np.isfinite(negative))),
                                       "max_abs_difference_from_mixture": float(np.nanmax(np.abs(negative - dense)))}
    out["bound_on_dropped_term_at_first_failure"] = float(stats.norm.cdf(-8.53))
    out["mixture_against_bounded"] = {"max_abs_difference": float(np.max(np.abs(dense - other))),
                                      "at_delta": float(grid[int(np.argmax(np.abs(dense - other)))])}
    out["value_at_zero"] = {"mixture": mixture(0.0), "frozen": frozen(0.0), "central_t": float(2 * stats.t.cdf(C, DF) - 1)}
    out["large_noncentrality"] = {str(d): {"mixture": mixture(d), "bounded": bounded(d)} for d in (100.0, 1e3, 1e6)}
    # convergence of the fixed quadrature: node count
    convergence = {}
    global _NODES, _WEIGHTS
    saved = (_NODES, _WEIGHTS)
    for nodes in (25, 50, 100, 200, 400):
        _NODES, _WEIGHTS = np.polynomial.legendre.leggauss(nodes)
        convergence[str(nodes)] = {str(d): mixture(d) for d in (0.0, 3.0, 8.55, 20.0, 40.0)}
    _NODES, _WEIGHTS = saved
    out["mixture_convergence_in_nodes_per_panel"] = convergence
    return out


def widening_is_separate() -> dict:
    """F-02 is a different matter: the same formula with a spread of batch means (lambda > 0)."""
    sys.path.insert(0, str(_HERE.parents[1] / "src"))
    from rnglab.prognosis.validity import coverage_probability
    rows = {}
    for delta in (1.0, 3.0, 5.0, 8.55, 12.0):
        rows[str(delta)] = {"lambda_0": float(coverage_probability(delta, 0.0)), "lambda_5": float(coverage_probability(delta, 5.0)),
                            "lambda_15": float(coverage_probability(delta, 15.0)), "mixture": mixture(delta)}
    return rows


def main() -> int:
    out = {"schema": "rnglab-v019-audit-f01-v1",
           "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                           "libc": " ".join(platform.libc_ver())},
           "critical_value": {"c": C, "c_exact_decimal": str(Decimal(C))},
           "failure_map": failure_map(), "against_reference": against_reference(), "properties": properties(),
           "widening_is_separate": widening_is_separate()}
    target = _HERE / "results/f01_oracle.json"
    target.write_text(json.dumps(out, indent=1, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({k: out[k] for k in ("failure_map",)}, indent=1)[:2500])
    print(json.dumps(out["against_reference"]["worst"], indent=1))
    print(json.dumps({k: v for k, v in out["properties"].items() if k != "mixture_convergence_in_nodes_per_panel"}, indent=1))
    print(json.dumps(out["properties"]["mixture_convergence_in_nodes_per_panel"]))
    print(json.dumps(out["widening_is_separate"]))
    print("reference converged everywhere:", out["against_reference"]["reference_converged_everywhere"],
          "largest last change:", out["against_reference"]["reference_last_relative_change_max"])
    for row in out["against_reference"]["points"]:
        cells = " ".join(f"{name}={'NaN' if value is None else format(value['abs_error'], '.1e')}"
                         f"/{'-' if value is None or value['rel_error'] is None else format(value['rel_error'], '.1e')}"
                         for name, value in row["candidates"].items())
        if row["delta"] not in (0.0, 8.55, 20.0, 43.6, 100.0, -20.0):
            continue
        print(f"delta={row['delta']:>8} q={row['reference_q'][:12]}{row['reference_q'][-6:]} lower={row['reference_lower'][:8]}{row['reference_lower'][-5:]} "
              f"step={row['step']} | {cells}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
