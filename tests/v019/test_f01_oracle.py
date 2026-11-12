"""F-01: the Stage D oracle probability, its failure, a 70-digit reference and stable equivalents.

The quantity is q(delta) = P(|T| <= c), T ~ noncentral t(15, delta). Reference values are read
from research/v019-audit/results/f01_oracle.json (written by f01_oracle.py with Python's decimal
module); a few are recomputed here to check that the stored file is what the code produces.
"""
import importlib.util as _ilu
import json
import math
import sys as _sys
import warnings
from decimal import Decimal, localcontext
from pathlib import Path as _Path

import numpy as np
import pytest
import scipy.special as sc
from scipy import stats

_spec = _ilu.spec_from_file_location("_audit_loader", _Path(__file__).with_name("_audit_loader.py"))
au = _sys.modules.get("_audit_loader") or _ilu.module_from_spec(_spec)
if "_audit_loader" not in _sys.modules:
    _sys.modules["_audit_loader"] = au
    _spec.loader.exec_module(au)

f01 = au.load("v019_audit_f01", au.REPO / "research/v019-audit/f01_oracle.py")
rc = au.recount()
STORED = json.loads((au.REPO / "research/v019-audit/results/f01_oracle.json").read_text())
POINTS = {row["delta"]: row for row in STORED["against_reference"]["points"]}


def quiet(function, *args):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return function(*args)


# --------------------------------------------------------------------------- the reference itself
def test_decimal_building_blocks():
    with localcontext() as context:
        context.prec = 60
        pi = f01._pi()
        assert str(pi)[:52] == "3.14159265358979323846264338327950288419716939937510"
        root_pi = pi.sqrt()
        for z in (0.0, 0.3, 1.0, 1.999, 2.0, 2.001, 3.5, 6.0, 12.0, 26.0, -0.7, -3.0):
            value = f01._erfc(Decimal(z), root_pi)
            assert float(value) == pytest.approx(math.erfc(z), rel=5e-15, abs=0)
        # the series and the continued fraction meet at their switch point
        below = f01._erfc(Decimal("2"), root_pi)
        above = f01._erfc(Decimal("2.0000000000000000000000001"), root_pi)
        assert abs(below - above) / below < Decimal("1e-23")
        # erfc(z) + erfc(-z) = 2 exactly
        assert f01._erfc(Decimal("1.3"), root_pi) + f01._erfc(Decimal("-1.3"), root_pi) == 2


def test_reference_agrees_with_the_central_t_and_is_converged():
    zero = f01.reference(0.0)
    assert zero["converged"] and zero["last_relative_change"] < 1e-45
    assert float(zero["q"]) == pytest.approx(2 * stats.t.cdf(f01.C, 15) - 1, abs=3e-16)
    assert float(zero["upper"] + zero["lower"]) == pytest.approx(1.0, abs=1e-15)  # symmetry of the central t
    assert STORED["against_reference"]["reference_converged_everywhere"]
    assert STORED["against_reference"]["reference_last_relative_change_max"] < 1e-45
    assert f"{zero['q']:.40E}" == POINTS[0.0]["reference_q"]  # the stored file is what the code computes


def test_reference_point_is_reproduced_and_even():
    for delta in (8.55, -10.0):
        again = f01.reference(delta)
        assert f"{again['q']:.40E}" == POINTS[delta]["reference_q"]
    assert POINTS[-10.0]["reference_q"] == POINTS[10.0]["reference_q"]
    assert POINTS[-20.0]["reference_q"] == POINTS[20.0]["reference_q"]


def test_reference_lower_term_respects_the_normal_bound():
    """T <= -c < 0 needs Z + delta < 0, so F(-c; delta) <= Phi(-delta) for delta >= 0."""
    for delta, row in POINTS.items():
        if delta >= 0:
            assert Decimal(row["reference_lower"]) <= Decimal(row["normal_bound_on_lower"]) * Decimal("1.0000000001") + Decimal("1e-320")
    assert float(Decimal(POINTS[8.55]["reference_lower"])) == pytest.approx(3.588e-24, rel=1e-3)
    assert stats.norm.cdf(-8.53) < 7.4e-18


# --------------------------------------------------------------------------- the frozen call
def test_frozen_call_fails_in_the_lower_term_only_with_a_library_no_result_error():
    failures = []
    for delta in np.round(np.arange(0.0, 45.01, 0.25), 2):
        upper, lower = quiet(f01.frozen_terms, float(delta))
        assert math.isfinite(upper)                       # the upper term never fails
        if not math.isfinite(lower):
            failures.append(float(delta))
    if not failures:
        pytest.skip("this SciPy no longer returns NaN for the lower term; F-01 does not occur here")
    assert min(failures) >= 8.5 and max(failures) <= 43.75
    with pytest.raises(sc.SpecialFunctionError, match="no result"):
        with sc.errstate(all="raise"):
            sc.nctdtr(15, failures[0], -f01.C)
    for delta in failures:                                # never for the mirrored sign
        assert all(math.isfinite(v) for v in quiet(f01.frozen_terms, -delta))
        assert math.isnan(quiet(f01.frozen, delta)) and math.isnan(quiet(rc.validity_probability_frozen, delta))


def test_stored_failure_map_is_consistent():
    m = STORED["failure_map"]
    assert m["upper_term_not_finite"] == 0 and m["lower_term_not_finite_for_negative_delta"] == 0
    assert m["nct_cdf_is_special_nctdtr"] and m["reflected_survival_not_finite"] == 0
    assert m["frozen_q_not_finite"] == m["lower_term_not_finite"]
    if m["lower_term_not_finite"]:
        assert m["lower_term_first_last"] == [8.53, 43.61]
        assert "no result" in m["library_report_for_lower_term"]["8.55"]


# --------------------------------------------------------------------------- stable equivalents
@pytest.mark.parametrize("name", ["bounded", "reflected", "mixture"])
def test_equivalents_are_finite_bounded_and_accurate(name):
    function = f01.CANDIDATES[name]
    for delta, row in POINTS.items():
        value = quiet(function, delta)
        q = Decimal(row["reference_q"])
        assert math.isfinite(value) and 0.0 <= value <= 1.0
        assert abs(Decimal(value) - q) < Decimal("3e-15")                 # absolute criterion, everywhere
        if q > Decimal("1e-250") and (name == "mixture" or delta >= 0):   # relative criterion away from underflow
            assert abs(Decimal(value) - q) / q < Decimal("2e-13")
    grid = np.round(np.arange(0.0, 60.01, 0.05), 2)
    values = np.array([quiet(function, float(d)) for d in grid])
    assert np.isfinite(values).all() and values.min() >= 0 and values.max() <= 0.95 + 1e-14
    assert np.all(np.diff(values) <= 1e-130)   # non-increasing in |delta| up to noise far below any score's resolution
    assert quiet(function, 1e6) == 0.0 and quiet(function, 0.0) == pytest.approx(0.95, abs=1e-14)


def test_bounded_changes_nothing_where_the_frozen_value_exists():
    grid = np.round(np.arange(0.0, 60.01, 0.05), 2)
    filled = 0
    for delta in grid:
        frozen, bounded = quiet(f01.frozen, float(delta)), f01.bounded(float(delta))
        if math.isfinite(frozen):
            assert bounded == frozen                       # bit for bit
        else:
            filled += 1
            upper, _ = quiet(f01.frozen_terms, float(delta))
            assert bounded == upper and upper <= stats.norm.cdf(-8.5) * 1e9  # the upper term alone; tiny
    assert f01.bounded(-10.0) == f01.bounded(10.0)
    if filled == 0:
        pytest.skip("no failure in this SciPy")


def test_mixture_is_exactly_even_and_monotone_and_converged_in_its_nodes():
    grid = np.round(np.arange(0.0, 60.01, 0.05), 2)
    values = np.array([f01.mixture(float(d)) for d in grid])
    assert np.all(np.diff(values) <= 0) and all(f01.mixture(-float(d)) == v for d, v in zip(grid[::40], values[::40]))
    table = STORED["properties"]["mixture_convergence_in_nodes_per_panel"]
    for delta in ("0.0", "3.0", "8.55", "20.0", "40.0"):
        column = [table[n][delta] for n in ("25", "50", "100", "200", "400")]
        assert max(column) / min(column) - 1 < 1e-12   # 25 nodes per panel already suffice
    assert rc.validity_probability_stable(12.0) == pytest.approx(f01.mixture(12.0), rel=1e-12)
    assert rc.validity_probability_stable(12.0) == pytest.approx(float(Decimal(POINTS[12.0]["reference_q"])), rel=1e-13)


def test_widening_is_a_different_question():
    """F-02: with a spread of batch means the probability is a different, larger number. F-01 is about
    evaluating the frozen formula (lambda = 0), not about which formula is the better model."""
    rows = STORED["widening_is_separate"]
    for delta, row in rows.items():
        assert row["lambda_0"] == pytest.approx(row["mixture"], rel=1e-9, abs=1e-15)
        assert row["lambda_0"] < row["lambda_5"] < row["lambda_15"]
    assert rows["3.0"]["lambda_15"] > 2 * rows["3.0"]["lambda_0"]


def test_failure_depends_on_the_last_bits_of_delta():
    """Whether the frozen call returns a number or NaN changes with a rounding-level change of delta.
    Two correct programs that compute delta in a different order of operations can therefore disagree
    on which runs have an oracle value, although both follow the frozen rule."""
    stored = STORED["rounding_sensitivity"]
    assert stored["points"] == 20000 and stored["range"] == [8.6, 43.5]
    grid = np.linspace(8.6, 43.5, 20000)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        here = np.isnan(sc.nctdtr(15, grid, -f01.C))
        neighbour = np.isnan(sc.nctdtr(15, np.nextafter(grid, np.inf), -f01.C))
    if not here.any():
        pytest.skip("this SciPy does not return NaN here")
    flips = np.flatnonzero(here != neighbour)
    assert stored["not_finite"] == int(here.sum()) and stored["flips_at_next_double"] == len(flips) > 100
    assert 0.01 < stored["share_flipping_at_1e-13"] < 0.10
    for index in flips[:: max(1, len(flips) // 25)]:
        a, b = float(grid[index]), float(np.nextafter(grid[index], np.inf))
        assert math.isnan(quiet(f01.frozen, a)) != math.isnan(quiet(f01.frozen, b))     # one has a value, the other has none
        assert abs(f01.bounded(a) - f01.bounded(b)) <= 1e-12 * f01.bounded(a) + 1e-300  # the filled rule does not jump
        assert abs(f01.mixture(a) - f01.mixture(b)) <= 1e-12 * f01.mixture(a) + 1e-300
