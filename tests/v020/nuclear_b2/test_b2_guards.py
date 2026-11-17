"""Regression tests for the second review of B2 (numerical correctness and statistical methodology).

Each test reproduces a confirmed finding and checks its correction:
- non-finite incomplete-beta values (NaN near n = 2**53) must raise, never steer a draw; the part limit
  2**50 keeps the sampler where Boost is finite; draws do not depend on the initial bracket;
- floating-point distribution functions decide the draw only up to their own error, so a draw can
  differ from exact midpoint inversion next to a step and nowhere else;
- moments from centred values (no cancellation), counts as exact integers (no int64 overflow);
- adaptive decimal precision tries the maximum before refusing;
- a non-finite reference fails D8 and is never counted as agreement;
- fixed-cell homogeneity and the per-interval photon test.

Draws only from the B2 test namespace (or none at all).
"""
import importlib.util
import math
import sys
from decimal import Context, Decimal, localcontext
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from scipy import special

from rnglab.nuclear import Channel, DecayScheme, Mode, Nuclide, state_probabilities
from rnglab.nuclear import binomial as binomial_module
from rnglab.nuclear.aggregate import AggregatePlan
from rnglab.nuclear.binomial import PART_MAX, UPPER_TAIL_SWITCH, binomial_inverse, split_counts
from rnglab.nuclear.comparison import allocation_check, consecutive_cells, homogeneity_check
from rnglab.nuclear.moments import PairMoments, categorical_pair_moments
from rnglab.nuclear.paths import MAX_DIGITS, bateman_decimal
from rnglab.seeds import namespace_seed

from .conftest import B2_TEST_NAMESPACE

REPO = Path(__file__).resolve().parents[3]


def _load(name: str, relative: str):
    spec = importlib.util.spec_from_file_location(name, REPO / relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


# --- 1. Non-finite distribution values, the part limit, bracket invariance -------------------------------

def test_part_limit_is_two_to_the_fifty():
    assert PART_MAX == 2 ** 50
    with pytest.raises(ValueError):
        binomial_inverse(2 ** 53 - 1, 0.5, 0.5)     # where Boost returns NaN; refused before any evaluation
    assert split_counts(2 ** 53 - 1, 8) == [PART_MAX] * 7 + [PART_MAX - 1]


def test_a_nan_from_the_lower_tail_raises(monkeypatch):
    real = special.betaincc

    def poisoned(a, b, x):
        out = np.array(real(a, b, x), dtype=np.float64)
        out[np.asarray(a) > np.asarray(b)] = np.nan   # NaN on one side of the centre, as Boost did
        return out
    monkeypatch.setattr(binomial_module.special, "betaincc", poisoned)
    with pytest.raises(FloatingPointError, match="F not finite"):
        binomial_inverse(10 ** 6 + 1, 0.5, 0.5)


def test_a_nan_from_the_upper_tail_raises(monkeypatch):
    monkeypatch.setattr(binomial_module.special, "betainc", lambda a, b, x: np.full(np.shape(a), np.nan))
    u = 1.0 - UPPER_TAIL_SWITCH / 4                # the far upper tail, consulted through S up to n = 1e9
    with pytest.raises(FloatingPointError, match="S not finite"):
        binomial_inverse(10 ** 6, 0.3, u)


def test_nan_inside_the_recurrence_table_raises(monkeypatch):
    monkeypatch.setattr(binomial_module.np, "cumprod", lambda x: np.full(np.shape(x), np.nan))
    with pytest.raises(FloatingPointError, match="recurrence table"):
        binomial_inverse(100, 0.3, 0.5)


@pytest.mark.parametrize("n", [10 ** 9 + 1, 10 ** 15 + 1, PART_MAX - 1])
def test_exact_symmetry_anchor(n):
    """Odd n, p = 1/2: F(m) = 1/2 exactly, so midpoint inversion gives m + 1 at u = 1/2 and m just below.
    Before the guard, n = 2**53 - 1 gave m + 3,843 at u = 1/2 with the default bracket."""
    m = (n - 1) // 2
    got = binomial_inverse(np.full(2, n), np.full(2, 0.5), np.array([0.5, 0.5 - 2.0 ** -53])).tolist()
    assert got == [m + 1, m]


def test_boost_is_finite_where_the_bisection_reaches_at_the_part_limit():
    for n in (PART_MAX - 1, PART_MAX):
        for p in (0.5, 0.3, 0.2):
            mean, sd = n * p, math.sqrt(n * p * (1 - p))
            ks = {int(math.floor(mean + sd * z)) for z in np.linspace(-9.0, 9.0, 13)}
            ks.update(range(int(mean) - 20, int(mean) + 21))
            ks = np.array(sorted(ks), dtype=np.float64)
            assert np.all(np.isfinite(special.betaincc(ks + 1.0, n - ks, p)))


def test_draws_do_not_depend_on_the_initial_bracket():
    us = np.array([0.0, 1e-9, 0.37, 0.5, 0.81, 1.0 - 1e-12])
    for n, p in ((PART_MAX, 0.5), (PART_MAX - 1, 0.1), (10 ** 12, 0.3)):
        draws = [binomial_inverse(np.full(us.size, n), np.full(us.size, p), us, half_width=w) for w in (1, 4, 64)]
        assert all(np.array_equal(d, draws[1]) for d in draws)


# --- 2. Floating-point F against exact inversion next to steps ----------------------------------------

def _exact_cdf(n: int, p: float) -> list[Fraction]:
    P = Fraction(p)
    total, out = Fraction(0), []
    for x in range(n + 1):
        total += math.comb(n, x) * P ** x * (1 - P) ** (n - x)
        out.append(total)
    return out


def _exact_draw(cdf: list[Fraction], k: int) -> int:
    v = Fraction(2 * k + 1, 2 ** 54)
    return next(x for x, f in enumerate(cdf) if f > v)


def test_draws_differ_from_exact_inversion_only_next_to_a_step():
    """The evaluated F errs by a few ulps; a draw may then fall on the other side of a step for the grid
    uniforms adjacent to it, never further away (review: n = 11, p = 1/2, u = 1/2 gives 5, exactly 6)."""
    assert int(binomial_inverse(11, 0.5, 0.5)) in (5, 6)
    adjacent = mismatches_far = 0
    for n in range(1, 41):
        for p in (0.5, 0.25, 0.1):
            cdf = _exact_cdf(n, p)
            for x in range(n):
                if not Fraction(1, 2 ** 54) < cdf[x] < 1 - Fraction(1, 2 ** 54):
                    continue
                k = math.ceil(cdf[x] * 2 ** 53 - Fraction(1, 2)) - 1       # u = k 2**-53 is the last cell below
                near = [k, k + 1]
                far = [k - 8, k + 9]
                for ks, is_near in ((near, True), (far, False)):
                    ks = [j for j in ks if 0 <= j < 2 ** 53]
                    got = binomial_inverse(np.full(len(ks), n), np.full(len(ks), p),
                                           np.array([j * 2.0 ** -53 for j in ks])).tolist()
                    for j, g in zip(ks, got):
                        want = _exact_draw(cdf, j)
                        if g != want:
                            assert is_near and abs(g - want) == 1, (n, p, x, j, g, want)
                            adjacent += 1
                        elif not is_near:
                            pass
                    if not is_near:
                        mismatches_far += sum(g != _exact_draw(cdf, j) for j, g in zip(ks, got))
    assert mismatches_far == 0
    assert adjacent > 0      # the floating-point CDF does decide some boundary cases differently


# --- 3. Moments -----------------------------------------------------------------------------------------

def test_rare_category_moments_keep_relative_accuracy():
    """Review: E[fg] - E[f] E[g] gave sigma_ab = -2048 and sigma_aa = 512 for the same pair at N = 2**62."""
    p = 2.0 ** -53
    m = categorical_pair_moments(np.array([p, 1 - p]), np.array([1.0, 2.0]), np.array([1.0, 2.0]), 2 ** 62)
    exact = 2 ** 62 * p * (1 - p)
    assert m.sigma_aa == pytest.approx(exact, rel=1e-12) and m.sigma_ab == pytest.approx(exact, rel=1e-12)


def test_moments_are_translation_invariant():
    rng = np.random.default_rng(namespace_seed(B2_TEST_NAMESPACE, "guards", "translation"))
    pi = rng.dirichlet(np.ones(6))
    f, g = rng.integers(0, 3, 6).astype(float), rng.integers(0, 3, 6).astype(float)
    a = categorical_pair_moments(pi, f, g, 10 ** 6)
    b = categorical_pair_moments(pi, f + 1e6, g - 3e5, 10 ** 6)
    for x, y in ((a.sigma_aa, b.sigma_aa), (a.sigma_bb, b.sigma_bb), (a.sigma_ab, b.sigma_ab), (a.mu_22, b.mu_22)):
        assert y == pytest.approx(x, rel=1e-6, abs=1e-6 * abs(a.sigma_aa))


def test_numpy_integer_counts_do_not_overflow():
    """Review: np.int64(1e12) made N (N - 1) wrap around and Var(s_ab) negative."""
    probabilities, values = np.array([0.5, 0.5]), np.array([0.0, 1.0])
    want = categorical_pair_moments(probabilities, values, values, 10 ** 12).sample_covariance_variance(400)
    got = categorical_pair_moments(probabilities, values, values, np.int64(10 ** 12)).sample_covariance_variance(400)
    sigma2 = 10 ** 12 * 0.25
    expected = (3 * sigma2 ** 2 - 10 ** 12 * 0.125 - sigma2 ** 2) / 400 + 2 * sigma2 ** 2 / (400 * 399)
    assert got == want > 0 and got == pytest.approx(expected, rel=1e-12)


def test_counts_and_replicates_must_be_integers():
    pi, v = np.array([0.5, 0.5]), np.array([0.0, 1.0])
    for bad in (-1, 2.5, True):
        with pytest.raises((TypeError, ValueError)):
            categorical_pair_moments(pi, v, v, bad)
    moments = categorical_pair_moments(pi, v, v, 10)
    with pytest.raises(ValueError):
        moments.sample_covariance_variance(1)
    with pytest.raises(ValueError):
        categorical_pair_moments(np.array([np.nan, 1.0]), v, v, 10)
    assert isinstance(moments, PairMoments)


# --- 4. Adaptive precision --------------------------------------------------------------------------------

def test_adaptive_precision_tries_the_maximum_before_refusing():
    """Review: 36 decay constants 2**-52 apart were refused although 700 digits suffice."""
    half_lives = [1 + i * 2.0 ** -52 for i in range(36)] + [math.inf]
    value, digits = bateman_decimal(half_lives, Decimal(1))
    reference, _ = bateman_decimal(half_lives, Decimal(1), min_digits=700)
    assert digits <= MAX_DIGITS
    with localcontext(Context(prec=60)):
        assert abs(value / reference - 1) < Decimal("1e-40")


def test_two_close_constants_and_a_tiny_interval():
    value, digits = bateman_decimal([1.0, math.nextafter(1.0, 2.0), math.inf], Decimal("1e-300"))
    assert value > 0 and digits <= MAX_DIGITS


# --- 5. Non-finite references -----------------------------------------------------------------------------

def _fast_scheme() -> DecayScheme:
    return DecayScheme("fast", [Nuclide("A", 150, 400, 1e-300, (Channel("B", Mode.BETA_MINUS, 1.0),)),
                                Nuclide("B", 151, 400, math.inf)])


def test_plan_diagnostics_never_count_a_nan_as_agreement():
    plan = AggregatePlan(_fast_scheme(), [0.0, 1e308], 10)
    assert plan.diagnostics["tree_comparisons"] == 0
    assert plan.diagnostics["tree_comparisons_unavailable"] == 1
    assert plan.interval_plans[1e308][0].law.as_float().tolist() == [0.0, 1.0]


def test_d8_fails_on_a_nonfinite_reference():
    verify = _load("b2_guards_run_verification", "research/v020-nuclear-b2-verify/run_verification.py")
    result = verify.law_checks({"fast": _fast_scheme()}, {"fast": {"times": [0.0, 1e308]}})
    # The Layer A matrix exponential of the whole scheme overflows, so both nuclides' references are non-finite.
    assert result["ok"] is False and result["nonfinite_references"] == ["fast A delta=1e+308", "fast B delta=1e+308"]
    assert result["max_abs_decimal_minus_tree"] == 0.0 and result["laws_checked"] == 2


def test_an_uncomputable_test_is_an_error_not_a_pass():
    verify = _load("b2_guards_run_verification2", "research/v020-nuclear-b2-verify/run_verification.py")
    result = verify.evaluate(homogeneity_check, [0, 0], [3, 4], cells=[[0], [1]])
    assert "error" in result and "p_value" not in result
    assert verify.errors_of({"x": result, "y": {"p_value": 0.5}}) == ["x"]


# --- 6. Fixed-cell homogeneity and per-interval photons ------------------------------------------------------

def test_fixed_cells_make_homogeneity_independent_of_the_margins():
    """Review: margin-based pooling made [100, 4] against itself undefined and [99, 5] a two-cell test."""
    for row in ([100, 4], [99, 5], [100, 0]):
        result = homogeneity_check(row, row, cells=[[0], [1]])
        assert result["dof"] == 1 and result["p_value"] == 1.0
    assert homogeneity_check([100, 0], [90, 0], cells=[[0], [1]])["empty_cells"] == 1
    with pytest.raises(ValueError):
        homogeneity_check([5, 1, 2], [5, 1, 0], cells=[[0], [1]])      # a count outside the cells


def test_consecutive_cells():
    assert consecutive_cells([3, 4, 5, 20, 1, 2], 10) == [[0, 1, 2], [3, 4, 5]]
    assert consecutive_cells([30, 1, 30], 10) == [[0], [1, 2]]
    assert consecutive_cells([1, 2], 10) == [[0, 1]]


def test_photon_allocation_over_time_is_tested():
    """Review: moving every photon to the earliest intervals (same total) passed every photon test."""
    decays = np.array([5000, 4000, 6000, 8000, 6000, 1600])
    cells = [[k] for k in range(6)]
    honest = np.array([2510, 1985, 3020, 3990, 2985, 790])
    moved = np.array([5000, 4000, 6000, 280, 0, 0])
    assert honest.sum() == moved.sum() == 15280
    assert allocation_check(honest, decays, 0.5, cells)["p_value"] > 0.05
    assert allocation_check(moved, decays, 0.5, cells)["p_value"] == 0.0
    with pytest.raises(ValueError):
        allocation_check(moved, decays, 0.5, [[0, 1, 2], [3, 4]])          # cells must partition the intervals


def test_photon_interval_cells_come_from_expectations_only(branching):
    verify = _load("b2_guards_run_verification3", "research/v020-nuclear-b2-verify/run_verification.py")
    case = {"initial": {"Xa-400": 100000}, "times": [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]}
    planned = verify.photon_interval_cells(branching, case)
    assert planned == verify.photon_interval_cells(branching, case)
    assert all(sorted(i for c in cells for i in c) == list(range(6)) for _, _, cells in planned)


def test_null_model_of_the_calibration_matches_exact_means(branching):
    calibration = _load("b2_guards_calibration", "research/v020-nuclear-b2-verify/calibration.py")
    times = [0.0, 0.5, 2.0, 8.0]
    model = calibration.NullModel(branching, times, 2000)
    rng = np.random.default_rng(namespace_seed(B2_TEST_NAMESPACE, "guards", "null-model"))
    draws = model.draw({"Xa-400": 2000}, rng, 4000)
    for it, t in enumerate(times[1:], start=1):
        p = state_probabilities(branching, "Xa-400", t)[0]
        mean = draws.populations[:, it, :].mean(axis=0)
        se = np.sqrt(2000 * p * (1 - p) / 4000)
        keep = se > 0
        assert np.all(np.abs(mean[keep] - 2000 * p[keep]) < 5 * se[keep])
    assert np.all(draws.populations.sum(axis=2) == 2000)
    assert np.all(draws.photons <= draws.decays)
