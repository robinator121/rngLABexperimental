"""B2 analytic moments (against exact enumeration), comparison checks and precursor groups."""
import itertools
import math

import numpy as np
import pytest
from scipy import stats

from rnglab.nuclear import state_probabilities
from rnglab.nuclear.comparison import (
    covariance_difference_check, homogeneity_check, plugin_covariance_variance, sample_covariance, z_check,
)
from rnglab.nuclear.decay import Mode
from rnglab.nuclear.moments import (
    categorical_pair_moments, state_count_moments, two_time_count_moments, two_time_law,
)
from rnglab.nuclear.precursors import PrecursorGroup, precursor_scheme


def test_pair_moments_and_sample_covariance_variance_by_exact_enumeration():
    """Two nuclei, three categories, three replicates: enumerate every outcome."""
    probs = np.array([0.5, 0.3, 0.2])
    f, g = np.array([1, 0, 1]), np.array([0, 1, 1])
    m = categorical_pair_moments(probs, f, g, count=2)
    replicate = []  # law of (A, B) for one replicate
    for outcome in itertools.product(range(3), repeat=2):
        weight = math.prod(probs[c] for c in outcome)
        replicate.append((weight, sum(f[c] for c in outcome), sum(g[c] for c in outcome)))
    ea = sum(w * a for w, a, _ in replicate)
    eb = sum(w * b for w, _, b in replicate)
    assert m.sigma_ab == pytest.approx(sum(w * (a - ea) * (b - eb) for w, a, b in replicate), abs=1e-15)
    assert m.mu_22 == pytest.approx(sum(w * (a - ea) ** 2 * (b - eb) ** 2 for w, a, b in replicate), abs=1e-15)
    r = 3
    mean_s = var_s = 0.0
    values = []
    for rows in itertools.product(replicate, repeat=r):
        w = math.prod(row[0] for row in rows)
        a = np.array([row[1] for row in rows], float)
        b = np.array([row[2] for row in rows], float)
        values.append((w, ((a - a.mean()) * (b - b.mean())).sum() / (r - 1)))
    mean_s = sum(w * s for w, s in values)
    var_s = sum(w * (s - mean_s) ** 2 for w, s in values)
    assert mean_s == pytest.approx(m.sigma_ab, abs=1e-14)
    assert var_s == pytest.approx(m.sample_covariance_variance(r), rel=1e-12)


def test_state_and_two_time_moments_are_consistent(branching):
    p2 = state_probabilities(branching, "Xa-400", 2.0)[0]
    m = state_count_moments(branching, "Xa-400", 1000, 2.0, 0, 1)
    assert m.sigma_ab == pytest.approx(-1000 * p2[0] * p2[1], rel=1e-12)
    assert m.sigma_aa == pytest.approx(1000 * p2[0] * (1 - p2[0]), rel=1e-12)
    joint = two_time_law(branching, "Xa-400", 2.0, 8.0)
    assert np.allclose(joint.sum(axis=1), p2, atol=1e-14)
    assert np.allclose(joint.sum(axis=0), state_probabilities(branching, "Xa-400", 8.0)[0], atol=1e-14)
    same = two_time_count_moments(branching, "Xa-400", 1000, 2.0, 8.0, 0, 0)
    p8 = state_probabilities(branching, "Xa-400", 8.0)[0]
    assert same.sigma_ab == pytest.approx(1000 * (p8[0] - p2[0] * p8[0]), rel=1e-10)


def test_z_and_homogeneity_checks():
    assert z_check(3.0, 1.0, 4.0)["z"] == 1.0
    result = homogeneity_check([30, 50, 20], [35, 40, 25])
    expected = stats.chi2_contingency(np.array([[30, 50, 20], [35, 40, 25]]), correction=False)
    assert result["p_value"] == pytest.approx(expected[1], rel=1e-12) and result["dof"] == 2
    pooled = homogeneity_check([500, 400, 2, 1, 0], [480, 420, 1, 3, 0])
    assert pooled["cells"] == [[0], [1, 2, 3]] and pooled["dof"] == 1
    with pytest.raises(ValueError):
        homogeneity_check([5, 0], [0, 0])
    with pytest.raises(ValueError):
        z_check(float("nan"), 0.0, 1.0)


def test_covariance_difference_and_plugin_variance():
    rng = np.random.default_rng(3)
    x = rng.normal(size=4000)
    y = 0.5 * x + rng.normal(size=4000)
    assert sample_covariance(x, y) == pytest.approx(np.cov(x, y, ddof=1)[0, 1], rel=1e-12)
    # For normal data Var(s_xy) = (sxx syy + sxy^2) / (R - 1); the plug-in estimate must be close.
    assert plugin_covariance_variance(x, y) == pytest.approx((1.0 * 1.25 + 0.25) / 3999, rel=0.1)
    assert covariance_difference_check(x, y, x, y)["p_value"] == 1.0


def test_precursor_groups(precursors):
    assert precursors.precursors() == pytest.approx({"G1": 0.3, "G2": 0.1, "G3": 1.0}, abs=1e-15)
    for nuclide in precursors.nuclides:
        assert nuclide.z >= 170
        for channel in nuclide.channels:
            assert channel.mode in (Mode.BETA_DELAYED_NEUTRON, Mode.BETA_MINUS)
    with pytest.raises(ValueError):
        PrecursorGroup("bad", 1.0, 0.0)
    with pytest.raises(ValueError):
        precursor_scheme([])


def test_scaled_chisquare_reduces_to_the_classical_test_for_normal_data():
    from rnglab.nuclear.comparison import scaled_chisquare_variance_check
    r, variance, s2 = 50, 2.0, 2.9
    result = scaled_chisquare_variance_check(s2, variance, 2 * variance ** 2 / (r - 1))
    statistic = (r - 1) * s2 / variance
    expected = 2 * min(stats.chi2.cdf(statistic, r - 1), stats.chi2.sf(statistic, r - 1))
    assert result["dof"] == pytest.approx(r - 1) and result["p_value"] == pytest.approx(expected, rel=1e-12)


def test_sum_moments_carry_the_covariance(branching):
    from rnglab.nuclear.moments import state_sum_moments, two_time_sum_moments
    pair = state_count_moments(branching, "Xa-400", 2000, 2.0, 0, 3)
    a = state_count_moments(branching, "Xa-400", 2000, 2.0, 0, 0).sigma_ab
    b = state_count_moments(branching, "Xa-400", 2000, 2.0, 3, 3).sigma_ab
    total = state_sum_moments(branching, "Xa-400", 2000, 2.0, 0, 3)
    assert total.sigma_ab == pytest.approx(a + b + 2 * pair.sigma_ab, rel=1e-12)
    two = two_time_count_moments(branching, "Xa-400", 2000, 2.0, 8.0, 0, 5)
    a8 = state_count_moments(branching, "Xa-400", 2000, 8.0, 5, 5).sigma_ab
    total2 = two_time_sum_moments(branching, "Xa-400", 2000, 2.0, 8.0, 0, 5)
    assert total2.sigma_ab == pytest.approx(a + a8 + 2 * two.sigma_ab, rel=1e-12)


def test_null_difference_check():
    from rnglab.nuclear.comparison import null_difference_check
    result = null_difference_check(5.0, 3.0, 2.0)
    assert result["z"] == pytest.approx(1.0) and result["p_value"] == pytest.approx(2 * stats.norm.sf(1.0))


def test_pooled_cells_follow_the_declared_rule():
    from rnglab.nuclear.comparison import pooled_cells
    assert pooled_cells([50, 30, 3, 3, 14], 5.0) == [[0], [1], [4], [2, 3]]
    assert pooled_cells([50, 30, 2, 1, 17], 5.0) == [[0], [1], [2, 3, 4]]
    assert pooled_cells([60, 38, 1, 1], 5.0) == [[0], [1, 2, 3]]
    assert pooled_cells([99, 1, 0], 5.0) == [[0, 1]]
    assert pooled_cells([0, 0], 5.0) == []
    with pytest.raises(ValueError):
        pooled_cells([1.0, float("nan")], 5.0)


def test_fixed_cells_chisquare_matches_scipy_and_refuses_bad_input():
    from rnglab.nuclear.comparison import fixed_cells_chisquare, summed_chisquare
    observed, probabilities = [48, 30, 2, 20], [0.5, 0.3, 0.02, 0.18]
    result = fixed_cells_chisquare(observed, probabilities, [[0], [1], [2, 3]])
    expected = stats.chisquare([48, 30, 22], [50, 30, 20])
    assert result["statistic"] == pytest.approx(expected.statistic, rel=1e-12) and result["dof"] == 2
    total = summed_chisquare([result, result])
    assert total["dof"] == 4 and total["p_value"] == pytest.approx(stats.chi2.sf(2 * expected.statistic, 4))
    for bad in (dict(cells=[[0, 1, 2, 3]]), dict(cells=[[0], [0, 1], [2, 3]]), dict(cells=[[0], [1]]),
                dict(observed=[0, 0, 0, 0]), dict(probabilities=[0.5, 0.3, float("nan"), 0.2])):
        arguments = dict(observed=observed, probabilities=probabilities, cells=[[0], [1], [2, 3]])
        arguments.update(bad)
        with pytest.raises(ValueError):
            fixed_cells_chisquare(**arguments)
    with pytest.raises(ValueError):
        summed_chisquare([])


def test_planned_path_cells_are_fixed_by_expected_counts(branching):
    from rnglab.nuclear.aggregate import AggregatePlan
    from rnglab.nuclear.comparison import planned_path_cells
    times = [0.0, 0.5, 2.0]
    plan = AggregatePlan(branching, times, 100000)
    expected = [state_probabilities(branching, "Xa-400", t)[0] * 100000 for t in times[:-1]]
    cells = planned_path_cells(plan, expected)
    assert set(cells[0]) == {0, 1}                    # the initial nuclide is tested in both intervals
    assert all(0 not in per for i, per in cells.items() if i != 0)   # nothing else exists at t = 0
    for per_interval in cells.values():
        for interval_cells in per_interval.values():
            assert len(interval_cells) >= 2
    again = planned_path_cells(plan, expected)
    assert again == cells


def test_variance_ratio_f_check():
    from rnglab.nuclear.comparison import variance_ratio_f_check
    r, variance = 60, 3.0
    result = variance_ratio_f_check(4.2, 2.8, variance, 2 * variance ** 2 / (r - 1))
    expected = 2 * min(stats.f.cdf(1.5, r - 1, r - 1), stats.f.sf(1.5, r - 1, r - 1))
    assert result["dof"] == pytest.approx(r - 1) and result["p_value"] == pytest.approx(expected, rel=1e-12)
    assert variance_ratio_f_check(2.0, 2.0, variance, 0.1)["p_value"] == pytest.approx(1.0, abs=1e-12)
    swapped = variance_ratio_f_check(2.8, 4.2, variance, 2 * variance ** 2 / (r - 1))
    assert swapped["p_value"] == pytest.approx(result["p_value"], rel=1e-12)
    with pytest.raises(ValueError):
        variance_ratio_f_check(0.0, 1.0, 1.0, 1.0)
