"""Layer A statistical checks: fixed conventions and hand-checked values."""
import math

import numpy as np
import pytest
from scipy import stats

from rnglab.nuclear.checks import binomial_check, holm, ks_check, multinomial_check, variance_ratio_check


def test_binomial_check_matches_scipy_and_handles_degenerate_probabilities():
    result = binomial_check(530, 1000, 0.5)
    assert result["p_value"] == pytest.approx(stats.binomtest(530, 1000, 0.5).pvalue, rel=1e-12)
    assert result["expected"] == 500.0 and result["z"] == pytest.approx(30 / math.sqrt(250))
    assert binomial_check(0, 10, 0.0)["p_value"] == 1.0 and binomial_check(1, 10, 0.0)["p_value"] == 0.0
    assert binomial_check(0, 0, 0.3)["p_value"] == 1.0
    with pytest.raises(ValueError):
        binomial_check(11, 10, 0.5)


def test_multinomial_check_pooling_rule():
    own_cell = multinomial_check([50, 30, 3, 3, 14], [0.5, 0.3, 0.03, 0.03, 0.14])
    assert own_cell["cells"] == [[0], [1], [4], [2, 3]]
    assert own_cell["dof"] == 3 and own_cell["cell_expected"] == pytest.approx([50, 30, 14, 6])
    into_smallest = multinomial_check([50, 30, 2, 1, 17], [0.5, 0.3, 0.02, 0.01, 0.17])
    assert into_smallest["cells"] == [[0], [1], [2, 3, 4]] and into_smallest["dof"] == 2
    merged = multinomial_check([60, 38, 1, 1], [0.6, 0.38, 0.01, 0.01])
    assert merged["cells"] == [[0], [1, 2, 3]] and merged["dof"] == 1


def test_multinomial_check_flags_impossible_outcomes_and_matches_scipy():
    assert multinomial_check([5, 5, 1], [0.5, 0.5, 0.0])["p_value"] == 0.0
    result = multinomial_check([48, 32, 20], [0.5, 0.3, 0.2])
    assert result["p_value"] == pytest.approx(stats.chisquare([48, 32, 20], [50, 30, 20]).pvalue, rel=1e-12)
    with pytest.raises(ValueError):
        multinomial_check([100, 1], [0.99, 0.01])
    with pytest.raises(ValueError):
        multinomial_check([1, 1], [0.5, 0.6])


def test_ks_and_variance_checks():
    sample = np.linspace(0.0005, 0.9995, 1000)
    assert ks_check(sample, lambda x: np.clip(x, 0, 1))["p_value"] > 0.99
    values = np.array([9.0, 11.0, 10.0, 10.0])
    result = variance_ratio_check(values, 2.0 / 3.0)
    assert result["ratio"] == pytest.approx(1.0) and result["dof"] == 3
    assert result["p_value"] == pytest.approx(2 * min(stats.chi2.cdf(3, 3), stats.chi2.sf(3, 3)))


def test_holm_wrapper():
    result = holm([0.001, 0.02, 0.04, 0.3], 0.05)
    assert result["reject"] == [True, False, False, False]
    assert result["adjusted_p_values"] == pytest.approx([0.004, 0.06, 0.08, 0.3])


def test_decay_invariants_hold_on_runs_and_catch_corruption(branching, competing, ns):
    import dataclasses

    from rnglab.nuclear import WordStream, labelled_source, simulate_decays
    from rnglab.nuclear.checks import decay_invariants

    for scheme, start in ((branching, "Xa-400"), (competing, "Pa-300")):
        for horizon in (0.0, 3.0, math.inf):
            run = simulate_decays(scheme, {start: 700}, horizon, WordStream(labelled_source(ns, "inv", start, horizon)))
            assert all(decay_invariants(run).values()), decay_invariants(run)
    run = simulate_decays(branching, {"Xa-400": 300}, 6.0, WordStream(labelled_source(ns, "inv", "corrupt")))

    def corrupted(column, change):
        events = {k: v.copy() for k, v in run.events.items()}
        change(events[column])
        return decay_invariants(dataclasses.replace(run, events=events))

    def bump(values):
        values[3] = values[3] + 1

    def late(values):
        values[0] = 7.0

    assert not corrupted("daughter", bump)["channel_defines_daughter_and_mode"]
    assert not corrupted("time", late)["times_within_horizon"]
    assert not corrupted("neutron", lambda v: v.__setitem__(0, not v[0]))["neutrons_only_from_beta_delayed_neutron"]


def test_non_finite_inputs_never_read_as_passing_checks():
    nan = float("nan")
    with pytest.raises(ValueError):
        multinomial_check([10, 10, 10], [nan, 0.5, 0.5])
    with pytest.raises(ValueError):
        variance_ratio_check([nan, 1.0], 1.0)
    with pytest.raises(ValueError):
        variance_ratio_check([1.0, 2.0], nan)
    with pytest.raises(ValueError):
        holm([nan, 0.5], 0.01)
    with pytest.raises(ValueError):
        ks_check(np.array([0.1, nan]), lambda x: x)
    with pytest.raises(ValueError):
        binomial_check(1, 10, nan)


@pytest.mark.parametrize("alpha", [float("nan"), 0.0, 1.0, -0.1, True])
def test_holm_refuses_an_invalid_alpha(alpha):
    with pytest.raises(ValueError):
        holm([0.001], alpha)
