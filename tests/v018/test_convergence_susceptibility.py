"""v0.18 split-R-hat and susceptibility estimators against explicit definitions."""
import math

import numpy as np
import pytest
from scipy.stats import norm, rankdata

from rnglab import convergence, susceptibility
from rnglab16.core import enumeration


def _rhat_reference(chains):
    """Plain loop transcription of rank-normalized split-R-hat (Blom scores) for one ensemble."""
    n = chains.shape[1] // 2
    halves = [c[:n] for c in chains] + [c[-n:] for c in chains]
    flat = np.concatenate(halves)
    z = norm.ppf((rankdata(flat) - 0.375) / (len(flat) + 0.25)).reshape(len(halves), n)
    w = np.mean([np.var(h, ddof=1) for h in z])
    b = n * np.var([h.mean() for h in z], ddof=1)
    return math.sqrt(((n - 1) / n * w + b / n) / w)


def test_split_rhat_matches_reference_and_separates_cases():
    rng = np.random.default_rng(11)
    mixed = rng.standard_normal((16, 1536))
    stuck = mixed + np.where(rng.random(16) < 0.5, -3.0, 3.0)[:, None]
    got = convergence.split_rhat(np.stack([mixed, stuck]))
    assert math.isclose(got[0], _rhat_reference(mixed), rel_tol=1e-12)
    assert math.isclose(got[1], _rhat_reference(stuck), rel_tol=1e-12)
    assert got[0] < 1.01 and got[1] > 1.5


def test_split_rhat_all_tied_is_nan():
    assert np.isnan(convergence.split_rhat(np.ones((4, 16)))[0])


def test_split_rhat_separated_constant_chains_is_infinite():
    # Codex review on PR #2: zero within-chain variance with distinct chain values is total non-mixing.
    assert np.isposinf(convergence.split_rhat(np.array([[-1.0] * 4, [1.0] * 4]))[0])


def test_split_rhat_is_invariant_to_sign_reversal():
    # Codex review on PR #2: Blom scores are antisymmetric in rank, so R-hat(x) == R-hat(-x).
    rng = np.random.default_rng(1)
    x = rng.standard_normal((4, 20)) + np.array([0, 0.3, 0.6, 0.9])[:, None]
    assert math.isclose(convergence.split_rhat(x)[0], convergence.split_rhat(-x)[0], rel_tol=1e-12)
    scores = norm.ppf((np.arange(1, 9) - 0.375) / (8 + 0.25))
    assert np.allclose(scores, -scores[::-1])


def test_split_rhat_odd_length_drops_middle():
    rng = np.random.default_rng(5)
    x = rng.standard_normal((4, 101))
    assert math.isclose(convergence.split_rhat(x)[0], convergence.split_rhat(np.delete(x, 50, axis=1))[0])


def test_plugin_symmetric_and_abs_definitions():
    m = np.array([[0.5, 0.5, -0.5, -0.5] * 8, [1.0] * 32])
    N, T = 16, 2.0
    assert np.allclose(susceptibility.plugin(m, N, T), [N / T * 0.25, 0.0])
    assert np.allclose(susceptibility.abs_plugin(m, N, T), [0.0, 0.0])
    s = susceptibility.symmetric(m, N, T)
    assert np.allclose(s["estimate"], [N / T * 0.25, N / T])
    assert s["batch"] == 2 and s["batches"] == 16


def test_ensemble_jackknife_matches_explicit_loop():
    rng = np.random.default_rng(9)
    m = np.clip(rng.normal(0, 0.6, (32, 64)) + np.where(rng.random(32) < 0.5, -0.4, 0.4)[:, None], -1, 1)
    N, T, K = 16, 1.5, 16
    got = susceptibility.ensemble_jackknife(m, N, T, K)
    for e in range(2):
        block = m[e * K:(e + 1) * K]
        A, B = (block ** 2).mean(axis=1), block.mean(axis=1)
        full = N / T * (A.mean() - B.mean() ** 2)
        loo = np.array([N / T * (np.delete(A, i).mean() - np.delete(B, i).mean() ** 2) for i in range(K)])
        corrected = K * full - (K - 1) * loo.mean()
        var = (K - 1) / K * np.sum((loo - loo.mean()) ** 2)
        assert math.isclose(got["estimate"][e], corrected, rel_tol=1e-12)
        assert math.isclose(got["variance"][e], var, rel_tol=1e-12)
    with pytest.raises(ValueError):
        susceptibility.ensemble_jackknife(m[:30], N, T, K)


def test_exact_references_from_enumeration():
    ref = susceptibility.exact_references(4, 1.5)
    exact = enumeration(4, 1.5)
    assert math.isclose(ref["chi"], 10.403226593514693, rel_tol=1e-12)
    expected_abs = 16 / 1.5 * (exact["susceptibility"] * 1.5 / 16 - exact["abs_magnetization"] ** 2)
    assert math.isclose(ref["chi_abs"], expected_abs, rel_tol=1e-12)
    assert ref["chi_abs"] < ref["chi"]
