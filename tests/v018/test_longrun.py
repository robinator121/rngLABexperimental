"""v0.18 long-run variance estimators against brute-force definitions and existing code."""
import math

import numpy as np
import pytest

from rnglab import longrun
from rnglab12.uncertainty import _geyer, summarize_series
from rnglab17.endpoints import batch_intervals


@pytest.fixture(scope="module")
def series():
    rng = np.random.default_rng(20261105)
    return np.cumsum(rng.standard_normal((3, 1024)), axis=1) * 0.05 + rng.standard_normal((3, 1024))


def test_batch_means_matches_existing_implementations(series):
    bm = longrun.batch_means(series, 64)
    old = batch_intervals(series, 0.0, 64)
    assert np.allclose(np.sqrt(bm["sigma2"] / bm["n_used"]), old["se"], rtol=1e-12)
    for row, sigma2 in zip(series, bm["sigma2"]):
        mcse = summarize_series(row, batch_sizes=(64,))["batch_means"]["64"]["mcse"]
        assert math.isclose(math.sqrt(sigma2 / 1024), mcse, rel_tol=1e-12)
    assert np.all(bm["df"] == 15)


def test_obm_matches_definition(series):
    b = 50
    got = longrun.overlapping_batch_means(series, b)["sigma2"]
    for row, value in zip(series, got):
        n = len(row)
        means = np.array([row[j:j + b].mean() for j in range(n - b + 1)])
        expected = n * b / ((n - b) * (n - b + 1)) * np.sum((means - row.mean()) ** 2)
        assert math.isclose(value, expected, rel_tol=1e-10)


@pytest.mark.parametrize("window", ["bartlett", "tukey-hanning"])
def test_spectral_variance_matches_definition(series, window):
    b = 40
    got = longrun.spectral_variance(series, b, window)
    for row, value in zip(series, got["sigma2"]):
        n = len(row)
        c = row - row.mean()
        gamma = np.array([np.dot(c[:n - k], c[k:]) / n for k in range(b)])
        k = np.arange(1, b)
        w = 1 - k / b if window == "bartlett" else (1 + np.cos(np.pi * k / b)) / 2
        assert math.isclose(value, gamma[0] + 2 * np.sum(w * gamma[1:]), rel_tol=1e-10)
    integral = 2 / 3 if window == "bartlett" else 0.75
    assert np.allclose(got["df"], 1024 / (b * integral))


def test_geyer_wraps_rnglab12_unchanged(series):
    got = longrun.geyer(series)
    for i, row in enumerate(series):
        old = _geyer(row, None)
        assert math.isclose(got["sigma2"][i], old["long_run_variance"], rel_tol=1e-15)
        assert math.isclose(got["tau"][i], old["tau"], rel_tol=1e-15)


def test_geyer_invalid_row_is_nan():
    alternating = np.tile([1.0, -1.0], 64)
    got = longrun.geyer(alternating)
    assert not got["valid"][0] and np.isnan(got["sigma2"][0])
    assert longrun.adaptive_status(128, got["tau"])[0] == "TAU_UNDEFINED"


def test_ar1_exact_variance_matches_double_sum():
    for phi in (0.0, 0.5, 0.9, -0.3):
        n = 37
        lags = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
        brute = (phi ** lags).sum() / n ** 2
        assert math.isclose(longrun.ar1_exact_var_mean(phi, n), brute, rel_tol=1e-12)


def test_ar1_series_is_stationary_with_unit_variance():
    from rnglab.seeds import replica_generators

    x = longrun.ar1_series(replica_generators(7, 2000), 64, 0.9)
    assert abs(x[:, 0].var() - 1) < 0.1 and abs(x[:, -1].var() - 1) < 0.1
    lag1 = np.mean(x[:, 1:] * x[:, :-1])
    assert abs(lag1 - 0.9) < 0.05


def test_interval_and_coverage_semantics():
    lo, hi = longrun.interval([0.0, 0.0, 0.0], [1.0, 0.0, np.nan], 100, [np.inf, 10.0, 10.0])
    assert math.isclose(hi[0], 1.959963984540054 * 0.1, rel_tol=1e-12)
    assert np.isnan(lo[1]) and np.isnan(hi[2])
    assert list(longrun.covers(lo, hi, 0.0)) == [True, False, False]


def test_diagnostic_rules():
    n = 8192
    assert list(longrun.adaptive_batch(n, [1.0, 20.0, 51.2, 700.0, np.nan])) == [16, 256, 512, 0, 0]
    assert list(longrun.adaptive_status(n, [1.0, 700.0])) == ["ADEQUATE", "RUN_TOO_SHORT"]
    assert list(longrun.batch_length_status(512, n, [51.2, 51.3])) == ["ADEQUATE", "WARN_BATCH_SHORT"]
    assert longrun.batch_length_status(1024, n, [1.0])[0] == "WARN_BATCH_SHORT"  # only 8 batches


def test_estimators_agree_on_iid_noise():
    rng = np.random.default_rng(3)
    x = rng.standard_normal((200, 4096))
    for method in ("BM", "OBM", "SV-B", "SV-TH"):
        assert abs(longrun.estimate(x, method, 64)["sigma2"].mean() - 1) < 0.03
    assert abs(np.nanmean(longrun.geyer(x)["sigma2"]) - 1) < 0.05


def test_rejects_bad_input():
    with pytest.raises(ValueError):
        longrun.batch_means(np.ones(10), 8)
    with pytest.raises(ValueError):
        longrun.spectral_variance(np.ones((2, 64)), 8, "parzen")
    with pytest.raises(ValueError):
        longrun.batch_means([[1.0, np.nan, 2.0, 3.0]], 1)
