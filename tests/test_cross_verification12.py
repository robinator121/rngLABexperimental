"""Alternative mathematical cross-verification checks for estimators."""
from itertools import product

import numpy as np
import pytest

from rnglab12.uncertainty import summarize_series


def test_geyer_matches_direct_odd_length_autocovariance():
    values = np.random.Generator(np.random.MT19937(933)).normal(size=513)
    centered = values - values.mean()
    n = len(values)
    # Direct inner products avoid the production FFT implementation entirely.
    covariances = np.array([np.dot(centered[:n-k], centered[k:]) / n for k in range(n)])
    positive = []
    for k in range(n // 2):
        pair = covariances[2*k] + covariances[2*k+1]
        if pair <= 0:
            break
        positive.append(pair)
    running_min = float('inf')
    monotone_sum = 0.
    for pair in positive:
        running_min = min(running_min, pair)
        monotone_sum += running_min
    variance = -covariances[0] + 2 * monotone_sum
    actual = summarize_series(values)['geyer']
    assert actual['long_run_variance'] == pytest.approx(variance, abs=1e-14)
    assert actual['tau'] == pytest.approx(variance / covariances[0], abs=1e-14)


@pytest.mark.parametrize('n,start,m,q', [(8, 2, 4, 2), (8, 1, 6, 2), (8, 2, 4, 4), (8, 0, 8, 2)])
def test_periodic_weighted_variance_and_fourth_moment_by_full_enumeration(n, start, m, q):
    # Independent symmetric two-point variables have variance=1, fourth moment=1.
    # Enumerating their full finite probability space is exact, not a simulation.
    native = np.array(list(product((-1., 1.), repeat=n)))
    faulty = native.copy()
    faulty[:, start:start+m] = np.tile(native[:, start:start+q], (1, m // q))
    means = faulty.mean(axis=1)
    theory_variance = (n-m) / n**2 + m**2 / (q * n**2)
    sum_weight_fourth = (n-m) / n**4 + q * (m / (q*n))**4
    theory_fourth = 3 * theory_variance**2 - 2 * sum_weight_fourth
    assert means.mean() == 0.
    assert np.mean(means**2) == pytest.approx(theory_variance, abs=1e-15)
    assert np.mean(means**4) == pytest.approx(theory_fourth, abs=1e-15)
