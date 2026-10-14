"""Finite-series uncertainty diagnostics; none establishes stationarity or mixing.

The Geyer initial-positive/initial-monotone sequence estimator is appropriate
under additional Markov-chain regularity (notably reversibility) assumptions.
A numeric ESS is a diagnostic, never a certificate for a corrupted trajectory.
"""
from __future__ import annotations

import math
from collections.abc import Sequence

import numpy as np
from scipy.stats import norm, t


_SCOPE = ('Approximate uncertainty only. Stationarity, mixing, reversibility for the '
          'Geyer sequence estimator, and batch-size adequacy are not established '
          'by these summaries; no uncertainty method repairs target bias.')


def _interval(mean: float, mcse: float, critical: float, reference: float | None) -> dict:
    bounds = [mean - critical * mcse, mean + critical * mcse]
    if not all(math.isfinite(value) for value in bounds):
        raise ValueError('Uncertainty interval overflowed')
    return {'interval_95': bounds,
            'covered': None if reference is None else bounds[0] <= reference <= bounds[1]}


def _geyer(values: np.ndarray, reference: float | None) -> dict:
    n = len(values)
    centered = values - values.mean()
    fft_n = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=fft_n)
    autocov = np.fft.irfft(spectrum * spectrum.conjugate(), n=fft_n)[:n] / n
    # Biased FFT autocovariances use a common n denominator. Pairing begins
    # with gamma_0 + gamma_1, as in Geyer's reversible-chain construction.
    pairs = autocov[:2 * (n // 2)].reshape(-1, 2).sum(axis=1)
    stop = next((j for j, value in enumerate(pairs) if value <= 0), len(pairs))
    positive = pairs[:stop]
    if not len(positive):
        return {'status': 'INVALID_LONG_RUN_VARIANCE', 'mcse': None, 'ess': None,
                'tau': None, 'interval_95': None, 'covered': None, 'positive_pairs': 0}
    monotone = np.minimum.accumulate(positive)
    lrv = float(-autocov[0] + 2 * monotone.sum())
    if not math.isfinite(lrv) or lrv <= 0 or not np.isfinite(autocov).all():
        return {'status': 'INVALID_LONG_RUN_VARIANCE', 'mcse': None, 'ess': None,
                'tau': None, 'interval_95': None, 'covered': None,
                'positive_pairs': len(positive)}
    tau = lrv / float(autocov[0])
    mcse = math.sqrt(lrv / n)
    return {'status': 'ESTIMATED', 'method': 'FFT Geyer initial-positive then initial-monotone pairs',
            'mcse': mcse, 'ess': n / tau, 'tau': tau, 'long_run_variance': lrv,
            'positive_pairs': len(positive), 'n': n,
            **_interval(float(values.mean()), mcse, float(norm.ppf(.975)), reference)}


def _batch(values: np.ndarray, size: int, reference: float | None) -> dict:
    count = len(values) // size
    used = count * size
    result = {'batch_size': size, 'batches': count, 'prefix_n': used,
              'discarded_tail': len(values) - used, 'mcse': None, 'ess': None,
              'interval_95': None, 'covered': None, 'prefix_mean': None}
    if count < 8:
        return {**result, 'status': 'INSUFFICIENT_BATCHES', 'minimum_batches': 8}
    prefix = values[:used]
    means = prefix.reshape(count, size).mean(axis=1)
    estimate = float(means.mean())
    variance = float(means.var(ddof=1))
    result['prefix_mean'] = estimate
    if variance <= 0 or not math.isfinite(variance):
        return {**result, 'status': 'DEGENERATE_BATCH_VARIANCE'}
    mcse = math.sqrt(variance / count)
    return {**result, 'status': 'ESTIMATED', 'method': 'nonoverlapping batch means',
            'mcse': mcse, 'ess': float(prefix.var(ddof=1)) / mcse**2,
            **_interval(estimate, mcse, float(t.ppf(.975, count - 1)), reference)}


def summarize_series(series: Sequence[float], reference: float | None = None,
                     batch_sizes: Sequence[int] = (32, 128)) -> dict:
    """Summarize one full scalar trace after a separately fixed burn-in rule.

    Preserve fixed batch size 32 as a historical baseline. Larger batches and
    the Geyer estimator are comparisons, not an automatically selected winner.
    Do not tune a method using a final holdout coverage outcome.
    """
    values = np.asarray(series, dtype=float)
    if values.ndim != 1 or len(values) < 32 or not np.isfinite(values).all():
        raise ValueError('A finite one-dimensional full trace with at least 32 values is required')
    if reference is not None and (isinstance(reference, bool) or not math.isfinite(reference)):
        raise ValueError('Reference must be finite or None')
    sizes = tuple(batch_sizes)
    if not sizes or len(set(sizes)) != len(sizes) or any(type(size) is not int or size < 2 for size in sizes):
        raise ValueError('Unique integer batch sizes >= 2 are required')
    with np.errstate(over='raise', invalid='raise'):
        try:
            mean = float(values.mean())
            variance = float(values.var(ddof=1))
        except FloatingPointError as error:
            raise ValueError('Trace moments overflowed') from error
    result = {'n': len(values), 'mean': mean, 'sample_variance': variance,
              'naive_mcse': math.sqrt(variance / len(values)),
              'reference': reference,
              'reference_status': 'KNOWN' if reference is not None else 'ABSTAIN_UNKNOWN_REFERENCE',
              'stationarity_verified': False, 'scope': _SCOPE}
    if variance == 0:
        return {**result, 'status': 'ABSTAIN_CONSTANT_TRACE',
                'geyer': {'status': 'ABSTAIN_CONSTANT_TRACE', 'mcse': None, 'ess': None,
                          'tau': None, 'interval_95': None, 'covered': None},
                'batch_means': {str(size): _batch(values, size, reference) for size in sizes}}
    with np.errstate(over='raise', invalid='raise', divide='raise'):
        try:
            geyer = _geyer(values, reference)
            batches = {str(size): _batch(values, size, reference) for size in sizes}
        except FloatingPointError as error:
            raise ValueError('Uncertainty calculation overflowed or became invalid') from error
    return {**result, 'status': 'ESTIMATED', 'geyer': geyer, 'batch_means': batches}
