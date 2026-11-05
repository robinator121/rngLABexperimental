"""RNGLab v0.18: long-run variance estimators and a batch-length diagnostic.

All functions take x with shape (n,) or (replicas, n) and return arrays with one
value per row. The target is sigma2_n = n Var(mean of x); an interval for the
mean is mean +- q * sqrt(sigma2_hat / n) with q from Student t (df given) or the
normal law (df = inf).

Estimators (definitions frozen in research/v018-estimators/protocol.json):

- BM(b): non-overlapping batch means, a = floor(n/b) batches over the prefix
  of length a*b; df = a - 1.
- OBM(b): overlapping batch means (Meketon & Schmeiser 1984); df = 1.5 n/b.
- SV(b, window): lag-window spectral variance at frequency zero with the
  Bartlett (1 - k/b) or Tukey-Hanning ((1 + cos(pi k/b))/2) window and biased
  FFT autocovariances; df = n / (b * integral of w^2), with integrals 2/3 and
  3/4.
- GEYER: rnglab12.uncertainty._geyer, initial positive then initial monotone
  sequence (Geyer 1992), applied row by row; normal critical value.

Consistency conditions for BM, OBM and SV are in Flegal & Jones (2010). None of
these estimators detects target bias or nonstationarity, and a numeric ESS is a
diagnostic, not a certificate.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.stats import norm, t as student_t

from rnglab12.uncertainty import _geyer

WINDOW_SQUARE_INTEGRAL = {"bartlett": 2.0 / 3.0, "tukey-hanning": 0.75}
ADEQUATE_BATCH_FACTOR = 10
MINIMUM_BATCHES = 16
MINIMUM_ADAPTIVE_BATCH = 16


def _rows(x) -> np.ndarray:
    a = np.asarray(x, dtype=float)
    if a.ndim == 1:
        a = a[None, :]
    if a.ndim != 2 or a.shape[1] < 2 or not np.isfinite(a).all():
        raise ValueError("finite series with shape (n,) or (replicas, n), n >= 2, required")
    return a


def _check_b(b: int, n: int, minimum_batches: int = 2) -> None:
    if type(b) is not int or b < 1 or n // b < minimum_batches:
        raise ValueError(f"batch length {b!r} needs an integer >= 1 and at least {minimum_batches} batches in n={n}")


def autocovariance(x, max_lag: int) -> np.ndarray:
    """Biased (divide by n) autocovariances at lags 0..max_lag, per row, by FFT."""
    a = _rows(x)
    n = a.shape[1]
    if type(max_lag) is not int or not 0 <= max_lag < n:
        raise ValueError("0 <= max_lag < n required")
    centered = a - a.mean(axis=1, keepdims=True)
    size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=size, axis=1)
    return np.fft.irfft(spectrum * spectrum.conj(), n=size, axis=1)[:, :max_lag + 1] / n


def batch_means(x, b: int) -> dict:
    a = _rows(x)
    n = a.shape[1]
    _check_b(b, n)
    count = n // b
    means = a[:, :count * b].reshape(len(a), count, b).mean(axis=2)
    sigma2 = b * means.var(axis=1, ddof=1)
    return {"sigma2": sigma2, "df": np.full(len(a), float(count - 1)), "n_used": count * b,
            "mean": means.mean(axis=1)}


def overlapping_batch_means(x, b: int) -> dict:
    a = _rows(x)
    n = a.shape[1]
    _check_b(b, n)
    total = np.concatenate([np.zeros((len(a), 1)), np.cumsum(a, axis=1)], axis=1)
    window_means = (total[:, b:] - total[:, :-b]) / b
    xbar = a.mean(axis=1, keepdims=True)
    sigma2 = n * b / ((n - b) * (n - b + 1)) * ((window_means - xbar) ** 2).sum(axis=1)
    return {"sigma2": sigma2, "df": np.full(len(a), 1.5 * n / b), "n_used": n, "mean": xbar[:, 0]}


def spectral_variance(x, b: int, window: str) -> dict:
    a = _rows(x)
    n = a.shape[1]
    _check_b(b, n)
    if window not in WINDOW_SQUARE_INTEGRAL:
        raise ValueError(f"window must be one of {sorted(WINDOW_SQUARE_INTEGRAL)}")
    gamma = autocovariance(a, b - 1)
    k = np.arange(1, b)
    weights = 1 - k / b if window == "bartlett" else (1 + np.cos(np.pi * k / b)) / 2
    sigma2 = gamma[:, 0] + 2 * (gamma[:, 1:] * weights).sum(axis=1)
    df = n / (b * WINDOW_SQUARE_INTEGRAL[window])
    return {"sigma2": sigma2, "df": np.full(len(a), df), "n_used": n, "mean": a.mean(axis=1)}


def geyer(x) -> dict:
    """Row-wise rnglab12 Geyer estimator; invalid rows give NaN sigma2 and tau."""
    a = _rows(x)
    sigma2 = np.full(len(a), np.nan)
    tau = np.full(len(a), np.nan)
    valid = np.zeros(len(a), dtype=bool)
    for i, row in enumerate(a):
        if row.var() == 0:
            continue
        result = _geyer(row, None)
        if result["status"] == "ESTIMATED":
            sigma2[i] = result["long_run_variance"]
            tau[i] = result["tau"]
            valid[i] = True
    return {"sigma2": sigma2, "tau": tau, "valid": valid, "df": np.full(len(a), np.inf),
            "n_used": a.shape[1], "mean": a.mean(axis=1)}


def estimate(x, method: str, b: int | None = None) -> dict:
    """Dispatch by protocol method ID: BM, OBM, SV-B, SV-TH, GEYER."""
    if method == "BM":
        return batch_means(x, b)
    if method == "OBM":
        return overlapping_batch_means(x, b)
    if method == "SV-B":
        return spectral_variance(x, b, "bartlett")
    if method == "SV-TH":
        return spectral_variance(x, b, "tukey-hanning")
    if method == "GEYER":
        return geyer(x)
    raise ValueError(f"unknown method {method!r}")


def interval(mean, sigma2, n_used: int, df, confidence: float = 0.95) -> tuple[np.ndarray, np.ndarray]:
    """Interval for the mean; NaN or non-positive sigma2 gives NaN bounds."""
    mean = np.asarray(mean, dtype=float)
    sigma2 = np.asarray(sigma2, dtype=float)
    df = np.broadcast_to(np.asarray(df, dtype=float), mean.shape)
    q = np.where(np.isinf(df), norm.ppf((1 + confidence) / 2),
                 student_t.ppf((1 + confidence) / 2, np.where(np.isinf(df), 1.0, df)))
    with np.errstate(invalid="ignore"):
        half = q * np.sqrt(np.where(sigma2 > 0, sigma2, np.nan) / n_used)
    return mean - half, mean + half


def covers(lo, hi, target: float) -> np.ndarray:
    """Coverage indicator; NaN bounds (numerical failure) count as not covering."""
    lo = np.asarray(lo)
    hi = np.asarray(hi)
    return np.isfinite(lo) & np.isfinite(hi) & (lo <= target) & (target <= hi)


def batch_length_status(b: int, n: int, tau_hat) -> np.ndarray:
    """ADEQUATE if b >= 10 tau_hat and n/b >= 16; WARN_BATCH_SHORT otherwise; TAU_UNDEFINED for NaN tau."""
    tau_hat = np.asarray(tau_hat, dtype=float)
    adequate = (b >= ADEQUATE_BATCH_FACTOR * tau_hat) & (n // b >= MINIMUM_BATCHES)
    return np.where(np.isnan(tau_hat), "TAU_UNDEFINED", np.where(adequate, "ADEQUATE", "WARN_BATCH_SHORT"))


def adaptive_batch(n: int, tau_hat) -> np.ndarray:
    """Smallest power of two in [16, n/16] with b >= 10 tau_hat; 0 means abstain."""
    tau_hat = np.asarray(tau_hat, dtype=float)
    candidates = [MINIMUM_ADAPTIVE_BATCH << k for k in range(64)
                  if (MINIMUM_ADAPTIVE_BATCH << k) <= n // MINIMUM_BATCHES]
    chosen = np.zeros(tau_hat.shape, dtype=np.int64)
    for b in reversed(candidates):
        chosen = np.where(b >= ADEQUATE_BATCH_FACTOR * tau_hat, b, chosen)
    return np.where(np.isnan(tau_hat), 0, chosen)


def adaptive_status(n: int, tau_hat) -> np.ndarray:
    tau_hat = np.asarray(tau_hat, dtype=float)
    chosen = adaptive_batch(n, tau_hat)
    return np.where(np.isnan(tau_hat), "TAU_UNDEFINED", np.where(chosen > 0, "ADEQUATE", "RUN_TOO_SHORT"))


def ar1_exact_var_mean(phi: float, n: int) -> float:
    """Var of the mean of a stationary unit-variance AR(1) of length n."""
    if not -1 < phi < 1 or type(n) is not int or n < 1:
        raise ValueError("|phi| < 1 and integer n >= 1 required")
    if phi == 0:
        return 1.0 / n
    return ((1 + phi) / (1 - phi) - 2 * phi * (1 - phi ** n) / (n * (1 - phi) ** 2)) / n


def ar1_series(generators, n: int, phi: float) -> np.ndarray:
    """Stationary unit-variance AR(1) rows, one per generator: x1 ~ N(0,1),
    x_t = phi x_{t-1} + sqrt(1 - phi^2) eps_t."""
    from scipy.signal import lfilter

    if not -1 < phi < 1 or type(n) is not int or n < 2:
        raise ValueError("|phi| < 1 and integer n >= 2 required")
    eps = np.stack([g.standard_normal(n) for g in generators])
    scale = math.sqrt(1 - phi * phi)
    eps[:, 0] /= scale
    return lfilter([scale], [1.0, -phi], eps, axis=1)
