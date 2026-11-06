"""RNGLab v0.18: zero-field Ising susceptibility estimators.

Estimand: chi = (N/T) Var_pi(m) = Var_pi(M)/(N T), with m = M/N the
magnetization per spin of an N-spin periodic lattice at zero field (the value
rnglab16.core.enumeration returns). Definitions are frozen in
research/v018-estimators/protocol.json:

- plugin (P): per chain, (N/T)(mean m^2 - (mean m)^2). Fails when a chain does
  not visit both signs of m; kept to reproduce that failure.
- symmetric (S): per chain, (N/T) mean m^2. Uses E_pi[m] = 0, which holds
  exactly here by global spin-flip symmetry; model-specific, not general.
  Interval: batch means with 16 batches on the series (N/T) m_t^2, t_15.
- ensemble_jackknife (G): K independent chains, pooled
  (N/T)(mean_i A_i - (mean_i B_i)^2) with A_i, B_i the chain means of m^2 and
  m; delete-one-chain jackknife bias correction and standard error, t_{K-1}.
  The formula does not use E[m] = 0, but it is valid only if the ensemble is
  distributed as pi; with flip-symmetric initial states and kernel this holds
  exactly for the sign of m, and not in general.
- abs_plugin: (N/T)(mean m^2 - (mean |m|)^2), a different estimand reported only
  for comparison.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import t as student_t

from .longrun import batch_means, interval

S_BATCHES = 16


def _m(m) -> np.ndarray:
    a = np.asarray(m, dtype=float)
    if a.ndim == 1:
        a = a[None]
    if a.ndim != 2 or a.shape[1] < S_BATCHES * 2 or not np.isfinite(a).all():
        raise ValueError("finite magnetization per spin with shape (chains, sweeps) required")
    return a


def plugin(m, N: int, T: float) -> np.ndarray:
    a = _m(m)
    return N / T * ((a * a).mean(axis=1) - a.mean(axis=1) ** 2)


def abs_plugin(m, N: int, T: float) -> np.ndarray:
    a = _m(m)
    return N / T * ((a * a).mean(axis=1) - np.abs(a).mean(axis=1) ** 2)


def symmetric(m, N: int, T: float, confidence: float = 0.95) -> dict:
    a = _m(m)
    series = N / T * a * a
    b = series.shape[1] // S_BATCHES
    bm = batch_means(series, b)
    lo, hi = interval(bm["mean"], bm["sigma2"], bm["n_used"], bm["df"], confidence)
    return {"estimate": bm["mean"], "variance": bm["sigma2"] / bm["n_used"], "lo": lo, "hi": hi,
            "batch": b, "batches": S_BATCHES, "series": series}


def ensemble_jackknife(m, N: int, T: float, K: int = 16, confidence: float = 0.95) -> dict:
    """Groups consecutive chains into ensembles of K; returns one estimate per ensemble."""
    a = _m(m)
    if type(K) is not int or K < 2 or len(a) % K:
        raise ValueError("chain count must be a positive multiple of K >= 2")
    A = (a * a).mean(axis=1).reshape(-1, K)
    B = a.mean(axis=1).reshape(-1, K)
    full = N / T * (A.mean(axis=1) - B.mean(axis=1) ** 2)
    loo_A = (A.sum(axis=1, keepdims=True) - A) / (K - 1)
    loo_B = (B.sum(axis=1, keepdims=True) - B) / (K - 1)
    leave_out = N / T * (loo_A - loo_B ** 2)
    loo_mean = leave_out.mean(axis=1)
    corrected = K * full - (K - 1) * loo_mean
    variance = (K - 1) / K * ((leave_out - loo_mean[:, None]) ** 2).sum(axis=1)
    q = student_t.ppf((1 + confidence) / 2, K - 1)
    half = q * np.sqrt(variance)
    return {"estimate": corrected, "uncorrected": full, "variance": variance,
            "lo": corrected - half, "hi": corrected + half, "K": K}


def exact_references(L: int, T: float) -> dict:
    """Exact finite-lattice values from rnglab16.core.enumeration (L <= 4)."""
    from rnglab16.core import enumeration

    exact = enumeration(L, T)
    N = L * L
    chi = exact["susceptibility"]
    return {"chi": chi, "chi_abs": chi - N * exact["abs_magnetization"] ** 2 / T,
            "energy_per_spin": exact["energy"], "N": N, "T": T, "L": L}
