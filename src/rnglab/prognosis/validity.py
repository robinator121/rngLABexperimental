"""Probability that a batch-means interval covers the reference when the mean drifts.

Model. The estimate at run length L is the mean of B batch means m_1..m_B of
equal size b = L/B. Under a severity path s(t) with bias slope beta and
long-run noise scale sigma_eff (per-sweep units, so that the standard error of
a mean over n sweeps is sigma_eff / sqrt(n)), the batch means are taken as
independent normal variables with

    E[m_j] = mu + beta * sbar_j,     Var[m_j] = sigma_eff^2 / b,

sbar_j being the mean severity inside batch j. Then, with nu = B - 1,

    mean of the m_j  ~  N(mu + beta * Sbar, sigma_eff^2 / L),
    nu * S^2 / (sigma_eff^2 / b)  ~  noncentral chi-square(nu, lambda),
    lambda = (b / sigma_eff^2) * beta^2 * sum_j (sbar_j - Sbar)^2,

and the two are independent. The interval mean +- c * S / sqrt(B) covers mu
with probability

    P(delta, lambda) = E_V[ Phi(c sqrt(V/nu) - delta) - Phi(-c sqrt(V/nu) - delta) ],
    V ~ noncentral chi-square(nu, lambda),   delta = beta * Sbar * sqrt(L) / sigma_eff.

lambda = 0 gives the noncentral-t formula that the Stage D pilot uses as its
oracle reference. lambda > 0 is the widening of the interval by drift between
batches, which that formula leaves out: the statistic is then doubly
noncentral. The normality and independence of batch means are assumptions of
this model, to be checked against data wherever it is applied.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy import special, stats

from .families import SeverityPath, block_means

FloatArray = NDArray[np.float64]
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(160)


def coverage_probability(delta: ArrayLike, lam: ArrayLike = 0.0, batches: int = 16, level: float = 0.95) -> FloatArray:
    """P(delta, lambda) by Gauss-Legendre quadrature over the noncentral chi-square mixing variable.

    The integration interval is the mixing variable's mean plus and minus a multiple of
    its standard deviation, widened towards the region where the integrand lives when
    |delta| is large. Absolute accuracy is about 1e-9 (tests compare with adaptive
    quadrature and simulation).
    """
    delta, lam = np.broadcast_arrays(np.abs(np.asarray(delta, dtype=np.float64)), np.asarray(lam, dtype=np.float64))
    nu = batches - 1
    c = float(special.stdtrit(nu, (1 + level) / 2))
    mean, sd = nu + lam, np.sqrt(2 * (nu + 2 * lam))
    low = np.maximum(mean - 12 * sd, 0.0)
    high = np.maximum(mean + 14 * sd, nu * ((delta + 9.0) / c)**2 * (delta > 4))  # reach c*sqrt(v/nu) ~ delta
    half, centre = (high - low) / 2, (high + low) / 2
    v = centre[..., None] + half[..., None] * _NODES
    density = np.where(lam[..., None] > 0, stats.ncx2.pdf(v, nu, np.maximum(lam[..., None], 1e-300)), stats.chi2.pdf(v, nu))
    root = c * np.sqrt(v / nu)
    inside = special.ndtr(root - delta[..., None]) - special.ndtr(-root - delta[..., None])
    return np.clip(np.sum(_WEIGHTS * inside * density, axis=-1) * half, 0.0, 1.0)


@dataclass(frozen=True)
class ValidityModel:
    """Maps a severity path to the probability that the interval at run length L covers the reference."""
    beta: float
    sigma_eff: float
    batches: int = 16
    level: float = 0.95

    def noncentralities(self, path: SeverityPath, length: int) -> tuple[FloatArray, FloatArray]:
        """(delta, lambda) of the path at run length `length`."""
        size = length // self.batches
        means = block_means(path, np.arange(self.batches + 1) * size)
        overall = means.mean(axis=-1)
        delta = self.beta * overall * np.sqrt(length) / self.sigma_eff
        lam = size / self.sigma_eff**2 * self.beta**2 * np.sum((means - overall[..., None])**2, axis=-1)
        return np.asarray(delta), np.asarray(lam)

    def probability(self, path: SeverityPath, length: int, *, drift_widening: bool = True) -> FloatArray:
        delta, lam = self.noncentralities(path, length)
        return coverage_probability(delta, lam if drift_widening else 0.0, self.batches, self.level)


@lru_cache(maxsize=8)
def _table(batches: int, level: float, delta_max: float, lam_max: float, points: int) -> tuple[FloatArray, FloatArray, FloatArray]:
    delta = np.linspace(0.0, delta_max, points)
    lam = np.linspace(0.0, np.sqrt(lam_max), points)**2  # denser near zero
    return delta, lam, coverage_probability(delta[:, None], lam[None, :], batches, level)


def coverage_probability_fast(delta: ArrayLike, lam: ArrayLike, batches: int = 16, level: float = 0.95,
                              delta_max: float = 14.0, lam_max: float = 2500.0, points: int = 281) -> FloatArray:
    """Bilinear interpolation of a precomputed table; exact evaluation outside the table.

    For forecasters that evaluate many hypotheses. The table error is below 5e-4 in
    probability (asserted in the tests). Arguments outside the table are evaluated exactly:
    a large lambda widens the interval enough to cover even a large delta.
    """
    delta, lam = np.broadcast_arrays(np.abs(np.asarray(delta, dtype=np.float64)), np.asarray(lam, dtype=np.float64))
    grid_d, grid_l, table = _table(batches, level, delta_max, lam_max, points)
    out = np.empty(delta.shape)
    outside = (lam > lam_max) | (delta > delta_max)
    inside = ~outside
    d = delta[inside]
    i = np.clip(np.searchsorted(grid_d, d, side="right") - 1, 0, points - 2)
    root = np.sqrt(lam[inside])
    step = np.sqrt(lam_max) / (points - 1)
    j = np.clip((root / step).astype(np.int64), 0, points - 2)
    wd = (d - grid_d[i]) / (grid_d[i + 1] - grid_d[i])
    wl = (lam[inside] - grid_l[j]) / (grid_l[j + 1] - grid_l[j])
    value = ((1 - wd) * (1 - wl) * table[i, j] + wd * (1 - wl) * table[i + 1, j]
             + (1 - wd) * wl * table[i, j + 1] + wd * wl * table[i + 1, j + 1])
    out[inside] = value
    if outside.any():
        out[outside] = coverage_probability(delta[outside], lam[outside], batches, level)
    return out
