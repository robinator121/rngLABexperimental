"""RNGLab v0.20 nuclear module, milestone B2: exact moments of population counts.

Draws no random numbers. For N nuclei started in one nuclide, each nucleus has
a categorical outcome: its state at time t, or the pair of its states at t1 and
t2 (Markov property: P(j at t1, k at t2) = p_j(t1) * P_kj(t2 - t1)). A count is a
sum of N independent indicator values, so its moments follow exactly from the
categorical law, with Y = f - E f and Z = g - E g:

    sigma_ab = N E[YZ],
    mu_22    = N E[Y^2 Z^2] + N (N - 1) (E[Y^2] E[Z^2] + 2 E[YZ]^2).

Every moment is computed from centred values: the uncentred form
E[fg] - E[f] E[g] cancels catastrophically for rare outcomes (review finding:
it gave a negative variance at N = 2**62, p = 2**-53). N is taken as a Python
integer, so N (N - 1) cannot overflow (a NumPy int64 count of 1e12 did).

and, for R independent replicates, the unbiased sample covariance s_ab has

    Var(s_ab) = (mu_22 - sigma_ab^2) / R + (sigma_aa sigma_bb + sigma_ab^2) / (R (R - 1)),

which is exact for any distribution with finite fourth moments.
"""
from __future__ import annotations

import operator
from dataclasses import dataclass

import numpy as np

from .bateman import state_probabilities
from .decay import DecayScheme


@dataclass(frozen=True)
class PairMoments:
    sigma_aa: float
    sigma_bb: float
    sigma_ab: float
    mu_22: float

    def sample_covariance_variance(self, replicates: int) -> float:
        r = _count(replicates, "replicates")
        if r < 2:
            raise ValueError("the sample covariance needs at least two replicates")
        return ((self.mu_22 - self.sigma_ab ** 2) / r
                + (self.sigma_aa * self.sigma_bb + self.sigma_ab ** 2) / (r * (r - 1)))


def _count(value, name: str) -> int:
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer, not a boolean")
    value = operator.index(value)  # NumPy integers become Python integers; floats are refused
    if value < 0:
        raise ValueError(f"{name} must be nonnegative")
    return value


def categorical_pair_moments(probabilities: np.ndarray, f: np.ndarray, g: np.ndarray, count: int) -> PairMoments:
    """Moments of A = sum f(c_i), B = sum g(c_i) over ``count`` iid categorical outcomes c_i."""
    n = _count(count, "count")
    pi = np.asarray(probabilities, dtype=np.float64)
    f = np.asarray(f, dtype=np.float64)
    g = np.asarray(g, dtype=np.float64)
    if not (np.all(np.isfinite(pi)) and np.all(pi >= 0.0) and np.all(np.isfinite(f)) and np.all(np.isfinite(g))):
        raise ValueError("probabilities and values must be finite, probabilities nonnegative")
    ef, eg = float(pi @ f), float(pi @ g)
    y, z = f - ef, g - eg
    eyy, ezz, eyz = float(pi @ (y * y)), float(pi @ (z * z)), float(pi @ (y * z))
    eyyzz = float(pi @ (y * y * z * z))
    pairs = float(n * (n - 1))  # exact integer product, rounded once
    return PairMoments(n * eyy, n * ezz, n * eyz, n * eyyzz + pairs * (eyy * ezz + 2.0 * eyz ** 2))


def single_time_law(scheme: DecayScheme, start: str, t: float) -> np.ndarray:
    """Per-nucleus state probabilities at t (matrix exponential, Layer A reference)."""
    return state_probabilities(scheme, start, t)[0]


def two_time_law(scheme: DecayScheme, start: str, t1: float, t2: float) -> np.ndarray:
    """P(state j at t1 and state k at t2), shape (nuclides, nuclides), for t1 < t2."""
    if not t2 > t1:
        raise ValueError("need t1 < t2")
    first = single_time_law(scheme, start, t1)
    size = len(scheme)
    joint = np.zeros((size, size))
    for j, nuclide in enumerate(scheme.nuclides):
        if first[j] > 0.0:
            joint[j] = first[j] * state_probabilities(scheme, nuclide.name, t2 - t1)[0]
    return joint


def state_count_moments(scheme: DecayScheme, start: str, count: int, t: float, j: int, k: int) -> PairMoments:
    """Moments of (N_j(t), N_k(t)) for ``count`` nuclei started in ``start``."""
    law = single_time_law(scheme, start, t)
    states = np.arange(len(scheme))
    return categorical_pair_moments(law, states == j, states == k, count)


def two_time_count_moments(scheme: DecayScheme, start: str, count: int, t1: float, t2: float, j: int,
                           k: int) -> PairMoments:
    """Moments of (N_j(t1), N_k(t2)) for ``count`` nuclei started in ``start``."""
    joint = two_time_law(scheme, start, t1, t2)
    size = len(scheme)
    first = np.repeat(np.arange(size), size)
    second = np.tile(np.arange(size), size)
    return categorical_pair_moments(joint.ravel(), first == j, second == k, count)


def state_sum_moments(scheme: DecayScheme, start: str, count: int, t: float, j: int, k: int) -> PairMoments:
    """Moments of N_j(t) + N_k(t) (paired with itself); its variance carries Cov(N_j, N_k)."""
    law = single_time_law(scheme, start, t)
    states = np.arange(len(scheme))
    f = (states == j).astype(np.float64) + (states == k).astype(np.float64)
    return categorical_pair_moments(law, f, f, count)


def two_time_sum_moments(scheme: DecayScheme, start: str, count: int, t1: float, t2: float, j: int,
                         k: int) -> PairMoments:
    """Moments of N_j(t1) + N_k(t2) (paired with itself); its variance carries Cov(N_j(t1), N_k(t2))."""
    joint = two_time_law(scheme, start, t1, t2)
    size = len(scheme)
    first = np.repeat(np.arange(size), size)
    second = np.tile(np.arange(size), size)
    f = (first == j).astype(np.float64) + (second == k).astype(np.float64)
    return categorical_pair_moments(joint.ravel(), f, f, count)
