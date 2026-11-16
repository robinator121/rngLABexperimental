"""RNGLab v0.20 nuclear module, Layer A: analytic references for decay.

Nothing here draws random numbers or imports the samplers. Two independent
analytic routes give the probability that one nucleus, started in a given
state at t = 0, is in state j at time t:

1. the matrix exponential of the decay generator, P(t) = expm(A t), where
   A[i, i] = -lambda_i and A[j, i] = sum of b_k * lambda_i over the channels
   i -> j (scipy.linalg.expm, scaling and squaring with Pade approximants);
2. the closed-form Bateman solution summed over every decay path, weighted by
   the product of branching ratios along the path.

Route 2 needs distinct decay constants along each path and loses accuracy
when several of them are close (the alternating sum cancels); ``bateman_chain``
refuses rates closer than ``DEGENERACY_RTOL`` and any evaluation whose
rounding-error bound exceeds ``ROUNDING_LIMIT``. Route 1 has no such restriction.
Both routes use the branching law ``Nuclide.channel_probabilities()``, the one
the samplers use.

Because nuclei are independent and a scheme is acyclic, the following laws
are exact for N nuclei all started in the same state:
- the counts per state at a fixed time are Multinomial(N, p(t));
- the number of decays of nuclide j within [0, t] is Binomial(N, q_j(t)),
  q_j = ``decay_probabilities``, because each nucleus decays out of j at most once;
- the decays of j in disjoint time bins, together with "no decay of j by the
  horizon", are Multinomial(N, bin probabilities);
- given D decays of j, the channel counts are Multinomial(D, branching ratios),
  and given C decays through a channel, its photon count is Binomial(C, yield).
"""
from __future__ import annotations

import math
import sys
from typing import Callable, Sequence

import numpy as np
from scipy.linalg import expm

from .decay import DecayScheme

DEGENERACY_RTOL = 1e-6
ROUNDING_LIMIT = 1e-10


def generator_matrix(scheme: DecayScheme) -> np.ndarray:
    """Decay generator A (column i: transitions out of nuclide i)."""
    size = len(scheme)
    matrix = np.zeros((size, size))
    for i, nuclide in enumerate(scheme.nuclides):
        if nuclide.stable:
            continue
        rate = math.log(2.0) / nuclide.half_life
        matrix[i, i] = -rate
        for channel, probability in zip(nuclide.channels, nuclide.channel_probabilities()):
            matrix[scheme.index[channel.daughter], i] += probability * rate
    return matrix


def _times(times: float | Sequence[float]) -> np.ndarray:
    times = np.atleast_1d(np.asarray(times, dtype=np.float64))
    if np.any(~np.isfinite(times)) or np.any(times < 0):
        raise ValueError("times must be finite and nonnegative")
    return times


def state_probabilities(scheme: DecayScheme, start: str, times: float | Sequence[float]) -> np.ndarray:
    """P(state j at t | start at 0) by matrix exponential; shape (len(times), nuclides)."""
    matrix, column = generator_matrix(scheme), scheme.index[start]
    return np.array([expm(matrix * t)[:, column] for t in _times(times)])


def decay_probabilities(scheme: DecayScheme, start: str, times: float | Sequence[float]) -> np.ndarray:
    """P(the nucleus has decayed out of nuclide j by t | start at 0); shape (len(times), nuclides).

    Uses the generator augmented with one counter per nuclide, dC_j/dt = lambda_j N_j.
    In an acyclic scheme the expected count equals the probability.
    """
    size = len(scheme)
    augmented = np.zeros((2 * size, 2 * size))
    augmented[:size, :size] = generator_matrix(scheme)
    for j, nuclide in enumerate(scheme.nuclides):
        if not nuclide.stable:
            augmented[size + j, j] = math.log(2.0) / nuclide.half_life
    column = scheme.index[start]
    return np.array([expm(augmented * t)[size:, column] for t in _times(times)])


def expected_populations(scheme: DecayScheme, initial: dict[str, int], times: float | Sequence[float]) -> np.ndarray:
    """Expected counts per nuclide at each time for a mixed initial population."""
    total = np.zeros((_times(times).size, len(scheme)))
    for name, count in initial.items():
        total += count * state_probabilities(scheme, name, times)
    return total


def bateman_chain(rates: Sequence[float], times: float | Sequence[float]) -> np.ndarray:
    """Closed-form Bateman probability of being at the last member of a linear chain.

    rates[q] is the decay constant of member q; the last may be 0 (stable).
    P_m(t) = (prod_{p<m} rates[p]) * sum_q exp(-rates[q] t) / prod_{r != q} (rates[r] - rates[q]).

    The result depends only on the products rate * t and on rate ratios, so it is
    evaluated on rates divided by the largest rate and on times multiplied by it;
    this keeps the products and differences away from underflow. The alternating
    sum still cancels badly when rates are close, and distinct rates alone do not
    prevent it. Every evaluation is therefore checked against a rounding-error
    bound, 4 (m + 1) eps * prefactor * sum_q |term_q| on the scaled problem, and
    refused (ValueError) if the bound exceeds ``ROUNDING_LIMIT``, if a scaled
    denominator is not a normal float (the bound would not hold), or if anything
    is not finite. Use the matrix exponential there.
    """
    rates = [float(r) for r in rates]
    if not rates or any(r < 0 or not math.isfinite(r) for r in rates) or any(r == 0 for r in rates[:-1]):
        raise ValueError("rates must be finite, nonnegative, and positive except possibly the last")
    scale = max(rates)
    for q in range(len(rates)):
        for r in range(q + 1, len(rates)):
            if abs(rates[q] - rates[r]) <= DEGENERACY_RTOL * scale:
                raise ValueError("rates along a path must be distinct for the closed form")
    times = _times(times)
    if scale == 0.0:  # a single stable member
        return np.ones_like(times)
    scaled = [r / scale for r in rates]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore", under="ignore"):  # checked below
        scaled_times = times * scale
        denominators = [math.prod(scaled[r] - rate for r in range(len(scaled)) if r != q)
                        for q, rate in enumerate(scaled)]
        prefactor = math.prod(scaled[:-1])
        terms = np.array([np.exp(-rate * scaled_times) / denominator
                          for rate, denominator in zip(scaled, denominators)])
        result = prefactor * terms.sum(axis=0)
        bound = 4.0 * len(scaled) * np.finfo(float).eps * prefactor * np.abs(terms).sum(axis=0)
    usable = (np.all(np.isfinite(scaled_times)) and all(abs(d) >= sys.float_info.min for d in denominators)
              and prefactor >= sys.float_info.min and np.all(np.isfinite(result)) and np.all(np.isfinite(bound)))
    if not usable or np.any(bound > ROUNDING_LIMIT):
        raise ValueError("closed form is ill-conditioned or not finite for these rates; use the matrix exponential")
    return result


def path_state_probabilities(scheme: DecayScheme, start: str, times: float | Sequence[float]) -> np.ndarray:
    """P(state j at t | start at 0) by the closed-form Bateman solution summed over paths."""
    times = _times(times)
    rates = [0.0 if n.stable else math.log(2.0) / n.half_life for n in scheme.nuclides]
    out = np.zeros((times.size, len(scheme)))
    for path, weight in scheme.paths(start):
        out[:, path[-1]] += weight * bateman_chain([rates[i] for i in path], times)
    return out


def survival_probability(half_life: float, times: float | Sequence[float]) -> np.ndarray:
    """P(a nucleus with this half-life has not decayed by t) = 2**(-t / half_life)."""
    return np.exp(-math.log(2.0) * _times(times) / half_life)


def truncated_exponential_cdf(rate: float, horizon: float) -> Callable[[np.ndarray], np.ndarray]:
    """CDF of an exponential(rate) waiting time conditioned on being at most ``horizon``."""
    if not rate > 0 or not horizon > 0:
        raise ValueError("rate and horizon must be positive")
    norm = -math.expm1(-rate * horizon) if math.isfinite(horizon) else 1.0

    def cdf(x):
        x = np.clip(np.asarray(x, dtype=np.float64), 0.0, horizon)
        return -np.expm1(-rate * x) / norm

    return cdf


def transient_equilibrium_ratio(parent_half_life: float, daughter_half_life: float) -> float:
    """Limit of daughter/parent activity, lambda_D / (lambda_D - lambda_P), for a longer-lived parent."""
    if not parent_half_life > daughter_half_life > 0:
        raise ValueError("needs a parent longer-lived than its daughter")
    parent, daughter = math.log(2.0) / parent_half_life, math.log(2.0) / daughter_half_life
    return daughter / (daughter - parent)
