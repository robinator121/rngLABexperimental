"""Declared degradation families: severity as a function of the sweep index.

A family is a law over severity paths s(0), ..., s(H-1). A forecast of a
future outcome is defined only relative to such a law (see
docs/v019/FORECASTING_METHODS.md, section 1): the path before issuance does
not determine the path after it.

Every family exposes `cumulative(n)`, the exact sum of s(t) over t < n, so that
the mean severity over any block of sweeps is a difference of two values and no
per-sweep array is needed.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

import numpy as np
from numpy.typing import ArrayLike, NDArray

FloatArray = NDArray[np.float64]


class SeverityPath(Protocol):
    """One severity path."""

    def cumulative(self, n: ArrayLike) -> FloatArray:
        """Sum of s(t) for t < n, elementwise in n."""
        ...


def block_means(path: SeverityPath, edges: ArrayLike) -> FloatArray:
    """Mean severity over each interval [edges[i], edges[i+1])."""
    edges = np.asarray(edges, dtype=np.float64)
    total = path.cumulative(edges)
    return np.diff(total, axis=-1) / np.diff(edges)


def _expand(parameter: ArrayLike, n: ArrayLike) -> FloatArray:
    """Reshape a parameter of shape P to P + (1,)*n.ndim, so that results have shape P + n.shape."""
    value = np.asarray(parameter, dtype=np.float64)
    return value.reshape(value.shape + (1,) * np.ndim(n))


def _ramp_sum(t0: FloatArray, t_cap: FloatArray, s_max: float, n: FloatArray) -> FloatArray:
    """Sum over t < n of min(s_max, s_max/t_cap * (t - t0)) for t >= t0; arguments broadcast."""
    m = np.maximum(n - 1 - t0, 0.0)        # positive offsets 1..m lie below n
    q = np.minimum(m, np.floor(t_cap))     # offsets still on the ramp
    return s_max / t_cap * q * (q + 1) / 2 + s_max * (m - q)


@dataclass(frozen=True)
class LinearRamp:
    """s(t) = min(s_max, r (t - t0)) for t >= t0, else 0, with r = s_max / t_cap (the Stage D profile).

    Fields may be arrays of a common shape P: `cumulative(n)` then has shape P + n.shape.
    """
    t0: ArrayLike
    t_cap: ArrayLike
    s_max: float

    def cumulative(self, n: ArrayLike) -> FloatArray:
        n = np.asarray(n, dtype=np.float64)
        return _ramp_sum(_expand(self.t0, n), _expand(self.t_cap, n), self.s_max, n)

    def severity(self, t: ArrayLike) -> FloatArray:
        t = np.asarray(t, dtype=np.float64)
        t0, t_cap = np.asarray(self.t0, dtype=np.float64), np.asarray(self.t_cap, dtype=np.float64)
        return np.where(t >= t0, np.minimum(self.s_max, self.s_max / t_cap * (t - t0)), 0.0)


@dataclass(frozen=True)
class StepHold:
    """s(t) = level for t >= t0, else 0. Not in the ramp family: an 'unseen' continuation."""
    t0: ArrayLike
    level: ArrayLike

    def cumulative(self, n: ArrayLike) -> FloatArray:
        n = np.asarray(n, dtype=np.float64)
        return _expand(self.level, n) * np.maximum(n - _expand(self.t0, n), 0.0)


@dataclass(frozen=True)
class RampThenHold:
    """A linear ramp that stops growing at sweep `freeze` and keeps the severity reached there.

    With freeze equal to the issuance sweep this path is identical to the ramp up to
    issuance and different after it: the constructive non-identifiable pair.
    """
    t0: ArrayLike
    t_cap: ArrayLike
    s_max: float
    freeze: ArrayLike

    def cumulative(self, n: ArrayLike) -> FloatArray:
        n = np.asarray(n, dtype=np.float64)
        t0, t_cap, freeze = _expand(self.t0, n), _expand(self.t_cap, n), _expand(self.freeze, n)
        level = np.where(freeze >= t0, np.minimum(self.s_max, self.s_max / t_cap * (freeze - t0)), 0.0)
        return _ramp_sum(t0, t_cap, self.s_max, np.minimum(n, freeze)) + level * np.maximum(n - freeze, 0.0)


@dataclass(frozen=True)
class SampledPath:
    """A path given by its values on blocks of `block` sweeps (random-walk and other simulated families)."""
    values: FloatArray  # shape (..., blocks): severity held constant inside each block
    block: int

    def cumulative(self, n: ArrayLike) -> FloatArray:
        n = np.asarray(n, dtype=np.float64)
        sums = np.concatenate([np.zeros(self.values.shape[:-1] + (1,)), np.cumsum(self.values, axis=-1) * self.block], axis=-1)
        index = n / self.block
        if not np.all(index == np.floor(index)):
            raise ValueError("a sampled path is defined on block boundaries only")
        return np.take(sums, index.astype(np.int64), axis=-1)


@dataclass(frozen=True)
class RampPrior:
    """The declared law of the ramp parameters: integer onset uniform on [t0_low, t0_high],
    t_cap log-uniform on [t_cap_low, t_cap_high), independent (Stage D protocol, section 2.2)."""
    t0_low: int = 2048
    t0_high: int = 8192
    t_cap_low: float = 2048.0
    t_cap_high: float = 163840.0
    s_max: float = 0.5

    def sample(self, rng: np.random.Generator, size: int) -> LinearRamp:
        t0 = rng.integers(self.t0_low, self.t0_high + 1, size=size)
        t_cap = np.exp(rng.uniform(np.log(self.t_cap_low), np.log(self.t_cap_high), size=size))
        return LinearRamp(t0, t_cap, self.s_max)

    def grid(self, onsets: int, rates: int) -> tuple[LinearRamp, FloatArray]:
        """A product grid of (t0, t_cap) with equal prior weight per cell, and those weights.

        Onsets are cell midpoints of an even partition of [t0_low, t0_high]; t_cap values are
        midpoints in log scale. Equal weights then follow from the declared uniform laws.
        """
        t0 = self.t0_low + (np.arange(onsets) + 0.5) * (self.t0_high - self.t0_low) / onsets
        log_cap = np.log(self.t_cap_low) + (np.arange(rates) + 0.5) * np.log(self.t_cap_high / self.t_cap_low) / rates
        a, b = np.meshgrid(t0, np.exp(log_cap), indexing="ij")
        return LinearRamp(a.ravel(), b.ravel(), self.s_max), np.full(onsets * rates, 1.0 / (onsets * rates))


def random_walk_paths(rng: np.random.Generator, size: int, blocks: int, block: int, *, t0_low: int, t0_high: int,
                      drift: float, volatility: float, s_max: float) -> SampledPath:
    """Severity that starts at 0 at a uniform onset and then follows a reflected random walk with drift,
    in steps of one block, clipped to [0, s_max]. Not in the ramp family."""
    onset_block = rng.integers(t0_low // block, t0_high // block + 1, size=size)
    steps = drift * block + volatility * np.sqrt(block) * rng.standard_normal((size, blocks))
    steps[np.arange(blocks)[None, :] < onset_block[:, None]] = 0.0
    values = np.empty((size, blocks))
    level = np.zeros(size)
    for j in range(blocks):
        level = np.clip(level + steps[:, j], 0.0, s_max)
        values[:, j] = level
    return SampledPath(values, block)
