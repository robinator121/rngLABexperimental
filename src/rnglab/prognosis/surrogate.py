"""A Gaussian surrogate of a monitored Monte Carlo run under a degrading input.

It reproduces the *structure* of the Stage D pilot and nothing of its physics:
a running estimate whose interval is a 16-batch batch-means interval, a bias
proportional to a severity path, a windowed monitoring statistic whose mean is
proportional to the mean severity in its window, a whole-horizon calibrated
alarm, validity at declared run lengths, and a sustained-exclusion failure
time. All noise is Gaussian and independent across blocks of sweeps, so every
quantity has a known law and forecasting methods can be checked against exact
answers. Results on the surrogate are statements about methods, not about the
path-integral sampler or about any random number generator.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from numpy.typing import NDArray
from scipy import special

from .families import SeverityPath, block_means

FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class SurrogateSpec:
    """Geometry and scales. Defaults follow the Stage D pilot's geometry; the scales are of the
    order measured for the LAG1 family in audit probe runs (docs/v019/FORECASTING_METHODS.md)."""
    horizon: int = 40960
    block: int = 16            # sweeps per simulated block
    stride: int = 256          # sweeps between checkpoints
    window: int = 512          # sweeps per monitoring window
    lengths: tuple[int, ...] = (10240, 20480, 40960)
    batches: int = 16
    level: float = 0.95
    beta: float = 0.06         # bias of the observable per unit severity
    sigma_eff: float = 0.613   # long-run noise scale: standard error of an n-sweep mean is sigma_eff/sqrt(n)
    kappa: float = 103.0       # mean of the monitoring statistic per unit window-mean severity
    s_max: float = 0.5
    sustain: int = 3           # consecutive excluded checkpoints that define failure

    def __post_init__(self) -> None:
        if self.window != 2 * self.stride or self.stride % (self.block * self.batches):
            raise ValueError("the window must be two strides and a stride a multiple of block*batches")
        if self.horizon % self.stride or any(length % self.stride for length in self.lengths):
            raise ValueError("horizon and lengths must be multiples of the stride")

    @property
    def blocks(self) -> int:
        return self.horizon // self.block

    @cached_property
    def ends(self) -> IntArray:
        """Exclusive end sweep of every checkpoint window."""
        return np.arange(self.window, self.horizon + 1, self.stride, dtype=np.int64)

    @property
    def censor_sweep(self) -> int:
        """Largest checkpoint at which a sustained exclusion can still be ascertained inside the horizon."""
        return int(self.ends[-self.sustain])

    @property
    def critical(self) -> float:
        return float(special.stdtrit(self.batches - 1, (1 + self.level) / 2))

    def window_mean_severity(self, path: SeverityPath) -> FloatArray:
        """Mean severity in each checkpoint window, shape (..., checkpoints)."""
        half = block_means(path, np.arange(0, self.horizon + 1, self.stride))
        return (half[..., :-1] + half[..., 1:]) / 2


@dataclass
class Trajectories:
    """Simulated runs. `path` is evaluation truth: forecasters receive `History` objects only."""
    spec: SurrogateSpec
    error_block: FloatArray   # (n, blocks): block means of observable minus reference
    statistic: FloatArray     # (n, checkpoints): monitoring statistic
    path: SeverityPath | None

    def __len__(self) -> int:
        return len(self.statistic)

    @cached_property
    def _cumulative(self) -> FloatArray:
        return np.concatenate([np.zeros((len(self), 1)), np.cumsum(self.error_block, axis=1)], axis=1)

    def interval(self, sweeps: int) -> tuple[FloatArray, FloatArray]:
        """(signed error of the estimate, interval half-width) of the batch-means interval over [0, sweeps)."""
        spec = self.spec
        size, rest = divmod(sweeps // spec.block, spec.batches)
        if sweeps % spec.block or rest or size == 0:
            raise ValueError("run length must be a positive multiple of block*batches sweeps")
        edges = np.arange(spec.batches + 1) * size
        means = np.diff(self._cumulative[:, edges], axis=1) / size
        return means.mean(axis=1), spec.critical * np.sqrt(means.var(axis=1, ddof=1) / spec.batches)

    def valid(self, sweeps: int) -> NDArray[np.bool_]:
        error, half = self.interval(sweeps)
        return np.abs(error) <= half

    @cached_property
    def excluded(self) -> NDArray[np.bool_]:
        """(n, checkpoints): the interval at each checkpoint end excludes the reference."""
        return np.stack([~self.valid(int(end)) for end in self.spec.ends], axis=1)

    @cached_property
    def failure(self) -> tuple[IntArray, NDArray[np.bool_]]:
        """(time, observed). Failure time = first checkpoint from which the interval excludes the
        reference at `sustain` consecutive checkpoints. Unobserved inside the horizon: time is the
        censoring sweep and observed is False (administrative right-censoring)."""
        spec, excluded = self.spec, self.excluded
        run = np.ones(excluded.shape[0:1] + (excluded.shape[1] - spec.sustain + 1,), dtype=bool)
        for k in range(spec.sustain):
            run &= excluded[:, k:k + run.shape[1]]
        observed = run.any(axis=1)
        first = np.argmax(run, axis=1)
        return np.where(observed, spec.ends[first], spec.censor_sweep).astype(np.int64), observed


def simulate(spec: SurrogateSpec, path: SeverityPath | None, size: int, rng: np.random.Generator) -> Trajectories:
    """Draw `size` runs. `path` holds one severity path per run (array-valued parameters), or None for healthy."""
    if path is None:
        severity = np.zeros((size, spec.blocks))
        half = np.zeros((size, spec.horizon // spec.stride))
    else:
        severity = np.broadcast_to(block_means(path, np.arange(spec.blocks + 1) * spec.block), (size, spec.blocks))
        half = np.broadcast_to(block_means(path, np.arange(0, spec.horizon + 1, spec.stride)),
                               (size, spec.horizon // spec.stride))
    error = spec.beta * severity + spec.sigma_eff / np.sqrt(spec.block) * rng.standard_normal((size, spec.blocks))
    # Independent half-window statistics; a window statistic is the normalised sum of its two halves,
    # which gives adjacent (overlapping) windows the correlation 1/2 they have for a linear statistic.
    pieces = spec.kappa * half / np.sqrt(2) + rng.standard_normal(half.shape)
    return Trajectories(spec, error, (pieces[:, :-1] + pieces[:, 1:]) / np.sqrt(2), path)


def calibrate_threshold(healthy: Trajectories, order: int = 950, required: int = 999) -> float:
    """Whole-horizon alarm threshold: the `order`-th order statistic of `required` healthy maxima."""
    if len(healthy) != required:
        raise ValueError(f"exactly {required} healthy calibration runs are required")
    return float(np.sort(healthy.statistic.max(axis=1))[order - 1])


def first_alarm(statistic: FloatArray, threshold: float) -> IntArray:
    """Index of the first checkpoint whose statistic strictly exceeds the threshold; -1 for none."""
    hit = statistic > threshold
    return np.where(hit.any(axis=1), np.argmax(hit, axis=1), -1).astype(np.int64)
