"""Bounded numerical-uniform diagnostics and externally calibrated decisions.

Scores are fixed statistics, not claimed analytical p-values for arbitrary
source processes. Calibrate the complete horizon and confirmation rule on fresh
healthy streams before supplying ``thresholds``. Source bits are not observed:
the monobit component uses an explicitly diagnostic 32-bit encoding of consumed
numerical uniforms.
"""
from __future__ import annotations

import math
from typing import Any
import numpy as np


RULES = ("monobit", "distribution", "serial", "temporal", "periodic", "cross", "cusum", "composite")
TEMPORAL_LAGS = (1, 2, 8, 32)
PERIODIC_LAGS = (128, 256)
COMPOSITE_SCALES = {"monobit": 1., "distribution": 1., "serial": 1.,
                    "temporal": 1., "periodic": 1., "cross": 1., "cusum": .5}
SCORE_DEFINITIONS = {
    "monobit": "max_j |2*ones_j-32*N|/sqrt(32*N), floor(2^32*u) diagnostic encoding",
    "distribution": "max_j max(0,(chi2_16bins_j-15)/sqrt(30))",
    "serial": "max_j |2*same_side_j-(N-1)|/sqrt(N-1), side=(u>=0.5)",
    "temporal": "max_j,lag in {1,2,8,32} sqrt(N-lag)*|Pearson(u[:-lag,j],u[lag:,j])|",
    "periodic": "max_j,lag in {128,256} sqrt(N-lag)*|Pearson(u[:-lag,j],u[lag:,j])|",
    "cross": "max_j<k sqrt(N)*|Pearson(u[:,j],u[:,k])|; zero for one column",
    "cusum": "max_j C_j at the current boundary; C=max(0,C+log(1.1) for same side, C+log(0.9) for flip), exact state across calls",
    "composite": "max fixed component scales: 1 for all components except 0.5 for CUSUM",
}


def _pearson_columns(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    """Return paired-column Pearson rho; constant-column pairs score zero.

    Constants are separately flagged and cause large distribution scores. A
    zero correlation score for a constant column is not a claim of independence.
    """
    a = left - left.mean(axis=0)
    b = right - right.mean(axis=0)
    numerator = np.sum(a * b, axis=0)
    denominator = np.sqrt(np.sum(a * a, axis=0) * np.sum(b * b, axis=0))
    rho = np.divide(numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0)
    return np.clip(rho, -1., 1.)


def window_scores(uniforms: np.ndarray) -> dict[str, float]:
    """Fixed scores on a chronological recent window, independent of thresholds."""
    # Columns are independent score channels; contiguous columns substantially
    # reduce the cost of repeated lag reductions for two/three-dimensional draws.
    u = np.asfortranarray(uniforms, dtype=np.float64)
    n, dimension = u.shape
    if not n:
        raise ValueError("empty monitor window")
    encoded = np.floor(u * 4294967296.).astype(np.uint32)
    ones = np.bitwise_count(encoded).sum(axis=0, dtype=np.int64)
    scores = {"monobit": float(np.max(np.abs(2 * ones - 32 * n)) / math.sqrt(32 * n))}

    bins = np.minimum((u * 16).astype(np.int64), 15)
    counts = np.zeros((dimension, 16), dtype=np.int64)
    for j in range(dimension):
        counts[j] = np.bincount(bins[:, j], minlength=16)
    expected = n / 16
    chisq = np.sum((counts - expected)**2 / expected, axis=1)
    scores["distribution"] = float(max(0., np.max((chisq - 15) / math.sqrt(30))))

    if n > 1:
        side = u >= .5
        same = np.sum(side[1:] == side[:-1], axis=0)
        scores["serial"] = float(np.max(np.abs(2 * same - (n - 1))) / math.sqrt(n - 1))
    else:
        scores["serial"] = 0.
    for name, lags in (("temporal", TEMPORAL_LAGS), ("periodic", PERIODIC_LAGS)):
        scores[name] = max((float(math.sqrt(n - lag) * np.max(np.abs(_pearson_columns(u[:-lag], u[lag:]))))
                            for lag in lags if n - lag >= 3), default=0.)
    cross = 0.
    if dimension > 1 and n >= 3:
        for j in range(dimension - 1):
            left = np.broadcast_to(u[:, j:j+1], (n, dimension - j - 1))
            cross = max(cross, float(math.sqrt(n) * np.max(np.abs(_pearson_columns(left, u[:, j+1:])))))
    scores["cross"] = cross
    return scores


class Monitor:
    """Recent-window scores plus exact sequential CUSUM and frozen decisions.

    Configuration keys:
    ``monitor_rows`` (or ``window_rows``): bounded recent window, default 4096;
    ``thresholds``: optional mapping from RULES to frozen upper score limits;
    ``confirmation_windows``: consecutive strict exceedances, default 2.

    ``update(u,time)`` records the score at the supplied proposal boundary. The
    consumed-row/draw counters are also returned, so the runner can explicitly
    verify its own time mapping. First excursions and confirmations are never
    backdated. Missing thresholds mean score collection without decisions.
    """
    def __init__(self, config: dict[str, Any], dimension: int):
        self.config = dict(config)
        self.dimension = int(dimension)
        self.capacity = int(config.get("monitor_rows", config.get("window_rows", 4096)))
        self.confirmation_windows = int(config.get("confirmation_windows", 2))
        if self.dimension < 1 or self.capacity < 2 or self.confirmation_windows < 1:
            raise ValueError("positive dimension/confirmation and monitor_rows>=2 required")
        self.thresholds = {str(k): float(v) for k, v in config.get("thresholds", {}).items()}
        if set(self.thresholds) - set(RULES):
            raise ValueError("unknown detector threshold")
        if any(math.isnan(v) or v < 0 for v in self.thresholds.values()):
            raise ValueError("thresholds must be nonnegative and not NaN")
        self._buffer = np.empty((self.capacity, self.dimension), dtype=np.float64)
        self._size = self._write = 0
        self.rows = self.looks = 0
        self.last_time: int | None = None
        self._last_side: np.ndarray | None = None
        self._cusum = np.zeros(self.dimension)
        self.first_excursions: dict[str, int | None] = {name: None for name in RULES}
        self.confirmed_alarms: dict[str, int | None] = {name: None for name in RULES}
        self._streaks = {name: 0 for name in RULES}

    def _append(self, u: np.ndarray) -> None:
        n = len(u)
        if n >= self.capacity:
            self._buffer[:] = u[-self.capacity:]
            self._write = 0
            self._size = self.capacity
            return
        first = min(n, self.capacity - self._write)
        self._buffer[self._write:self._write + first] = u[:first]
        if first < n:
            self._buffer[:n - first] = u[first:]
        self._write = (self._write + n) % self.capacity
        self._size = min(self.capacity, self._size + n)

    def _ordered(self) -> np.ndarray:
        if self._size < self.capacity:
            return self._buffer[:self._size]
        if self._write == 0:
            return self._buffer
        return np.concatenate((self._buffer[self._write:], self._buffer[:self._write]), axis=0)

    def _update_cusum(self, u: np.ndarray) -> None:
        side = u >= .5
        if self._last_side is None:
            same = side[1:] == side[:-1]
        else:
            same = np.concatenate((side[:1] == self._last_side, side[1:] == side[:-1]), axis=0)
        if len(same):
            increments = np.where(same, math.log(1.1), math.log(.9))
            unreflected = self._cusum + np.cumsum(increments, axis=0)
            running_floor = np.minimum.accumulate(np.minimum(unreflected, 0.), axis=0)
            self._cusum = (unreflected - running_floor)[-1]
        self._last_side = side[-1].copy()

    def update(self, uniforms: np.ndarray, time: int) -> dict[str, Any]:
        u = np.asarray(uniforms, dtype=np.float64)
        if u.ndim != 2 or u.shape[1] != self.dimension or not len(u):
            raise ValueError("nonempty uniform matrix with configured dimension required")
        if not np.isfinite(u).all() or np.any(u < 0) or np.any(u >= 1):
            raise ValueError("uniform variates must be finite in [0,1)")
        if not isinstance(time, (int, np.integer)) or time < 0 or (self.last_time is not None and time <= self.last_time):
            raise ValueError("proposal boundary must be an increasing nonnegative integer")
        self._update_cusum(u)
        self._append(u)
        self.rows += len(u)
        self.looks += 1
        self.last_time = int(time)
        recent = self._ordered()
        scores = window_scores(recent)
        scores["cusum"] = float(np.max(self._cusum))
        scores["composite"] = max(scores[name] * factor for name, factor in COMPOSITE_SCALES.items())
        excursions, confirmations = {}, {}
        for name, cutoff in self.thresholds.items():
            above = scores[name] > cutoff
            self._streaks[name] = self._streaks[name] + 1 if above else 0
            excursions[name] = bool(above)
            if above and self.first_excursions[name] is None:
                self.first_excursions[name] = int(time)
            confirmed = self._streaks[name] >= self.confirmation_windows
            confirmations[name] = bool(confirmed)
            if confirmed and self.confirmed_alarms[name] is None:
                self.confirmed_alarms[name] = int(time)
        return {"time": int(time), "look": self.looks, "window_rows": len(recent),
                "window_start_row": self.rows - len(recent), "window_end_row": self.rows,
                "consumed_uniforms": self.rows * self.dimension,
                "encoding": "diagnostic uint32=floor(2^32*u); original source bits unavailable",
                "degenerate_columns": np.flatnonzero(np.ptp(recent, axis=0) == 0).tolist(),
                **scores, "excursions": excursions, "confirmations": confirmations,
                "first_excursions": dict(self.first_excursions), "confirmed_alarms": dict(self.confirmed_alarms)}

    def state_dict(self) -> dict[str, Any]:
        """JSON-compatible exact state for checkpoint/replay verification."""
        return {"config": self.config, "dimension": self.dimension, "rows": self.rows,
                "looks": self.looks, "last_time": self.last_time, "buffer": self._ordered().tolist(),
                "last_side": None if self._last_side is None else self._last_side.tolist(),
                "cusum": self._cusum.tolist(), "first_excursions": dict(self.first_excursions),
                "confirmed_alarms": dict(self.confirmed_alarms), "streaks": dict(self._streaks)}

    @classmethod
    def from_state_dict(cls, state: dict[str, Any]) -> "Monitor":
        monitor = cls(state["config"], state["dimension"])
        buffer = np.asarray(state["buffer"], dtype=np.float64).reshape(-1, monitor.dimension)
        if len(buffer):
            monitor._append(buffer)
        monitor.rows, monitor.looks, monitor.last_time = state["rows"], state["looks"], state["last_time"]
        monitor._last_side = None if state["last_side"] is None else np.asarray(state["last_side"], dtype=bool)
        monitor._cusum = np.asarray(state["cusum"], dtype=np.float64)
        monitor.first_excursions = dict(state["first_excursions"])
        monitor.confirmed_alarms = dict(state["confirmed_alarms"])
        monitor._streaks = dict(state["streaks"])
        return monitor
