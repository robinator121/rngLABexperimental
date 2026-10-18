"""Frozen rank calibration at the independent trajectory/chain level.

For exchangeable null scores, (1 + #{calibration >= test})/(n + 1) is a
conservative marginal p-value, including ties. It does not certify the null
model, independence, target correctness, or conditional calibration. The score
and protocol must be fixed before evaluating held-out results.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from typing import Mapping, Sequence


def _canonical(value: object) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


@dataclass(frozen=True)
class FrozenRankCalibration:
    scores: tuple[float, ...]
    protocol_json: str
    alpha: float = 0.05

    def __post_init__(self) -> None:
        if not 0.0 < self.alpha < 1.0 or not math.isfinite(self.alpha):
            raise ValueError("alpha must be finite and strictly between zero and one")
        if not isinstance(self.scores, tuple) or not self.scores:
            raise ValueError("nonempty immutable calibration scores required")
        if not all(math.isfinite(s) for s in self.scores):
            raise ValueError("finite calibration scores required")
        if len(self.scores) + 1 < math.ceil(1.0 / self.alpha):
            raise ValueError("too few calibration scores to resolve the requested alpha")
        if _canonical(json.loads(self.protocol_json)) != self.protocol_json:
            raise ValueError("protocol_json must use canonical JSON")

    @classmethod
    def fit(cls, scores: Sequence[float], protocol: Mapping[str, object],
            alpha: float = 0.05) -> "FrozenRankCalibration":
        return cls(tuple(float(s) for s in scores), _canonical(dict(protocol)), alpha)

    @property
    def digest(self) -> str:
        return hashlib.sha256(_canonical({"scores": self.scores,
            "protocol": json.loads(self.protocol_json), "alpha": self.alpha}).encode()).hexdigest()

    def p_value(self, score: float, protocol: Mapping[str, object] | None = None) -> float:
        if not math.isfinite(score):
            raise ValueError("finite held-out score required")
        if protocol is not None and _canonical(dict(protocol)) != self.protocol_json:
            raise ValueError("held-out protocol differs from frozen calibration protocol")
        return (1 + sum(s >= score for s in self.scores)) / (len(self.scores) + 1)

    def alarm(self, score: float, protocol: Mapping[str, object] | None = None,
              alpha: float | None = None) -> bool:
        level = self.alpha if alpha is None else alpha
        if not math.isfinite(level) or not 0.0 < level < 1.0:
            raise ValueError("alarm alpha must be finite and strictly between zero and one")
        return self.p_value(score, protocol) <= level
