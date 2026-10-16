"""NIST SP 800-90B Online Health Tests.

Implements:
1. Repetition Count Test (RCT) - Section 4.4.1
2. Adaptive Proportion Test (APT) - Section 4.4.2

Designed to continuously monitor entropy sources for catastrophic stuck-bit
or distribution-collapse faults with mathematically defined false alarm rates.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


@dataclass(frozen=True)
class HealthTestResult:
    passed: bool
    test_name: str
    sample_index: int
    current_value: int
    cutoff: int
    details: str


class RepetitionCountTest:
    """NIST SP 800-90B Section 4.4.1 Repetition Count Test (RCT).

    Detects when an entropy source outputs identical consecutive samples,
    indicating hardware lockup or stuck-bit failure.
    Cutoff C = 1 + ceil(-log2(alpha) / H)
    """

    def __init__(self, min_entropy_per_sample: float = 1.0, alpha: float = 2**-20):
        if min_entropy_per_sample <= 0:
            raise ValueError("Min-entropy must be strictly positive")
        self.H = min_entropy_per_sample
        self.alpha = alpha
        # SP 800-90B formula
        self.cutoff = 1 + math.ceil((-math.log2(alpha)) / self.H)

        self._last_sample: int | None = None
        self._rep_count: int = 0
        self._total_samples: int = 0
        self._failures: int = 0

    def step(self, sample: int) -> HealthTestResult:
        """Feed a single sample and return health evaluation."""
        self._total_samples += 1
        if self._last_sample is None:
            self._last_sample = sample
            self._rep_count = 1
            return HealthTestResult(
                passed=True,
                test_name="RCT",
                sample_index=self._total_samples,
                current_value=1,
                cutoff=self.cutoff,
                details="Initial sample",
            )

        if sample == self._last_sample:
            self._rep_count += 1
        else:
            self._last_sample = sample
            self._rep_count = 1

        passed = self._rep_count < self.cutoff
        if not passed:
            self._failures += 1

        return HealthTestResult(
            passed=passed,
            test_name="RCT",
            sample_index=self._total_samples,
            current_value=self._rep_count,
            cutoff=self.cutoff,
            details=f"Repetition count {self._rep_count} {'<' if passed else '>='} cutoff {self.cutoff}",
        )


class AdaptiveProportionTest:
    """NIST SP 800-90B Section 4.4.2 Adaptive Proportion Test (APT).

    Detects distribution collapse or severe bias by counting the frequency
    of a baseline sample within a sliding/fixed window W.
    """

    def __init__(
        self,
        window_size: int = 512,
        min_entropy_per_sample: float = 1.0,
        alpha: float = 2**-20,
    ):
        if window_size < 64:
            raise ValueError("Window size W must be at least 64")
        self.W = window_size
        self.H = min_entropy_per_sample
        self.alpha = alpha

        # Binomial upper tail approximation for cutoff C
        p = 2**(-self.H)
        # Rational approximation to normal quantile for small alpha (Abramowitz & Stegun 26.2.23)
        t = math.sqrt(-2.0 * math.log(max(1e-15, self.alpha)))
        z = t - ((0.010328 * t + 0.802853) * t + 2.515517) / (((0.001308 * t + 0.189269) * t + 1.432788) * t + 1.0)
        self.cutoff = 1 + math.ceil((self.W - 1) * p + z * math.sqrt((self.W - 1) * p * (1 - p)))

        self._target_sample: int | None = None
        self._count: int = 0
        self._window_pos: int = 0
        self._total_samples: int = 0
        self._failures: int = 0

    def step(self, sample: int) -> HealthTestResult:
        """Feed single sample and evaluate against window cutoff."""
        self._total_samples += 1

        if self._window_pos == 0:
            self._target_sample = sample
            self._count = 1
            self._window_pos = 1
            return HealthTestResult(
                passed=True,
                test_name="APT",
                sample_index=self._total_samples,
                current_value=1,
                cutoff=self.cutoff,
                details="Window start",
            )

        if sample == self._target_sample:
            self._count += 1
        self._window_pos += 1

        passed = self._count < self.cutoff
        if not passed:
            self._failures += 1

        result = HealthTestResult(
            passed=passed,
            test_name="APT",
            sample_index=self._total_samples,
            current_value=self._count,
            cutoff=self.cutoff,
            details=f"Count {self._count} in window pos {self._window_pos}/{self.W} {'<' if passed else '>='} cutoff {self.cutoff}",
        )

        # Reset window after W samples
        if self._window_pos >= self.W:
            self._window_pos = 0

        return result
