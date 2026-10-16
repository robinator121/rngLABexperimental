"""Monte Carlo Pi Estimation baseline application for RNGLab v0.13+."""
from __future__ import annotations

from dataclasses import dataclass
import math

from rnglab13.sources.base import RandomSource


@dataclass(frozen=True)
class PiMCResult:
    total_samples: int
    inside_count: int
    pi_estimate: float
    true_pi: float
    absolute_error: float
    theoretical_std: float
    z_score: float


def run_pi_monte_carlo(source: RandomSource, samples: int = 10000) -> PiMCResult:
    """Run Monte Carlo pi estimation consuming random coordinates from source."""
    inside = 0
    for _ in range(samples):
        # Coordinates in [-1, 1]
        x = source.random_float() * 2.0 - 1.0
        y = source.random_float() * 2.0 - 1.0
        if x * x + y * y <= 1.0:
            inside += 1

    p_hat = inside / samples
    pi_est = 4.0 * p_hat
    p_true = math.pi / 4.0
    theo_std = 4.0 * math.sqrt(p_true * (1.0 - p_true) / samples)
    abs_err = abs(pi_est - math.pi)
    z = (pi_est - math.pi) / theo_std if theo_std > 0 else 0.0

    return PiMCResult(
        total_samples=samples,
        inside_count=inside,
        pi_estimate=round(pi_est, 6),
        true_pi=round(math.pi, 6),
        absolute_error=round(abs_err, 6),
        theoretical_std=round(theo_std, 6),
        z_score=round(z, 3),
    )
