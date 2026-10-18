"""RNGLab v0.14: Classical Monte Carlo Pi Estimation with Channel Isolation.

Supports independent stochastic channel assignment:
- Coordinate proposal channel (x, y generation)
- Resampling / seed initialization channel
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Optional, Union

from rnglab14.sources.base import RandomSource, StochasticChannels


@dataclass(frozen=True)
class PiMCResult:
    total_samples: int
    inside_circle: int
    estimated_pi: float
    true_pi: float
    error: float
    relative_error: float
    standard_error: float
    z_score: float
    confidence_interval_95: tuple[float, float]
    provenance_hash: str


def run_pi_monte_carlo(
    source_or_channels: Union[RandomSource, StochasticChannels],
    samples: int = 10000,
) -> PiMCResult:
    """Estimate pi via hit-or-miss Monte Carlo inside the unit quadrant."""
    if isinstance(source_or_channels, StochasticChannels):
        src = source_or_channels.proposal
    else:
        src = source_or_channels

    inside = 0
    for _ in range(samples):
        x = src.random_float()
        y = src.random_float()
        if x * x + y * y <= 1.0:
            inside += 1

    p_hat = inside / samples
    pi_hat = 4.0 * p_hat
    true_pi = math.pi
    err = pi_hat - true_pi
    rel_err = abs(err) / true_pi

    # Variance of indicator is p*(1-p); variance of 4*p is 16*p*(1-p)/N
    var_pi = 16.0 * p_hat * (1.0 - p_hat) / samples if samples > 1 else 1.0
    se = math.sqrt(var_pi)
    z = err / se if se > 0 else 0.0
    ci = (pi_hat - 1.95996 * se, pi_hat + 1.95996 * se)

    return PiMCResult(
        total_samples=samples,
        inside_circle=inside,
        estimated_pi=round(pi_hat, 6),
        true_pi=round(true_pi, 6),
        error=round(err, 6),
        relative_error=round(rel_err, 6),
        standard_error=round(se, 6),
        z_score=round(z, 4),
        confidence_interval_95=(round(ci[0], 6), round(ci[1], 6)),
        provenance_hash=src.stream_digest,
    )
