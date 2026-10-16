"""MATH-06: Recovery, Hysteresis, and Estimator Memory.

Separates and models:
1. Random-Source Recovery (instantaneous upon hardware/stream restore).
2. Markov-Chain Relaxation Memory (exponential decay governed by spectral gap gamma):
   ||pi_t - pi_0||_TV <= C * (1 - gamma)^(t - T_restore) = C * exp(-(t - T_restore) / tau_relax)
3. Cumulative Estimator Memory (algebraic persistent decay):
   Bias(t) = Delta_mu * (tau_defect / t)  for cumulative averages
4. Regenerative / Reset Estimators (complete memory erasure after burn-in).
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import random
from typing import Sequence


@dataclass(frozen=True)
class RecoveryTrajectoryResult:
    time_step: int
    true_mean: float
    is_source_healthy: bool
    source_bias: float
    instantaneous_sample_mean: float
    cumulative_estimator_mean: float
    cumulative_estimator_bias: float
    reset_estimator_mean: float | None
    reset_estimator_bias: float | None
    theoretical_cumulative_bias: float


def simulate_recovery_trajectory(
    t_healthy_pre: int = 200,
    t_defect: int = 300,
    t_healthy_post: int = 500,
    defect_bias: float = 0.30,
    seed: int = 42,
) -> list[RecoveryTrajectoryResult]:
    """Simulate complete recovery trajectory showing algebraic cumulative memory vs reset estimator."""
    rng = random.Random(seed)
    total_steps = t_healthy_pre + t_defect + t_healthy_post
    t_onset = t_healthy_pre
    t_restore = t_healthy_pre + t_defect

    trajectory = []
    cum_sum = 0.0
    reset_sum = 0.0
    reset_count = 0
    true_target = 0.5

    for t in range(1, total_steps + 1):
        # Determine source state
        in_defect = (t_onset < t <= t_restore)
        bias = defect_bias if in_defect else 0.0
        p_val = 0.5 + bias

        sample = 1.0 if rng.random() < p_val else 0.0
        cum_sum += sample
        cum_est = cum_sum / t
        cum_bias = cum_est - true_target

        # Theoretical cumulative bias
        if t <= t_onset:
            theo_cum_bias = 0.0
        elif t <= t_restore:
            theo_cum_bias = defect_bias * (t - t_onset) / t
        else:
            theo_cum_bias = defect_bias * t_defect / t

        # Reset estimator (starts accumulating only after restore)
        if t > t_restore:
            reset_sum += sample
            reset_count += 1
            reset_est = reset_sum / reset_count
            reset_bias = reset_est - true_target
        else:
            reset_est = None
            reset_bias = None

        trajectory.append(RecoveryTrajectoryResult(
            time_step=t,
            true_mean=true_target,
            is_source_healthy=(not in_defect),
            source_bias=round(bias, 4),
            instantaneous_sample_mean=round(sample, 4),
            cumulative_estimator_mean=round(cum_est, 6),
            cumulative_estimator_bias=round(cum_bias, 6),
            reset_estimator_mean=round(reset_est, 6) if reset_est is not None else None,
            reset_estimator_bias=round(reset_bias, 6) if reset_bias is not None else None,
            theoretical_cumulative_bias=round(theo_cum_bias, 6),
        ))

    return trajectory


def predict_recovery_horizon(
    tau_defect: int,
    defect_bias: float,
    acceptable_tolerance: float = 0.01,
) -> int:
    """Predict number of healthy samples required for cumulative bias to decay below tolerance:

    tolerance >= defect_bias * (tau_defect / (T_restore + H))
    => H >= tau_defect * (defect_bias / tolerance) - tau_defect
    """
    if acceptable_tolerance <= 0:
        raise ValueError("Tolerance must be positive")
    required_total_samples = int(math.ceil(tau_defect * abs(defect_bias) / acceptable_tolerance))
    return max(0, required_total_samples - tau_defect)
