"""MATH-03: Degradation Trajectory Families for Forecasting Evaluation.

Implements reproducible trajectories:
1. Constant bias
2. Linear drift
3. Accelerating (quadratic) drift
4. Gradual autocorrelation increase
5. Stepwise failure
6. Intermittent / periodic faults
7. Fault recovery
8. Stochastic drift (Ornstein-Uhlenbeck / Random Walk)
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import math
import random
from typing import Any, Callable, Sequence


class TrajectoryFamily(str, Enum):
    CONSTANT_BIAS = "constant_bias"
    LINEAR_DRIFT = "linear_drift"
    ACCELERATING_DRIFT = "accelerating_drift"
    GRADUAL_AUTOCORRELATION = "gradual_autocorrelation"
    STEPWISE_FAILURE = "stepwise_failure"
    INTERMITTENT_FAULT = "intermittent_fault"
    FAULT_RECOVERY = "fault_recovery"
    STOCHASTIC_DRIFT = "stochastic_drift"


@dataclass(frozen=True)
class DegradationPoint:
    t: int
    parameter: float
    is_fault_active: bool
    description: str


class TrajectoryGenerator:
    """Generates ground-truth parameter evolution theta_t across time."""

    def __init__(self, family: TrajectoryFamily, seed: int = 42, **params: Any):
        self.family = family
        self.seed = seed
        self.params = params
        self.rng = random.Random(seed)

    def generate_path(self, total_steps: int) -> list[DegradationPoint]:
        path = []
        family = self.family

        if family == TrajectoryFamily.CONSTANT_BIAS:
            bias = float(self.params.get("bias", 0.15))
            for t in range(1, total_steps + 1):
                path.append(DegradationPoint(t, bias, True, "constant_bias"))

        elif family == TrajectoryFamily.LINEAR_DRIFT:
            start_val = float(self.params.get("start_val", 0.0))
            drift_rate = float(self.params.get("drift_rate", 0.0005))
            for t in range(1, total_steps + 1):
                val = start_val + drift_rate * t
                path.append(DegradationPoint(t, val, val > 0.05, "linear_drift"))

        elif family == TrajectoryFamily.ACCELERATING_DRIFT:
            start_val = float(self.params.get("start_val", 0.0))
            coeff = float(self.params.get("coeff", 1e-6))
            for t in range(1, total_steps + 1):
                val = start_val + coeff * (t ** 2)
                path.append(DegradationPoint(t, val, val > 0.05, "accelerating_drift"))

        elif family == TrajectoryFamily.GRADUAL_AUTOCORRELATION:
            start_rho = float(self.params.get("start_rho", 0.0))
            drift_rate = float(self.params.get("drift_rate", 0.0004))
            for t in range(1, total_steps + 1):
                rho = min(0.95, start_rho + drift_rate * t)
                path.append(DegradationPoint(t, rho, rho > 0.1, "gradual_autocorrelation"))

        elif family == TrajectoryFamily.STEPWISE_FAILURE:
            step_time = int(self.params.get("step_time", total_steps // 2))
            baseline = float(self.params.get("baseline", 0.0))
            failure_level = float(self.params.get("failure_level", 0.40))
            for t in range(1, total_steps + 1):
                val = failure_level if t >= step_time else baseline
                path.append(DegradationPoint(t, val, t >= step_time, "stepwise_failure"))

        elif family == TrajectoryFamily.INTERMITTENT_FAULT:
            period = int(self.params.get("period", 200))
            duty_cycle = float(self.params.get("duty_cycle", 0.3))
            amplitude = float(self.params.get("amplitude", 0.35))
            for t in range(1, total_steps + 1):
                in_pulse = ((t % period) / period) < duty_cycle
                val = amplitude if in_pulse else 0.0
                path.append(DegradationPoint(t, val, in_pulse, "intermittent_fault"))

        elif family == TrajectoryFamily.FAULT_RECOVERY:
            t_onset = int(self.params.get("t_onset", total_steps // 4))
            t_restore = int(self.params.get("t_restore", 3 * total_steps // 4))
            severity = float(self.params.get("severity", 0.30))
            for t in range(1, total_steps + 1):
                in_fault = (t_onset <= t <= t_restore)
                val = severity if in_fault else 0.0
                path.append(DegradationPoint(t, val, in_fault, "fault_recovery"))

        elif family == TrajectoryFamily.STOCHASTIC_DRIFT:
            # Random walk with positive drift
            cur = float(self.params.get("start_val", 0.0))
            drift = float(self.params.get("drift", 0.0003))
            vol = float(self.params.get("volatility", 0.005))
            for t in range(1, total_steps + 1):
                cur += drift + self.rng.gauss(0.0, vol)
                path.append(DegradationPoint(t, cur, cur > 0.05, "stochastic_drift"))

        return path
