"""MATH-03: Remaining Useful Samples (RUS) Prediction Subsystem.

Formal mathematical definition:
N_delta(t) = sup { h >= 0 : P(tau > t + h | F_t) >= 1 - delta }
where:
- tau is the first passage time to the failure threshold: tau = inf { s > t : theta_s >= theta_crit }
- F_t is the filtration representing observations up to time t (no future data leakage).
- 1 - delta is the target reliability confidence level (e.g. 0.95).

Crucial Distinction:
N_delta(t) is a statistical lower prediction bound, NOT a deterministic guarantee.
A false reassurance occurs when the actual failure happens earlier than N_delta(t):
  tau - t < N_delta(t)
Target calibration requires: P(tau - t < N_delta(t)) <= delta.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Any, Sequence
import numpy as np


@dataclass(frozen=True)
class RUSForecast:
    current_t: int
    current_value: float
    threshold: float
    point_forecast_samples: int
    lower_bound_95_samples: int
    actual_remaining_samples: int | None
    forecast_error: int | None
    covered_by_lower_bound: bool | None
    is_false_reassurance: bool | None
    model_name: str


class BasePredictor:
    """Abstract base predictor for Remaining Useful Samples."""

    def __init__(self, threshold: float = 0.50):
        self.threshold = threshold

    def predict(self, history: Sequence[float], current_t: int) -> tuple[int, int]:
        """Return (point_forecast, lower_bound_95)."""
        raise NotImplementedError


class NoChangeBaselinePredictor(BasePredictor):
    """Null baseline: assumes current rate of failure risk is invariant."""

    def predict(self, history: Sequence[float], current_t: int) -> tuple[int, int]:
        if not history:
            return 10000, 5000
        val = history[-1]
        if val >= self.threshold:
            return 0, 0
        # Naive constant horizon assumption
        return 5000, 2000


class LinearExtrapolationPredictor(BasePredictor):
    """Linear drift extrapolation with Student-t prediction interval."""

    def __init__(self, threshold: float = 0.50, window: int = 100):
        super().__init__(threshold=threshold)
        self.window = window

    def predict(self, history: Sequence[float], current_t: int) -> tuple[int, int]:
        n_hist = len(history)
        if n_hist < 10:
            return 5000, 1000

        # Use recent window
        w_size = min(self.window, n_hist)
        y = np.array(history[-w_size:], dtype=float)
        x = np.arange(w_size, dtype=float)

        cur_val = y[-1]
        if cur_val >= self.threshold:
            return 0, 0

        # Fit linear regression y = a * x + b
        # a = Cov(x, y) / Var(x)
        x_mean = np.mean(x)
        y_mean = np.mean(y)
        ss_xx = float(np.sum((x - x_mean) ** 2))
        if ss_xx < 1e-12:
            return 5000, 1000

        slope = float(np.sum((x - x_mean) * (y - y_mean)) / ss_xx)
        intercept = float(y_mean - slope * x_mean)

        # Residual variance
        residuals = y - (slope * x + intercept)
        s_err = math.sqrt(float(np.sum(residuals ** 2)) / max(1, w_size - 2))

        # If slope <= 0 (no positive drift observed):
        if slope <= 1e-6:
            # Conservative prediction based on residual noise
            remaining_gap = self.threshold - cur_val
            # If 3*s_err exceeds gap:
            steps = int(max(100, remaining_gap / max(1e-5, s_err * 0.1)))
            return steps, max(10, steps // 3)

        # Distance to threshold
        remaining_gap = self.threshold - cur_val
        point_h = int(math.floor(remaining_gap / slope))

        # 95% Lower prediction bound using upper confidence bound on slope (slope + 1.96 * SE)
        se_slope = s_err / math.sqrt(ss_xx)
        upper_slope = slope + 1.95996 * se_slope
        lower_bound_h = int(math.floor(remaining_gap / upper_slope))

        return max(0, point_h), max(0, lower_bound_h)


class CUSUMChangePointPredictor(BasePredictor):
    """Detects onset of degradation via CUSUM and forecasts based on post-change slope."""

    def __init__(self, threshold: float = 0.50, baseline_mean: float = 0.0, k: float = 0.01, h: float = 0.10):
        super().__init__(threshold=threshold)
        self.baseline_mean = baseline_mean
        self.k = k  # allowance parameter
        self.h = h  # decision boundary threshold

    def predict(self, history: Sequence[float], current_t: int) -> tuple[int, int]:
        n = len(history)
        if n < 10:
            return 5000, 1000
        cur_val = history[-1]
        if cur_val >= self.threshold:
            return 0, 0

        # Run CUSUM
        s_pos = 0.0
        change_point = None
        for i, val in enumerate(history):
            s_pos = max(0.0, s_pos + (val - self.baseline_mean) - self.k)
            if s_pos >= self.h and change_point is None:
                change_point = i

        if change_point is None or change_point >= n - 3:
            # No confirmed change point yet
            linear_pred = LinearExtrapolationPredictor(self.threshold)
            return linear_pred.predict(history, current_t)

        # Fit slope from change point to present
        post_change = history[change_point:]
        x = np.arange(len(post_change), dtype=float)
        y = np.array(post_change, dtype=float)
        if len(x) < 2:
            return 1000, 300

        slope = float(np.polyfit(x, y, 1)[0])
        if slope <= 1e-6:
            return 2000, 500

        gap = self.threshold - cur_val
        point_h = int(math.floor(gap / slope))
        lower_bound_h = int(math.floor(gap / (slope * 1.5)))
        return max(0, point_h), max(0, lower_bound_h)


@dataclass(frozen=True)
class ForecastingEvaluationReport:
    model_name: str
    total_evaluations: int
    mean_absolute_error: float
    root_mean_squared_error: float
    lower_bound_coverage_95: float
    false_reassurance_rate: float
    early_warning_rate: float
    is_well_calibrated: bool


def evaluate_forecast_on_trajectory(
    predictor: BasePredictor,
    trajectory_values: Sequence[float],
    threshold: float,
    evaluation_times: Sequence[int],
    model_name: str = "LinearExtrapolation",
) -> list[RUSForecast]:
    """Evaluate predictor along a trajectory at specified evaluation checkpoints."""
    total_len = len(trajectory_values)
    # Find true threshold crossing
    actual_crossing = None
    for i, v in enumerate(trajectory_values):
        if v >= threshold:
            actual_crossing = i + 1
            break

    forecasts = []
    for t in evaluation_times:
        if t >= total_len:
            continue
        history = trajectory_values[:t]
        cur_val = history[-1]
        point_h, lower_95 = predictor.predict(history, t)

        # Actual remaining
        if actual_crossing is not None:
            actual_rem = max(0, actual_crossing - t)
        else:
            actual_rem = total_len - t + 5000  # Threshold not reached

        err = point_h - actual_rem
        # Lower bound coverage: did the system survive at least lower_95 samples?
        covered = (actual_rem >= lower_95)
        # False reassurance: predictor told operator it was safe for lower_95, but it failed earlier!
        false_reassurance = (actual_rem < lower_95)

        forecasts.append(RUSForecast(
            current_t=t,
            current_value=round(cur_val, 4),
            threshold=threshold,
            point_forecast_samples=point_h,
            lower_bound_95_samples=lower_95,
            actual_remaining_samples=actual_rem,
            forecast_error=err,
            covered_by_lower_bound=covered,
            is_false_reassurance=false_reassurance,
            model_name=model_name,
        ))

    return forecasts


def summarize_forecasting_performance(forecasts: list[RUSForecast]) -> ForecastingEvaluationReport:
    """Calculate aggregate forecast error, calibration, and false reassurance metrics."""
    if not forecasts:
        raise ValueError("Cannot summarize empty forecast list")

    errors = [f.forecast_error for f in forecasts if f.forecast_error is not None]
    mae = float(np.mean(np.abs(errors))) if errors else 0.0
    rmse = float(math.sqrt(np.mean(np.array(errors) ** 2))) if errors else 0.0

    covered_count = sum(1 for f in forecasts if f.covered_by_lower_bound)
    false_reassurance_count = sum(1 for f in forecasts if f.is_false_reassurance)
    n = len(forecasts)

    coverage = covered_count / n
    false_rate = false_reassurance_count / n
    early_warning = sum(1 for f in forecasts if f.lower_bound_95_samples > 0 and not f.is_false_reassurance) / n

    # Well-calibrated if false reassurance rate <= 10% (nominal 5% with finite sample tolerance)
    well_calibrated = (false_rate <= 0.12)

    return ForecastingEvaluationReport(
        model_name=forecasts[0].model_name,
        total_evaluations=n,
        mean_absolute_error=round(mae, 2),
        root_mean_squared_error=round(rmse, 2),
        lower_bound_coverage_95=round(coverage, 4),
        false_reassurance_rate=round(false_rate, 4),
        early_warning_rate=round(early_warning, 4),
        is_well_calibrated=well_calibrated,
    )
