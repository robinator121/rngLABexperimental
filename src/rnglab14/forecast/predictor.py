"""Research horizon predictors and honest pre-event evaluation.

Historical lower_bound_95 field names are retained for compatibility. Heuristic
slope bounds are NOT calibrated95% predictions. Censoring can leave outcomes
unresolved; abstention is not a successful forecast; no automatic certification.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
import math
from typing import List, Optional, Sequence, Tuple

import numpy as np
from scipy.stats import t as student_t


@dataclass(frozen=True)
class RUSForecast:
    current_t: int
    current_value: float
    threshold: float
    point_forecast_samples: Optional[int]
    lower_bound_95_samples: Optional[int]
    actual_remaining_samples: Optional[int]
    is_censored: bool
    abstained: bool
    forecast_error: Optional[int]
    covered_by_lower_bound: Optional[bool]
    is_false_reassurance: Optional[bool]
    model_name: str
    status_note: str
    observed_followup_samples: int = 0
    trajectory_id: Optional[str] = None
    bound_assurance: str = "HEURISTIC_NOT_CALIBRATED"


class BasePredictor(ABC):
    bound_assurance = "HEURISTIC_NOT_CALIBRATED"
    """Abstract base predictor for Remaining Useful Samples."""

    def __init__(self, threshold: float = 0.50, max_horizon: int = 2000):
        if not math.isfinite(threshold) or isinstance(max_horizon,bool) or not isinstance(max_horizon,int) or max_horizon<0:
            raise ValueError("finite threshold and nonnegative integer horizon required")
        self.threshold = threshold
        self.max_horizon = max_horizon

    @abstractmethod
    def predict(self, history: Sequence[float], current_t: int) -> Tuple[Optional[int], Optional[int], bool]:
        """Return (point_forecast, lower_bound_95, abstained)."""
        raise NotImplementedError


class NoChangeBaselinePredictor(BasePredictor):
    """Null baseline: assumes current degradation state remains invariant."""

    def predict(self, history: Sequence[float], current_t: int) -> Tuple[Optional[int], Optional[int], bool]:
        if not history:
            return None, None, True
        val = history[-1]
        if val >= self.threshold:
            return 0, 0, False
        # If no change is assumed and current value is below threshold,
        # failure will never occur under the model assumption.
        # Clamped to max_horizon.
        return self.max_horizon, self.max_horizon // 2, False


class LinearExtrapolationPredictor(BasePredictor):
    """Linear slope heuristic with t slope screen; NOT a first-passage prediction interval."""

    def __init__(self, threshold: float = 0.50, window: int = 50, max_horizon: int = 2000):
        super().__init__(threshold=threshold, max_horizon=max_horizon)
        if isinstance(window,bool) or not isinstance(window,int) or window<10:
            raise ValueError("window must be an integer >=10")
        self.window = window

    def predict(self, history: Sequence[float], current_t: int) -> Tuple[Optional[int], Optional[int], bool]:
        n_hist = len(history)
        if n_hist < 10:
            return None, None, True  # Insufficient data to fit trend

        w_size = min(self.window, n_hist)
        y = np.array(history[-w_size:], dtype=float)
        x = np.arange(w_size, dtype=float)

        cur_val = y[-1]
        if cur_val >= self.threshold:
            return 0, 0, False

        # Fit OLS
        x_mean = np.mean(x)
        y_mean = np.mean(y)
        ss_xx = float(np.sum((x - x_mean) ** 2))
        if ss_xx < 1e-12:
            return None, None, True

        slope = float(np.sum((x - x_mean) * (y - y_mean)) / ss_xx)
        intercept = float(y_mean - slope * x_mean)

        residuals = y - (slope * x + intercept)
        dof = max(1, w_size - 2)
        s_err = math.sqrt(float(np.sum(residuals ** 2)) / dof)
        se_slope = s_err / math.sqrt(ss_xx)

        if se_slope > 1e-12:
            t_stat = slope / se_slope
        else:
            t_stat = float("inf") if slope > 1e-8 else 0.0

        # Uncertainty abstention: If slope is not statistically significantly positive (t < 1.96),
        # trend is indistinguishable from noise -> ABSTAIN!
        critical = float(student_t.ppf(.975, dof))
        if t_stat < critical:
            return None, None, True

        gap = self.threshold - cur_val
        point_h = int(math.floor(gap / slope))
        upper_slope = slope + critical * se_slope
        lower_bound_h = int(math.floor(gap / upper_slope)) if upper_slope > 0 else 0

        # Clamp to max horizon
        point_h = min(self.max_horizon, max(0, point_h))
        lower_bound_h = min(self.max_horizon, max(0, lower_bound_h))

        return point_h, lower_bound_h, False


class CUSUMChangePointPredictor(BasePredictor):
    """CUSUM change-point detector with post-change slope estimation and uncertainty abstention."""

    def __init__(
        self,
        threshold: float = 0.50,
        baseline_mean: float = 0.0,
        k: float = 0.01,
        h: float = 0.05,
        max_horizon: int = 2000,
    ):
        super().__init__(threshold=threshold, max_horizon=max_horizon)
        if not all(math.isfinite(v) for v in [baseline_mean,k,h]) or k<0 or h<=0:
            raise ValueError("invalid finite CUSUM parameters")
        self.baseline_mean = baseline_mean
        self.k = k
        self.h = h

    def predict(self, history: Sequence[float], current_t: int) -> Tuple[Optional[int], Optional[int], bool]:
        n = len(history)
        if n < 10:
            return None, None, True
        cur_val = history[-1]
        if cur_val >= self.threshold:
            return 0, 0, False

        # Run CUSUM
        s_pos = 0.0
        change_point = None
        for i, val in enumerate(history):
            s_pos = max(0.0, s_pos + (val - self.baseline_mean) - self.k)
            if s_pos >= self.h and change_point is None:
                change_point = i

        if change_point is None or change_point >= n - 5:
            # No change point confirmed: state remains in nominal baseline
            # Under nominal baseline (mean=0.0 < threshold), no imminent failure
            return self.max_horizon, self.max_horizon // 2, False

        # Fit post-change slope
        post_change = history[change_point:]
        m_len = len(post_change)
        x = np.arange(m_len, dtype=float)
        y = np.array(post_change, dtype=float)

        x_mean = np.mean(x)
        y_mean = np.mean(y)
        ss_xx = float(np.sum((x - x_mean) ** 2))
        if ss_xx < 1e-12:
            return None, None, True

        slope = float(np.sum((x - x_mean) * (y - y_mean)) / ss_xx)
        if slope <= 1e-5:
            # Change point was a transient shock that flattened out -> Abstain or claim safe
            return self.max_horizon, self.max_horizon // 2, False

        gap = self.threshold - cur_val
        point_h = min(self.max_horizon, max(0, int(math.floor(gap / slope))))
        lower_bound_h = min(self.max_horizon, max(0, int(math.floor(gap / (slope * 1.5)))))

        return point_h, lower_bound_h, False


class LipschitzStructuralBoundPredictor(BasePredictor):
    """Mathematical lower bound predictor under known Lipschitz drift rate bound L (MATH-04)."""

    def __init__(self, threshold: float = 0.50, lipschitz_L: float = 0.001, error_margin: float = 0.01):
        super().__init__(threshold=threshold)
        if not math.isfinite(lipschitz_L) or lipschitz_L<=0 or not math.isfinite(error_margin) or error_margin<0:
            raise ValueError("positive finite L and nonnegative finite state margin required")
        self.L = lipschitz_L
        self.eps = error_margin
        self.bound_assurance = "CONDITIONAL_ON_VALID_FUTURE_L_AND_STATE_MARGIN"

    def predict(self, history: Sequence[float], current_t: int) -> Tuple[Optional[int], Optional[int], bool]:
        if not history:
            return None, None, True
        cur_val = history[-1]
        if cur_val >= self.threshold:
            return 0, 0, False

        # N_delta(t) = floor((threshold - cur_val - eps) / L)
        remaining_gap = self.threshold - cur_val - self.eps
        if remaining_gap <= 0:
            return 0, 0, False

        bound_steps = int(math.floor(remaining_gap / self.L))
        point_h = int(math.floor((self.threshold - cur_val) / self.L))
        return point_h, bound_steps, False


class FixedHorizonPredictor(BasePredictor):
    """Noninformative, fixed useful-horizon baseline, not an assurance claim."""
    def __init__(self, horizon=50, threshold=.5):
        super().__init__(threshold=threshold, max_horizon=horizon)
    def predict(self, history, current_t):
        return (self.max_horizon,self.max_horizon,False) if len(history) else (None,None,True)


@dataclass(frozen=True)
class ForecastingEvaluationReport:
    model_name: str
    total_evaluations: int
    uncensored_events_count: int  # Historical name: selected event LANDMARKS
    censored_events_count: int
    abstention_count: int
    abstention_rate: float
    event_mae: Optional[float]
    event_rmse: Optional[float]
    lower_bound_coverage_95: Optional[float]  # Descriptive, resolved landmarks only
    false_reassurance_count: int
    false_reassurance_rate: Optional[float]
    is_well_calibrated: Optional[bool]
    unresolved_bound_count: int = 0
    resolved_bound_count: int = 0
    emitted_count: int = 0
    unique_trajectory_ids_count: Optional[int] = None
    calibration_status: str = "NOT_ESTABLISHED; repeated landmarks are not independent trials"


def evaluate_forecast_on_trajectory(predictor, trajectory_values, threshold, evaluation_times,
                                    model_name="Predictor", trajectory_id=None):
    """Evaluate only pre-first-event landmarks.

    Bound semantics: remaining event time R>=h; a censored followup c supports
    h<=c, otherwise is UNRESOLVED. R is undefined for censored trajectories.
    This is not a model of survival probabilities and supplies no Brier score.
    """
    values=list(trajectory_values)
    if not values or not all(math.isfinite(float(v)) for v in values):
        raise ValueError("trajectory must be nonempty and finite")
    if not math.isfinite(threshold) or predictor.threshold!=threshold:
        raise ValueError("evaluation and predictor thresholds must agree")
    total_len=len(values)
    crossing=next((i+1 for i,v in enumerate(values) if v>=threshold),None)
    forecasts=[]
    for at in evaluation_times:
        if isinstance(at,bool) or not isinstance(at,int) or not 1<=at<=total_len:
            raise ValueError("landmark must be an integer in the observed trajectory")
        if crossing is not None and crossing<=at:
            continue  # Detection after event is not a forecast.
        point,lower,abstain=predictor.predict(values[:at],at)
        if not isinstance(abstain,(bool,np.bool_)):
            raise ValueError("abstention must be boolean")
        if abstain and (point is not None or lower is not None):
            raise ValueError("abstentions must not emit point/bound values")
        if not abstain and any(isinstance(v,(bool,np.bool_)) or not isinstance(v,(int,np.integer)) or v<0 for v in [point,lower]):
            raise ValueError("emitted forecasts require finite nonnegative integer horizons")
        if not abstain and lower>point:
            raise ValueError("lower bound exceeds point horizon")
        followup=total_len-at
        remaining=None if crossing is None else crossing-at
        covered=unsafe=error=None
        if abstain:
            note="ABSTAINED_UNCERTAINTY"
        elif remaining is not None:
            error=point-remaining;covered=remaining>=lower;unsafe=not covered
            note="EVENT_BOUND_SUPPORTED" if covered else "EVENT_BOUND_VIOLATED"
        elif lower<=followup:
            covered=True;unsafe=False;note="CENSORED_BOUND_SUPPORTED"
        else:
            note="CENSORED_UNRESOLVED"
        forecasts.append(RUSForecast(at,float(values[at-1]),threshold,point,lower,remaining,
            crossing is None,bool(abstain),error,covered,unsafe,model_name,note,
            followup,trajectory_id,predictor.bound_assurance))
    return forecasts


def summarize_forecasting_performance(forecasts):
    if not forecasts: raise ValueError("no eligible pre-event forecasts")
    names={f.model_name for f in forecasts}
    if len(names)!=1:raise ValueError("summarize one predictor at a time")
    total=len(forecasts);abstained=sum(f.abstained for f in forecasts)
    emitted=[f for f in forecasts if not f.abstained]
    event=[f for f in emitted if f.forecast_error is not None]
    resolved=[f for f in emitted if f.covered_by_lower_bound is not None]
    errors=[f.forecast_error for f in event]
    unsafe=sum(bool(f.is_false_reassurance) for f in resolved)
    ids={f.trajectory_id for f in forecasts}
    return ForecastingEvaluationReport(
        forecasts[0].model_name,total,len(event),sum(f.is_censored for f in forecasts),
        abstained,abstained/total,
        float(np.mean(np.abs(errors))) if errors else None,
        float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
        sum(bool(f.covered_by_lower_bound) for f in resolved)/len(resolved) if resolved else None,
        unsafe,unsafe/len(resolved) if resolved else None,None,
        len(emitted)-len(resolved),len(resolved),len(emitted),
        len(ids) if None not in ids else None)
