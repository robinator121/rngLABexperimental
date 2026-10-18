"""RNGLab v0.14: Dual Source-Level and Application-Level Integrity Monitor.

Implements the 4-quadrant operational surveillance matrix:
- Q1 [Nominal]: Source Healthy & Application Acceptable
- Q2 [Early Warning]: Source Unhealthy & Application Acceptable (Pre-emptive intervention window)
- Q3 [Synchronized]: Source Unhealthy & Application Failed (Standard alarm trigger)
- Q4 [Observable Masking / Covert Hazard]: Source Apparently Healthy & Application Failed
  (Occurs when subtle correlations or unmeasured defect channels escape source health tests
   or when zero-variance observable masking hides physical corruption).
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
import math
from typing import Any, Dict, List, Optional, Sequence
from rnglab14.monitoring.calibration import FrozenRankCalibration


def evaluate_alarm_trace(flags: Sequence[bool], fault_onset: Optional[int],
                         failure_onset: Optional[int]) -> dict:
    """Evaluate one complete held-out trace against independent fault truth.

    Pre-fault alarms are recorded but never credited as true detections. The
    active-step hazardous fraction differs from a cumulative first-alert rule.
    Repeated steps here are not independent experimental units.
    """
    if not flags or not all(type(x) is bool for x in flags):
        raise ValueError("nonempty sequence of boolean alarm flags required")
    n=len(flags)
    if (fault_onset is None)!=(failure_onset is None):
        raise ValueError("provide both independent fault and failure onset, or neither for healthy")
    if fault_onset is not None and (type(fault_onset) is not int or type(failure_onset) is not int
                                   or not 0<=fault_onset<=failure_onset<n):
        raise ValueError("fault and failure onset must be ordered indices inside the trace")
    hits=[j for j,flag in enumerate(flags) if flag]
    first=None if not hits else hits[0]
    eligible=[] if fault_onset is None else [j for j in hits if j>=fault_onset]
    detection=None if not eligible else eligible[0]
    return {"first_alarm":first,"first_detection_from_fault_onset":detection,
        "pre_fault_or_healthy_false_alarm":bool(hits) if fault_onset is None else any(flags[:fault_onset]),
        "missed_failure":None if failure_onset is None else detection is None,
        "detection_delay_conditional_on_detection":None if detection is None else detection-failure_onset,
        "failure_steps_before_first_true_detection":None if failure_onset is None else
            max(0,(n if detection is None else detection)-failure_onset),
        "active_failure_steps_without_alarm":None if failure_onset is None else
            sum(not x for x in flags[failure_onset:]),
        "active_failure_step_hazard_fraction":None if failure_onset is None else
            sum(not x for x in flags[failure_onset:])/(n-failure_onset)}


class MonitoringQuadrant(Enum):
    Q1_NOMINAL = "Q1_NOMINAL"                               # Source OK, App OK
    Q2_EARLY_WARNING = "Q2_EARLY_WARNING"                   # Source Unhealthy, App OK
    Q3_SYNCHRONIZED_FAILURE = "Q3_SYNCHRONIZED_FAILURE"     # Source Unhealthy, App Failed
    Q4_HAZARDOUS_GAP = "Q4_HAZARDOUS_GAP"                   # Source OK, App Failed (Hazardous!)


@dataclass(frozen=True)
class DualMonitorSnapshot:
    step_t: int
    source_healthy: bool
    source_alarm_active: bool
    source_metric_value: float
    app_acceptable: bool
    app_failure_active: bool
    app_metric_value: float
    quadrant: MonitoringQuadrant
    detection_gap_steps: Optional[int]


@dataclass(frozen=True)
class DualMonitoringReport:
    total_snapshots: int
    q1_nominal_count: int
    q2_early_warning_count: int
    q3_synchronized_count: int
    q4_hazardous_gap_count: int
    hazardous_gap_rate: float
    early_warning_rate: float
    mean_early_warning_lead_time_steps: Optional[float]
    mean_detection_lag_steps: Optional[float]
    combined_monitoring_superior: Optional[bool]
    source_alarm_observed: bool
    application_failure_observed: bool
    comparison_status: str
    calibration_status: str


class DualIntegrityMonitor:
    """Synchronizes source-level health tests with application-level observable checks."""

    def __init__(
        self,
        source_alarm_threshold: float = 0.05,
        app_failure_threshold: float = 0.10,
        source_calibration: Optional[FrozenRankCalibration] = None,
        app_calibration: Optional[FrozenRankCalibration] = None,
        family_alpha: float = .05,
        horizon: Optional[int] = None,
    ):
        if (not math.isfinite(source_alarm_threshold) or source_alarm_threshold <= 0
                or not math.isfinite(app_failure_threshold) or app_failure_threshold <= 0):
            raise ValueError("finite positive thresholds required")
        if (source_calibration is None) != (app_calibration is None):
            raise ValueError("both channel calibrations are required")
        if not math.isfinite(family_alpha) or not 0<family_alpha<1:
            raise ValueError("family_alpha must be finite and between zero and one")
        if source_calibration is not None:
            if type(horizon) is not int or horizon < 1:
                raise ValueError("fixed positive horizon required for calibrated monitoring")
            import json
            if any(json.loads(c.protocol_json).get("horizon") != horizon
                   for c in (source_calibration,app_calibration)):
                raise ValueError("calibration horizon differs from monitoring horizon")
        self.source_calibration=source_calibration
        self.app_calibration=app_calibration
        self.family_alpha=family_alpha
        self.horizon=horizon
        self.source_thresh = source_alarm_threshold
        self.app_thresh = app_failure_threshold
        self.history: List[DualMonitorSnapshot] = []
        self._first_source_alarm_t: Optional[int] = None
        self._first_app_failure_t: Optional[int] = None

    def record_step(
        self,
        step_t: int,
        source_defect_magnitude: float,
        app_observable_error: float,
    ) -> DualMonitorSnapshot:
        if type(step_t) is not int or step_t < 0 or (self.history and step_t <= self.history[-1].step_t):
            raise ValueError("strictly increasing nonnegative integer step times required")
        if not all(math.isfinite(v) and v >= 0 for v in (source_defect_magnitude, app_observable_error)):
            raise ValueError("finite nonnegative monitoring magnitudes required")
        if self.source_calibration is not None:
            if len(self.history)>=self.horizon:
                raise ValueError("frozen calibrated monitoring horizon exhausted")
            # Whole-trajectory calibration maxima bound any threshold crossing
            # over this fixed horizon; Bonferroni allocates alpha/2 per channel.
            source_unhealthy=self.source_calibration.alarm(source_defect_magnitude,alpha=self.family_alpha/2)
            app_failed=self.app_calibration.alarm(app_observable_error,alpha=self.family_alpha/2)
        else:
            source_unhealthy = (source_defect_magnitude >= self.source_thresh)
            app_failed = (app_observable_error >= self.app_thresh)

        if source_unhealthy and self._first_source_alarm_t is None:
            self._first_source_alarm_t = step_t

        if app_failed and self._first_app_failure_t is None:
            self._first_app_failure_t = step_t

        # Determine quadrant
        if not source_unhealthy and not app_failed:
            quad = MonitoringQuadrant.Q1_NOMINAL
        elif source_unhealthy and not app_failed:
            quad = MonitoringQuadrant.Q2_EARLY_WARNING
        elif source_unhealthy and app_failed:
            quad = MonitoringQuadrant.Q3_SYNCHRONIZED_FAILURE
        else:
            quad = MonitoringQuadrant.Q4_HAZARDOUS_GAP

        # Detection gap calculation
        gap = None
        if self._first_app_failure_t is not None and self._first_source_alarm_t is not None:
            gap = self._first_app_failure_t - self._first_source_alarm_t

        snapshot = DualMonitorSnapshot(
            step_t=step_t,
            source_healthy=not source_unhealthy,
            source_alarm_active=source_unhealthy,
            source_metric_value=round(source_defect_magnitude, 4),
            app_acceptable=not app_failed,
            app_failure_active=app_failed,
            app_metric_value=round(app_observable_error, 4),
            quadrant=quad,
            detection_gap_steps=gap,
        )
        self.history.append(snapshot)
        return snapshot

    def summarize(self) -> DualMonitoringReport:
        if not self.history:
            raise ValueError("No monitoring history recorded")

        n = len(self.history)
        q1 = sum(1 for s in self.history if s.quadrant == MonitoringQuadrant.Q1_NOMINAL)
        q2 = sum(1 for s in self.history if s.quadrant == MonitoringQuadrant.Q2_EARLY_WARNING)
        q3 = sum(1 for s in self.history if s.quadrant == MonitoringQuadrant.Q3_SYNCHRONIZED_FAILURE)
        q4 = sum(1 for s in self.history if s.quadrant == MonitoringQuadrant.Q4_HAZARDOUS_GAP)

        lead_time = None
        lag_time = None

        if self._first_source_alarm_t is not None and self._first_app_failure_t is not None:
            diff = self._first_app_failure_t - self._first_source_alarm_t
            lead_time = 0.0
            lag_time = 0.0
            if diff > 0:
                lead_time = float(diff)  # Source gave warning ahead of application failure
            elif diff < 0:
                lag_time = float(-diff)  # Source was late; application failed first!
        hazard_rate = q4 / n
        early_rate = q2 / n

        return DualMonitoringReport(
            total_snapshots=n,
            q1_nominal_count=q1,
            q2_early_warning_count=q2,
            q3_synchronized_count=q3,
            q4_hazardous_gap_count=q4,
            hazardous_gap_rate=round(hazard_rate, 4),
            early_warning_rate=round(early_rate, 4),
            mean_early_warning_lead_time_steps=None if lead_time is None else round(lead_time, 2),
            mean_detection_lag_steps=None if lag_time is None else round(lag_time, 2),
            combined_monitoring_superior=None,
            source_alarm_observed=self._first_source_alarm_t is not None,
            application_failure_observed=self._first_app_failure_t is not None,
            comparison_status="DESCRIPTIVE_SINGLE_TRAJECTORY_NOT_A_SUPERIORITY_TEST",
            calibration_status=("UNCALIBRATED_THRESHOLDS" if self.source_calibration is None else
                "FROZEN_EXCHANGEABLE_TRAJECTORY_RANKS_BONFERRONI_FINITE_HORIZON"),
        )
