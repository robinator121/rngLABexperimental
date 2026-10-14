"""Bounded prospective monitoring with frozen, separate calibration histories.

Finite-window scores reuse RNGLab v0.3. Whole-history order calibration reuses
v0.7. Neither scores nor the prototype's default settings claim a measured FAR.
"""
from __future__ import annotations

from collections import deque
from copy import deepcopy
import hashlib
import json
import math
from typing import Any, Iterable

import numpy as np

from rnglab.v03.monitor import window_scores
from rnglab.v07.calibration import order_certificate, proportion
from .statistics import uniform_array, _alpha, _positive_int


RULES = ("marginal", "cross", "cusum", "likelihood", "combined")


def _digest(value: dict[str, Any]) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def _signed_likelihood(side: np.ndarray, delta: float) -> tuple[np.ndarray, np.ndarray]:
    return (np.where(side, math.log1p(2*delta), math.log1p(-2*delta)),
            np.where(side, math.log1p(-2*delta), math.log1p(2*delta)))


class StreamingMonitor:
    """Chunk-invariant fixed boundaries, bounded multiscale windows, and state.

    New CUSUM/likelihood scores target threshold-bit bias under an IID
    Bernoulli(.5) null and fixed alternatives .5+/-delta. They do not test every
    possible dependence defect. Alarm decisions need healthy *whole-history*
    calibration records under identical config, including evaluation horizon.
    """
    def __init__(self, window: int = 256, step: int = 64, channels: int = 1,
                 scales: Iterable[int] | None = None, calibration: dict[str, Any] | None = None,
                 cusum_delta: float = .05, horizon: int = 4096):
        self.window = _positive_int(window, "window", 8)
        self.step = _positive_int(step, "step")
        self.channels = _positive_int(channels, "channels")
        if self.channels > 32:
            raise ValueError("at most 32 channels are supported by this bounded monitor")
        self.scales = tuple(sorted(set(_positive_int(s, "scale", 8)
                                      for s in (scales if scales is not None else (self.window,)))))
        if not self.scales or any(s > 65536 for s in self.scales):
            raise ValueError("nonempty scales in [8,65536] required")
        if max(self.scales)*self.channels > 1_000_000:
            raise ValueError("window storage is limited to 1000000 numerical uniform values")
        if self.window not in self.scales:
            raise ValueError("window must appear among scales")
        self.horizon = _positive_int(horizon, "horizon", max(self.scales))
        if (isinstance(cusum_delta, bool) or not isinstance(cusum_delta, (int, float, np.number))
                or not math.isfinite(float(cusum_delta)) or not 0 < cusum_delta < .5):
            raise ValueError("cusum_delta must be finite in (0,.5)")
        self.delta = float(cusum_delta)
        self._buffer: deque[np.ndarray] = deque(maxlen=max(self.scales))
        self.rows = 0
        self._cusum_plus = np.zeros(self.channels)
        self._cusum_minus = np.zeros(self.channels)
        self._log_plus = np.zeros(self.channels)
        self._log_minus = np.zeros(self.channels)
        self.max_scores = {rule: 0. for rule in RULES}
        self.first_alerts: dict[str, int | None] = {rule: None for rule in RULES}
        self.calibration = deepcopy(calibration)
        if self.calibration is not None:
            body = {k: v for k, v in self.calibration.items() if k != "sha256"}
            if self.calibration.get("sha256") != _digest(body):
                raise ValueError("calibration digest mismatch")
            if self.calibration.get("monitor_config") != self.config():
                raise ValueError("calibration must match frozen monitor settings and horizon")
            if self.calibration.get("status") != "FROZEN_CALIBRATION":
                raise ValueError("insufficient calibration cannot authorize alarm decisions")
            thresholds = self.calibration.get("thresholds")
            rules = self.calibration.get("rules")
            if (not isinstance(thresholds, dict) or not thresholds or not isinstance(rules, list)
                    or len(set(rules)) != len(rules) or set(rules) != set(thresholds)
                    or set(rules)-set(RULES)):
                raise ValueError("calibration must contain explicit known frozen rules")
            if any(isinstance(v, bool) or not isinstance(v, (int, float))
                   or not math.isfinite(v) or v < 0 for v in thresholds.values()):
                raise ValueError("calibration cutoffs must be finite nonnegative numbers")

    def config(self) -> dict[str, Any]:
        return {"version": 1, "window": self.window, "step": self.step,
                "channels": self.channels, "scales": list(self.scales),
                "cusum_delta": self.delta, "horizon": self.horizon}

    def update(self, samples: Any) -> list[dict[str, Any]]:
        u = uniform_array(samples)
        if u.shape[1] != self.channels:
            raise ValueError("sample channel count differs from frozen monitor config")
        # Preflight an entire chunk so an over-horizon update is atomic.
        if self.rows + len(u) > self.horizon:
            raise ValueError("monitor horizon exhausted; no calibrated lifetime extension is permitted")
        events = []
        for row in u:
            inc_plus, inc_minus = _signed_likelihood(row >= .5, self.delta)
            self._cusum_plus = np.maximum(0., self._cusum_plus+inc_plus)
            self._cusum_minus = np.maximum(0., self._cusum_minus+inc_minus)
            self._log_plus += inc_plus
            self._log_minus += inc_minus
            self.rows += 1
            self._buffer.append(row.copy())
            if self.rows < max(self.scales) or self.rows % self.step:
                continue
            recent = np.asarray(self._buffer)
            per_scale = {str(scale): window_scores(recent[-scale:]) for scale in self.scales}
            marginal = max(v[name] for v in per_scale.values()
                           for name in ("monobit", "distribution", "serial", "temporal", "periodic"))
            cross = max(v["cross"] for v in per_scale.values())
            cusum = float(max(self._cusum_plus.max(), self._cusum_minus.max()))
            likelihood = float(max(0., np.max(np.logaddexp(self._log_plus, self._log_minus)-math.log(2))))
            scores = {"marginal": marginal, "cross": cross, "cusum": cusum, "likelihood": likelihood}
            scores["combined"] = max(scores.values())
            for rule, value in scores.items():
                self.max_scores[rule] = max(self.max_scores[rule], float(value))
            alerts = {}
            thresholds = {} if self.calibration is None else self.calibration["thresholds"]
            for rule, cutoff in thresholds.items():
                alerts[rule] = scores[rule] > cutoff
                if alerts[rule] and self.first_alerts[rule] is None:
                    self.first_alerts[rule] = self.rows
            events.append({"row": self.rows, "window_rows": len(recent), "scores": scores,
                "per_scale": per_scale, "alerts": alerts, "first_alerts": dict(self.first_alerts),
                "decision_status": "CALIBRATED_RULE" if self.calibration else "SCORES_ONLY_UNCALIBRATED",
                "calibration_sha256": None if self.calibration is None else self.calibration["sha256"]})
        return events

    def calibration_record(self, dataset_id: str, role: str, healthy: bool) -> dict[str, Any]:
        if not isinstance(dataset_id, str) or not dataset_id.strip() or role != "calibration" or healthy is not True:
            raise ValueError("named, explicitly healthy calibration datasets are required")
        if self.rows != self.horizon:
            raise ValueError("calibration records require the complete frozen monitoring horizon")
        return {"dataset_id": dataset_id, "role": role, "healthy": healthy,
                "monitor_config": self.config(), "max_scores": dict(self.max_scores)}

    def state_dict(self) -> dict[str, Any]:
        return {"version": 1, "config": self.config(), "calibration": deepcopy(self.calibration),
                "rows": self.rows, "buffer": [row.tolist() for row in self._buffer],
                "cusum_plus": self._cusum_plus.tolist(), "cusum_minus": self._cusum_minus.tolist(),
                "log_plus": self._log_plus.tolist(), "log_minus": self._log_minus.tolist(),
                "max_scores": dict(self.max_scores), "first_alerts": dict(self.first_alerts)}

    @classmethod
    def from_state_dict(cls, state: dict[str, Any]) -> "StreamingMonitor":
        if state.get("version") != 1:
            raise ValueError("unknown monitor state version")
        cfg = dict(state["config"])
        if cfg.pop("version", None) != 1:
            raise ValueError("unknown monitor config version")
        result = cls(**cfg, calibration=state["calibration"])
        rows = state["rows"]
        if isinstance(rows, bool) or not isinstance(rows, int) or not 0 <= rows <= result.horizon:
            raise ValueError("invalid row checkpoint")
        if rows:
            buffer = uniform_array(state["buffer"])
            if buffer.shape != (min(rows, max(result.scales)), result.channels):
                raise ValueError("invalid monitor checkpoint buffer")
            result._buffer.extend(row.copy() for row in buffer)
        elif state["buffer"]:
            raise ValueError("nonempty buffer with zero rows")
        for name in ("cusum_plus", "cusum_minus", "log_plus", "log_minus"):
            vector = np.asarray(state[name], float)
            if vector.shape != (result.channels,) or not np.isfinite(vector).all():
                raise ValueError("invalid sequential-score checkpoint")
            setattr(result, "_"+name, vector)
        for rule in RULES:
            value = state["max_scores"][rule]
            alert = state["first_alerts"][rule]
            if not math.isfinite(value) or value < 0 or (alert is not None and (not isinstance(alert, int) or not 1 <= alert <= rows)):
                raise ValueError("invalid score/alert checkpoint")
        result.rows = rows
        result.max_scores = dict(state["max_scores"])
        result.first_alerts = dict(state["first_alerts"])
        return result


def fit_calibration(records: Iterable[dict[str, Any]], alpha: float = .01,
                    rules: Iterable[str] = ("combined",), family_size: int | None = None) -> dict[str, Any]:
    """Freeze upper order-statistic cutoffs on separate healthy complete histories.

    Rank ceil((n+1)*(1-alpha/family_size)) gives a marginal predictive bound
    under exchangeable IID calibration/new healthy histories. Overlapping looks
    are already included in each history maximum. A population-tolerance bound
    is separately recorded; neither is a measured or exactly matched FAR.
    Too few histories returns INSUFFICIENT_CALIBRATION, never a default cutoff.
    """
    alpha = _alpha(alpha)
    selected = tuple(rules)
    if not selected or len(set(selected)) != len(selected) or set(selected)-set(RULES):
        raise ValueError("unique known rules are required")
    family = _positive_int(len(selected) if family_size is None else family_size, "family_size")
    if family < len(selected):
        raise ValueError("family_size must include all simultaneously calibrated rules")
    rows = list(records)
    if not rows or len(rows)>100000:
        raise ValueError("bounded nonempty calibration records are required")
    seen = set()
    config = rows[0].get("monitor_config")
    if not isinstance(config, dict):
        raise ValueError("frozen monitor config is required")
    check = dict(config)
    if check.pop("version", None) != 1:
        raise ValueError("unknown calibration config version")
    # Reuse monitor configuration validation before trusting supplied records.
    if StreamingMonitor(**check).config() != config:
        raise ValueError("noncanonical monitor config")
    values = {name: [] for name in selected}
    for row in rows:
        name = row.get("dataset_id")
        if not isinstance(name, str) or not name.strip() or name in seen:
            raise ValueError("unique dataset IDs are required")
        if row.get("role") != "calibration" or row.get("healthy") is not True or row.get("monitor_config") != config:
            raise ValueError("only separate healthy complete calibration histories with identical config are allowed")
        seen.add(name)
        for rule in selected:
            value = row["max_scores"][rule]
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value < 0:
                raise ValueError("finite nonnegative maximum scores are required")
            values[rule].append(value)
    n = len(rows)
    rank = math.ceil((n+1)*(1-alpha/family))
    sufficient = rank <= n
    thresholds = {name: float(np.sort(value)[rank-1]) for name, value in values.items()} if sufficient else {}
    result = {"version": 1, "status": "FROZEN_CALIBRATION" if sufficient else "INSUFFICIENT_CALIBRATION",
              "monitor_config": config, "dataset_ids": sorted(seen), "training_histories": n,
              "alpha_per_history_family": alpha, "family_size": family, "rules": list(selected),
              "one_based_rank": rank, "thresholds": thresholds,
              "false_alarm_status": "PROSPECTIVE_BOUND_ONLY; matched realized FAR unmeasured",
              "certificate": order_certificate(n, rank, target=alpha/family, family_size=family) if sufficient else None,
              "assumptions": "IID independent complete healthy histories and frozen settings/horizon; no development/confirmation reuse"}
    result["sha256"] = _digest(result)
    return result


def detection_metrics(records: Iterable[dict[str, Any]], confidence: float = .95) -> dict[str, Any]:
    """Exact counts, pointwise binomial intervals, and observed/censored latency.

    Records: end (inclusive), onset (first defective row or None), alerts (row
    list), severity, optional failure (application failure row). An alert at the
    same row as failure is not prior warning. Misses are right-censored, not
    assigned invented detection times or included in mean observed latency.
    """
    confidence = _alpha(confidence)
    rows = list(records)
    if not rows:
        raise ValueError("at least one detector history is required")
    healthy, healthy_alarm, defect, detected, failures_without_prior = 0, 0, 0, 0, 0
    false_alarms, histories_with_false_alarm = 0, 0
    detail = []
    for row in rows:
        end = _positive_int(row["end"], "end")
        onset, failure = row.get("onset"), row.get("failure")
        if onset is not None and (_positive_int(onset, "onset") > end):
            raise ValueError("onset cannot exceed observed horizon")
        if failure is not None and (_positive_int(failure, "failure") > end):
            raise ValueError("failure cannot exceed observed horizon")
        alerts = sorted(set(_positive_int(v, "alert") for v in row.get("alerts", [])))
        if alerts and alerts[-1] > end:
            raise ValueError("alerts must lie inside observed history")
        false = [t for t in alerts if onset is None or t < onset]
        false_alarms += len(false)
        histories_with_false_alarm += int(bool(false))
        if onset is None:
            healthy += 1
            healthy_alarm += int(bool(alerts))
            first = latency = censored_time = None
        else:
            defect += 1
            first = next((t for t in alerts if t >= onset), None)
            detected += int(first is not None)
            latency = None if first is None else first-onset
            censored_time = end-onset if first is None else None
        no_warning = failure is not None and not any(t < failure for t in alerts)
        failures_without_prior += int(no_warning)
        severity = row.get("severity")
        if isinstance(severity, float) and not math.isfinite(severity):
            raise ValueError("severity cannot be nonfinite")
        detail.append({"dataset_id": row.get("dataset_id"), "severity": severity,
            "first_post_onset_alert": first, "detection_latency": latency,
            "right_censored": onset is not None and first is None, "censored_followup": censored_time,
            "false_alarm_events": len(false), "application_failure_without_prior_alert": no_warning})
    return {"histories": len(rows), "defect_histories": defect, "detections": detected,
            "missed_events": defect-detected, "false_alarm_events": false_alarms,
            "histories_with_false_alarm": histories_with_false_alarm,
            "detection_probability": proportion(detected, defect, confidence) if defect else None,
            "healthy_history_false_alarm_probability": proportion(healthy_alarm, healthy, confidence) if healthy else None,
            "application_failures_without_prior_alert": failures_without_prior, "events": detail,
            "interval_scope": "Pointwise exact binomial intervals across independent histories; not simultaneous or time-uniform"}
