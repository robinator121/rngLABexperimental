"""Application-aware event times and forecast evaluation on independent runs.

Event times are observed decision boundaries. A sustained crossing is dated at
the confirming observation, never backdated. This module computes estimands;
it neither fits a model nor authenticates the calibration references supplied
by an experiment operator.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import islice
import math
from typing import Any, Iterable

from scipy.stats import binomtest


UNITS = frozenset({"raw_bits", "conditioned_bits", "generated_uniforms", "drbg_bytes",
                   "generator_calls", "simulation_steps", "effective_samples"})
VALIDATED_FORECAST_STATUSES = frozenset({"VALIDATED", "VALIDATED_PROTOTYPE"})
SOURCE_CALIBRATION_STATUSES = frozenset({"FROZEN_CALIBRATION", "CALIBRATED_RULE"})
FORECAST_STATUSES = VALIDATED_FORECAST_STATUSES | frozenset({"NO_PREDICTION", "INSUFFICIENT_EVIDENCE",
                                                          "EXPLORATORY", "EXPLORATORY_CALIBRATED"})


def _integer(value: int, name: str, minimum: int = 0) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _finite(value: float, name: str, minimum: float | None = None) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite numeric data")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return float(value)


def _identity(value: str, name: str) -> None:
    if not isinstance(value, str) or not value.strip() or len(value) > 512:
        raise ValueError(f"a bounded nonempty {name} is required")


@dataclass(frozen=True)
class ApplicationObservation:
    time: int
    error: float
    source_alarm: bool = False
    source_calibration_id: str | None = None
    source_calibration_status: str = "UNCALIBRATED"

    def __post_init__(self):
        _integer(self.time, "observation time")
        _finite(self.error, "absolute application error", 0.)
        if type(self.source_alarm) is not bool:
            raise ValueError("source_alarm must be boolean")
        if self.source_calibration_id is not None:
            _identity(self.source_calibration_id, "source calibration ID")
        _identity(self.source_calibration_status, "source calibration status")


@dataclass(frozen=True)
class RiskForecast:
    issued_at: int
    horizon: int
    probability: float | None = None
    status: str = "NO_PREDICTION"
    model_id: str | None = None
    calibration_id: str | None = None
    validation_id: str | None = None
    unit: str = "simulation_steps"
    event_definition_id: str | None = None
    reason: str | None = None
    risk_upper: float | None = None
    risk_scope: str = "application_horizon"
    adversarial_model_id: str | None = None

    def __post_init__(self):
        _integer(self.issued_at, "forecast issue time")
        _integer(self.horizon, "forecast horizon", 1)
        if self.status not in FORECAST_STATUSES or self.unit not in UNITS:
            raise ValueError("known forecast status and explicit operation unit are required")
        for field in ("model_id", "calibration_id", "validation_id", "event_definition_id", "adversarial_model_id"):
            value = getattr(self, field)
            if value is not None:
                _identity(value, field)
        if self.probability is not None and not 0 <= _finite(self.probability, "forecast probability") <= 1:
            raise ValueError("forecast probability must lie in [0,1]")
        if self.risk_upper is not None:
            if not 0 <= _finite(self.risk_upper, "forecast risk upper bound") <= 1:
                raise ValueError("risk upper bound must lie in [0,1]")
            if self.probability is None or self.risk_upper < self.probability:
                raise ValueError("upper risk bound requires a point probability no greater than the bound")
        if self.status in VALIDATED_FORECAST_STATUSES and self.probability is None:
            raise ValueError("a validated risk forecast requires an actual probability")
        if self.status in {"NO_PREDICTION", "INSUFFICIENT_EVIDENCE"}:
            if self.probability is not None:
                raise ValueError("an abstaining forecast cannot simultaneously assert a probability")
            _identity(self.reason, "abstention reason")
        if self.risk_scope not in {"application_horizon", "residual_tail_proxy", "unvalidated_model"}:
            raise ValueError("unknown risk estimand scope")


@dataclass(frozen=True)
class ResearchRun:
    run_id: str
    independent_group_id: str
    application: str
    unit: str
    threshold: float
    end: int
    observations: tuple[ApplicationObservation, ...]
    forecasts: tuple[RiskForecast, ...] = ()
    failure_consecutive: int = 2
    recovery_consecutive: int = 2
    onset: int | None = None
    source_restore: int | None = None
    healthy_control: bool = False
    degradation_kind: str = "stochastic"
    adversarial_model_id: str | None = None
    pair_id: str | None = None
    scenario_id: str = "default"
    event_definition_id: str = "absolute_error_sustained/v1"

    def __post_init__(self):
        for field in ("run_id", "independent_group_id", "application", "scenario_id", "event_definition_id"):
            _identity(getattr(self, field), field)
        if self.unit not in UNITS:
            raise ValueError("an explicit supported operation unit is required")
        _finite(self.threshold, "absolute error threshold", 0.)
        _integer(self.end, "observed end")
        _integer(self.failure_consecutive, "failure_consecutive", 1)
        _integer(self.recovery_consecutive, "recovery_consecutive", 1)
        if type(self.healthy_control) is not bool:
            raise ValueError("healthy_control must be boolean")
        if self.degradation_kind not in {"stochastic", "adversarial", "unknown"}:
            raise ValueError("degradation_kind must distinguish stochastic/adversarial/unknown")
        if self.adversarial_model_id is not None:
            _identity(self.adversarial_model_id, "adversarial model ID")
        if self.pair_id is not None:
            _identity(self.pair_id, "matched pair ID")
        for field in ("onset", "source_restore"):
            value = getattr(self, field)
            if value is not None and _integer(value, field) > self.end:
                raise ValueError(f"{field} cannot exceed observed end")
        if self.healthy_control and self.onset is not None:
            raise ValueError("a healthy negative control cannot declare a defect onset")
        if self.source_restore is not None and self.onset is not None and self.source_restore < self.onset:
            raise ValueError("source restoration cannot precede defect onset")
        if not isinstance(self.observations, tuple) or not self.observations or len(self.observations) > 1_000_000:
            raise ValueError("freeze a nonempty bounded tuple of observations")
        if any(not isinstance(o, ApplicationObservation) for o in self.observations):
            raise ValueError("typed ApplicationObservation records are required")
        times = [o.time for o in self.observations]
        if any(b <= a for a, b in zip(times, times[1:])) or times[-1] != self.end:
            raise ValueError("observation times must strictly increase and end at the declared observed horizon")
        if not isinstance(self.forecasts, tuple) or len(self.forecasts) > 100_000:
            raise ValueError("freeze a bounded tuple of forecasts")
        seen = set()
        for f in self.forecasts:
            if not isinstance(f, RiskForecast) or f.issued_at > self.end:
                raise ValueError("typed forecasts must be issued inside the observed horizon")
            key = (f.issued_at, f.horizon)
            if key in seen:
                raise ValueError("one frozen forecast per issue time and horizon is required")
            seen.add(key)


def _sustained_time(observations, length, predicate, start=None):
    streak = 0
    for observation in observations:
        if start is not None and observation.time < start:
            continue
        streak = streak + 1 if predicate(observation) else 0
        if streak >= length:
            return observation.time
    return None


def _forecast_eligibility(run: ResearchRun, forecast: RiskForecast) -> tuple[bool, str]:
    if forecast.status not in VALIDATED_FORECAST_STATUSES:
        return False, forecast.reason or "Forecast has no validated application-horizon status"
    if not forecast.model_id or not forecast.calibration_id or not forecast.validation_id:
        return False, "Missing frozen model/calibration/independent-validation provenance"
    if forecast.unit != run.unit or forecast.event_definition_id != run.event_definition_id:
        return False, "Forecast unit or failure event definition differs from evaluated application"
    if forecast.risk_scope != "application_horizon":
        return False, "A residual-tail proxy or unvalidated model is not application failure probability"
    if run.degradation_kind == "unknown":
        return False, "NO_PREDICTION: degradation population is unspecified"
    if run.degradation_kind == "adversarial":
        if not run.adversarial_model_id or forecast.adversarial_model_id != run.adversarial_model_id:
            return False, "NO_PREDICTION: arbitrary adversarial manipulation has no assumed predictable trend"
    return True, "Eligible under declared matching validation population; provenance is not independently authenticated by this module"


def evaluate_run(run: ResearchRun, warning_probability: float = .5) -> dict[str, Any]:
    """Event times, signed warning lead, censoring and fully resolved horizons.

    t_A is the first sustained error>threshold confirmation; t_S is the first
    source alarm with frozen calibration provenance. t_P is a validated risk
    warning's issue time. Positive lead=t_A-t_warning is a prior warning; an
    equal-time alert is not prior. No future event time is inferred at censoring.
    """
    if not isinstance(run, ResearchRun):
        raise ValueError("a typed ResearchRun is required")
    if not 0 <= _finite(warning_probability, "warning probability") <= 1:
        raise ValueError("warning_probability must lie in [0,1]")
    failure = _sustained_time(run.observations, run.failure_consecutive, lambda o: o.error > run.threshold)
    source_valid = [o.time for o in run.observations if o.source_alarm and o.source_calibration_id
                    and o.source_calibration_status in SOURCE_CALIBRATION_STATUSES]
    valid_source_times = set(source_valid)
    source_invalid = [o.time for o in run.observations if o.source_alarm and o.time not in valid_source_times]
    source = min(source_valid, default=None)
    source_post_onset = (min((time for time in source_valid if time >= run.onset), default=None)
                         if run.onset is not None else None)
    warnings, forecasts = [], []
    for forecast in sorted(run.forecasts, key=lambda f: (f.issued_at, f.horizon)):
        eligible, reason = _forecast_eligibility(run, forecast)
        prior_failure = failure is not None and failure <= forecast.issued_at
        warning = eligible and not prior_failure and forecast.probability >= warning_probability
        if warning:
            warnings.append(forecast.issued_at)
        horizon_end = forecast.issued_at + forecast.horizon
        if prior_failure:
            outcome_status, outcome = "NOT_AT_RISK", None
        elif failure is not None and forecast.issued_at < failure <= horizon_end:
            outcome_status, outcome = "OBSERVED_EVENT", 1
        elif run.end >= horizon_end:
            outcome_status, outcome = "FULL_HORIZON_NO_EVENT", 0
        else:
            outcome_status, outcome = "RIGHT_CENSORED", None
        forecasts.append({"issued_at": forecast.issued_at, "horizon": forecast.horizon,
            "horizon_end": horizon_end, "probability": forecast.probability,
            "risk_upper": forecast.risk_upper, "eligible": eligible, "eligibility_reason": reason,
            "warning": bool(warning), "forecast_status": forecast.status,
            "model_id": forecast.model_id, "calibration_id": forecast.calibration_id,
            "validation_id": forecast.validation_id, "outcome_status": outcome_status,
            "outcome": outcome, "brier": (forecast.probability-outcome)**2 if eligible and outcome is not None else None})
    prediction = min(warnings, default=None)
    recovery = None
    if failure is not None and run.source_restore is not None:
        # Healthy errors before failure cannot constitute application recovery.
        recovery = _sustained_time(run.observations, run.recovery_consecutive,
                                  lambda o: o.error <= run.threshold,
                                  start=max(run.source_restore, failure))
    return {"run_id": run.run_id, "independent_group_id": run.independent_group_id,
        "pair_id": run.pair_id, "scenario_id": run.scenario_id, "application": run.application,
        "unit": run.unit, "event_definition_id": run.event_definition_id,
        "threshold": run.threshold, "end": run.end, "onset": run.onset,
        "healthy_control": run.healthy_control, "degradation_kind": run.degradation_kind,
        "t_A": failure, "t_S": source, "t_P": prediction, "t_restore": run.source_restore, "t_R": recovery,
        "application_failure_status": "OBSERVED" if failure is not None else "RIGHT_CENSORED",
        "application_failure_followup": run.end,
        "source_alarm_status": "OBSERVED_CALIBRATED" if source is not None else "RIGHT_CENSORED_OR_UNAVAILABLE",
        "prediction_status": "VALIDATED_WARNING_ISSUED" if prediction is not None else "NO_VALIDATED_WARNING",
        "source_warning_lead": failure-source if failure is not None and source is not None else None,
        "prediction_warning_lead": failure-prediction if failure is not None and prediction is not None else None,
        "source_detection_delay_from_onset": source-run.onset if source is not None and run.onset is not None else None,
        "first_post_onset_source_alarm": source_post_onset,
        "post_onset_source_detection_delay": source_post_onset-run.onset if source_post_onset is not None else None,
        "failure_after_onset": failure-run.onset if failure is not None and run.onset is not None else None,
        "pre_onset_application_failure": failure is not None and run.onset is not None and failure < run.onset,
        "pre_onset_source_alarm": source is not None and run.onset is not None and source < run.onset,
        "failure_without_prior_source_alert": failure is not None and (source is None or source >= failure),
        "failure_without_post_onset_source_alert":
            (source_post_onset is None or source_post_onset >= failure)
            if failure is not None and run.onset is not None and failure >= run.onset else None,
        "failure_without_prior_prediction_warning": failure is not None and (prediction is None or prediction >= failure),
        "recovery_status": "OBSERVED" if recovery is not None else
                           ("RIGHT_CENSORED" if failure is not None and run.source_restore is not None else "NOT_APPLICABLE"),
        "recovery_delay_after_source_restoration": recovery-run.source_restore if recovery is not None else None,
        "source_recovery_is_application_recovery": False,
        "ignored_uncalibrated_source_alarm_times": source_invalid, "forecasts": forecasts,
        "scientific_scope": "Observed fixed-schedule event estimands; calibrated provenance declarations require independent audit"}


def _proportion(count, n, confidence):
    if not n:
        return None
    ci = binomtest(count, n).proportion_ci(confidence_level=confidence, method="exact")
    return {"count": count, "n": n, "rate": count/n, "low": float(ci.low), "high": float(ci.high),
            "confidence": confidence, "scope": "pointwise across independent trajectories"}


def aggregate_runs(runs: Iterable[ResearchRun], warning_probability: float = .5, *,
                   forecast_issued_at: int | None = None, forecast_horizon: int | None = None,
                   confidence: float = .95) -> dict[str, Any]:
    """One independent trajectory per denominator; no window-level replication.

    Forecast scoring requires a prespecified common issue time and horizon.
    Multiple correlated runs within one independent group must be analyzed by
    scenario separately or using paired_differences; they cannot be pooled.
    Right-censored horizons stay excluded, not relabeled as nonfailures. Brier
    on resolved outcomes is descriptive when censoring can be informative.
    """
    rows = list(islice(runs, 100_001))
    if not rows or len(rows) > 100_000:
        raise ValueError("a bounded nonempty collection of independent runs is required")
    if not 0 < _finite(confidence, "confidence") < 1:
        raise ValueError("confidence must lie in (0,1)")
    if any(not isinstance(r, ResearchRun) for r in rows):
        raise ValueError("typed ResearchRun records required")
    ids = [r.run_id for r in rows]
    groups = [r.independent_group_id for r in rows]
    if len(set(ids)) != len(ids) or len(set(groups)) != len(groups):
        raise ValueError("pseudoreplication: one run per independent group; separate scenarios or use paired_differences")
    pair_ids = [r.pair_id for r in rows if r.pair_id is not None]
    if len(set(pair_ids)) != len(pair_ids):
        raise ValueError("pseudoreplication: the same matched pair cannot be relabeled as different independent groups")
    definitions = {(r.application, r.unit, r.threshold, r.end, r.failure_consecutive, r.recovery_consecutive,
                    r.event_definition_id, r.scenario_id, tuple(o.time for o in r.observations)) for r in rows}
    if len(definitions) != 1:
        raise ValueError("aggregate only a frozen common application, event definition, units, and scenario")
    if (forecast_issued_at is None) != (forecast_horizon is None):
        raise ValueError("prespecify both forecast issue time and horizon")
    if forecast_issued_at is not None:
        _integer(forecast_issued_at, "forecast_issued_at")
        _integer(forecast_horizon, "forecast_horizon", 1)
    evaluated = [evaluate_run(r, warning_probability) for r in rows]
    failures = [r for r in evaluated if r["t_A"] is not None]
    healthy = [r for r in evaluated if r["healthy_control"]]
    selected, missing = [], 0
    if forecast_issued_at is not None:
        for run in evaluated:
            found = [f for f in run["forecasts"] if f["issued_at"] == forecast_issued_at and f["horizon"] == forecast_horizon]
            if found:
                selected.append(found[0])
            else:
                missing += 1
    eligible = [f for f in selected if f["eligible"]]
    resolved = [f for f in eligible if f["outcome"] is not None]
    return {"independent_trajectories": len(rows), "windows_as_independent_units": False,
        "application_failure_events": len(failures),
        "application_right_censored": len(rows)-len(failures),
        "pre_onset_application_failures": sum(r["pre_onset_application_failure"] for r in evaluated),
        "pre_onset_source_alarms": sum(r["pre_onset_source_alarm"] for r in evaluated),
        "failures_without_post_onset_source_alert": sum(r["failure_without_post_onset_source_alert"] is True for r in evaluated),
        "failures_without_prior_source_alert": sum(r["failure_without_prior_source_alert"] for r in failures),
        "failures_without_prior_prediction_warning": sum(r["failure_without_prior_prediction_warning"] for r in failures),
        "source_prior_warning_among_observed_failures": _proportion(sum(not r["failure_without_prior_source_alert"] for r in failures), len(failures), confidence),
        "prediction_prior_warning_among_observed_failures": _proportion(sum(not r["failure_without_prior_prediction_warning"] for r in failures), len(failures), confidence),
        "healthy_negative_controls": len(healthy),
        "healthy_source_alarm_probability": _proportion(sum(r["t_S"] is not None for r in healthy), len(healthy), confidence),
        "healthy_prediction_warning_probability": _proportion(sum(r["t_P"] is not None for r in healthy), len(healthy), confidence),
        "healthy_application_threshold_crossings": sum(r["t_A"] is not None for r in healthy),
        "forecast_scoring": {"status": "PRESPECIFIED_ENDPOINT" if forecast_issued_at is not None else "NOT_REQUESTED",
            "issued_at": forecast_issued_at, "horizon": forecast_horizon,
            "selected_trajectories": len(selected), "missing_forecast_trajectories": missing,
            "eligible_trajectories": len(eligible), "abstained_or_ineligible": len(selected)-len(eligible),
            "resolved_trajectories": len(resolved),
            "right_censored_trajectories": sum(f["outcome_status"] == "RIGHT_CENSORED" for f in eligible),
            "not_at_risk_trajectories": sum(f["outcome_status"] == "NOT_AT_RISK" for f in eligible),
            "mean_brier": sum(f["brier"] for f in resolved)/len(resolved) if resolved else None,
            "calibration_observed_frequency": sum(f["outcome"] for f in resolved)/len(resolved) if resolved else None,
            "mean_predicted_probability": sum(f["probability"] for f in resolved)/len(resolved) if resolved else None,
            "limitations": "Resolved-outcome descriptive score only; informative censoring can bias comparisons. No empirical calibration claim or confidence inferred."},
        "events": evaluated,
        "limitations": "Independent groups are operator declarations. Matched/clustered trajectories must not be pooled; event-grid spacing and sustained rules affect timing."}


def paired_differences(runs: Iterable[ResearchRun], scenario_a: str, scenario_b: str, *,
                       metric: str = "failure_without_prior_source_alert") -> dict[str, Any]:
    """One matched difference per independent group; no pooled-window CI.

    Supports binary failure/alert outcomes. Pair matching IDs, application,
    horizon, error rules, and units must agree. No paired mean is a confirmation
    conclusion; this deterministic helper only prepares independent contrasts.
    """
    allowed = {"failure_without_prior_source_alert", "failure_without_prior_prediction_warning",
               "pre_onset_application_failure"}
    if metric not in allowed or scenario_a == scenario_b:
        raise ValueError("distinct scenarios and a supported paired binary metric required")
    groups = {}
    run_ids = set()
    bounded_runs = list(islice(runs, 200_001))
    if len(bounded_runs)>200_000:
        raise ValueError("at most 100000 independent matched pairs are supported")
    for run in bounded_runs:
        if not isinstance(run, ResearchRun) or run.scenario_id not in {scenario_a, scenario_b}:
            raise ValueError("typed runs from exactly the two declared scenarios required")
        if run.run_id in run_ids:
            raise ValueError("duplicate run ID")
        run_ids.add(run.run_id)
        pair = groups.setdefault(run.independent_group_id, {})
        if run.scenario_id in pair:
            raise ValueError("one run per scenario per independent matched group required")
        pair[run.scenario_id] = run
    if not groups or len(groups)>100_000:
        raise ValueError("bounded nonempty matched groups required")
    differences = []
    seen_pairs = set()
    for group_id, pair in sorted(groups.items()):
        if set(pair) != {scenario_a, scenario_b}:
            raise ValueError("unmatched groups cannot be silently omitted")
        a, b = pair[scenario_a], pair[scenario_b]
        signature = lambda r: (r.application, r.unit, r.threshold, r.end, r.failure_consecutive,
                               r.recovery_consecutive, r.event_definition_id, tuple(o.time for o in r.observations))
        if not a.pair_id or a.pair_id != b.pair_id or signature(a) != signature(b):
            raise ValueError("matched pair identity and observation/event definitions must agree")
        if a.pair_id in seen_pairs:
            raise ValueError("the same matched pair cannot be relabeled as independent groups")
        seen_pairs.add(a.pair_id)
        ea, eb = evaluate_run(a), evaluate_run(b)
        differences.append({"independent_group_id": group_id, "pair_id": a.pair_id,
                            "a": int(ea[metric]), "b": int(eb[metric]),
                            "difference_b_minus_a": int(eb[metric])-int(ea[metric])})
    return {"metric": metric, "scenario_a": scenario_a, "scenario_b": scenario_b,
            "independent_pairs": len(differences), "pairs": differences,
            "mean_difference_b_minus_a": sum(v["difference_b_minus_a"] for v in differences)/len(differences),
            "inference_status": "DESCRIPTIVE_ONLY: no confirmation test or sample size frozen"}
