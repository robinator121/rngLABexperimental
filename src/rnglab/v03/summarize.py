"""Evidence summaries from checksum-verified V0.3 trial checkpoints.

No thresholds or experiments are fitted here. Statistics never pool physical
workloads or correlated scenario components. Trial ensembles, rather than
within-chain IID assumptions, supply uncertainty. Repeated complete streams
retain descriptive outcomes but suppress independent-replication intervals.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import time

import numpy as np
from scipy import stats

from rnglab.v02.runner import atomic_json, canonical
from rnglab.v03.runner import applicable, boundary_events, payloads, settings


VERSION = "v03-checkpoint-summary-1"
TRIAL_FIELDS = ("workload", "scenario", "trial", "seed", "seed_namespace", "purpose", "kind",
    "status", "calibration_key", "model_target", "mean", "signed_error", "absolute_error",
    "samples", "burn_in", "horizon", "proposals", "dimension", "emitted_uniforms",
    "initialization_uniforms", "total_simulation_uniforms", "elapsed_seconds", "peak_rss_mib",
    "fault_onset", "analysis_onset", "administrative_control_epoch", "source_recovery", "nonstationary_input", "initial_state", "sampler_backend",
    "stream_sha256", "initialization_uniform_sha256", "initial_key_sha256", "final_key_sha256",
    "initial_key_reused_within_group", "stream_reused_within_group", "any_pre_onset_application",
    "any_pre_onset_detector", "globally_at_risk", "ess", "tau", "mcse", "batch_mcse", "mixing_warning",
    "base_uniform_draws", "decision_uniform_draws", "initial_uniform_draws", "normal_random_draws",
    "inverse_normal_transforms", "raw_pcg64_outputs", "seed_material_bits", "reseed_events",
    "endpoint_recent_mean", "endpoint_recent_error", "endpoint_recent_score", "endpoint_cumulative_score",
    "horizon_max_cumulative_error", "horizon_max_recent_error", "eligible_max_cumulative_error",
    "eligible_max_recent_error", "horizon_max_cumulative_score", "horizon_max_recent_score",
    "auxiliary", "source_provenance", "events", "analysis_events")
EVENT_FIELDS = ("workload", "scenario", "trial", "family", "rule", "fault_onset", "source_recovery",
    "horizon", "first_excursion", "confirmation", "post_onset_excursion", "post_onset_confirmation",
    "pre_onset_confirmation", "post_recovery_within_band", "selected_event", "globally_at_risk",
    "analysis_onset", "analysis_selected_event", "analysis_pre_onset_confirmation", "event_origin")
WINDOW_FIELDS = ("workload", "scenario", "trial", "start", "time", "retained", "phase", "window_mean",
    "cumulative_mean", "recent_mean", "cumulative_error", "recent_error", "acceptance",
    "cumulative_score", "recent_score", "cumulative_tolerance", "recent_tolerance")
MONITOR_FIELDS = ("workload", "scenario", "trial", "time", "monobit", "distribution", "serial",
                  "temporal", "periodic", "cross", "cusum", "composite")


def _finite(value):
    return value is not None and math.isfinite(float(value))


def _json_value(value):
    return canonical(value) if isinstance(value, (dict, list, tuple)) else value


def _write_table(path, rows, fields=None):
    fields = fields or list(dict.fromkeys(key for row in rows for key in row))
    if not fields:
        fields = ["no_records"]
    temporary = path.with_name(path.name + f".{os.getpid()}.tmp")
    with temporary.open("w", newline="") as file:
        writer = csv.DictWriter(file, fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: _json_value(row.get(key)) for key in fields})
    os.replace(temporary, path)


def mean_interval(values, alpha=.05, eligible=True):
    """Across-independent-trial Student-t interval, with explicit suppression."""
    a = np.asarray(values, dtype=float)
    if not eligible or len(a) < 2 or not np.isfinite(a).all():
        return None
    mean = float(a.mean())
    radius = float(stats.t.ppf(1-alpha/2, len(a)-1) * a.std(ddof=1)/np.sqrt(len(a)))
    return [mean-radius, mean+radius]


def binomial_interval(successes, denominator, alpha=.05, eligible=True):
    """Exact Clopper-Pearson interval for independent trial event indicators."""
    if denominator == 0 or not eligible:
        return None
    if not 0 <= successes <= denominator:
        raise ValueError("invalid event counts")
    low = 0. if successes == 0 else float(stats.beta.ppf(alpha/2, successes, denominator-successes+1))
    high = 1. if successes == denominator else float(stats.beta.ppf(1-alpha/2, successes+1, denominator-successes))
    return [low, high]


def _proportion(successes, denominator, eligible=True, alpha=.05):
    return {"successes": int(successes), "denominator": int(denominator),
            "rate": successes/denominator if denominator else None,
            "ci95": binomial_interval(successes, denominator, alpha, eligible),
            "interval": "Clopper-Pearson pointwise" if eligible else "suppressed: repeated/unverified streams"}


def _event_time(event, onset):
    return event.get("confirmation") if onset is None else event.get("post_onset_confirmation")


def _pre_event(event, onset):
    confirmation = event.get("confirmation")
    return bool(event.get("pre_onset_confirmation", False) or
                (onset is not None and confirmation is not None and confirmation <= onset))


def _analysis_events(row):
    return row.get("analysis_events") or row["events"]


def _clock(row):
    return row.get("analysis_onset", row["fault_onset"])


def _administrative_events(payload, marker, calibration):
    """Replay a predeclared control epoch without changing checkpoint evidence."""
    result = payload["result"]
    if not 0 <= marker < result["horizon"]:
        raise ValueError("administrative analysis onset must lie within observation horizon")
    stratum = calibration["strata"].get(result["calibration_key"])
    if stratum is None or stratum["target"] != result["model_target"]:
        raise ValueError("administrative replay has no matching calibrated target/stratum")
    k = stratum["settings"]["confirmation_windows"]
    monitors = payload["monitor"]
    windows = [w for w in payload["windows"] if w["retained"] >= stratum["settings"]["min_measurement"]]
    events = {"detectors": {}, "applications": {}}
    for rule, threshold in stratum["monitor_cutoffs"].items():
        events["detectors"][rule] = boundary_events([m[rule] for m in monitors], [m["time"] for m in monitors],
                                                   threshold, k, onset=marker)
    for criterion, threshold in stratum["application_cutoffs"].items():
        # Stored scores are required: the original calibrated window processing
        # is evidence, and analysis must not silently invent an alternative scale.
        scores = [w[criterion+"_score"] for w in windows]
        events["applications"][criterion] = boundary_events(scores, [w["time"] for w in windows],
                                                           threshold, k, onset=marker)
    return events


def _independence(records):
    hashes = [r.get("stream_sha256") for r in records]
    seeds = [r.get("seed") for r in records]
    duplicated_streams = len(hashes)-len(set(hashes))
    duplicated_seeds = len(seeds)-len(set(seeds))
    eligible = len(records) >= 1 and all(hashes) and all(seeds) and not duplicated_streams and not duplicated_seeds
    return {"independent_intervals_eligible": bool(eligible),
            "duplicate_full_streams": duplicated_streams,
            "duplicate_seed_identities": duplicated_seeds,
            "unique_full_streams": len(set(h for h in hashes if h)),
            "independence_scope": "independent trial streams assumed from disjoint identities; no within-chain IID claim",
            "suppression_reason": None if eligible else "no trials or missing/repeated complete stream/seed identity"}


def _distribution(values):
    values = np.asarray([v for v in values if _finite(v)], dtype=float)
    if not len(values):
        return {"observed": 0, "mean": None, "median": None, "q25": None, "q75": None}
    return {"observed": len(values), "mean": float(values.mean()), "median": float(np.median(values)),
            "q25": float(np.quantile(values, .25)), "q75": float(np.quantile(values, .75))}


def kaplan_meier(delays, observed, followups):
    """Descriptive product-limit curve and RMST under administrative censoring.

    RMST truncates at the smallest follow-up for common support. Delay means
    among observed events are reported separately and must not replace RMST.
    """
    if not len(followups):
        return {"n": 0, "median_delay": None, "restricted_mean_delay": None, "restriction": None,
                "terminal_survival": None, "curve": []}
    times = np.asarray([d if o else f for d, o, f in zip(delays, observed, followups)], dtype=float)
    flags = np.asarray(observed, dtype=bool)
    limit = float(min(followups))
    survival, previous, area = 1., 0., 0.
    median = None
    curve = []
    for current in sorted(set(times.tolist())):
        if current < 0:
            raise ValueError("negative event/follow-up delay")
        if current <= limit:
            area += survival*(current-previous)
            previous = current
        risk = int(np.count_nonzero(times >= current))
        events = int(np.count_nonzero((times == current) & flags))
        censored = int(np.count_nonzero((times == current) & ~flags))
        survival *= 1-events/risk
        if median is None and survival <= .5:
            median = float(current)
        curve.append({"delay": float(current), "at_risk": risk, "events": events,
                      "censored": censored, "survival": survival})
    if previous < limit:
        # Only curve values up to common support contribute to this integral.
        at_limit = 1.
        for point in curve:
            if point["delay"] <= limit:
                at_limit = point["survival"]
        area += at_limit*(limit-previous)
    return {"n": len(followups), "median_delay": median, "restricted_mean_delay": area,
            "restriction": limit, "terminal_survival": survival, "curve": curve,
            "scope": "descriptive KM; fixed administrative horizons, no proportional-hazards model"}


def _paired_bootstrap(scenario, control, repetitions, seed):
    rng = np.random.default_rng(seed)
    n = len(scenario)
    estimates = np.empty((repetitions, 4))
    for start in range(0, repetitions, 256):
        count = min(256, repetitions-start)
        indices = rng.integers(0, n, size=(count, n))
        s, c = scenario[indices], control[indices]
        estimates[start:start+count, 0] = (s-c).mean(axis=1)
        estimates[start:start+count, 1] = (np.abs(s)-np.abs(c)).mean(axis=1)
        estimates[start:start+count, 2] = (s*s-c*c).mean(axis=1)
        estimates[start:start+count, 3] = np.sqrt((s*s).mean(axis=1))-np.sqrt((c*c).mean(axis=1))
    return np.quantile(estimates, [.025, .975], axis=0)


def _group_bootstrap(errors, repetitions, seed):
    """Independent-trial percentile intervals for skewed error distributions."""
    rng = np.random.default_rng(seed)
    estimates = np.empty((repetitions, 4))
    n = len(errors)
    for start in range(0, repetitions, 256):
        count = min(256, repetitions-start)
        values = errors[rng.integers(0, n, size=(count, n))]
        estimates[start:start+count, 0] = values.mean(axis=1)
        estimates[start:start+count, 1] = np.abs(values).mean(axis=1)
        estimates[start:start+count, 2] = np.sqrt((values*values).mean(axis=1))
        estimates[start:start+count, 3] = values.var(axis=1, ddof=1)
    return np.quantile(estimates, [.025, .975], axis=0)


def _group_summary(records, family_size, alpha, bootstrap_replicates):
    first = records[0]
    targets = {r["model_target"] for r in records}
    if len(targets) != 1:
        raise ValueError("mixed exact targets within workload/scenario group")
    means = np.asarray([r["mean"] for r in records])
    target = float(first["model_target"])
    errors = means-target
    n = len(records)
    independence = _independence(records)
    eligible = independence["independent_intervals_eligible"]
    bias = float(errors.mean())
    population_variance = float(np.var(means))
    mse = float(np.mean(errors*errors))
    keys = [r["initial_key_sha256"] for r in records if r.get("initial_key_sha256")]
    key_duplicates = len(keys)-len(set(keys))
    durations = [r["elapsed_seconds"] for r in records if _finite(r.get("elapsed_seconds"))]
    result = {"workload": first["workload"], "scenario": first["scenario"], "trials": n,
        "target": target, "mean": float(means.mean()), "bias": bias, "mae": float(np.mean(np.abs(errors))),
        "mse": mse, "rmse": float(np.sqrt(mse)), "bias_squared": bias*bias,
        "sampling_variance_population": population_variance,
        "sampling_variance_unbiased": float(np.var(means, ddof=1)) if n > 1 else None,
        "mse_decomposition_residual": mse-bias*bias-population_variance,
        "decomposition_scope": "empirical MSE = squared ensemble sample bias + n-denominator variance; n-1 variance reported separately",
        "mean_ci95": mean_interval(means, alpha, eligible),
        "mean_ci_familywise": mean_interval(means, alpha/family_size, eligible),
        "familywise_mean_groups": family_size, "family_alpha": alpha,
        "initial_key_duplicates": key_duplicates,
        "initial_key_uncertainty_scope": "conditional on deliberately common key if keys repeat; no seed-diversity/security conclusion" if key_duplicates else None,
        "nonstationary_trials": sum(bool(r.get("nonstationary_input")) for r in records),
        "initial_states": dict(Counter(r.get("initial_state", "unspecified") for r in records)),
        "elapsed_seconds_total": float(sum(durations)), "trial_seconds": _distribution(durations),
        "emitted_uniforms_total": sum(r.get("emitted_uniforms", 0) for r in records),
        "initialization_uniforms_total": sum(r.get("initialization_uniforms", 0) for r in records),
        "total_simulation_uniforms": sum(r.get("total_simulation_uniforms", r.get("emitted_uniforms", 0)) for r in records),
        "ess_descriptive": _distribution([r.get("ess") for r in records]),
        "mixing_warning_trials": sum(bool(r.get("mixing_warning")) for r in records),
        "peak_rss_mib": max((r["peak_rss_mib"] for r in records if _finite(r.get("peak_rss_mib"))), default=None),
        **independence}
    bootstrap = _group_bootstrap(errors, bootstrap_replicates,
        [40404, int.from_bytes(hashlib.sha256((first["workload"]+"/"+first["scenario"]).encode()).digest()[:4], "big")]) if eligible and n >= 2 else None
    result["bootstrap_ci95"] = {name: bootstrap[:, index].tolist() if bootstrap is not None else None
        for index, name in enumerate(("bias", "mae", "rmse", "sampling_variance_unbiased"))}
    result["bootstrap_scope"] = "pointwise descriptive percentile intervals across independent trials; not simultaneous"
    for field in ("endpoint_recent_error", "endpoint_recent_score", "endpoint_cumulative_score",
        "horizon_max_cumulative_error", "horizon_max_recent_error", "eligible_max_cumulative_error",
        "eligible_max_recent_error", "horizon_max_cumulative_score", "horizon_max_recent_score"):
        result[field] = _distribution([r.get(field) for r in records])
    for field in ("base_uniform_draws", "decision_uniform_draws", "initial_uniform_draws", "normal_random_draws",
                  "inverse_normal_transforms", "raw_pcg64_outputs", "seed_material_bits", "reseed_events"):
        present = [r[field] for r in records if _finite(r.get(field))]
        result[field+"_total"] = int(sum(present)) if present else None
        result[field+"_reported_trials"] = len(present)
    result["auxiliary"] = {}
    for name in sorted(set().union(*(r["auxiliary"].keys() for r in records))):
        values = [r["auxiliary"][name]["mean"] for r in records if name in r["auxiliary"]]
        entry = {"trials": len(values), "mean": float(np.mean(values)),
                 "ci95": mean_interval(values, alpha, eligible and len(values) == n)}
        targets = {r["auxiliary"][name].get("target") for r in records if name in r["auxiliary"]}
        if len(targets) == 1 and None not in targets:
            entry["target"] = float(next(iter(targets)))
            entry["bias"] = entry["mean"]-entry["target"]
        result["auxiliary"][name] = entry
    return result


def _compact_survival(delays, observed, followups):
    result = kaplan_meier(delays, observed, followups)
    result["curve_points"] = len(result.pop("curve"))
    result["curve_data"] = "regenerate from event_results.csv with kaplan_meier; curves are not duplicated across detector/application pairs"
    return result


def _event_rates(records, control, alpha):
    output = []
    first = records[0]
    for family in ("detectors", "applications"):
        rules = sorted(set().union(*(_analysis_events(r).get(family, {}).keys() for r in records)))
        for rule in rules:
            available = [r for r in records if rule in _analysis_events(r).get(family, {})]
            for scope in ("unconditional", "at_risk_selected_rule", "at_risk_all_rules"):
                rows = (available if scope == "unconditional" else
                        [r for r in available if not _pre_event(_analysis_events(r)[family][rule], _clock(r))]
                        if scope == "at_risk_selected_rule" else [r for r in available if r["globally_at_risk"]])
                events = [_event_time(_analysis_events(r)[family][rule], _clock(r)) for r in rows]
                pre = sum(_pre_event(_analysis_events(r)[family][rule], _clock(r)) for r in rows)
                independence = _independence(rows)
                output.append({"workload": first["workload"], "scenario": first["scenario"],
                    "family": family, "rule": rule, "scope": scope, "available_trials": len(available),
                    "excluded_pre_onset_by_scope": len(available)-len(rows),
                    "confirmed": _proportion(sum(t is not None for t in events), len(rows),
                                              independence["independent_intervals_eligible"], alpha),
                    "pre_onset_confirmations": pre,
                    "whole_horizon_confirmations": sum(r["events"][family][rule].get("confirmation") is not None for r in rows),
                    "censored": sum(t is None for t in events),
                    "healthy_control": bool(first["scenario"] == control and all(r["fault_onset"] is None for r in records)),
                    "administrative_control_epoch": any(r.get("administrative_control_epoch", False) for r in records),
                    "clock": "first confirmation when analysis onset absent; reset post-onset confirmation for actual fault or explicitly administrative control epoch",
                    **independence})
    return output


def _operational(records, eligible, alpha, control):
    """Unconditional and globally at-risk conditional event summaries."""
    first = records[0]
    output, recoveries = [], []
    detector_rules = sorted(set().union(*(_analysis_events(r).get("detectors", {}).keys() for r in records)))
    application_rules = sorted(set().union(*(_analysis_events(r).get("applications", {}).keys() for r in records)))
    for detector in detector_rules:
        for application in application_rules:
            complete = [r for r in records if detector in _analysis_events(r).get("detectors", {}) and application in _analysis_events(r).get("applications", {})]
            for scope in ("unconditional", "at_risk_selected_pair", "at_risk_all_rules"):
                rows = (complete if scope == "unconditional" else
                    [r for r in complete if not _pre_event(_analysis_events(r)["detectors"][detector], _clock(r))
                     and not _pre_event(_analysis_events(r)["applications"][application], _clock(r))]
                    if scope == "at_risk_selected_pair" else [r for r in complete if r["globally_at_risk"]])
                categories = Counter()
                deltas, d_delays, f_delays, followups, observed_d, observed_f = [], [], [], [], [], []
                pre_app = pre_det = 0
                for row in rows:
                    onset = _clock(row)
                    start = onset if onset is not None else 0
                    d = _event_time(_analysis_events(row)["detectors"][detector], onset)
                    f = _event_time(_analysis_events(row)["applications"][application], onset)
                    pre_app += _pre_event(_analysis_events(row)["applications"][application], onset)
                    pre_det += _pre_event(_analysis_events(row)["detectors"][detector], onset)
                    followups.append(max(0, row["horizon"]-start))
                    observed_d.append(d is not None); observed_f.append(f is not None)
                    d_delays.append(d-start if d is not None else None)
                    f_delays.append(f-start if f is not None else None)
                    if d is not None and f is not None:
                        delta = f-d
                        deltas.append(delta)
                        category = "detector_before_failure" if delta > 0 else "simultaneous" if delta == 0 else "failure_before_detector"
                    else:
                        category = "detected_no_observed_failure" if d is not None else "failure_no_observed_detection" if f is not None else "neither_observed"
                    categories[category] += 1
                n = len(rows)
                interval_eligible = eligible and _independence(rows)["independent_intervals_eligible"]
                output.append({"workload": first["workload"], "scenario": first["scenario"],
                    "detector": detector, "application": application, "scope": scope,
                    "available_trials": len(complete), "denominator": n,
                    "excluded_pre_onset_by_scope": len(complete)-n if scope != "unconditional" else 0,
                    "excluded_pre_onset_any_rule": len(complete)-n if scope == "at_risk_all_rules" else None,
                    "selected_pair_pre_onset_application": pre_app, "selected_pair_pre_onset_detector": pre_det,
                    "healthy_control": bool(first["scenario"] == control and all(r["fault_onset"] is None for r in records)),
                    "administrative_control_epoch": any(r.get("administrative_control_epoch", False) for r in records),
                    "event_clock": "first confirmation when analysis onset absent; reset post-onset confirmation for actual fault or explicitly administrative control epoch",
                    "detected": _proportion(sum(observed_d), n, interval_eligible, alpha),
                    "application_failed": _proportion(sum(observed_f), n, interval_eligible, alpha),
                    "detector_before_failure": _proportion(categories["detector_before_failure"], n, interval_eligible, alpha),
                    "both_observed": _proportion(len(deltas), n, interval_eligible, alpha),
                    "detector_censored": n-sum(observed_d), "application_censored": n-sum(observed_f),
                    "categories": dict(categories), "delta_failure_minus_detector": _distribution(deltas),
                    "delta_ci95_observed_pairs": mean_interval(deltas, alpha, interval_eligible),
                    "observed_detector_delays": _distribution(d_delays), "observed_failure_delays": _distribution(f_delays),
                    "detector_survival": _compact_survival(d_delays, observed_d, followups),
                    "application_survival": _compact_survival(f_delays, observed_f, followups),
                    "inference_scope": "conditional on observed pairs for delta; administrative censoring retained for event probabilities/KM"})
    for criterion in application_rules:
        available = [r for r in records if criterion in _analysis_events(r).get("applications", {})]
        for scope in ("unconditional", "at_risk_selected_application", "at_risk_all_rules"):
            rows = (available if scope == "unconditional" else
                [r for r in available if not _pre_event(_analysis_events(r)["applications"][criterion], _clock(r))]
                if scope == "at_risk_selected_application" else [r for r in available if r["globally_at_risk"]])
            defined, times = [], []
            failed_without_source_recovery = never_failed = invalid_recovery_order = 0
            for row in rows:
                event = _analysis_events(row)["applications"][criterion]
                failure = _event_time(event, _clock(row))
                recovery = row["source_recovery"]
                if failure is None:
                    never_failed += 1
                    continue
                if recovery is None:
                    failed_without_source_recovery += 1
                    continue
                # Recovery of an application failure requires that failure was
                # already defined before source recovery, not a later new failure.
                if failure > recovery:
                    invalid_recovery_order += 1
                    continue
                regain = event.get("post_recovery_within_band")
                if regain is not None and regain < max(failure, recovery):
                    invalid_recovery_order += 1
                    continue
                defined.append(row)
                times.append(regain-recovery if regain is not None else None)
            interval_eligible = eligible and _independence(defined)["independent_intervals_eligible"]
            recoveries.append({"workload": first["workload"], "scenario": first["scenario"],
                "application": criterion, "scope": scope, "available_trials": len(rows),
                "defined_prior_failures": len(defined), "never_failed": never_failed,
                "failed_without_source_recovery": failed_without_source_recovery,
                "not_prior_to_recovery_or_invalid_order": invalid_recovery_order,
                "recovered": _proportion(sum(t is not None for t in times), len(defined), interval_eligible, alpha),
                "censored_recoveries": sum(t is None for t in times),
                "observed_recovery_delay": _distribution(times),
                "estimand": "estimator within-band confirmation after source recovery, conditional on prior application failure; not microscopic-chain recovery"})
    return output, recoveries


def summarize(output, control="iid", alpha=.05, bootstrap_replicates=2000, write_windows=True):
    """Summarize a snapshot of completed checkpoints and return the saved dict."""
    start = time.perf_counter()
    output = Path(output)
    if not 0 < alpha < 1 or bootstrap_replicates < 100:
        raise ValueError("alpha in (0,1) and at least100 bootstrap replicates required")
    manifest_path = output / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    config = manifest.get("identity", {}).get("config", {})
    workload_settings = {w["name"]: settings(config, w) for w in config.get("workloads", [])}
    administrative_marker = config.get("analysis_onset")
    administrative_calibration = None
    administrative_provenance = None
    if administrative_marker is not None:
        if isinstance(administrative_marker, bool) or not isinstance(administrative_marker, int) or administrative_marker < 0:
            raise ValueError("analysis_onset must be a nonnegative integer proposal boundary")
        if not config.get("calibration_file"):
            raise ValueError("administrative control replay requires a frozen calibration artifact")
        calibration_path = Path(config["calibration_file"])
        if not calibration_path.is_absolute():
            calibration_path = Path(__file__).parents[3] / calibration_path
        administrative_calibration = json.loads(calibration_path.read_text())
        calibration_hash = hashlib.sha256(canonical(administrative_calibration).encode()).hexdigest()
        if calibration_hash != manifest.get("identity", {}).get("calibration_sha256"):
            raise ValueError("administrative calibration checksum differs from experiment manifest")
        if administrative_calibration.get("engine_sha256") != manifest.get("identity", {}).get("engine_sha256"):
            raise ValueError("administrative calibration engine differs from experiment manifest")
        administrative_provenance = {"analysis_onset": administrative_marker,
            "calibration_file": str(calibration_path.resolve()), "calibration_sha256": calibration_hash,
            "scope": "counterfactual administrative control epoch; original fault_onset/events/checkpoints remain unchanged"}
    groups, seen, stream_members = defaultdict(list), set(), defaultdict(list)
    warnings = {}
    def warning(code, identity):
        entry = warnings.setdefault(code, {"count": 0, "examples": []})
        entry["count"] += 1
        if len(entry["examples"]) < 5:
            entry["examples"].append(identity)
    window_file = monitor_file = None
    temporary_window = output / f"window_results.csv.gz.{os.getpid()}.tmp"
    temporary_monitor = output / f"monitor_results.csv.gz.{os.getpid()}.tmp"
    event_rows = []
    digest = hashlib.sha256()
    if write_windows:
        window_file = gzip.open(temporary_window, "wt", newline="", compresslevel=3)
        monitor_file = gzip.open(temporary_monitor, "wt", newline="", compresslevel=3)
        window_writer = csv.DictWriter(window_file, WINDOW_FIELDS, extrasaction="ignore")
        monitor_writer = csv.DictWriter(monitor_file, MONITOR_FIELDS, extrasaction="ignore")
        window_writer.writeheader(); monitor_writer.writeheader()
    try:
        for payload in payloads(output):
            result = payload["result"]
            identity = (result["workload"], result["scenario"], result["trial"])
            if identity in seen:
                raise ValueError("duplicate checkpoint trial identity")
            seen.add(identity)
            if result.get("status") != "complete":
                warning("noncomplete_checkpoint_omitted", identity)
                continue
            digest.update(hashlib.sha256(canonical(payload).encode()).digest())
            row = {key: result.get(key) for key in TRIAL_FIELDS}
            row["seed"] = canonical(result.get("seed")) if result.get("seed") is not None else None
            row["initialization_uniforms"] = result.get("initialization_uniforms", 0)
            row["total_simulation_uniforms"] = result.get("total_simulation_uniforms", result.get("emitted_uniforms", 0))
            row["auxiliary"] = {key: {k: value[k] for k in ("mean", "target", "bias") if k in value}
                                for key, value in result.get("auxiliary", {}).items()}
            details = result.get("source_details", {})
            for key in ("initial_key_sha256", "final_key_sha256", "seed_material_bits"):
                row[key] = details.get(key)
            row["reseed_events"] = len(details["reseed_events"]) if "reseed_events" in details else None
            row.update({key: details.get("latent_draws", {}).get(key) for key in (
                "base_uniform_draws", "decision_uniform_draws", "initial_uniform_draws", "normal_random_draws",
                "inverse_normal_transforms", "raw_pcg64_outputs")})
            for key in ("ess", "tau", "mcse", "batch_mcse", "mixing_warning"):
                row[key] = result.get("statistics", {}).get(key)
            row["events"] = payload.get("events", {"detectors": {}, "applications": {}})
            onset = result.get("fault_onset")
            row["analysis_onset"] = onset
            row["administrative_control_epoch"] = bool(administrative_marker is not None and result["scenario"] == control)
            if row["administrative_control_epoch"]:
                if onset is not None:
                    raise ValueError("administrative control epoch cannot replace an actual fault onset")
                row["analysis_onset"] = administrative_marker
                row["analysis_events"] = _administrative_events(payload, administrative_marker, administrative_calibration)
            else:
                row["analysis_events"] = None
            effective_events = _analysis_events(row)
            row["any_pre_onset_application"] = any(_pre_event(e, _clock(row)) for e in effective_events.get("applications", {}).values())
            row["any_pre_onset_detector"] = any(_pre_event(e, _clock(row)) for e in effective_events.get("detectors", {}).values())
            row["globally_at_risk"] = not row["any_pre_onset_application"] and not row["any_pre_onset_detector"]
            for family in ("detectors", "applications"):
                for rule, event in row["events"].get(family, {}).items():
                    confirmation = event.get("confirmation")
                    expected_pre = bool(onset is not None and confirmation is not None and confirmation <= onset)
                    if bool(event.get("pre_onset_confirmation", False)) != expected_pre:
                        warning("pre_onset_flag_inconsistent", identity)
                    post = event.get("post_onset_confirmation")
                    if onset is not None and post is not None and post <= onset:
                        raise ValueError("post-onset confirmation is not after fault onset")
                    selected = _event_time(event, onset)
                    if selected is not None and selected > result["horizon"]:
                        raise ValueError("event outside observation horizon")
                    event_rows.append({**{key: row.get(key) for key in ("workload", "scenario", "trial", "fault_onset", "source_recovery", "horizon", "globally_at_risk")},
                                       "family": family, "rule": rule, **event, "selected_event": selected,
                                       "analysis_onset": _clock(row),
                                       "analysis_selected_event": _event_time(effective_events[family][rule], _clock(row)),
                                       "analysis_pre_onset_confirmation": _pre_event(effective_events[family][rule], _clock(row)),
                                       "event_origin": "administrative_control_epoch" if row["administrative_control_epoch"] else "checkpoint"})
            windows = payload.get("windows", [])
            minimum = workload_settings.get(row["workload"], {}).get("min_measurement", 0)
            eligible = [w for w in windows if w.get("retained", 0) >= minimum and w.get("cumulative_error") is not None]
            measurement = [w for w in windows if w.get("cumulative_error") is not None]
            last = measurement[-1] if measurement else {}
            for criterion in ("cumulative", "recent"):
                row[f"horizon_max_{criterion}_error"] = max((w[f"{criterion}_error"] for w in measurement if _finite(w.get(f"{criterion}_error"))), default=None)
                row[f"eligible_max_{criterion}_error"] = max((w[f"{criterion}_error"] for w in eligible if _finite(w.get(f"{criterion}_error"))), default=None)
                row[f"horizon_max_{criterion}_score"] = max((w[f"{criterion}_score"] for w in eligible if _finite(w.get(f"{criterion}_score"))), default=None)
                row[f"endpoint_{criterion}_score"] = last.get(f"{criterion}_score")
            row["endpoint_recent_mean"] = last.get("recent_mean")
            row["endpoint_recent_error"] = last.get("recent_error")
            if write_windows:
                prefix = {key: row[key] for key in ("workload", "scenario", "trial")}
                for window in windows:
                    window_writer.writerow({**prefix, **window})
                for monitor in payload.get("monitor", []):
                    monitor_writer.writerow({**prefix, **monitor})
            source_spec = row.get("source_provenance", {}).get("spec", {})
            if source_spec.get("provider_spec", {}).get("provider") == "triggered" and not row["reseed_events"]:
                warning("triggered_provider_has_no_reseed_events", identity)
            if not row.get("stream_sha256"):
                warning("missing_stream_identity", identity)
            else:
                stream_members[row["stream_sha256"]].append(identity)
            groups[identity[:2]].append(row)
    finally:
        if window_file is not None:
            window_file.close(); monitor_file.close()
    if not groups:
        raise ValueError("no completed trial checkpoints")
    if write_windows:
        os.replace(temporary_window, output / "window_results.csv.gz")
        os.replace(temporary_monitor, output / "monitor_results.csv.gz")
    for records in groups.values():
        keys = Counter(r.get("initial_key_sha256") for r in records if r.get("initial_key_sha256"))
        streams = Counter(r.get("stream_sha256") for r in records if r.get("stream_sha256"))
        for row in records:
            row["initial_key_reused_within_group"] = bool(row.get("initial_key_sha256") and keys[row["initial_key_sha256"]] > 1)
            row["stream_reused_within_group"] = bool(row.get("stream_sha256") and streams[row["stream_sha256"]] > 1)
    family_size = len(groups)
    summaries, operational, recoveries, event_rates = [], [], [], []
    for identity, records in sorted(groups.items()):
        group = _group_summary(records, family_size, alpha, bootstrap_replicates)
        summaries.append(group)
        operations, recovery = _operational(records, group["independent_intervals_eligible"], alpha, control)
        operational.extend(operations); recoveries.extend(recovery)
        event_rates.extend(_event_rates(records, control, alpha))
    paired = []
    comparisons = [(key, records) for key, records in sorted(groups.items())
                   if key[1] != control and (key[0], control) in groups]
    paired_family = max(1, 2*len(comparisons))
    for (workload, scenario), records in comparisons:
        controls = {r["trial"]: r for r in groups[(workload, control)]}
        matches = [(r, controls[r["trial"]]) for r in records if r["trial"] in controls and r["seed"] == controls[r["trial"]]["seed"]]
        if not matches:
            warning("no_matching_control_seeds", (workload, scenario))
            continue
        s = np.array([r["mean"]-r["model_target"] for r, c in matches])
        c = np.array([r["mean"]-r["model_target"] for q, r in matches])
        eligible = _independence([r for r, _ in matches])["independent_intervals_eligible"] and _independence([r for _, r in matches])["independent_intervals_eligible"]
        signed, absolute = s-c, np.abs(s)-np.abs(c)
        bootstrap = _paired_bootstrap(s, c, bootstrap_replicates,
            [30303, int.from_bytes(hashlib.sha256((workload+"/"+scenario).encode()).digest()[:4], "big")]) if eligible and len(matches) >= 2 else None
        paired.append({"workload": workload, "scenario": scenario, "control": control,
            "matched_trials": len(matches), "unmatched_scenario_trials": len(records)-len(matches),
            "independent_intervals_eligible": bool(eligible),
            "signed_mean_difference": float(signed.mean()), "signed_mean_ci95": mean_interval(signed, alpha, eligible),
            "signed_mean_ci_familywise": mean_interval(signed, alpha/paired_family, eligible),
            "mae_difference": float(absolute.mean()), "mae_difference_ci95": mean_interval(absolute, alpha, eligible),
            "mae_difference_ci_familywise": mean_interval(absolute, alpha/paired_family, eligible),
            "mse_difference": float(np.mean(s*s-c*c)), "rmse_difference": float(np.sqrt(np.mean(s*s))-np.sqrt(np.mean(c*c))),
            "bootstrap_ci95": {name: bootstrap[:, index].tolist() if bootstrap is not None else None
                               for index, name in enumerate(("signed_mean_difference", "mae_difference", "mse_difference", "rmse_difference"))},
            "standardized_mean_difference_control_sd": float(signed.mean()/c.std(ddof=1)) if len(c)>1 and c.std(ddof=1)>0 else None,
            "familywise_primary_contrasts": paired_family,
            "pairing": "same workload/trial and exact seed; paired resampling; scenarios not independent pooled evidence"})
    all_records = [r for records in groups.values() for r in records]
    raw = sorted(all_records, key=lambda r: (r["workload"], r["scenario"], r["trial"]))
    duplicates = [{"stream_sha256": key, "members": [list(m) for m in members],
                   "within_group_repetition": len(set(m[:2] for m in members)) < len(members),
                   "interpretation": "repeat within group suppresses independent intervals; across scenarios may be planned pairing"}
                  for key, members in stream_members.items() if len(members)>1]
    expected = (config.get("trials", 0)*sum(applicable(w, scenario)
                for w in config.get("workloads", []) for scenario in config.get("scenarios", []))) or None
    executions = []
    for path in sorted(output.glob("execution-*.json")):
        executions.append(json.loads(path.read_text()))
    metadata = {"version": VERSION, "created_unix": time.time(), "output": str(output.resolve()),
        "summary_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "administrative_control_epoch": administrative_provenance,
        "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest() if manifest_path.exists() else None,
        "checkpoint_snapshot_sha256": digest.hexdigest(), "completed_trial_executions": len(raw),
        "configured_trial_executions": expected, "snapshot_complete": len(raw)==expected if expected is not None else None,
        "distinct_seed_blocks": len({(r["workload"], r["seed"]) for r in raw}),
        "distinct_complete_uniform_streams": len(stream_members),
        "failed_attempt_records": len(list((output / "failures").glob("*.json"))), "execution_metadata": executions,
        "groups": family_size, "control": control, "alpha": alpha, "bootstrap_replicates": bootstrap_replicates,
        "trial_rss_reported": sum(_finite(r.get("peak_rss_mib")) for r in raw),
        "summarization_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        "elapsed_seconds": time.perf_counter()-start,
        "scope": ["Trial executions, independent seed blocks and distinct streams are different counts.",
            "Mean t intervals use independent-trial ensemble asymptotics; familywise Bonferroni correction across primary groups.",
            "Paired primary signed/absolute-error contrasts have separate Bonferroni t intervals; bootstrap intervals are pointwise descriptive.",
            "Event proportions use exact pointwise binomial intervals; no joint all-detector false-alarm claim.",
            "Selected-rule/pair at-risk scopes exclude their own pre-onset confirmations; stricter all-rules sensitivity excludes any confirmed pre-onset event. Unconditional counts remain visible.",
            "Administrative control epochs replay frozen rules only when predeclared analysis_onset is present; actual checkpoint events/fault_onset remain original.",
            "Source recovery and estimator within-band recovery are distinct; recent and cumulative criteria are reported separately.",
            "Within-chain ESS/MCSE are descriptive, especially under nonstationarity; independent trial ensembles supply inference.",
            "Fixed-source repetitions do not count as independent replications; conditional means with repeated initial keys do not establish seed diversity."]}
    result = {"metadata": metadata, "groups": summaries, "paired": paired, "operational": operational,
              "event_rates": event_rates,
              "recovery": recoveries, "duplicates": duplicates, "pipeline_warnings": warnings}
    _write_table(output / "trial_results.csv", raw, TRIAL_FIELDS)
    _write_table(output / "event_results.csv", event_rows, EVENT_FIELDS)
    _write_table(output / "event_rate_results.csv", event_rates)
    _write_table(output / "group_results.csv", summaries)
    _write_table(output / "paired_results.csv", paired)
    _write_table(output / "operational_results.csv", operational)
    _write_table(output / "recovery_results.csv", recoveries)
    atomic_json(output / "summary.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--control", default="iid")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--no-windows", action="store_true")
    args = parser.parse_args()
    result = summarize(args.output, args.control, bootstrap_replicates=args.bootstrap, write_windows=not args.no_windows)
    print(json.dumps({"trials": result["metadata"]["completed_trial_executions"],
                      "groups": len(result["groups"]), "warnings": result["pipeline_warnings"]}))


if __name__ == "__main__":
    main()
