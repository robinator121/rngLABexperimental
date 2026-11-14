"""Stage D pilot: calibration, issuance, healthy and fault summaries, implementation checks, stop criteria.

Each function writes one output file of protocol.json `outputs` and returns it. Thresholds
and rules are the frozen ones (constants in pilot_common). Where the protocol makes a
quantity unavailable, the output says so; nothing is imputed, and no denominator shrinks.

Calls prescribed by the protocol are used as prescribed: h1lib.threshold and
h1lib.application_constants for calibration, C_analyze.cp for every binomial interval, and
the measurement layer's issuance, onset_classification, eligibility and issuance_endpoint.
"""
from __future__ import annotations

import datetime
import math
import sys
from pathlib import Path
from typing import Any, Mapping

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))
import pilot_common as pc  # noqa: E402
import pilot_generate as pg  # noqa: E402


def load_run(context: pc.Context, condition: str, replica: int) -> tuple[dict, dict | None]:
    """A planned run's record and arrays; (failure record, None) for a run without a trace.

    Raises if the run has no verified record: analysis never proceeds on a store that is not complete.
    """
    state = pg.classify(context, condition, replica, verify_hash=False)
    if state["state"] == "FAILED_RECORD":
        return pc.read_json(context.record_path(condition, replica)), None
    if state["state"] != "COMPLETE":
        raise pc.PilotError(f"run {condition}/{replica:05d} is {state['state']}: {state['reason']}")
    return context.measure.load_run(context.root, condition, replica)  # verifies file and array hashes


def cp(context: pc.Context, k: int, n: int) -> list[float] | None:
    return context.c_analyze.cp(k, n, pc.OTHER_LEVEL)


def _write(context: pc.Context, name: str, body: Any) -> Any:
    data = pc.canonical_jsonl(body) if name.endswith(".jsonl") else pc.canonical_json(body)
    pc.atomic_write(context.result(name), data)
    return body


def _require_science_plan(context: pc.Context) -> None:
    if context.plan.kind == "smoke":
        raise pc.PilotError("the smoke is for engineering and timing only: no threshold, model or criterion "
                            "is computed from it (protocol.json smoke.purpose)")


# ---------------------------------------------------------------------------
# Calibration (protocol.json operational_alarm.threshold, C_prime, signal_schema.application_statistics, oracle.sigma_eff)
# ---------------------------------------------------------------------------
def calibrate(context: pc.Context) -> dict:
    _require_science_plan(context)
    plan, measure, h1lib = context.plan, context.measure, context.h1lib
    required = plan.calibration_required
    if plan.count("calibration") != required:
        raise pc.PilotError("the plan's calibration condition does not have the required number of runs")
    maxima: dict[str, list[float | None]] = {"full": [], "prime": []}
    records, hashes, mcse = [], {}, {length: [] for length in pc.LENGTHS}
    application_complete = True
    for replica in range(required):
        record, arrays = load_run(context, "calibration", replica)
        hashes[f"calibration/{replica:05d}"] = pc.sha256_file(context.record_path("calibration", replica))
        checkpoints = record.get("conventional") or []
        for label, prime in (("full", False), ("prime", True)):
            try:
                value = measure.issuance(checkpoints, 0.0, prime=prime)["whole_horizon_max"] if arrays is not None else None
            except ValueError:
                value = None
            maxima[label].append(value)
        application = record.get("application") or []
        if arrays is None or len(application) != len(measure.APPLICATION_ENDS) or any(
                pc.finite_or_none(c.get("mean_D")) is None or pc.finite_or_none(c.get("uphill_fraction")) is None
                for c in application):
            application_complete = False
        records.append(record)
        for length in pc.LENGTHS:
            mcse[length].append(pc.finite_or_none(((record.get("endpoints") or {}).get(str(length)) or {}).get("mcse")))
    out: dict[str, Any] = {"schema": "rnglab-v019-d-pilot-calibration-v1", "namespace": plan.namespace,
                           "protocol_sha256": pc.protocol_sha256(), "calibration_runs": required,
                           "order_statistic": plan.threshold_order, "input_record_sha256": hashes, "complete": {}}
    for label in ("full", "prime"):
        values = maxima[label]
        complete = all(v is not None for v in values)
        tau = h1lib.threshold(values, order_statistic=plan.threshold_order) if complete else None
        out[f"maxima_{label}"] = values
        out[f"sorted_maxima_{label}"] = sorted(values) if complete else None
        out[f"threshold_{label}"] = tau
        out[f"ties_at_threshold_{label}"] = sum(1 for v in values if v == tau) if complete else None
        out["complete"][f"maxima_{label}"] = complete
    constants = None
    if application_complete:
        constants = h1lib.application_constants(records)
        if not all(math.isfinite(v) for v in constants.values()) or constants["sigma_D"] <= 0 or constants["sigma_f"] <= 0:
            constants = None
    out["application_constants"] = constants
    out["application_checkpoints_pooled"] = required * len(measure.APPLICATION_ENDS) if application_complete else None
    out["complete"]["application_constants"] = constants is not None
    out["sigma_eff"] = {}
    for length in pc.LENGTHS:
        values = mcse[length]
        sigma = math.sqrt(length) * float(np.mean(values)) if all(v is not None for v in values) else None
        if sigma is not None and not (math.isfinite(sigma) and sigma > 0):
            sigma = None
        out["sigma_eff"][str(length)] = sigma
        out["complete"][f"sigma_eff_{length}"] = sigma is not None
    out["status"] = "COMPLETE" if out["complete"]["maxima_full"] else "INSUFFICIENT_INFORMATION"
    return _write(context, "calibration.json", out)


def read_calibration(context: pc.Context) -> tuple[dict, str]:
    path = context.result("calibration.json")
    if not path.is_file():
        raise pc.PilotError("results/calibration.json does not exist: calibrate first")
    return pc.read_json(path), pc.sha256_file(path)


# ---------------------------------------------------------------------------
# Issuance (protocol.json operational_alarm.issuance, onset_classification, eligibility; endpoint.issuance)
# ---------------------------------------------------------------------------
def issuance_row(context: pc.Context, record: Mapping[str, Any], arrays: Mapping[str, np.ndarray] | None,
                 calibration: Mapping[str, Any], threshold_sha256: str) -> dict:
    measure = context.measure
    t0 = (record.get("profile") or {}).get("t0")
    row: dict[str, Any] = {"condition": record["condition"], "family": record["family"], "replica": record["replica"],
                           "split": context.plan.split(record["condition"], record["replica"]), "t0": t0,
                           "threshold_sha256": threshold_sha256, "run_status": record.get("status")}
    endpoints = record.get("endpoints") or {}
    row["valid"] = {str(length): bool((endpoints.get(str(length)) or {}).get("valid", False)) for length in pc.LENGTHS}
    for label, prime in (("full", False), ("prime", True)):
        tau = calibration[f"threshold_{label}"]
        if arrays is None or tau is None:
            issued = {"status": "UNKNOWN", "alarm_sweep": None, "S_at_alarm": None, "whole_horizon_max": None,
                      "whole_horizon_alarm": None}
        else:
            issued = measure.issuance(record["conventional"], tau, prime=prime)
        known = issued["status"] != "UNKNOWN"
        alarm = issued["alarm_sweep"]
        block = {"status": issued["status"], "alarm_sweep": alarm, "S_at_alarm": issued["S_at_alarm"],
                 "whole_horizon_max": issued.get("whole_horizon_max"), "whole_horizon_alarm": issued.get("whole_horizon_alarm"),
                 **measure.onset_classification(alarm, t0, known=known),
                 "eligibility": {str(length): measure.eligibility(alarm, length, known=known) for length in pc.LENGTHS}}
        if known and arrays is not None:
            at_issuance = measure.issuance_endpoint(arrays["CV"], alarm)
            endpoint = at_issuance["endpoint"]
            block["issuance_endpoint"] = {"n": at_issuance["n"], "already_invalid": at_issuance["already_invalid"],
                                          "valid": None if endpoint is None else endpoint["valid"],
                                          "invalid_reason": None if endpoint is None else endpoint["invalid_reason"]}
        else:
            block["issuance_endpoint"] = {"n": None, "already_invalid": None, "valid": None, "invalid_reason": None}
        row[label] = block
    return row


def issue(context: pc.Context) -> list[dict]:
    """results/issuance.jsonl for every healthy-validation and fault run, in plan order."""
    _require_science_plan(context)
    calibration, digest = read_calibration(context)
    rows = []
    for condition, _, replica in context.plan.runs():
        if condition == "calibration":
            continue
        record, arrays = load_run(context, condition, replica)
        rows.append(issuance_row(context, record, arrays, calibration, digest))
    return _write(context, "issuance.jsonl", rows)


def read_issuance(context: pc.Context) -> list[dict]:
    path = context.result("issuance.jsonl")
    if not path.is_file():
        raise pc.PilotError("results/issuance.jsonl does not exist: run the issuance stage first")
    import json
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]


# ---------------------------------------------------------------------------
# Healthy runs (protocol.json stop_criteria.1, operational_alarm.V2_D, endpoint.issuance.short_interval_check)
# ---------------------------------------------------------------------------
def healthy_summary(context: pc.Context) -> dict:
    _require_science_plan(context)
    plan, measure = context.plan, context.measure
    planned = plan.count("healthy")
    rows = [r for r in read_issuance(context) if r["condition"] == "healthy"]
    if len(rows) != planned:
        raise pc.PilotError("issuance.jsonl does not hold every planned healthy run")
    out: dict[str, Any] = {"schema": "rnglab-v019-d-pilot-healthy-v1", "runs_planned": planned,
                           "runs_with_failure_record": sum(1 for r in rows if r["run_status"] == "FAILED")}
    out["coverage"] = {}
    for length in pc.LENGTHS:
        k = sum(1 for r in rows if r["valid"][str(length)])  # anything not valid counts as invalid; denominator fixed
        out["coverage"][str(length)] = {"covered": k, "runs": planned, "rate": k / planned, "cp95": cp(context, k, planned)}
    for label in ("full", "prime"):
        alarms = sum(1 for r in rows if r[label]["status"] == "ALARM")
        unknown = sum(1 for r in rows if r[label]["status"] == "UNKNOWN")
        interval = cp(context, alarms, planned)
        block = {"alarms": alarms, "no_alarm": sum(1 for r in rows if r[label]["status"] == "NO_ALARM"),
                 "unknown": unknown, "runs": planned, "rate": alarms / planned, "cp95": interval}
        if label == "full":
            block["result"] = ("INSUFFICIENT_INFORMATION" if unknown else
                               "PASS" if interval[1] <= pc.V2D_MAXIMUM else "FAIL")
        out["V2_D" if label == "full" else "C_prime_alarms_descriptive"] = block
    # short intervals at every battery checkpoint, all healthy runs; then conditional on the alarm
    ends = [int(e) for e in measure.BATTERY_ENDS]
    covered = {end: 0 for end in ends}
    available = 0
    for r in rows:
        record, arrays = load_run(context, "healthy", r["replica"])
        short = record.get("short_intervals") if arrays is not None else None
        if not short or [int(e["L"]) for e in short] != ends:
            continue
        available += 1
        for entry in short:
            covered[int(entry["L"])] += bool(entry["valid"])
    out["short_interval_fixed_checkpoints"] = {
        "runs_planned": planned, "runs_with_short_intervals": available,
        "per_checkpoint": [{"n": end, "covered": covered[end], "runs": planned, "cp95": cp(context, covered[end], planned)}
                           for end in ends]}
    for label in ("full", "prime"):
        alarmed = [r for r in rows if r[label]["status"] == "ALARM"]
        per = {}
        for r in alarmed:
            n = r[label]["issuance_endpoint"]["n"]
            cell = per.setdefault(n, {"n": n, "runs": 0, "covered": 0})
            cell["runs"] += 1
            cell["covered"] += bool(r[label]["issuance_endpoint"]["valid"])
        pooled = sum(c["covered"] for c in per.values())
        out[f"short_interval_given_alarm_{label}"] = {
            "runs_with_alarm": len(alarmed), "pooled_covered": pooled, "pooled_cp95": cp(context, pooled, len(alarmed)),
            "per_issuance_n": [{**per[n], "cp95": cp(context, per[n]["covered"], per[n]["runs"])} for n in sorted(per)]}
    return _write(context, "healthy.json", out)


# ---------------------------------------------------------------------------
# Fault runs, descriptive (PROTOCOL.md section 4.2)
# ---------------------------------------------------------------------------
def _five(values: list[float]) -> dict | None:
    if not values:
        return None
    q = np.quantile(np.asarray(values, dtype=np.float64), [0, 0.25, 0.5, 0.75, 1], method="linear")
    return {"count": len(values), "min": float(q[0]), "q25": float(q[1]), "median": float(q[2]), "q75": float(q[3]),
            "max": float(q[4])}


def fault_summary(context: pc.Context) -> dict:
    _require_science_plan(context)
    rows = read_issuance(context)
    out: dict[str, Any] = {"schema": "rnglab-v019-d-pilot-faults-v1", "families": {}}
    for condition, family, planned in context.plan.conditions:
        if family == "healthy":
            continue
        family_rows = [r for r in rows if r["condition"] == condition]
        if len(family_rows) != planned:
            raise pc.PilotError(f"issuance.jsonl does not hold every planned {condition} run")
        block = {}
        for population in ("all", "fitting", "evaluation"):
            chosen = [r for r in family_rows if population == "all" or r["split"] == population]
            cell: dict[str, Any] = {"runs": len(chosen), "alarm_class": {}, "eligibility": {}, "validity_pattern": {},
                                    "issuance_to_length": {}}
            for r in chosen:
                label = r["full"]["onset_classification"] or "unclassified"
                cell["alarm_class"][label] = cell["alarm_class"].get(label, 0) + 1
                pattern = "".join("1" if r["valid"][str(length)] else "0" for length in pc.LENGTHS)
                cell["validity_pattern"][pattern] = cell["validity_pattern"].get(pattern, 0) + 1
            cell["alarm_minus_t0"] = _five([r["full"]["alarm_minus_t0"] for r in chosen
                                            if r["full"]["alarm_minus_t0"] is not None])
            for length in pc.LENGTHS:
                reasons: dict[str, int] = {}
                for r in chosen:
                    e = r["full"]["eligibility"][str(length)]
                    key = "eligible" if e["eligible"] else str(e["reason"])
                    reasons[key] = reasons.get(key, 0) + 1
                cell["eligibility"][str(length)] = reasons
                transition = {"valid_at_issuance_valid_at_L": 0, "valid_at_issuance_invalid_at_L": 0,
                              "invalid_at_issuance_valid_at_L": 0, "invalid_at_issuance_invalid_at_L": 0}
                for r in chosen:
                    if r["full"]["eligibility"][str(length)]["eligible"]:
                        first = "valid" if r["full"]["issuance_endpoint"]["valid"] else "invalid"
                        second = "valid" if r["valid"][str(length)] else "invalid"
                        transition[f"{first}_at_issuance_{second}_at_L"] += 1
                cell["issuance_to_length"][str(length)] = transition
            for a, b in ((0, 1), (0, 2), (1, 2)):
                key = f"{pc.LENGTHS[a]}:{pc.LENGTHS[b]}"
                counts = {"11": 0, "10": 0, "01": 0, "00": 0}
                for r in chosen:
                    counts[("1" if r["valid"][str(pc.LENGTHS[a])] else "0") + ("1" if r["valid"][str(pc.LENGTHS[b])] else "0")] += 1
                cell.setdefault("pairwise_validity", {})[key] = counts
            block[population] = cell
        block["per_run_diagnostics"] = f"results/runs/{condition}/<replica>.json, keys endpoints and severity path in the trace"
        out["families"][family] = block
    return _write(context, "faults.json", out)


# ---------------------------------------------------------------------------
# Implementation checks (protocol.json implementation_checks)
# ---------------------------------------------------------------------------
def implementation_checks(context: pc.Context, replay: Mapping[str, Any] | None) -> dict:
    """Aggregate I-1 to I-3 from the stored records and I-4 from the replay result."""
    failing: dict[str, list[str]] = {"I-1": [], "I-2": [], "I-3": [], "I-3-step": []}
    counted = {"I-1": 0, "I-2": 0, "I-3": 0, "I-3-step": 0}
    step_expected = []
    for condition, family, replica in context.plan.runs():
        record, arrays = load_run(context, condition, replica)
        run_id = f"{condition}/{replica:05d}"
        if family != "healthy" and replica == 0:
            step_expected.append(run_id)
        if arrays is None:
            continue  # a run with a failure record is handled under failed_runs, not as an I-check failure
        names = ("I-1", "I-2") + (("I-3",) if family != "healthy" else ()) + (
            ("I-3-step",) if family != "healthy" and replica == 0 else ())
        for name in names:
            counted[name] += 1
            if not (record.get(name) or {}).get("passed", False):
                failing[name].append(run_id)
    checks = {name: {"runs_checked": counted[name], "failing_runs": failing[name], "passed": not failing[name]}
              for name in failing}
    checks["I-3-step"]["expected_runs"] = step_expected
    checks["I-3-step"]["passed"] = checks["I-3-step"]["passed"] and counted["I-3-step"] == len(step_expected)
    checks["I-4"] = dict(replay) if replay is not None else {"passed": False, "reason": "replay has not been run"}
    passed = all(c["passed"] for c in checks.values())
    return _write(context, "implementation_checks.json",
                  {"schema": "rnglab-v019-d-pilot-implementation-checks-v1", "checks": checks, "passed": passed,
                   "on_failure": "interpretation is blocked (protocol.json implementation_checks.on_failure)"})


# ---------------------------------------------------------------------------
# Stop criteria (protocol.json stop_criteria)
# ---------------------------------------------------------------------------
def criterion_1(healthy: Mapping[str, Any]) -> dict:
    low, high = pc.COVERAGE_BAND
    per, labels = {}, []
    for length in pc.LENGTHS:
        cell = healthy["coverage"][str(length)]
        lower, upper = cell["cp95"]
        label = ("PASS" if lower >= low and upper <= high else
                 "OUTSIDE_BAND" if upper < low or lower > high else "STRADDLES_BAND")
        per[str(length)] = {**cell, "label": label}
        labels.append(label)
    result = "PASS" if all(v == "PASS" for v in labels) else "FAIL" if "OUTSIDE_BAND" in labels else "INSUFFICIENT_INFORMATION"
    return {"per_length": per, "result": result,
            "classification": ("MISCALIBRATION" if "OUTSIDE_BAND" in labels else
                               "INSUFFICIENT_PRECISION" if "STRADDLES_BAND" in labels else None),
            "action": "stop D as designed" if "OUTSIDE_BAND" in labels else None}


def _vote(votes: Mapping[str, bool | None], needed: int) -> str:
    passed = sum(1 for v in votes.values() if v is True)
    unresolved = sum(1 for v in votes.values() if v is None)
    return "PASS" if passed >= needed else "INSUFFICIENT_INFORMATION" if passed + unresolved >= needed else "FAIL"


def criteria_2_3(rows: list[dict]) -> dict:
    """Criteria 2 and 3 on all planned fault runs of each family with the primary alarm below 20,480."""
    length = str(pc.CRITERION_LENGTH)
    per, votes_2, votes_3 = {}, {}, {}
    for family in pc.FAMILIES:
        family_rows = [r for r in rows if r["family"] == family]
        unknown = sum(1 for r in family_rows if r["full"]["status"] == "UNKNOWN")
        eligible = [r for r in family_rows if r["full"]["eligibility"][length]["eligible"]]
        n = len(eligible)
        valid = sum(1 for r in eligible if r["valid"][length])
        already = sum(1 for r in eligible if r["full"]["issuance_endpoint"]["already_invalid"])
        median = float(np.median([pc.CRITERION_LENGTH - r["full"]["alarm_sweep"] for r in eligible])) if n else None
        sufficient = n >= pc.MINIMUM_ELIGIBLE and unknown == 0
        p, a = (valid / n, already / n) if n else (None, None)
        votes_2[family] = (pc.BALANCE_BAND[0] <= p <= pc.BALANCE_BAND[1] and a < pc.ALREADY_INVALID_MAXIMUM) if sufficient else None
        votes_3[family] = (median >= pc.MINIMUM_REMAINING_MEDIAN) if sufficient else None
        per[family] = {"planned": len(family_rows), "unknown_issuance": unknown, "eligible": n, "valid": valid,
                       "fraction_valid": p, "already_invalid": already, "fraction_already_invalid": a,
                       "median_remaining": median, "criterion_2": votes_2[family], "criterion_3": votes_3[family]}
    return {"per_family": per, "criterion_2": _vote(votes_2, 2), "criterion_3": _vote(votes_3, 2)}


def criterion_4(scores: Mapping[str, Any]) -> dict:
    length = str(pc.CRITERION_LENGTH)
    per, votes = {}, {}
    for family in pc.FAMILIES:
        cell = ((scores.get("families") or {}).get(family) or {}).get(length)
        if cell is None:
            votes[family], per[family] = None, {"reason": "no scores"}
            continue
        n = cell["counts"]["eligible"]
        model, baseline = cell["models"].get("F-B-lite") or {}, cell["models"].get("F0+") or {}
        interval = model.get("skill_interval_criterion_4")
        ok = (n >= pc.MINIMUM_ELIGIBLE_EVALUATION and baseline.get("brier") not in (None, 0.0)
              and model.get("skill") is not None and interval is not None)
        votes[family] = bool(model["skill"] >= pc.MINIMUM_SKILL and interval[0] > 0) if ok else None
        per[family] = {"eligible_evaluation": n, "brier_F0plus": baseline.get("brier"), "skill": model.get("skill"),
                       "interval": interval, "level": pc.CRITERION_4_LEVEL, "pass": votes[family]}
    return {"per_family": per, "result": _vote(votes, 1)}


def _three_months_before(date: datetime.date) -> datetime.date:
    """Three calendar months earlier, clamped to the last day of the target month."""
    month, year = date.month - 3, date.year
    if month <= 0:
        month, year = month + 12, year - 1
    last = (datetime.date(year + (month == 12), month % 12 + 1, 1) - datetime.timedelta(days=1)).day
    return datetime.date(year, month, min(date.day, last))


def criterion_5(tdk_deadline: str | None, engineering_complete_on: str | None, results_available_on: str | None) -> dict:
    out = {"tdk_deadline": tdk_deadline, "engineering_complete_on": engineering_complete_on,
           "results_available_on": results_available_on, "engineering_deadline": pc.ENGINEERING_DEADLINE,
           "timezone": "Europe/Budapest", "results_deadline": None}
    if not (tdk_deadline and engineering_complete_on and results_available_on):
        return {**out, "result": "INSUFFICIENT_INFORMATION"}
    deadline = _three_months_before(datetime.date.fromisoformat(tdk_deadline))
    out["results_deadline"] = deadline.isoformat()
    on_time = (datetime.date.fromisoformat(engineering_complete_on) <= datetime.date.fromisoformat(pc.ENGINEERING_DEADLINE)
               and datetime.date.fromisoformat(results_available_on) <= deadline)
    return {**out, "result": "PASS" if on_time else "FAIL"}


# Which stop criteria use the runs of a condition. Calibration runs give the thresholds of criteria 2 to 4;
# healthy runs give criterion 1 and the gate V2-D, on which criteria 2 to 4 depend; fault runs give criteria 2 to 4.
CRITERIA_USING = {"calibration": ("2", "3", "4"), "healthy": ("1", "2", "3", "4"), "fault": ("2", "3", "4")}


def unresolved_infrastructure(context: pc.Context) -> dict:
    """Runs whose infrastructure retries were exhausted, and the criteria they keep from passing.

    Such a run has a failure record with phase 'infrastructure'. It is counted like any failed run
    in the analysis (protocol.json failed_runs). Its status says that it has no result after the
    allowed attempts, not what the cause was, so it is not read as evidence either way: a criterion
    that uses it cannot be PASS (the author's decision of 2026-11-13, confirmed 2026-11-14).
    """
    runs = []
    for condition, family, replica in context.plan.runs():
        if pg.classify(context, condition, replica, verify_hash=False).get("failure_phase") == pg.INFRASTRUCTURE_PHASE:
            runs.append({"condition": condition, "replica": replica})
    blocked = sorted({key for run in runs for key in CRITERIA_USING["fault" if run["condition"].startswith("fault/")
                                                                    else run["condition"]]})
    return {"status": pg.UNRESOLVED if runs else None, "status_meaning": pg.INTERRUPTED_RUN_RULES["unresolved_status"],
            "runs": runs,
            "per_condition": {condition: sum(1 for run in runs if run["condition"] == condition)
                              for condition, _, _ in context.plan.conditions},
            "criteria_that_cannot_pass": blocked, "max_infrastructure_retries": pg.MAX_INFRASTRUCTURE_RETRIES}


def criteria(context: pc.Context, *, engineering_complete_on: str | None = None,
             results_available_on: str | None = None) -> dict:
    """results/criteria.json with the precedence rules: a failed I-check blocks interpretation, and a
    failed or unresolved V2-D makes criteria 2 to 4 uninterpretable. The two dates of criterion 5 are
    supplied by the author; the protocol does not say where else they would come from."""
    _require_science_plan(context)
    healthy = pc.read_json(context.result("healthy.json"))
    checks = pc.read_json(context.result("implementation_checks.json"))
    scores = pc.read_json(context.result("scores.json")) if context.result("scores.json").is_file() else {}
    manifest_path = context.root / "pre_generation_manifest.json"
    deadline = pc.read_json(manifest_path).get("tdk_deadline") if manifest_path.is_file() else None
    one, two_three, four = criterion_1(healthy), criteria_2_3([r for r in read_issuance(context) if r["family"] != "healthy"]), criterion_4(scores)
    five = criterion_5(deadline, engineering_complete_on, results_available_on)
    gate = healthy["V2_D"]["result"]
    raw = {"1": one["result"], "2": two_three["criterion_2"], "3": two_three["criterion_3"], "4": four["result"],
           "5": five["result"]}
    final = dict(raw)
    if gate != "PASS":
        for key in ("2", "3", "4"):
            final[key] = "UNINTERPRETABLE"
    unresolved = unresolved_infrastructure(context)
    for key in unresolved["criteria_that_cannot_pass"]:
        if final[key] == "PASS":
            final[key] = "INSUFFICIENT_INFORMATION"
    out = {"schema": "rnglab-v019-d-pilot-criteria-v1", "implementation_checks_passed": checks["passed"],
           "unresolved_infrastructure": unresolved,
           "interpretation_blocked": not checks["passed"], "V2_D": healthy["V2_D"],
           "criterion_1": one, "criteria_2_3": two_three, "criterion_4": four, "criterion_5": five,
           "results_before_precedence": raw, "results": final,
           "counts_as_failed_for_the_stop_rule": {k: v in ("FAIL", "INSUFFICIENT_INFORMATION", "UNINTERPRETABLE")
                                                  for k, v in final.items()},
           "note": "Each criterion carries its own action (PROTOCOL.md section 5); no single verdict is formed."}
    return _write(context, "criteria.json", out)
