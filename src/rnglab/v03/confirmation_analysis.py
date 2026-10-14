"""Exact frozen V0.3 confirmation criteria, with separate independent stages.

The protocol is read-only. All required components of each compound claim must
pass at its Bonferroni allocation; replication requires both predeclared stages.
No exploratory trial, altered cutoff, replacement claim, or optional stop enters
this analysis. Administrative control epochs never alter checkpoint evidence.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy import stats

from rnglab.v02.runner import atomic_json, canonical, digest
from rnglab.v03.runner import boundary_events, payloads


EXPECTED_PROTOCOL_SHA256 = "2c2a6d5e834cef7b50e68fbb4e28badec35781c068efaa9e627878a4d8fa52a4"


def mean_evidence(values, reference=0., alpha=.05/3):
    """Protocol's central 1-alpha Student-t interval across independent trials."""
    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or not np.isfinite(x).all():
        raise ValueError("finite one-dimensional trial means required")
    n = len(x)
    if n < 2:
        return {"n": n, "mean": float(x.mean()) if n else None, "reference": reference,
                "difference": float(x.mean()-reference) if n else None, "ci": None,
                "difference_ci": None, "standard_error": None, "alpha": alpha}
    mean = float(x.mean());se = float(x.std(ddof=1)/np.sqrt(n))
    radius = float(stats.t.ppf(1-alpha/2, n-1)*se)
    return {"n": n, "mean": mean, "reference": float(reference), "difference": mean-reference,
            "ci": [mean-radius, mean+radius],
            "difference_ci": [mean-reference-radius, mean-reference+radius],
            "standard_error": se, "alpha": alpha,
            "confidence_level": 1-alpha, "df": n-1,
            "method": "independent-trial Student-t ensemble; paired vectors when supplied"}


def _tail_p(value, standard_error, df, alternative):
    if standard_error == 0:
        holds = value < 0 if alternative == "less" else value > 0
        return 0. if holds else 1.
    t = value/standard_error
    return float(stats.t.cdf(t, df) if alternative == "less" else stats.t.sf(t, df))


def equivalence(values, reference, margin, alpha=.05/3):
    """Strict containment of the conservative central protocol interval."""
    if margin <= 0:
        raise ValueError("positive frozen equivalence margin required")
    result = mean_evidence(values, reference, alpha)
    interval = result["difference_ci"]
    result.update({"margin": margin, "pass": bool(interval is not None and interval[0] > -margin and interval[1] < margin)})
    if interval is None:
        result.update(p_tost=None, p_central_interval_rule=None)
    else:
        delta, se, df = result["difference"], result["standard_error"], result["df"]
        p_lower = _tail_p(delta+margin, se, df, "greater")
        p_upper = _tail_p(delta-margin, se, df, "less")
        result["p_tost"] = max(p_lower, p_upper)
        result["p_central_interval_rule"] = min(1., 2*max(p_lower, p_upper))
    result["interpretation"] = "equivalence within declared margin; not exact unbiasedness"
    return result


def signed_effect(values, direction, minimum, alpha=.05/3):
    if direction not in {-1, 1} or minimum < 0:
        raise ValueError("direction +/-1 and nonnegative minimum required")
    x = np.asarray(values)*direction
    result = mean_evidence(x, minimum, alpha)
    result.update(direction=direction, minimum=minimum,
        **{"pass": bool(result["difference_ci"] is not None and result["difference_ci"][0] > 0)})
    result["p_one_sided"] = (_tail_p(result["difference"], result["standard_error"], result["df"], "greater")
                             if result["difference_ci"] is not None else None)
    result["p_central_interval_rule"] = min(1., 2*result["p_one_sided"]) if result["p_one_sided"] is not None else None
    return result


def cp_lower(successes, n, tail):
    return 0. if successes == 0 else float(stats.beta.ppf(tail, successes, n-successes+1))


def cp_upper(successes, n, tail):
    return 1. if successes == n else float(stats.beta.ppf(1-tail, successes+1, n-successes))


def paired_risk_difference(fault, control, minimum=.02, alpha=.05/3):
    """Paired multinomial discordances; exact one-sided CP union bound.

    L(p_fault-only; alpha/2)-U(p_control-only; alpha/2) is the frozen primary
    lower bound. The separately reported wider central interval splits alpha
    into four tails; it must not silently replace the primary lower bound.
    """
    a, b = np.asarray(fault, dtype=bool), np.asarray(control, dtype=bool)
    if a.ndim != 1 or a.shape != b.shape:
        raise ValueError("one-to-one matched event vectors required")
    n = len(a)
    if not n:
        return {"n": 0, "pass": False, "reason": "zero common eligible pairs", "lower_risk_difference": None,
                "p_bound_margin": None, "p_mcnemar_greater": None}
    fault_only = int(np.sum(a & ~b));control_only = int(np.sum(~a & b))
    both = int(np.sum(a & b));neither = n-fault_only-control_only-both
    lower = cp_lower(fault_only, n, alpha/2)-cp_upper(control_only, n, alpha/2)
    upper = cp_upper(fault_only, n, alpha/4)-cp_lower(control_only, n, alpha/4)
    central_lower = cp_lower(fault_only, n, alpha/4)-cp_upper(control_only, n, alpha/4)
    discordant = fault_only+control_only
    mcnemar = float(stats.binomtest(fault_only, discordant, .5, alternative="greater").pvalue) if discordant else 1.
    # Invert the same conservative lower-bound family to retain a p value for
    # the declared >minimum claim, rather than only testing >zero.
    if cp_lower(fault_only, n, .5)-cp_upper(control_only, n, .5) <= minimum:
        p_bound = 1.
    else:
        lo, hi = 0., 1.
        for _ in range(60):
            midpoint = (lo+hi)/2
            bound = cp_lower(fault_only, n, midpoint/2)-cp_upper(control_only, n, midpoint/2)
            if bound > minimum:
                hi = midpoint
            else:
                lo = midpoint
        p_bound = hi
    return {"n": n, "fault_bad": int(a.sum()), "control_bad": int(b.sum()),
        "fault_only": fault_only, "control_only": control_only, "both_bad": both, "neither_bad": neither,
        "risk_difference": float(a.mean()-b.mean()), "lower_risk_difference": lower,
        "minimum_excess": minimum, "marginal_one_sided_tail": alpha/2,
        "central_conservative_interval": [central_lower, upper],
        "p_bound_margin": p_bound, "p_mcnemar_greater": mcnemar,
        "pass": bool(lower > minimum and mcnemar < alpha),
        "method": "exact CP marginal union lower bound plus one-sided exact paired McNemar",
        "zero_count_note": "zero discordance retains nonzero conservative uncertainty; no degenerate bootstrap-only bound"}


def recovery_joint_probability(successes, n, threshold=.5, alpha=.05/3):
    if n == 0:
        return {"n": 0, "successes": successes, "lower": None, "p_greater": None,
                "pass": False, "reason": "zero jointly eligible prior failures"}
    lower = cp_lower(successes, n, alpha)
    p = float(stats.binomtest(successes, n, threshold, alternative="greater").pvalue)
    return {"n": n, "successes": successes, "rate": successes/n, "lower": lower,
            "threshold": threshold, "p_greater": p, "alpha": alpha,
            "pass": bool(lower > threshold), "method": "one-sided exact Clopper-Pearson/binomial"}


def reset_estimates(payload, recovery, budgets):
    """Exact reset averages from complete, supplied-draw proposal windows."""
    windows = payload["windows"]
    result = {}
    for budget in budgets:
        chosen = [w for w in windows if w["start"] >= recovery and w["time"] <= recovery+budget]
        cursor = recovery;total = 0.
        for window in chosen:
            if window["start"] != cursor or window["time"] <= cursor:
                raise ValueError("reset window grid has gaps or overlap")
            total += window["window_mean"]*(window["time"]-window["start"])
            cursor = window["time"]
        if cursor != recovery+budget:
            raise ValueError("reset budget is not exactly available on saved window grid")
        result[budget] = total/budget
    at_recovery = next((w["cumulative_mean"] for w in windows if w["time"] == recovery), None)
    if at_recovery is None:
        raise ValueError("cumulative checkpoint at source recovery is missing")
    return result, at_recovery


def administrative_events(payload, onset, calibration):
    stratum = calibration["strata"][payload["result"]["calibration_key"]]
    k = stratum["settings"]["confirmation_windows"]
    monitor = payload["monitor"]
    windows = [w for w in payload["windows"] if w["retained"] >= stratum["settings"]["min_measurement"]]
    return {"detectors": {rule: boundary_events([m[rule] for m in monitor], [m["time"] for m in monitor], cutoff, k, onset)
             for rule, cutoff in stratum["monitor_cutoffs"].items()},
            "applications": {rule: boundary_events([w[rule+"_score"] for w in windows], [w["time"] for w in windows], cutoff, k, onset)
             for rule, cutoff in stratum["application_cutoffs"].items()}}


def _pre(event, onset):
    return bool(event.get("pre_onset_confirmation") or
                (event.get("confirmation") is not None and event["confirmation"] <= onset))


def _bad(events, criterion, detector):
    failure = events["applications"][criterion]["post_onset_confirmation"]
    alarm = events["detectors"][detector]["post_onset_confirmation"]
    return failure is not None and (alarm is None or failure < alarm)


def analyze_records(records, protocol, stage_label):
    """Apply frozen criteria to already validated, independent thin records."""
    alpha = protocol["per_claim_alpha"]
    claims = {x["id"]: x for x in protocol["primary_claims"]}
    c1, c2, c3 = (claims[x] for x in ("C1", "C2", "C3"))
    controls = records[(c1["workload"], "iid")]
    predictions = {p["rho"]: p for p in c1["predictions"]}
    conditions = {}
    healthy_energy = c1["predictions"][0]["healthy_energy"]
    conditions["healthy_position_accuracy"] = equivalence([r["mean"] for r in controls], 0., c1["position_prediction_tolerance"], alpha)
    conditions["healthy_energy_accuracy"] = equivalence([r["energy"] for r in controls], healthy_energy, c1["energy_accuracy_tolerance"], alpha)
    for scenario in c1["scenarios"]:
        rows = records[(c1["workload"], scenario)]
        rho = rows[0]["rho"]
        prediction = predictions[rho]
        direction = 1 if prediction["position_cumulative_prediction"] > 0 else -1
        conditions[scenario+"_prediction_equivalence"] = equivalence([r["mean"] for r in rows], prediction["position_cumulative_prediction"], c1["position_prediction_tolerance"], alpha)
        conditions[scenario+"_paired_position_effect"] = signed_effect([r["mean"]-c["mean"] for r, c in zip(rows, controls)], direction, c1["minimum_signed_position_effect"], alpha)
        conditions[scenario+"_same_chain_energy_accuracy"] = equivalence([r["energy"] for r in rows], healthy_energy, c1["energy_accuracy_tolerance"], alpha)
    claim1 = {"id": "C1", "pass": all(c["pass"] for c in conditions.values()), "conditions": conditions,
              "reference_numerical_allowances": c1.get("numerical_allowances"),
              "allowance_scope": "reported metadata; frozen equivalence margins are not expanded"}

    faults = records[(c2["workload"], c2["scenario"])]
    baseline = records[(c2["workload"], "iid")]
    common = []
    pre_counts = Counter()
    onset = protocol["budgets"]["onset"]
    for fault, control in zip(faults, baseline):
        fe, ce = fault["events"], control["analysis_events"]
        fp = _pre(fe["applications"][c2["application"]], onset) or _pre(fe["detectors"][c2["detector"]], onset)
        cp = _pre(ce["applications"][c2["application"]], onset) or _pre(ce["detectors"][c2["detector"]], onset)
        pre_counts["fault_pre"] += fp;pre_counts["control_pre"] += cp
        pre_counts["pre_flag_mismatch"] += fp != cp
        if not fp and not cp:
            common.append((_bad(fe, c2["application"], c2["detector"]), _bad(ce, c2["application"], c2["detector"])))
    risk = paired_risk_difference([x[0] for x in common], [x[1] for x in common], c2["minimum_excess_risk"], alpha)
    claim2 = {"id": "C2", "pass": risk["pass"], "total_pairs": len(faults), "common_at_risk_pairs": len(common),
              "excluded_pairs": len(faults)-len(common), "pre_onset_counts": dict(pre_counts), "risk": risk,
              "unconditional_fault_bad": sum(_bad(r["events"], c2["application"], c2["detector"]) for r in faults),
              "unconditional_control_bad": sum(_bad(r["analysis_events"], c2["application"], c2["detector"]) for r in baseline)}

    burst = records[(c3["workload"], c3["scenario"])]
    baseline = records[(c3["workload"], "iid")]
    recovery = c3["source_recovery"]
    jointly_eligible, joint_successes, pre_flagged, sensitivity_n, sensitivity_successes = 0, 0, 0, 0, 0
    categories = Counter()
    for row in burst:
        events = row["events"]["applications"]
        f = [events[k]["post_onset_confirmation"] for k in ("recent", "cumulative")]
        if not all(t is not None and t < recovery for t in f):
            continue
        jointly_eligible += 1
        recent, cumulative = (events[k]["post_recovery_within_band"] for k in ("recent", "cumulative"))
        success = recent is not None and recent <= recovery+16384 and cumulative is None
        joint_successes += success
        flagged = any(_pre(events[k], onset) for k in ("recent", "cumulative"))
        pre_flagged += flagged
        if not flagged:
            sensitivity_n += 1;sensitivity_successes += success
        categories["recent_censored" if recent is None else "recent_observed"] += 1
        categories["cumulative_censored" if cumulative is None else "cumulative_observed"] += 1
    joint = recovery_joint_probability(joint_successes, jointly_eligible, .5, alpha)
    recovery_conditions = {"joint_recovery_event": joint}
    for budget, tolerance in zip(c3["reset_budgets"], c3["reset_tolerances"]):
        recovery_conditions[f"reset_{budget}_accuracy"] = equivalence([r["reset"][budget] for r in burst], 0., tolerance, alpha)
        recovery_conditions[f"reset_{budget}_paired_difference"] = equivalence([r["reset"][budget]-c["reset"][budget] for r, c in zip(burst, baseline)], 0., tolerance, alpha)
    residuals = [r["decomposition_residual"] for r in burst]
    identity = {"n": len(residuals), "maximum_absolute_residual": max(map(abs, residuals)),
                "tolerance": c3["decomposition_tolerance"], "failed_trials": sum(abs(x)>c3["decomposition_tolerance"] for x in residuals),
                "pass": all(abs(x)<=c3["decomposition_tolerance"] for x in residuals)}
    recovery_conditions["cumulative_decomposition"] = identity
    claim3 = {"id": "C3", "pass": all(c["pass"] for c in recovery_conditions.values()),
              "total_trials": len(burst), "jointly_eligible_prior_failures": jointly_eligible,
              "selected_pre_onset_flagged_among_eligible": pre_flagged,
              "sensitivity_excluding_selected_pre_onset": recovery_joint_probability(sensitivity_successes, sensitivity_n, .5, alpha),
              "censoring_counts": dict(categories), "conditions": recovery_conditions,
              "interpretation": c3["interpretation"]}
    return {"stage": stage_label, "alpha_per_compound_claim": alpha, "claims": [claim1, claim2, claim3],
            "all_three_pass": all(x["pass"] for x in (claim1, claim2, claim3)),
            "zero_eligible_denominator": jointly_eligible == 0 or len(common) == 0}


def load_stage(directory, stage_plan, protocol, root):
    directory = Path(directory)
    manifest = json.loads((directory/"manifest.json").read_text())
    identity = manifest["identity"];config = identity["config"]
    config_path = root/stage_plan["config"]
    if hashlib.sha256(config_path.read_bytes()).hexdigest() != stage_plan["config_sha256"]:
        raise ValueError("frozen stage configuration changed")
    if json.loads(config_path.read_text()) != config:
        raise ValueError("checkpoint manifest differs from frozen configuration")
    if config["seed_namespace"] != stage_plan["namespace"] or config["trials"] != stage_plan["trials_per_cell"]:
        raise ValueError("stage namespace/trial budget differs from frozen protocol")
    calibration_path = Path(config["calibration_file"])
    if not calibration_path.is_absolute():
        calibration_path = root/calibration_path
    calibration = json.loads(calibration_path.read_text())
    if digest(calibration) != identity["calibration_sha256"] or calibration["engine_sha256"] != identity["engine_sha256"]:
        raise ValueError("frozen calibration identity differs from checkpoint manifest")
    scenarios = {s["name"]:s for s in config["scenarios"]}
    if set(scenarios) != set(protocol["scenarios"]):
        raise ValueError("a frozen scenario was dropped/replaced")
    expected = {(w["name"], scenario, trial) for w in config["workloads"] for scenario in scenarios
                for trial in range(stage_plan["trials_per_cell"])}
    records = defaultdict(list);seen=set();hashes=defaultdict(set);seeds=defaultdict(set)
    checkpoint_hash = hashlib.sha256()
    c3 = next(c for c in protocol["primary_claims"] if c["id"] == "C3")
    budgets=protocol["budgets"]
    for payload in payloads(directory):
        result=payload["result"];ident=(result["workload"],result["scenario"],result["trial"])
        if ident in seen or ident not in expected or result["status"] != "complete":
            raise ValueError("duplicate, unexpected, or incomplete trial")
        seen.add(ident);group=ident[:2]
        if result["seed_namespace"] != stage_plan["namespace"] or result["seed"][1] != stage_plan["namespace"]:
            raise ValueError("nonindependent stage seed identity")
        seed=tuple(result["seed"]);stream=result["stream_sha256"]
        if seed in seeds[group] or stream in hashes[group]:
            raise ValueError("repeated trial stream/seed cannot supply independent confirmation")
        seeds[group].add(seed);hashes[group].add(stream)
        if result["samples"] != budgets["retained_samples"] or result["burn_in"] != budgets["burn_in"] or result["horizon"] != budgets["raw_horizon"]:
            raise ValueError("checkpoint observation budget differs from frozen protocol")
        if ident[1] == "iid" and result["fault_onset"] is not None:
            raise ValueError("healthy control has been falsely assigned a fault onset")
        if ident[1] != "iid" and result["fault_onset"] != budgets["onset"]:
            raise ValueError("fault onset differs from frozen protocol")
        row={"trial":ident[2],"seed":seed,"mean":result["mean"],
             "energy":result["auxiliary"]["energy"]["mean"],"events":payload["events"],
             "rho":scenarios[ident[1]]["stream"].get("rho",0.),"stream_sha256":stream}
        if ident[1]=="iid":
            row["analysis_events"]=administrative_events(payload,budgets["onset"],calibration)
        if ident[0]==c3["workload"]:
            reset,at_recovery=reset_estimates(payload,c3["source_recovery"],c3["reset_budgets"])
            row["reset"]=reset
            row["decomposition_residual"]=result["mean"]-.5*at_recovery-.5*reset[32768]
        records[group].append(row)
        checkpoint_hash.update(bytes.fromhex(digest(payload)))
    if seen != expected or len(seen) != stage_plan["executions"]:
        raise ValueError("stage not complete: exactly all frozen trials required; no optional stopping")
    for (workload,scenario),rows in records.items():
        rows.sort(key=lambda r:r["trial"])
        control=sorted(records[(workload,"iid")],key=lambda r:r["trial"])
        if any(r["trial"]!=c["trial"] or r["seed"]!=c["seed"] for r,c in zip(rows,control)):
            raise ValueError("scenario/control identifiers are not paired one-to-one")
    return records,{"directory":str(directory.resolve()),"namespace":stage_plan["namespace"],
        "complete_trials":len(seen),"independent_seed_blocks":len({(w,r["seed"]) for (w,s),rows in records.items() for r in rows}),
        "checkpoint_snapshot_sha256":checkpoint_hash.hexdigest(),"manifest_sha256":hashlib.sha256((directory/"manifest.json").read_bytes()).hexdigest(),
        "calibration_sha256":identity["calibration_sha256"],"failed_attempt_records":len(list((directory/"failures").glob("*.json")))}


def analyze(protocol_path, output_path, expected_protocol_sha256=EXPECTED_PROTOCOL_SHA256):
    protocol_path=Path(protocol_path).resolve();output_path=Path(output_path)
    protocol_bytes=protocol_path.read_bytes();actual_hash=hashlib.sha256(protocol_bytes).hexdigest()
    if actual_hash != expected_protocol_sha256:
        raise ValueError("frozen protocol checksum changed")
    protocol=json.loads(protocol_bytes);root=protocol_path.parents[2]
    artifact_checks=[]
    for relative,expected_hash in protocol["artifact_sha256"].items():
        actual=hashlib.sha256((root/relative).read_bytes()).hexdigest()
        artifact_checks.append({"path":relative,"sha256":actual,"matches":actual==expected_hash})
    if not all(a["matches"] for a in artifact_checks):
        raise ValueError("a frozen protocol artifact changed: "+", ".join(a["path"] for a in artifact_checks if not a["matches"]))
    stages=[]
    for stage_plan in protocol["stage_plan"]:
        records,provenance=load_stage(root/stage_plan["output"],stage_plan,protocol,root)
        stage=analyze_records(records,protocol,stage_plan["label"])
        stage["provenance"]=provenance;stages.append(stage)
    if len({s["provenance"]["namespace"] for s in stages}) != len(stages):
        raise ValueError("independent stages reuse a namespace")
    replicated=[]
    for ident in ("C1","C2","C3"):
        passes={s["stage"]:next(c["pass"] for c in s["claims"] if c["id"]==ident) for s in stages}
        replicated.append({"id":ident,"stage_passes":passes,"independently_replicated":all(passes.values()),
                           "interpretation":"passing in both predeclared independent stages required; no pooled rescue"})
    result={"version":"v03-frozen-confirmation-analysis-1","protocol_sha256":actual_hash,
        "analysis_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"created_unix":time.time(),
        "artifact_checks":artifact_checks,"stages":stages,"replicated_claims":replicated,
        "secondary_scope":"predeclared secondary scenarios retained as executed evidence; they never replace a failed primary",
        "statistical_limits":["Central mean intervals use normal/independent-ensemble asymptotics, not distribution-free coverage.",
            "IUT compound claims require every frozen condition; Bonferroni allocation is over three claims.",
            "Control administrative timing never changes its actual fault_onset or checkpoint events.",
            "C3 primary eligibility retains flagged pre-onset cases as frozen; exclusion is separately labeled sensitivity.",
            "No observed alarm/recovery time is imputed; cumulative censoring is an explicit joint endpoint condition.",
            "Fixed reset means assess finite-budget mean accuracy, not instantaneous full-state stationarity."]}
    output_path.parent.mkdir(parents=True,exist_ok=True);atomic_json(output_path,result)
    lines=["# Frozen independent confirmation", "", f"Protocol: `{actual_hash}`.", "",
           "Every required component must pass; a claim is independently replicated only when both stages pass.", ""]
    for claim in replicated:
        labels=", ".join(f"{key}: {'passed' if value else 'not established'}" for key,value in claim["stage_passes"].items())
        lines.append(f"- {claim['id']}: {labels}; independently replicated: **{claim['independently_replicated']}**.")
    lines += ["", "Full component estimates, confidence intervals, p-values, paired counts, censoring and artifact checks are saved in the companion JSON.",
              "No secondary positive finding substitutes for a failed primary. No new threshold is applied retrospectively to exploratory data.", ""]
    output_path.with_suffix(".md").write_text("\n".join(lines))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol",type=Path,default=Path("docs/v0.3/CONFIRMATION_PROTOCOL.json"))
    parser.add_argument("--output",type=Path,default=Path("results/v0.3/confirmatory-analysis/CONFIRMATION_ANALYSIS.json"))
    args=parser.parse_args();result=analyze(args.protocol,args.output)
    print(json.dumps(result["replicated_claims"],indent=2))


if __name__=="__main__":main()
