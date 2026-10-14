"""V0.4 verified-checkpoint descriptives and frozen-claim inference.

This analysis does not alter generators, models, thresholds or old data.
Matched classical/quantum labels do not create additional independent samples.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path
import resource
import time

import numpy as np
from scipy import stats

from rnglab.v02.runner import atomic_json, canonical, digest
from rnglab.v03.runner import boundary_events, payloads
from rnglab.v03.confirmation_analysis import mean_evidence, cp_lower, cp_upper, paired_risk_difference
from rnglab.v04.runner import ROOT, applicable, engine_hash, KERNEL_SHA256, COMPILER_FLAGS

EXPECTED_PROTOCOL_SHA256="7742820d6b46f9f99284d5c4aae193a9e468b8eac944d27a72f85e619e637f3c"


def equivalence(values, reference, margin, allowance=0., alpha=.05/4):
    """Central CI plus outward numerical allowance, inside an unchanged margin."""
    if margin <= 0 or allowance < 0:
        raise ValueError("positive margin and nonnegative numerical allowance required")
    evidence = mean_evidence(values, reference, alpha)
    raw = evidence["difference_ci"]
    widened = [raw[0]-allowance, raw[1]+allowance] if raw is not None else None
    evidence.update(margin=margin, numerical_allowance=allowance, widened_difference_ci=widened,
                    **{"pass": bool(widened is not None and widened[0] > -margin and widened[1] < margin)})
    if raw is not None:
        delta, se, df = evidence["difference"], evidence["standard_error"], evidence["df"]
        effective_margin = margin-allowance
        if effective_margin <= 0:
            p = 1.
        elif se == 0:
            p = 0. if abs(delta) < effective_margin else 1.
        else:
            p = min(1., 2*max(float(stats.t.sf((delta+effective_margin)/se,df)),
                               float(stats.t.cdf((delta-effective_margin)/se,df))))
        evidence["p_central_interval_rule"] = p
    else:
        evidence["p_central_interval_rule"] = None
    return evidence


def signed_effect(values, direction, minimum, alpha=.05/4):
    if direction not in {-1,1} or minimum<0:
        raise ValueError("direction +/-1 and nonnegative minimum required")
    x=np.asarray(values,dtype=float)*direction
    evidence=mean_evidence(x,minimum,alpha)
    evidence.update(direction=direction,minimum=minimum,
                    **{"pass":bool(evidence["difference_ci"] is not None and evidence["difference_ci"][0]>0)})
    if evidence["difference_ci"] is None:
        evidence["p_central_interval_rule"]=None
    elif evidence["standard_error"]==0:
        evidence["p_central_interval_rule"]=0. if evidence["difference"]>0 else 1.
    else:
        evidence["p_central_interval_rule"]=min(1.,2*float(stats.t.sf(evidence["difference"]/evidence["standard_error"],evidence["df"])))
    return evidence


def paired_detection_difference(joint, marginal, minimum=.25, alpha=.05/4):
    result=paired_risk_difference(joint,marginal,minimum,alpha)
    result.update(first_event="joint_detector",second_event="marginal_detector")
    return result


def window_estimate(payload, start, stop, metric="x", minimum_rows=8192):
    """Exact saved-block estimate with explicit available-budget censoring."""
    horizon=payload["result"]["horizon"]
    if start is None:
        return {"available_rows":0,"estimate":None,"censored":True,"reason":"no causal restart event"}
    stop=min(stop,horizon);available=max(0,stop-start)
    if available<minimum_rows:
        return {"available_rows":available,"estimate":None,"censored":True,"reason":"less than frozen minimum usable budget"}
    selected=[w for w in payload["windows"] if w["start"]>=start and w["time"]<=stop]
    cursor=start;total=0.
    for w in selected:
        if w["start"]!=cursor or w["time"]<=cursor:
            raise ValueError("saved estimator blocks contain gaps/overlap")
        total+=w["window_means"][metric]*(w["time"]-w["start"]);cursor=w["time"]
    if cursor!=stop:
        raise ValueError("requested estimator endpoints are not represented by saved blocks")
    return {"available_rows":available,"estimate":total/available,"censored":False,"metric":metric,
            "start":start,"stop":stop,"scope":"finite-budget mean; not instantaneous state-equilibrium certificate"}


def causal_detector_reentry(monitors, rule, cutoff, consecutive=2):
    """First alarm then K healthy looks, without consulting true source end."""
    high=low=0;alarm=None
    for row in monitors:
        if alarm is None:
            high=high+1 if row[rule]>cutoff else 0
            if high>=consecutive:alarm=row["time"]
        else:
            low=low+1 if row[rule]<=cutoff else 0
            if low>=consecutive:
                return {"alarm":alarm,"reentry":row["time"],"right_censored":False,
                        "scope":"causal detector policy; no oracle source-recovery input"}
    return {"alarm":alarm,"reentry":None,"right_censored":True,
            "scope":"causal detector policy; no oracle source-recovery input"}


def administrative_events(payload,onset,calibration):
    stratum=calibration["strata"][payload["result"]["calibration_key"]]
    k=stratum["settings"]["confirmation_windows"]
    monitor=payload["monitor"]
    windows=[w for w in payload["windows"] if w["retained"]>=stratum["settings"]["min_measurement"]]
    return {"detectors":{key:boundary_events([m[key] for m in monitor],[m["time"] for m in monitor],cut,k,onset)
                for key,cut in stratum["monitor_cutoffs"].items()},
            "applications":{key:boundary_events([w[key+"_score"] for w in windows],[w["time"] for w in windows],cut,k,onset)
                for key,cut in stratum["application_cutoffs"].items()}}


def _pre(event,onset):
    return bool(event.get("pre_onset_confirmation") or
                (onset is not None and event.get("confirmation") is not None and event["confirmation"]<=onset))


def _post_confirmed(event,onset,horizon):
    time_index=event.get("post_onset_confirmation")
    return time_index is not None and onset<time_index<=horizon


def _cp(successes,n,alpha=.05,eligible=True):
    return {"successes":successes,"denominator":n,"rate":successes/n if n else None,
            "ci":([cp_lower(successes,n,alpha/2),cp_upper(successes,n,alpha/2)] if n and eligible else None)}


def _engine_checks(manifest):
    identity=manifest["identity"]
    return {"engine_matches_current":identity.get("engine_sha256")==engine_hash(),
            "kernel_matches_current":identity.get("kernel_sha256")==KERNEL_SHA256,
            "compiler_flags_match_current":identity.get("compiler_flags")==COMPILER_FLAGS}


def load_records(output, require_complete=False, control_epoch=None):
    """Stream verified envelopes into bounded thin records, retaining all nulls."""
    output=Path(output);manifest=json.loads((output/"manifest.json").read_text());config=manifest["identity"]["config"]
    engine_checks=_engine_checks(manifest)
    if require_complete and not all(engine_checks.values()):
        raise ValueError("saved executed engine/kernel/compiler provenance differs from frozen current core")
    workload_specs={w["name"]:w for w in config["workloads"]}
    expected={(w["name"],s["name"],t) for w in config["workloads"] for s in config["scenarios"] if applicable(w,s)
              for t in range(config["trials"])}
    calibration=None
    if config.get("calibration_file"):
        path=Path(config["calibration_file"]);path=path if path.is_absolute() else ROOT/path
        calibration=json.loads(path.read_text())
        if digest(calibration)!=manifest["identity"]["calibration_sha256"]:
            raise ValueError("analysis calibration checksum differs from completed experiment")
        if calibration["engine_sha256"]!=manifest["identity"].get("engine_sha256"):
            raise ValueError("calibration engine differs from saved executed engine")
    groups=defaultdict(list);seen=set();snapshot=hashlib.sha256()
    onset=control_epoch if control_epoch is not None else config.get("analysis_onset")
    for payload in payloads(output):
        r=payload["result"];identity=(r["workload"],r["scenario"],r["trial"])
        if identity in seen or identity not in expected or r["status"]!="complete":
            raise ValueError("duplicate, unexpected or incomplete checkpoint identity")
        group=workload_specs[r["workload"]].get("seed_group",r["workload"])
        expected_seed=[config["master_seed"],config["seed_namespace"],
                       int.from_bytes(hashlib.sha256(group.encode()).digest()[:4],"big"),r["trial"]]
        if r["seed"]!=expected_seed or r["seed_namespace"]!=config["seed_namespace"] or r.get("seed_group",r["workload"])!=group:
            raise ValueError("checkpoint seed does not match frozen deterministic trial identity")
        seen.add(identity);snapshot.update(bytes.fromhex(digest(payload)))
        metrics={name:{k:value[k] for k in ("mean","target","error","mcse","batch_mcse","ess","tau","mixing_warning","intended_target_inside_mcse95") if k in value}
                 for name,value in r["metrics"].items()}
        events=payload.get("events",{"detectors":{},"applications":{}})
        analysis_events=administrative_events(payload,onset,calibration) if r["scenario"]=="iid" and onset is not None and calibration else events
        row={"workload":r["workload"],"scenario":r["scenario"],"trial":r["trial"],"seed":r["seed"],
            "seed_group":r.get("seed_group",r["workload"]),"namespace":r["seed_namespace"],"metrics":metrics,
            "targets":r["targets"],"model_spec":r["model_spec"],"samples":r["samples"],"burn_in":r["burn_in"],
            "horizon":r["horizon"],"fault_onset":r["fault_onset"],"source_recovery":r["source_recovery"],
            "stream_sha256":r["stream_sha256"],"position_sha256":r["position_sha256"],"acceptance_sha256":r["acceptance_sha256"],
            "elapsed_seconds":r["elapsed_seconds"],"emitted_uniforms":r["emitted_uniforms"],
            "source_provenance":r["source_provenance"],"events":events,"analysis_events":analysis_events,
            "analysis_onset":onset if r["scenario"]=="iid" else r["fault_onset"],
            "windows":payload["windows"],"monitor":payload["monitor"]}
        row["uniform_diagnostics"]={k:r[k] for k in ("uniform_means","uniform_variances","uniform_cross_correlation") if k in r}
        # Windows are short block-level records, not full trajectories. They
        # preserve variable-budget post-detector estimates without new simulation.
        if calibration:
            stratum=calibration["strata"][r["calibration_key"]]
            policy=causal_detector_reentry(payload["monitor"],"normal_cross",stratum["monitor_cutoffs"]["normal_cross"],
                                           stratum["settings"]["confirmation_windows"])
            policy["premature_relative_to_oracle"]=(policy["reentry"] is not None and r["source_recovery"] is not None and policy["reentry"]<r["source_recovery"])
            row["detector_restart"]=policy
        groups[identity[:2]].append(row)
    if require_complete and seen!=expected:
        raise ValueError("frozen stage incomplete; no optional stopping or silent missing data")
    for rows in groups.values():rows.sort(key=lambda r:r["trial"])
    return groups,{"manifest":manifest,"config":config,"calibration":calibration,
        "executed_engine_checks":engine_checks,
        "completed":len(seen),"expected":len(expected),"complete":seen==expected,
        "checkpoint_snapshot_sha256":snapshot.hexdigest(),
        "independent_seed_identities":len({(r["seed_group"],tuple(r["seed"])) for rows in groups.values() for r in rows}),
        "failed_attempt_records":len(list((output/"failures").glob("*.json")))}


def paired_rows(groups,workload,fault,control="iid"):
    a=groups[(workload,fault)];b=groups[(workload,control)]
    if len({r["trial"] for r in a})!=len(a) or len({r["trial"] for r in b})!=len(b):
        raise ValueError("duplicate paired trial identity")
    lookup={r["trial"]:r for r in b};pairs=[]
    for row in a:
        other=lookup.get(row["trial"])
        if other is None or row["seed"]!=other["seed"]:
            raise ValueError("missing/mismatched paired trial identity")
        pairs.append((row,other))
    if len(a)!=len(b):raise ValueError("paired cell lengths differ")
    return pairs


def _raw_latent_counters(source):
    inputs=source.get("pre_relationship_inputs",[source])
    return {name:sum(item.get(name,0) for item in inputs) for name in
            ("base_raw64_words","raw_fault_initial_bernoulli_draws","decision_uniform_draws")}


def _trial_csv_row(r):
    return {"workload":r["workload"],"scenario":r["scenario"],"trial":r["trial"],"seed_group":r["seed_group"],
        "seed":canonical(r["seed"]),"namespace":r["namespace"],"metrics":canonical(r["metrics"]),
        "stream_sha256":r["stream_sha256"],"position_sha256":r["position_sha256"],"acceptance_sha256":r["acceptance_sha256"],
        "source_provenance":canonical(r["source_provenance"]),"uniform_diagnostics":canonical(r["uniform_diagnostics"]),
        "events":canonical(r["events"]),"analysis_events":canonical(r["analysis_events"]),
        "elapsed_seconds":r["elapsed_seconds"]}


def _post_initial_key_summary(rows):
    post_keys={(key["key_sha256"],key["nonce_epoch"]) for r in rows
               for key in r["source_provenance"].get("experimental_public_keys",[]) if key["epoch"]>0}
    repeated=[r["source_provenance"].get("repeated_key_nonce_pairs",0) for r in rows]
    collapse=bool(post_keys and len(post_keys)==1 and len(rows)>1)
    return {"trials_with_repeated_key_nonce_blocks":sum(n>0 for n in repeated),
        "maximum_repeated_key_nonce_pairs_per_trial":max(repeated,default=0),
        "distinct_post_initial_key_nonce_pairs_across_trials":len(post_keys),
        "post_initial_key_collapse_across_trials":collapse,
        "interval_scope":("Independent initial healthy-key trials conditional on a fixed shared degraded suffix; full-streamSHAuniqueness does not establish independent post-initial source replicates."
                          if collapse else "Independent trial-key ensemble conditional on the explicitly simulated source/failure model; no cryptographic security inference.")}


def describe(output,write_csv=True,*,records=None):
    output=Path(output);groups,metadata=records if records is not None else load_records(output)
    summaries=[];event_rates=[];contrasts=[]
    scenario_specs={s["name"]:s for s in metadata["config"]["scenarios"]}
    for (workload,scenario),rows in sorted(groups.items()):
        streams=[r["stream_sha256"] for r in rows];eligible=len(set(streams))==len(rows)
        group={"workload":workload,"scenario":scenario,"trials":len(rows),"unique_streams":len(set(streams)),
            "independent_intervals_eligible":eligible,"metrics":{},"elapsed_seconds":sum(r["elapsed_seconds"] for r in rows),
            "emitted_uniforms":sum(r["emitted_uniforms"] for r in rows)}
        group.update(is_control=scenario_specs[scenario].get("is_control",scenario=="iid"),
                     fault_onset=rows[0]["fault_onset"],analysis_onset=rows[0]["analysis_onset"],
                     source_recovery=rows[0]["source_recovery"],
                     input_generators=sorted({r["source_provenance"]["generator"] for r in rows}))
        group["draw_accounting"]={name:sum(r["source_provenance"].get(name,0) for r in rows) for name in
            ("base_uniforms","cipher_bytes","raw_source_bits_consumed_for_generation","raw_input_bits_total_for_generation",
             "initial_raw_bits_consumed","post_initial_raw_bits_consumed","clean_initial_material_bits",
             "output_raw64_words","discarded_word_bits","clean_initial_raw64_words")}
        group["draw_accounting"].update({name:sum(_raw_latent_counters(r["source_provenance"].get("raw_source",{}))[name] for r in rows)
            for name in ("base_raw64_words","raw_fault_initial_bernoulli_draws","decision_uniform_draws")})
        group["separate_observational_probe_bits"]=sum(r["source_provenance"].get("observational_raw_probe",{}).get("bits",0) for r in rows)
        group["generation_seconds"]=sum(r["source_provenance"].get("elapsed_seconds",0.) for r in rows)
        group["stream_reuse_diagnostics"]=_post_initial_key_summary(rows)
        for name in rows[0]["metrics"]:
            values=np.array([r["metrics"][name]["mean"] for r in rows]);target=rows[0]["targets"].get(name)
            entry={"mean":float(values.mean()),"trial_variance":float(values.var(ddof=1)) if len(values)>1 else None,
                "ci95":mean_evidence(values,0.,.05)["ci"] if eligible else None,
                "ci_scope":group["stream_reuse_diagnostics"]["interval_scope"],
                "target":target,"mcse_target_coverage":_cp(sum(bool(r["metrics"][name].get("intended_target_inside_mcse95")) for r in rows),
                    sum("intended_target_inside_mcse95" in r["metrics"][name] for r in rows),eligible=eligible)}
            if target is not None:
                error=values-target;entry.update(bias=float(error.mean()),mae=float(np.abs(error).mean()),rmse=float(np.sqrt((error*error).mean())))
            group["metrics"][name]=entry
        summaries.append(group)
        for family in ("detectors","applications"):
            for rule in rows[0]["analysis_events"].get(family,{}):
                for scope in ("unconditional","at_risk_selected_rule"):
                    selected=rows if scope=="unconditional" else [r for r in rows if not _pre(r["analysis_events"][family][rule],r["analysis_onset"])]
                    events=[r["analysis_events"][family][rule].get("post_onset_confirmation") if r["analysis_onset"] is not None else r["analysis_events"][family][rule].get("confirmation") for r in selected]
                    event_rates.append({"workload":workload,"scenario":scenario,"family":family,"rule":rule,"scope":scope,
                        "confirmed":_cp(sum(t is not None for t in events),len(selected),eligible=eligible),
                        "right_censored":sum(t is None for t in events),"excluded_pre":len(rows)-len(selected),
                        "observed_delays":[t-(r["analysis_onset"] or 0) for r,t in zip(selected,events) if t is not None],
                        "source_operational_scope":"frozen input-family calibrated reference; pipeline-strategy rates descriptive unless separately validated"})
        if scenario!="iid" and (workload,"iid") in groups:
            pairs=paired_rows(groups,workload,scenario)
            contrast={"workload":workload,"scenario":scenario,"matched_pairs":len(pairs),"metrics":{}}
            for name in rows[0]["metrics"]:
                differences=[a["metrics"][name]["mean"]-b["metrics"][name]["mean"] for a,b in pairs]
                contrast["metrics"][name]=mean_evidence(differences,0.,.05) if eligible else {"mean":float(np.mean(differences)),"ci":None}
            contrasts.append(contrast)
    result={"version":"v04-descriptives-1","groups":summaries,"paired_contrasts":contrasts,"event_rates":event_rates,
        "metadata":{k:v for k,v in metadata.items() if k not in {"manifest","config","calibration"}},
        "scope":"descriptive pilot/calibration/study snapshot; does not establish a frozen confirmatory claim"}
    atomic_json(output/"descriptive_summary.json",result)
    if write_csv and groups:
        with (output/"trial_results.csv").open("w",newline="") as file:
            writer=None
            for rows in groups.values():
                for r in rows:
                    row=_trial_csv_row(r)
                    if writer is None:writer=csv.DictWriter(file,list(row));writer.writeheader()
                    writer.writerow(row)
    return result


def _means(rows,metric):
    return [r["metrics"][metric]["mean"] for r in rows]


def _paired_means(groups,workload,fault,control,metric):
    return [a["metrics"][metric]["mean"]-b["metrics"][metric]["mean"] for a,b in paired_rows(groups,workload,fault,control)]


def _identity(groups,first,second,scenarios):
    comparisons=[]
    for scenario in scenarios:
        a=groups[(first,scenario)];lookup={r["trial"]:r for r in groups[(second,scenario)]}
        if len(a)!=len(lookup):raise ValueError("matched model cell lengths differ")
        for row in a:
            other=lookup.get(row["trial"])
            if other is None or row["seed"]!=other["seed"]:
                raise ValueError("matched classical/quantum seeds differ")
            comparisons.append(all(row[k]==other[k] for k in ("stream_sha256","position_sha256","acceptance_sha256")))
    return {"comparisons":len(comparisons),"failures":sum(not x for x in comparisons),"pass":all(comparisons),
            "evidence_scope":"matched paths, not additional independent replications"}


def _pipeline_identity(groups,workload,pairs):
    failures=[];count=0
    for fault,control in pairs:
        for a,b in paired_rows(groups,workload,fault,control):
            count+=1
            checks={k:a[k]==b[k] for k in ("stream_sha256","position_sha256","acceptance_sha256")}
            for label,row in (("fault",a),("control",b)):
                source=row["source_provenance"]
                checks[label+"_generation_raw_zero"]=source["raw_source_bits_consumed_for_generation"]==0
                checks[label+"_post_initial_raw_zero"]=source["post_initial_raw_bits_consumed"]==0
                checks[label+"_probe_separate"]=source.get("observational_raw_probe",{}).get("feeds_generator") is False
            if not all(checks.values()):failures.append({"scenario":fault,"trial":a["trial"],"checks":checks})
    return {"comparisons":count,"failures":failures,"pass":not failures,
            "scope":"causal absence of unused source input, not a cryptographic security certificate"}


def _primary_claims(groups,metadata,protocol,stage_label):
    alpha=protocol["per_claim_alpha"];output=[]
    for claim in protocol["primary_claims"]:
        if stage_label not in claim["stages"]:continue
        ident=claim["id"];w=claim["workload"];conditions={};details={}
        if ident=="H1":
            healthy=groups[(w,claim["control"])]
            conditions["healthy_x_accuracy"]=equivalence(_means(healthy,"x"),0.,claim["healthy_x_margin"],claim["target_allowance"],alpha)
            conditions["healthy_energy_accuracy"]=equivalence(_means(healthy,"energy"),protocol["intended_energy"],claim["energy_accuracy_margin"],claim["target_allowance"],alpha)
            for scenario in claim["scenarios"]:
                rows=groups[(w,scenario)];prediction=claim["predictions"][scenario]["cumulative_x"]
                conditions[scenario+"_prediction_equivalence"]=equivalence(_means(rows,"x"),prediction,claim["prediction_x_margin"],claim["prediction_x_allowance"],alpha)
                direction=1 if prediction>0 else -1
                conditions[scenario+"_paired_signed_x"]=signed_effect(_paired_means(groups,w,scenario,claim["control"],"x"),direction,claim["minimum_signed_paired_x"],alpha)
                conditions[scenario+"_same_chain_energy_accuracy"]=equivalence(_means(rows,"energy"),protocol["intended_energy"],claim["energy_accuracy_margin"],claim["target_allowance"],alpha)
            conditions["matched_classical_quantum_paths"]=_identity(groups,w,"class_quartic",[s["name"] for s in metadata["config"]["scenarios"]])
        elif ident=="H2":
            scenario=claim["scenario"];rows=groups[(w,scenario)]
            conditions["paired_energy_increase"]=signed_effect(_paired_means(groups,w,scenario,claim["control"],"energy"),1,claim["minimum_paired_energy_increase"],alpha)
            conditions["unfitted_word_energy_prediction"]=equivalence(_means(rows,"energy"),claim["stationary_energy_prediction"],claim["prediction_energy_margin"],claim["prediction_energy_allowance"],alpha)
            for fault,control in claim["restoration_pairs"]:
                for metric,target,margin in (("x",0.,claim["restoration_x_margin"]),("energy",protocol["intended_energy"],claim["restoration_energy_margin"])):
                    conditions[fault+"_"+metric+"_target_accuracy"]=equivalence(_means(groups[(w,fault)],metric),target,margin,claim["target_allowance"],alpha)
                    conditions[fault+"_"+metric+"_paired_healthy_equivalence"]=equivalence(_paired_means(groups,w,fault,control,metric),0.,margin,0.,alpha)
            conditions["no_reseed_exact_negative_controls"]=_pipeline_identity(groups,w,claim["no_reseed_identity_pairs"])
            details["reference_scope"]=claim["independence_approximation"]
            details["source_operations_scope"]=claim["fault_clock"]
        elif ident=="H3":
            scenario=claim["scenario"];pairs=paired_rows(groups,w,scenario,claim["control"])
            recovery=claim["recovery"];end=recovery+claim["reset_samples"]
            residuals={"x":[],"energy":[]};reset_details={}
            for metric,target,margin,field in (("x",0.,claim["reset_x_margin"],"cumulative_mean"),
                                               ("energy",protocol["intended_energy"],claim["reset_energy_margin"],"energy_cumulative_mean")):
                values=[];differences=[]
                for fault,control in pairs:
                    a=window_estimate({"result":fault,"windows":fault["windows"]},recovery,end,metric)
                    b=window_estimate({"result":control,"windows":control["windows"]},recovery,end,metric)
                    if a["censored"] or b["censored"]:raise ValueError("fixed oracle budget missing; no imputation")
                    values.append(a["estimate"]);differences.append(a["estimate"]-b["estimate"])
                    prefix=next(row[field] for row in fault["windows"] if row["time"]==recovery)
                    residuals[metric].append(fault["metrics"][metric]["mean"]-.5*prefix-.5*a["estimate"])
                conditions["oracle_reset_"+metric+"_target_accuracy"]=equivalence(values,target,margin,claim["target_allowance"],alpha)
                conditions["oracle_reset_"+metric+"_paired_healthy_equivalence"]=equivalence(differences,0.,margin,0.,alpha)
                reset_details[metric]={"mean":float(np.mean(values)),"rows_per_trial":claim["reset_samples"]}
            conditions["old_cumulative_negative_x"]=signed_effect(_means(groups[(w,scenario)],"x"),-1,claim["minimum_negative_cumulative_x"],alpha)
            conditions["exact_x_and_energy_half_identity"]={"maximum_absolute_residuals":{k:max(map(abs,v)) for k,v in residuals.items()},
                "trials":len(pairs),"tolerance":claim["decomposition_tolerance"],
                "pass":all(abs(value)<=claim["decomposition_tolerance"] for values in residuals.values() for value in values)}
            details["reset_means"]=reset_details;details["interpretation"]=claim["interpretation"]
        elif ident=="H4":
            rows=groups[(w,claim["scenario"])];eligible=[];eligible_times=[];pre_joint=pre_marginal=0
            for row in rows:
                events=row["analysis_events"]["detectors"]
                j,m=events[claim["detector"]],events[claim["comparison"]]
                pj,pm=_pre(j,protocol["budgets"]["onset"]),_pre(m,protocol["budgets"]["onset"])
                pre_joint+=pj;pre_marginal+=pm
                if not pj and not pm:
                    eligible.append(tuple(_post_confirmed(event,protocol["budgets"]["onset"],protocol["budgets"]["horizon"]) for event in (j,m)))
                    eligible_times.append({"trial":row["trial"],"joint_confirmation":j["post_onset_confirmation"],
                                           "marginal_confirmation":m["post_onset_confirmation"]})
            conditions["paired_detection_advantage"]=paired_detection_difference([r[0] for r in eligible],[r[1] for r in eligible],claim["minimum_detection_probability_excess"],alpha)
            details.update(total_trials=len(rows),common_no_prealarm_trials=len(eligible),excluded_pre=len(rows)-len(eligible),
                joint_pre_confirmations=pre_joint,marginal_pre_confirmations=pre_marginal,
                unconditional_joint_confirmed=sum(_post_confirmed(r["analysis_events"]["detectors"][claim["detector"]],protocol["budgets"]["onset"],protocol["budgets"]["horizon"]) for r in rows),
                unconditional_marginal_confirmed=sum(_post_confirmed(r["analysis_events"]["detectors"][claim["comparison"]],protocol["budgets"]["onset"],protocol["budgets"]["horizon"]) for r in rows),
                common_eligible_confirmation_times=eligible_times,
                false_alarm_scope=claim["false_alarm_scope"])
        else:raise ValueError("unexpected frozen primary claim")
        output.append({"id":ident,"name":claim["name"],"pass":all(c["pass"] for c in conditions.values()),
                       "conditions":conditions,"details":details,"per_claim_alpha":alpha})
    return output


def _healthy_detector_rates(groups,protocol):
    result=[];onset=protocol["budgets"]["onset"]
    for (workload,scenario),rows in groups.items():
        if scenario!="iid":continue
        common=[r for r in rows if all(name in r["analysis_events"]["detectors"] and
                     not _pre(r["analysis_events"]["detectors"][name],onset) for name in ("normal_cross","marginal"))]
        for rule in rows[0]["events"].get("detectors",{}):
            whole=[r["events"]["detectors"][rule] for r in rows]
            selected=[r for r in rows if not _pre(r["analysis_events"]["detectors"][rule],onset)]
            result.append({"workload":workload,"detector":rule,
                "unconditional_whole_horizon":_cp(sum(e["confirmation"] is not None for e in whole),len(rows)),
                "same_epoch_preflag_free":_cp(sum(r["analysis_events"]["detectors"][rule]["post_onset_confirmation"] is not None for r in selected),len(selected)),
                "same_epoch_common_joint_marginal_preflag_free":(_cp(sum(_post_confirmed(r["analysis_events"]["detectors"][rule],onset,protocol["budgets"]["horizon"]) for r in common),len(common))
                    if rule in {"normal_cross","marginal"} else None),
                "administrative_epoch":onset,"actual_fault_onset":None,"preflag_excluded":len(rows)-len(selected),
                "scope":"conditional epoch false-alarm rate does not automatically inherit unconditional1%certificate"})
    return result


def _batch_interval(row,start,stop,metric,target):
    if start is None:return {"blocks":0,"ci95":None,"target_inside":None}
    values=[w["window_means"][metric] for w in row["windows"] if w["start"]>=start and w["time"]<=stop]
    if len(values)<4:return {"blocks":len(values),"ci95":None,"target_inside":None}
    interval=mean_evidence(values,0.,.05)["ci"]
    return {"blocks":len(values),"ci95":interval,"target_inside":interval[0]<=target<=interval[1],
            "scope":"coarse2048-block mean-t approximation; not recomputed whole-chain ESS or guaranteed stationary uncertainty"}


def _secondary_recovery(groups,metadata,protocol):
    if ("q_quartic","burst_p030") not in groups:return []
    result=[];recovery=protocol["budgets"]["source_end"]
    for fault,healthy in paired_rows(groups,"q_quartic","burst_p030","iid"):
        policy=fault.get("detector_restart",{"alarm":None,"reentry":None,"right_censored":True})
        prior={name:event.get("post_onset_confirmation") is not None and event["post_onset_confirmation"]<recovery
               for name,event in fault["events"].get("applications",{}).items()}
        configurations=[("causal",discard,None,policy["reentry"]+discard if policy["reentry"] is not None else None,
                         fault["horizon"]) for discard in (0,2048,8192)]
        configurations += [("oracle",discard,budget,recovery+discard,min(recovery+discard+budget,fault["horizon"]))
                           for discard in (0,2048,8192) for budget in (8192,32768)]
        configurations.append(("final_rolling",0,8192,fault["horizon"]-8192,fault["horizon"]))
        for policy_name,discard,requested,start,stop in configurations:
            metrics={}
            for name,target in (("x",0.),("energy",protocol["intended_energy"])):
                a=window_estimate({"result":fault,"windows":fault["windows"]},start,stop,name)
                b=window_estimate({"result":healthy,"windows":healthy["windows"]},start,stop,name)
                metrics[name]={"fault":a,"healthy_same_fault_driven_budget":b,
                    "paired_difference":a["estimate"]-b["estimate"] if not a["censored"] and not b["censored"] else None,"target":target,
                    "coarse_batch_target_containment":_batch_interval(fault,start,stop,name,target) if not a["censored"] else {"blocks":0,"ci95":None,"target_inside":None}}
            result.append({"trial":fault["trial"],"policy_name":policy_name,"detector_policy":policy,
                "discard":discard,"requested_rows":requested,"start":start,"stop":stop,"prior_post_onset_failures":prior,"metrics":metrics,
                "achieved_rows":metrics["x"]["fault"]["available_rows"],
                "requested_budget_right_censored":requested is not None and metrics["x"]["fault"]["available_rows"]<requested,
                "minimum_usable_budget_right_censored":metrics["x"]["fault"]["censored"],
                "budget_scope":"variable achieved budget limited by horizon; requested32k plus discard may be right-censored as32k despite an observed shorter>=8192 mean",
                "scope":"secondary finite-budget estimate; observed means conditional on availability and selected prior failure where stated; no instantaneous equilibration claim"})
    return result


def _secondary_recovery_summaries(rows):
    grouped=defaultdict(list);targets={}
    strata=("all_trials","cumulative_prior_failure","recent_prior_failure","both_position_prior_failures")
    for row in rows:
        targets[(row["stage"],row["policy_name"],row["discard"],row["requested_rows"])]= {name:entry["target"] for name,entry in row["metrics"].items()}
        for stratum in strata:
            grouped[(row["stage"],row["policy_name"],row["discard"],row["requested_rows"],stratum)]
            prior=row["prior_post_onset_failures"]
            selected=(stratum=="all_trials" or
                (stratum=="cumulative_prior_failure" and prior.get("cumulative",False)) or
                (stratum=="recent_prior_failure" and prior.get("recent",False)) or
                (stratum=="both_position_prior_failures" and prior.get("cumulative",False) and prior.get("recent",False)))
            if selected:grouped[(row["stage"],row["policy_name"],row["discard"],row["requested_rows"],stratum)].append(row)
    output=[]
    for (stage,policy,discard,requested,stratum),selected in grouped.items():
        metric_summary={}
        for name in ("x","energy"):
            observed=[r["metrics"][name] for r in selected if not r["metrics"][name]["fault"]["censored"]]
            target=targets[(stage,policy,discard,requested)][name]
            available=[r["fault"]["available_rows"] for r in observed]
            containment=[r["coarse_batch_target_containment"]["target_inside"] for r in observed
                         if r["coarse_batch_target_containment"]["target_inside"] is not None]
            metric_summary[name]={"observed":len(observed),"right_censored":len(selected)-len(observed),
                "available_rows_min":min(available) if available else None,"available_rows_max":max(available) if available else None,
                "fault_mean_evidence":mean_evidence([r["fault"]["estimate"] for r in observed],target,.05),
                "paired_same_budget_evidence":mean_evidence([r["paired_difference"] for r in observed],0.,.05),
                "coarse_batch_target_containment":_cp(sum(containment),len(containment)),
                "scope":"descriptive pointwise95% ensemble interval conditional on availability/stratum; not another primary test"}
        output.append({"stage":stage,"policy":policy,"discard":discard,"requested_rows":requested,"stratum":stratum,
            "stratum_trials":len(selected),"metrics":metric_summary,
            "requested_budget_right_censored_trials":sum(r["requested_budget_right_censored"] for r in selected),
            "budget_scope":"mean evidence uses explicitly achieved horizon-limited rows, never imputed requested-budget rows",
            "premature_detector_reentry":sum(bool(r["detector_policy"].get("premature_relative_to_oracle")) for r in selected) if policy=="causal" else None})
    return output


def frozen_analysis(protocol_path,output_dir):
    tick=time.perf_counter()
    protocol_path=Path(protocol_path).resolve();output_dir=Path(output_dir)
    raw=protocol_path.read_bytes();protocol_hash=hashlib.sha256(raw).hexdigest()
    if protocol_hash!=EXPECTED_PROTOCOL_SHA256:raise ValueError("frozen V04protocol checksum changed")
    protocol=json.loads(raw);project=protocol_path.parents[2];checks=[]
    for name,expected in protocol["artifact_sha256"].items():
        actual=hashlib.sha256((project/name).read_bytes()).hexdigest();checks.append({"path":name,"matches":actual==expected,"sha256":actual})
    if not all(c["matches"] for c in checks):raise ValueError("frozen scientific artifact changed")
    # Refuse partial labels: all four planned stages must be complete first.
    for stage in protocol["stage_plan"]:
        directory=project/stage["output"]
        if not (directory/"manifest.json").exists() or len(list((directory/"trials").glob("*.json.gz")))!=stage["executions"]:
            raise ValueError("definitive analysis requires all complete frozen stages; no partial success labels")
    stages=[];adaptive=[];prior_seeds=set();prior_streams=set()
    for stage in protocol["stage_plan"]:
        config_path=project/stage["config"]
        if hashlib.sha256(config_path.read_bytes()).hexdigest()!=stage["config_sha256"]:raise ValueError("frozen config changed")
        groups,meta=load_records(project/stage["output"],require_complete=True,
                                 control_epoch=protocol["budgets"]["onset"] if stage["label"].startswith("main_") else None)
        if meta["config"]!=json.loads(config_path.read_text()):raise ValueError("manifest differs from frozen stage config")
        if meta["config"]["seed_namespace"]!=stage["namespace"] or meta["config"]["trials"]!=stage["trials_per_cell"]:
            raise ValueError("fresh namespace/budget differs from frozen plan")
        if meta["completed"]!=stage["executions"] or meta["independent_seed_identities"]!=stage["independent_seed_groups"]:
            raise ValueError("execution/independent-identity count differs from frozen plan")
        for rows in groups.values():
            if len({tuple(r["seed"]) for r in rows})!=len(rows) or len({r["stream_sha256"] for r in rows})!=len(rows):
                raise ValueError("repeated seed/stream within frozen cell; no independent inference")
            if any(r["namespace"]!=stage["namespace"] for r in rows):raise ValueError("wrong checkpoint namespace")
        seeds={tuple(r["seed"]) for rows in groups.values() for r in rows}
        streams={r["stream_sha256"] for rows in groups.values() for r in rows}
        if seeds & prior_seeds or streams & prior_streams:raise ValueError("fresh stages reuse seeds or exact full streams")
        prior_seeds.update(seeds);prior_streams.update(streams)
        primary=_primary_claims(groups,meta,protocol,stage["label"])
        stages.append({"stage":stage["label"],"namespace":stage["namespace"],"claims":primary,
            "provenance":{k:v for k,v in meta.items() if k not in {"manifest","config","calibration"}},
            "healthy_detector_rates":_healthy_detector_rates(groups,protocol) if stage["label"].startswith("main_") else []})
        for row in _secondary_recovery(groups,meta,protocol):adaptive.append({"stage":stage["label"],**row})
        # Descriptive summaries preserve all secondary cells, not replacements.
        describe(project/stage["output"],write_csv=True,records=(groups,meta))
        del groups,meta,seeds,streams,rows
    replicated=[]
    for claim in protocol["primary_claims"]:
        passes={s["stage"]:next(c["pass"] for c in s["claims"] if c["id"]==claim["id"]) for s in stages if s["stage"] in claim["stages"]}
        replicated.append({"id":claim["id"],"stage_passes":passes,"independently_replicated":len(passes)==2 and all(passes.values()),
                           "scope":"every frozen component in both applicable stages; no pooling or replacements"})
    result={"version":"v04-frozen-analysis-1","protocol_sha256":protocol_hash,
        "analysis_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"created_unix":time.time(),
        "analysis_elapsed_seconds":time.perf_counter()-tick,
        "peak_analysis_rss_mib_linux":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        "artifact_checks":checks,"stages":stages,"replicated_claims":replicated,
        "total_executions":sum(s["provenance"]["completed"] for s in stages),
        "independent_seed_groups":sum(s["provenance"]["independent_seed_identities"] for s in stages),
        "secondary_adaptive_recovery":adaptive,
        "secondary_recovery_summaries":_secondary_recovery_summaries(adaptive),
        "limitations":["Mean intervals use independent-trial ensemble approximations; numerical allowances widen intervals, not equivalence margins.",
            "Matched classical/quantum and no-reseed counterfactuals reuse inputs and do not create additional independent evidence.",
            "Source pipeline operational counts againstPCG calibration are descriptive, not per-strategy validated1%certificates.",
            "Raw-source clocks remain distinct from application proposal clocks; absent source onset is not a healthy designation.",
            "Detector restart is causal; observed conditional estimates and budget censoring are not unconditional policy reliability.",
            "No microscopic instantaneous-equilibrium, hardware entropy, NISTcompliance or cryptographic-security claim is established."]}
    output_dir.mkdir(parents=True,exist_ok=True);atomic_json(output_dir/"CONFIRMATION_ANALYSIS.json",result)
    lines=["# V0.4 frozen independent confirmation","",f"Protocol: `{protocol_hash}`.","",
        f"Executed {result['total_executions']:,} frozen trials using {result['independent_seed_groups']:,} independent seed groups with intentional paired configurations.","",
        "Each compound claim requires all frozen conditions in both applicable stages.",""]
    for row in replicated:
        labels=", ".join(f"{k}: {'passed' if v else 'not established'}" for k,v in row["stage_passes"].items())
        lines.append(f"- {row['id']}: {labels}; independently replicated: **{row['independently_replicated']}**.")
    lines += ["","Component estimates, intervals, numerical padding, discordant paired counts, nulls and censoring are retained in the companion JSON.",
              "Secondary results never replace a failed primary. Detector-driven estimates have explicitly variable budgets and conditional availability.",""]
    (output_dir/"CONFIRMATION_ANALYSIS.md").write_text("\n".join(lines))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--describe",type=Path)
    parser.add_argument("--protocol",type=Path,default=Path("docs/v0.4/CONFIRMATION_PROTOCOL.json"))
    parser.add_argument("--output",type=Path,default=Path("results/v0.4/confirmatory-analysis"))
    args=parser.parse_args()
    result=describe(args.describe) if args.describe else frozen_analysis(args.protocol,args.output)
    print(json.dumps(result.get("metadata",result.get("replicated_claims",{})),indent=2))


if __name__=="__main__":main()
