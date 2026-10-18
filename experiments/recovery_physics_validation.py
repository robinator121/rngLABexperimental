"""Bounded frozen-protocol verification, never writes historical experiment paths.

First call --freeze SPEC.json, then --run SPEC.json --out NEW_DIRECTORY.
The second call refuses existing outputs; preserved protocol/results use SHA-256.
Run with one BLAS thread. No provider, service, campaign, or network calls.
"""
from __future__ import annotations

import argparse
import dataclasses
import datetime as dt
import hashlib
import json
import math
from pathlib import Path
import resource
import sys
import time

import numpy as np
from scipy import stats

from rnglab14.models.pimc_harmonic import (
    PathIntegralHarmonicOscillator, exact_finite_bead_energy, harmonic_energy_estimators,
)
from rnglab14.models.vmc_harmonic import HarmonicOscillatorVMC
from rnglab14.monitoring.calibration import FrozenRankCalibration
from rnglab14.monitoring.dual_monitor import DualIntegrityMonitor, evaluate_alarm_trace
from rnglab14.sources.base import StochasticChannels
from rnglab14.sources.reference import PCG64Source


SPEC = {
    "schema": "rnglab.physics.recovery.v3", "classification": "prospectively_frozen_recovery_validation",
    "claim_scope": "synthetic numerical/reference controls, not new physical discoveries",
    "cpu_seconds": 120, "memory_address_space_bytes": 2_000_000_000,
    "pimc": {"beta": 2., "beads": [4, 8, 16, 32], "iid_seed": 810271,
             "iid_replications": 128, "paths_per_replication": 512,
             "markov_beads": [4, 16], "markov_replicates": 32,
             "sweeps": 3000, "burn_in": 500, "delta": .4, "markov_seed_base": 820000},
    "exact_reference_mode": {"beta":2.,"beads":[4,16],"path_counts":[512,2048],
        "replications":128,"seed_base":10100000,"confidence":.95,
        "protocol":"EXACT_FINITE_BEAD_GAUSSIAN_REFERENCE_V1",
        "source":"PCG64; fresh seeds disjoint from v1 controls",
        "primary_endpoint":"nominal95% interval coverage for TD and CV; retain failures"},
    "vmc": {"steps": 2000, "burn_in": 500, "alpha": .5, "delta": 1.,
            "calibration_count": 199, "healthy_holdout_count": 128, "fault_holdout_count": 64,
            "alpha_alarm": .05, "calibration_seed_base": 910000,
            "healthy_seed_base": 920000, "fault_seed_base": 930000,
            "fault": "proposal_uniform_shift_plus_0.08_clipped_to_[0,1); acceptance healthy"},
    "quadrants": {"horizon": 64, "failure_onset": 24, "alpha_alarm": .05,
        "calibration_count": 199, "healthy_holdout_count": 128, "per_fault_family_count": 64,
        "calibration_seed": 941001, "heldout_seed": 954001,
        "null": "independent absolute N(0,1) source and application metric trajectories",
        "families": ["source_early", "app_only", "source_late", "weak_both"],
        "combined_rule": "source rank p<=alpha/2 OR application rank p<=alpha/2",
        "unit": "one independently simulated whole fixed-horizon trajectory",
        "endpoints": "active-step hazard; detection after known fault onset; separate pre-fault false alarms",
        "limits": "synthetic metrics with explicit external failure truth; not calibrated physical models"},
    "primary_checks": ["PIMC means agree with finite-bead Gaussian quadratic reference",
        "old TD coefficient fails independent action derivative",
        "report heldout coverage, false positives, missed failures and conditional detection delay",
        "do not tune thresholds, horizons, seeds or endpoints after viewing outcomes"],
    "success_policy": "formula agreement is deterministic; heldout rates are reported even when poor",
    "stop": "CPU or address-space limit; any nonfinite sample; output collision; protocol mismatch",
}


def canonical(x):
    return json.dumps(x, sort_keys=True, indent=2, allow_nan=False).encode()


def proportion(k, n):
    if n == 0:
        return {"events": k, "independent_units": n, "rate": None, "exact_95_interval": None}
    ci = stats.binomtest(k, n).proportion_ci(confidence_level=.95, method="exact")
    return {"events": int(k), "independent_units": int(n), "rate": k/n,
            "exact_95_interval": [float(ci.low), float(ci.high)]}


def channels(seed, shifted=False):
    proposal = PCG64Source(seed=seed)
    if shifted:
        class ShiftedProposal:
            def random_float(self):
                return min(np.nextafter(1., 0.), proposal.random_float()+.08)
        prop = ShiftedProposal()
    else:
        prop = proposal
    return StochasticChannels(prop, PCG64Source(seed=seed+1), PCG64Source(seed=seed+2))


def pimc(spec):
    r = np.random.default_rng(spec["iid_seed"])
    beta = spec["beta"]
    rows = []
    for m in spec["beads"]:
        identity = np.eye(m); d = np.roll(identity, -1, axis=0)-identity; lap = d.T@d
        cov = np.linalg.inv(beta/m*identity+m/beta*lap)
        chol = np.linalg.cholesky(cov)
        b = identity/(2*m)-m*lap/(2*beta**2)
        c = identity/m-np.ones((m,m))/(2*m*m)
        b_old = identity/(2*m)-lap/(2*beta)
        exact_reference = float(m/(2*beta)+np.trace(b@cov))
        row = {"beads": m, "exact_gaussian_energy": exact_reference,
               "closed_form_finite_bead_energy": exact_finite_bead_energy(m,beta),
               "continuum_energy": .5/math.tanh(beta/2),
               "old_td_expectation": float(m/(2*beta)+np.trace(b_old@cov)),
               "corrected_td_variance": float(2*np.trace(b@cov@b@cov)),
               "centroid_virial_variance": float(2*np.trace(c@cov@c@cov)),
               "independent_replications": spec["iid_replications"],
               "independent_paths_per_replication": spec["paths_per_replication"]}
        summary = {name: {"means": [], "covers": 0} for name in ["old_td", "corrected_td", "centroid_virial"]}
        for _ in range(spec["iid_replications"]):
            paths = r.normal(size=(spec["paths_per_replication"], m))@chol.T
            td, cv = harmonic_energy_estimators(paths,beta)
            springs = np.sum((np.roll(paths,-1,axis=1)-paths)**2, axis=1)
            old = m/(2*beta)-springs/(2*beta)+np.mean(.5*paths*paths,axis=1)
            for name, values in [("old_td",old),("corrected_td",td),("centroid_virial",cv)]:
                mean = float(values.mean()); se = float(values.std(ddof=1)/math.sqrt(len(values)))
                critical = float(stats.t.ppf(.975,len(values)-1))
                summary[name]["means"].append(mean)
                summary[name]["covers"] += bool(abs(mean-exact_reference)<=critical*se)
        for name, s in summary.items():
            row[name] = {"mean_of_replication_means": float(np.mean(s["means"])),
                "between_replication_se": float(np.std(s["means"],ddof=1)/math.sqrt(len(s["means"]))),
                "iid_t_interval_coverage": proportion(s["covers"],spec["iid_replications"])}
        rows.append(row)
    markov = []
    for m in spec["markov_beads"]:
        results = []
        for j in range(spec["markov_replicates"]):
            result = PathIntegralHarmonicOscillator(m, 1/beta, spec["delta"]).run_simulation(
                channels(spec["markov_seed_base"]+10000*m+10*j), spec["sweeps"],spec["burn_in"])
            results.append(dataclasses.asdict(result))
        coverage = {}
        for estimator, field in [("centroid_virial", "virial_uncertainty"), ("thermodynamic", "thermodynamic_uncertainty")]:
            for method in ["geyer", "batch_32", "batch_128"]:
                checks = [r[field]["geyer"] if method=="geyer" else r[field]["batch_means"][method[6:]] for r in results]
                eligible = [v for v in checks if v.get("covered") is not None]
                coverage[f"{estimator}_{method}"] = {**proportion(sum(bool(x["covered"]) for x in eligible),len(eligible)),
                    "abstentions": len(checks)-len(eligible)}
        markov.append({"beads":m,"coverage":coverage,"complete_results":results,
            "limits":"32 chains per bead count; fixed burn-in, no established equilibration; undercoverage is preserved"})
    return {"independent_exact_path_controls": rows, "serial_metropolis_controls":markov}


def vmc(spec):
    model = HarmonicOscillatorVMC(spec["alpha"],spec["delta"])
    def run(seed, calibration=None, shifted=False):
        return model.run_simulation(channels(seed,shifted),spec["steps"],spec["burn_in"],calibration)
    cal_results = [run(spec["calibration_seed_base"]+10*j) for j in range(spec["calibration_count"])]
    # Use full-precision trace mean, not the rounded display bias.
    cal_scores = [abs(r.x2_uncertainty["mean"]-r.exact_x2) for r in cal_results]
    calibration = FrozenRankCalibration.fit(cal_scores, model.monitoring_protocol(spec["steps"],spec["burn_in"]),spec["alpha_alarm"])
    healthy = [run(spec["healthy_seed_base"]+10*j,calibration) for j in range(spec["healthy_holdout_count"])]
    fault = [run(spec["fault_seed_base"]+10*j,calibration,True) for j in range(spec["fault_holdout_count"])]
    def brief(r):
        return {"energy":r.energy_mean,"x2":r.x2_mean,"x2_full_precision":r.x2_uncertainty["mean"],
                "legacy_exploratory_flag":r.exploratory_masking_flag,"calibrated_alarm":r.observable_masking_detected,
                "p_value":r.masking_p_value}
    return {"calibration_digest":calibration.digest,"calibration_scores":cal_scores,
        "protocol_json":calibration.protocol_json,
        "healthy_alarm_rate":proportion(sum(r.observable_masking_detected for r in healthy),len(healthy)),
        "healthy_legacy_flag_rate":proportion(sum(r.exploratory_masking_flag for r in healthy),len(healthy)),
        "fault_alarm_rate":proportion(sum(r.observable_masking_detected for r in fault),len(fault)),
        "healthy_results":[brief(r) for r in healthy],"fault_results":[brief(r) for r in fault],
        "limits":"one fixed sampler protocol and one asymmetric-proposal alternative; validity requires exchangeable null whole chains"}


def exact_reference_mode(spec):
    rows=[]
    for beads in spec["beads"]:
        model=PathIntegralHarmonicOscillator(beads,1/spec["beta"])
        for n in spec["path_counts"]:
            results=[model.run_exact_reference(PCG64Source(seed=spec["seed_base"]+100000*beads+10*n+j),
                independent_paths=n,confidence=spec["confidence"]) for j in range(spec["replications"])]
            row={"beads":beads,"independent_paths_per_replication":n,"replications":len(results),
                "protocol":spec["protocol"],"complete_results":results}
            for field in ["thermodynamic","centroid_virial"]:
                means=np.array([r[field]["mean"] for r in results])
                row[field]={"coverage":proportion(sum(r[field]["covered"] for r in results),len(results)),
                    "mean_of_replication_means":float(means.mean()),
                    "between_replication_se":float(means.std(ddof=1)/math.sqrt(len(means))),
                    "finite_bead_reference":results[0]["finite_bead_reference"]}
            rows.append(row)
    return {"controls":rows,"limits":"one healthy PRNG family, beta2, beads4/16 and path counts512/2048; finite t intervals remain approximate"}


def quadrants(spec):
    cal_rng=np.random.default_rng(spec["calibration_seed"])
    n=spec["calibration_count"]; h=spec["horizon"]; alpha=spec["alpha_alarm"]
    cal=np.abs(cal_rng.normal(size=(n,h,2)))
    protocol={"horizon":h,"null":spec["null"],"score":"maximum over fixed-horizon trajectory"}
    source=FrozenRankCalibration.fit(np.max(cal[:,:,0],axis=1),{**protocol,"channel":"source"},alpha)
    app=FrozenRankCalibration.fit(np.max(cal[:,:,1],axis=1),{**protocol,"channel":"application"},alpha)
    rng=np.random.default_rng(spec["heldout_seed"])
    cohorts={}
    methods=["source_only","application_only","combined"]
    for family in ["healthy"]+spec["families"]:
        count=spec["healthy_holdout_count"] if family=="healthy" else spec["per_fault_family_count"]
        per_path=[]
        for _ in range(count):
            x=np.abs(rng.normal(size=(h,2))); onset=None if family=="healthy" else spec["failure_onset"]
            if family=="source_early":x[onset-8:,0]+=5.;x[onset+6:,1]+=5.
            elif family=="app_only":x[onset:,1]+=5.
            elif family=="source_late":x[onset+8:,0]+=5.;x[onset:,1]+=5.
            elif family=="weak_both":x[onset:,0]+=.8;x[onset:,1]+=1.2
            ps=np.array([source.p_value(float(v)) for v in x[:,0]])
            pa=np.array([app.p_value(float(v)) for v in x[:,1]])
            alarms={"source_only":ps<=alpha,"application_only":pa<=alpha,"combined":(ps<=alpha/2)|(pa<=alpha/2)}
            # Exercise the actual monitor's integrated calibration path.
            monitor=DualIntegrityMonitor(source_calibration=source,app_calibration=app,
                family_alpha=alpha,horizon=h)
            actual_combined=[]
            for t in range(h):
                snapshot=monitor.record_step(t,float(x[t,0]),float(x[t,1]))
                actual_combined.append(snapshot.source_alarm_active or snapshot.app_failure_active)
            if not np.array_equal(actual_combined,alarms["combined"]):
                raise AssertionError("integrated monitor differs from independent rank/Bonferroni computation")
            fault_onset=None if onset is None else onset-8 if family=="source_early" else onset
            row={"true_failure_onset":onset,"true_fault_onset":fault_onset}
            for method, flags in alarms.items():
                row[method]=evaluate_alarm_trace(flags.tolist(),fault_onset,onset)
            per_path.append(row)
        summary={}
        for method in methods:
            observed=[r[method] for r in per_path]
            delays=[r["detection_delay_conditional_on_detection"] for r in observed if r["detection_delay_conditional_on_detection"] is not None]
            summary[method]={"any_alarm_rate":proportion(sum(r["first_alarm"] is not None for r in observed),count),
                "pre_fault_or_healthy_false_alarm_rate":proportion(sum(r["pre_fault_or_healthy_false_alarm"] for r in observed),count),
                "missed_failure_rate":None if family=="healthy" else proportion(sum(r["missed_failure"] for r in observed),count),
                "median_delay_conditional_on_true_detection":None if not delays else float(np.median(delays)),
                "fraction_failure_steps_before_first_true_detection":None if family=="healthy" else float(np.mean([r["failure_steps_before_first_true_detection"] for r in observed]))/(h-spec["failure_onset"]),
                "active_failure_step_hazard_rate":None if family=="healthy" else float(np.mean([r["active_failure_step_hazard_fraction"] for r in observed])),
                "trajectories_with_any_active_hazard":None if family=="healthy" else proportion(sum(r["active_failure_steps_without_alarm"]>0 for r in observed),count)}
        cohorts[family]={"summary":summary,"paths":per_path}
    return {"matched_finite_horizon_null_target":alpha,"source_calibration_digest":source.digest,
        "application_calibration_digest":app.digest,"calibration_max_scores":cal.max(axis=1).tolist(),
        "cohorts":cohorts,"superiority_demonstrated":False,
        "limits":"conditional heldout binomial intervals given one frozen calibration; repeated times are not independent trials; source/app test metrics are synthetic; Bonferroni target matched, realized FPRs need not equal"}


def main():
    p=argparse.ArgumentParser(); modes=p.add_mutually_exclusive_group(required=True)
    modes.add_argument("--freeze",type=Path);modes.add_argument("--run",type=Path);p.add_argument("--out",type=Path)
    p.add_argument("--only",choices=["pimc","vmc","quadrants","exact_reference_mode"],nargs="+")
    a=p.parse_args()
    if a.freeze:
        with a.freeze.open("xb") as f:f.write(canonical({"frozen_at":dt.datetime.now(dt.timezone.utc).isoformat(),"specification":SPEC}))
        print(hashlib.sha256(a.freeze.read_bytes()).hexdigest());return
    document=json.loads(a.run.read_bytes());spec=document["specification"]
    if spec!=SPEC:raise ValueError("specification differs from frozen implementation contract")
    if a.out is None:raise ValueError("--out required")
    a.out.mkdir(parents=True,exist_ok=False)
    resource.setrlimit(resource.RLIMIT_CPU,(spec["cpu_seconds"],spec["cpu_seconds"]+1))
    resource.setrlimit(resource.RLIMIT_AS,(spec["memory_address_space_bytes"],spec["memory_address_space_bytes"]))
    started=time.monotonic();results={};timings={}
    source_root=Path(__file__).resolve().parents[1]
    def source_hashes():
        return {str(p.relative_to(source_root)):hashlib.sha256(p.read_bytes()).hexdigest()
            for root in [source_root/"src/rnglab14",source_root/"src/rnglab12"] for p in root.rglob("*.py")}
    source_before=source_hashes()
    for name,fn in [("pimc",pimc),("vmc",vmc),("quadrants",quadrants),("exact_reference_mode",exact_reference_mode)]:
        if a.only is not None and name not in a.only:continue
        phase=time.monotonic();results[name]=fn(spec[name]);timings[name]=time.monotonic()-phase
        print(name,"completed",round(timings[name],3),flush=True)
        with (a.out/(name+".json")).open("xb") as f:f.write(canonical(results[name]))
    manifest={"spec_sha256":hashlib.sha256(a.run.read_bytes()).hexdigest(),
        "script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "source_sha256_before":source_before,"source_sha256_after":source_hashes(),
        "phases_executed":list(results),
        "python":sys.version,"numpy":np.__version__,"timings_seconds":timings,
        "wall_seconds":time.monotonic()-started,"peak_rss_kib":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
        "artifacts":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in a.out.glob("*.json")},
        "evidence_status":"prospectively frozen bounded recovery controls; not preregistered external confirmation"}
    manifest["source_unchanged_during_execution"]=(manifest["source_sha256_before"]==manifest["source_sha256_after"])
    with (a.out/"manifest.json").open("xb") as f:f.write(canonical(manifest))


if __name__=="__main__":main()
