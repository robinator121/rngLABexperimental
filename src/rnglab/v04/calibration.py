"""Fresh independent scale/threshold fitting for V0.4 operating rules.

Artifacts are separate for PCG64 and public deterministic-key ChaCha20 cohorts.
Each criterion receives a simultaneous tolerance-confidence allocation within
the predeclared global family. No fault or holdout outcomes select cutoffs.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import time
import numpy as np

from rnglab.v02.runner import atomic_json, digest
from rnglab.v03.calibration import horizon_score
from rnglab.v03.conservative_calibration import tolerance_rank
from rnglab.v04.runner import DETECTORS, calibration_key, engine_hash, settings
from rnglab.v04.references import quadrature_targets

APPLICATIONS = ("cumulative", "recent", "energy_cumulative", "energy_recent")


def _eligible(payload, minimum):
    return [row for row in payload["windows"] if row["retained"] >= minimum]


def _load(directory):
    from rnglab.v03.runner import payloads
    directory = Path(directory)
    manifest_path = directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    identity, config = manifest["identity"], manifest["identity"]["config"]
    if identity["engine_sha256"] != engine_hash():
        raise ValueError("training engine differs from current frozen V0.4 engine")
    family = config.get("input_family")
    if family not in {"pcg64", "chacha20"}:
        raise ValueError("explicit PCG64/ChaCha20 input_family required")
    if config["seed_namespace"] < 50000 or config.get("calibration_file") or config.get("_calibration"):
        raise ValueError("fresh unfitted V0.4 calibration namespace required")
    if len(config["scenarios"]) != 1:
        raise ValueError("one independent healthy scenario per family stage required")
    scenario = config["scenarios"][0]
    if scenario.get("base_generator", "pcg64") != family or scenario.get("copula") or scenario.get("pipeline") or any(scenario.get(k) is not None for k in ("onset", "source_recovery")):
        raise ValueError("training requires unmodified healthy uniforms of the declared family")
    workloads = {w["name"]: w for w in config["workloads"]}
    expected = {(name, i) for name in workloads for i in range(config["trials"])}
    seen, hashes, seeds, groups, checksums = set(), set(), set(), {}, []
    stratum_hashes, stratum_seeds = {}, {}
    for payload in payloads(directory):
        r = payload["result"]
        identity_pair = (r["workload"], r["trial"])
        if identity_pair not in expected or identity_pair in seen or r["status"] != "complete" or r["scenario"] != scenario["name"]:
            raise ValueError("unexpected or duplicate training identity")
        seen.add(identity_pair)
        if r["seed_namespace"] != config["seed_namespace"] or r["seed"][1] != config["seed_namespace"]:
            raise ValueError("training namespace differs from recorded seed")
        if r["source_provenance"].get("generator") != family:
            raise ValueError("actual training generator differs from declared input_family")
        key = calibration_key(config, workloads[r["workload"]])
        if r["calibration_key"] != key or r["model_target"] != 0.:
            raise ValueError("training calibration key or intended position target mismatch")
        recorded_spec = {k:v for k,v in workloads[r["workload"]].items() if k not in {"name", "seed_group"}}
        if r["model_spec"] != recorded_spec:
            raise ValueError("recorded model differs from manifest model")
        # The tolerance proof requires independent trials within each stratum,
        # not independent metrics/strata. Different kernels may intentionally
        # share trial seeds; keep that pairing visible rather than inventing
        # additional independent source identities.
        local_hashes = stratum_hashes.setdefault(key, set())
        local_seeds = stratum_seeds.setdefault(key, set())
        if r["stream_sha256"] in local_hashes or tuple(r["seed"]) in local_seeds:
            raise ValueError("reused streams/seeds within a calibration stratum")
        local_hashes.add(r["stream_sha256"]);local_seeds.add(tuple(r["seed"]))
        hashes.add(r["stream_sha256"]);seeds.add(tuple(r["seed"]))
        if r.get("fault_onset") is not None or r.get("source_recovery") is not None:
            raise ValueError("fault phase in healthy training data")
        groups.setdefault(key, []).append(payload);checksums.append(digest(payload))
    if seen != expected:
        raise ValueError("incomplete healthy training cohort")
    for group in groups.values():
        if len({p["result"]["workload"] for p in group}) != 1:
            raise ValueError("same-kernel model aliases cannot inflate independent calibration counts")
        group.sort(key=lambda p: p["result"]["trial"])
    return {"directory": str(directory.resolve()), "manifest": manifest, "config": config, "family": family,
        "groups": groups, "hashes": hashes, "seeds": seeds, "manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "payload_digest": digest(sorted(checksums))}


def build(scale_dir, threshold_dir, out_path, global_metrics=32, alpha=.01, confidence=.95):
    """Freeze four application and four monitor cutoffs from independent A/B."""
    output = Path(out_path)
    if output.exists():
        raise ValueError("frozen calibration exists; output must be a new path")
    a, b = _load(scale_dir), _load(threshold_dir)
    if a["family"] != b["family"] or a["manifest"]["identity"]["engine_sha256"] != b["manifest"]["identity"]["engine_sha256"]:
        raise ValueError("scale/threshold family or engine mismatch")
    if a["config"]["seed_namespace"] == b["config"]["seed_namespace"] or a["hashes"] & b["hashes"] or a["seeds"] & b["seeds"]:
        raise ValueError("scale and threshold cohorts must use disjoint namespaces/seeds/streams")
    if set(a["groups"]) != set(b["groups"]):
        raise ValueError("scale/threshold kernel, family, or settings strata differ")
    actual_metrics = len(a["groups"]) * (len(APPLICATIONS) + len(DETECTORS))
    if not isinstance(global_metrics, int) or global_metrics < actual_metrics:
        raise ValueError("global multiplicity family cannot omit artifact criteria")
    strata = {}
    selections = []
    for key, ag in a["groups"].items():
        bg = b["groups"][key]
        if len(ag) < 2:
            raise ValueError("at least two independent scale trials required")
        s = settings(a["config"])
        if settings(b["config"]) != s or s["confirmation_windows"] != 2:
            raise ValueError("fixed two-look confirmation/settings mismatch")
        spec = ag[0]["result"]["model_spec"]
        energy_target = float(ag[0]["result"]["targets"]["energy"])
        independent_energy = quadrature_targets(float(spec["alpha"]), float(spec.get("coupling", 0.)), float(spec.get("trial_quartic", 0.)))["energy"]
        if not np.isfinite(energy_target) or abs(energy_target - independent_energy) > 5e-11:
            raise ValueError("energy target disagrees with independent variational reference")
        if any(p["result"]["targets"]["energy"] != energy_target for p in ag + bg):
            raise ValueError("energy targets differ across healthy stages")
        expected_times = [end for end in range(s["window"], s["burn_in"] + s["samples"] + 1, s["window"]) if end - s["burn_in"] >= s["min_measurement"]]
        all_times = list(range(s["window"], s["burn_in"] + s["samples"] + 1, s["window"]))
        if len(expected_times) < 2:
            raise ValueError("insufficient eligible application looks")
        for p in ag + bg:
            if [r["time"] for r in _eligible(p, s["min_measurement"])] != expected_times or [m["time"] for m in p["monitor"]] != all_times:
                raise ValueError("training checkpoint grids differ from declared settings")
            if any(r["retained"] != max(0, r["time"] - s["burn_in"]) for r in p["windows"]):
                raise ValueError("training retained-sample clock mismatch")
        order = tolerance_rank(len(bg), alpha, confidence, global_metrics)
        selections.append(order)
        scale, app_cuts, monitor_cuts, statistics = {}, {}, {}, {"application": {}, "monitor": {}}
        for criterion in APPLICATIONS:
            target = energy_target if criterion.startswith("energy_") else 0.
            floor = max(1e-12, 128 * np.finfo(float).eps * max(1., abs(target)))
            values = np.array([[row[criterion + "_mean"] - target for row in _eligible(p, s["min_measurement"])] for p in ag], dtype=float)
            if not np.isfinite(values).all():
                raise ValueError("nonfinite scale estimates")
            sd = np.maximum(values.std(axis=0, ddof=1), floor)
            scale[criterion] = {str(end): float(value) for end, value in zip(expected_times, sd)}
            horizon = []
            for p in bg:
                rows = _eligible(p, s["min_measurement"])
                errors = np.array([row[criterion + "_error"] for row in rows], dtype=float)
                means = np.array([row[criterion + "_mean"] for row in rows], dtype=float)
                if not np.isfinite(errors).all() or np.any(errors < 0) or not np.allclose(errors, abs(means - target), rtol=1e-10, atol=floor):
                    raise ValueError("training application error/target mismatch")
                horizon.append(horizon_score(errors / sd, 2))
            cutoff = float(np.partition(horizon, order["rank"] - 1)[order["rank"] - 1])
            app_cuts[criterion] = cutoff
            statistics["application"][criterion] = {**order, "cutoff": cutoff, "target": target, "scale_floor": floor}
        for detector in DETECTORS:
            horizon = [horizon_score([m[detector] for m in p["monitor"]], 2) for p in bg]
            cutoff = float(np.partition(horizon, order["rank"] - 1)[order["rank"] - 1])
            monitor_cuts[detector] = cutoff
            statistics["monitor"][detector] = {**order, "cutoff": cutoff}
        strata[key] = {"scale": scale, "application_cutoffs": app_cuts, "monitor_cutoffs": monitor_cuts,
            "input_family": a["family"], "sampling_spec": {k:spec.get(k, default) for k, default in [("alpha", None), ("trial_quartic", 0.), ("coupling", 0.), ("proposal", 1.5)]},
            "settings": s, "position_target": 0., "energy_target": energy_target,
            "energy_interpretation": "Quantum trial-family variational observable; for matched classical paths this same energy is a counterfactual diagnostic, not classical thermodynamic energy.",
            "trainingcounts": {"scale_trials": len(ag), "threshold_trials": len(bg)}, "order_statistics": statistics,
            "provenance": {"scale_manifest_sha256": a["manifest_sha256"], "threshold_manifest_sha256": b["manifest_sha256"], "scale_payload_digest": a["payload_digest"], "threshold_payload_digest": b["payload_digest"]}}
    # A complete global certificate requires every external artifact allocated
    # in this family to follow the same budget. Dependence between metrics is fine.
    worst_failure = max(order["per_metric_coverage_failure_bound"] for order in selections)
    artifact = {"version": "v04-independent-simultaneous-tolerance-calibration-1",
        "engine_sha256": a["manifest"]["identity"]["engine_sha256"], "input_family": a["family"],
        "fitting_namespaces": [a["config"]["seed_namespace"], b["config"]["seed_namespace"]], "strata": strata,
        "created_unix": time.time(), "calibration_builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "reference_check_source_sha256": hashlib.sha256((Path(__file__).parent / "references.py").read_bytes()).hexdigest(),
        "unique_fitting_seed_identities": {"scale": len(a["seeds"]), "threshold": len(b["seeds"])},
        "healthy_only": True, "fault_or_holdout_values_used_to_select_thresholds": False,
        "certificate": {"target_realized_probability_upper": alpha, "requested_confidence": confidence,
            "global_metrics": global_metrics, "metrics_in_this_artifact": actual_metrics,
            "global_confidence_lower_bound_if_all_allocated_metrics_follow_policy": max(0., 1 - global_metrics * worst_failure),
            "this_artifact_confidence_lower_bound": max(0., 1 - actual_metrics * worst_failure),
            "scope": "Each individual confirmed whole-horizon event probability; not an alpha bound on any event over all rules.",
            "assumptions": "Conditionally IID healthy trial scores from fixed distributions given independent A scales and frozen settings. Finite exchangeability alone supports only marginal prediction. Separate generator artifacts must all use the allocated global budget.",
            "predictive_distinction": "Each fitted order has a recorded marginal predictive false-alarm bound; neither that average bound nor the realized upper target is a claim of exact equality."},
        "decision": "Strict score>cutoff for two consecutive eligible looks; confirmation times never backdated.",
        "scale_stage": {"directory": a["directory"], "manifest_sha256": a["manifest_sha256"]},
        "threshold_stage": {"directory": b["directory"], "manifest_sha256": b["manifest_sha256"]},
        "source_scope": "Independent PCG64/ChaCha20 numerical cohorts. Source-pipeline controls must be evaluated separately; this certificate is not entropy estimation or cryptographic security proof."}
    output.parent.mkdir(parents=True, exist_ok=True);atomic_json(output, artifact)
    return artifact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scale", required=True, type=Path);parser.add_argument("--threshold", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path);parser.add_argument("--global-metrics", type=int, default=32)
    parser.add_argument("--alpha", type=float, default=.01);parser.add_argument("--confidence", type=float, default=.95)
    args = parser.parse_args();artifact = build(args.scale, args.threshold, args.output, args.global_metrics, args.alpha, args.confidence)
    print(json.dumps({"family": artifact["input_family"], "strata": len(artifact["strata"]), "output": str(args.output), "certificate": artifact["certificate"]}, indent=2))


if __name__ == "__main__":
    main()
