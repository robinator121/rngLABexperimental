"""Independent scale learning and finite-horizon predictive calibration.

Stage A learns a fixed score scaling. Independently generated stage B sets each
whole-horizon cutoff. A later independent stage D validates the frozen pipeline.
No fault observations enter this module. The predictive bound is marginal over
the stage-B training sample and a future exchangeable healthy trial, conditional
on stage A; it is not a conditional guarantee for every realized fitted cutoff.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np

from rnglab.v02.runner import atomic_json, digest
from rnglab.v03.monitor import RULES, SCORE_DEFINITIONS, COMPOSITE_SCALES
from rnglab.v03.physics import make_model
from rnglab.v03.runner import calibration_key, engine_hash, model_spec, payloads, settings


def horizon_score(values, consecutive: int) -> float:
    """A strict cutoff is crossed K times consecutively iff this score exceeds it."""
    a = np.asarray(values, dtype=float)
    if a.ndim != 1 or not np.isfinite(a).all() or consecutive < 1 or len(a) < consecutive:
        raise ValueError("finite one-dimensional scores with at least K looks required")
    windows = np.lib.stride_tricks.sliding_window_view(a, consecutive)
    return float(np.max(np.min(windows, axis=1)))


def predictive_cutoff(scores, alpha: float = .01) -> dict:
    """Finite-sample upper order statistic, with strict future-score exceedance.

    For exchangeable trials, P(future score > cutoff) <= (n+1-rank)/(n+1)
    <= alpha. Ties make the bound conservative. n>=1/alpha-1 is necessary
    for a finite cutoff. Scores are not modified or interpolated.
    """
    a = np.asarray(scores, dtype=float)
    if not 0 < alpha < 1 or a.ndim != 1 or not len(a) or not np.isfinite(a).all():
        raise ValueError("finite nonempty scores and alpha in (0,1) required")
    n = len(a)
    rank = int(math.ceil((1 - alpha) * (n + 1)))
    if rank > n:
        raise ValueError(f"{n} threshold trials insufficient for alpha={alpha}; require at least {math.ceil(1/alpha - 1)}")
    return {"cutoff": float(np.partition(a, rank - 1)[rank - 1]), "rank": rank,
            "training_trials": n, "predictive_false_alarm_bound": (n + 1 - rank) / (n + 1),
            "decision": "strict score > cutoff"}


def _load_stage(directory: Path) -> dict:
    path = directory / "manifest.json"
    manifest = json.loads(path.read_text())
    identity = manifest["identity"]
    config = identity["config"]
    if identity["engine_sha256"] != engine_hash():
        raise ValueError("calibration source engine mismatch")
    scenarios = config["scenarios"]
    if len(scenarios) != 1:
        raise ValueError("one independently generated healthy scenario per calibration stage required")
    scenario = scenarios[0]
    stream = scenario["stream"]
    kind = stream.get("kind", stream.get("process", "iid"))
    if kind not in {"iid", "clean"} or any(scenario.get(k) is not None for k in ("onset", "source_recovery", "end")):
        raise ValueError("calibration sources must be unmodified healthy IID uniforms")
    if config.get("calibration_file") or config.get("_calibration"):
        raise ValueError("calibration training must not already use fitted cutoffs")
    rows = list(payloads(directory))
    workload_by_name = {w["name"]: w for w in config["workloads"]}
    expected = {(w, trial) for w in workload_by_name for trial in range(config["trials"])}
    seen, groups, hashes, seeds = set(), {}, set(), set()
    for p in rows:
        r = p["result"]
        ident = (r["workload"], r["trial"])
        if ident in seen or ident not in expected or r["scenario"] != scenario["name"] or r["status"] != "complete":
            raise ValueError("duplicate or unexpected calibration trial")
        seen.add(ident)
        w = workload_by_name[r["workload"]]
        key = calibration_key(config, w)
        if r["calibration_key"] != key or r["model_spec"] != model_spec(w) or r["seed_namespace"] != config["seed_namespace"]:
            raise ValueError("calibration model/settings/namespace mismatch")
        if r["stream_sha256"] in hashes or tuple(r["seed"]) in seeds:
            raise ValueError("repeated calibration streams/seeds are not independent evidence")
        hashes.add(r["stream_sha256"]); seeds.add(tuple(r["seed"]))
        groups.setdefault(key, []).append(p)
    if seen != expected:
        raise ValueError("incomplete calibration dataset")
    for key, group in groups.items():
        names = {p["result"]["workload"] for p in group}
        if len(names) != 1:
            raise ValueError("duplicate calibration strata under different workload names")
        group.sort(key=lambda p: p["result"]["trial"])
    return {"directory": str(directory.resolve()), "manifest": manifest, "config": config,
            "groups": groups, "hashes": hashes, "seeds": seeds,
            "manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "checkpoint_payload_sha256": digest(sorted(digest(p) for p in rows))}


def _eligible(payload, minimum: int) -> list:
    return [w for w in payload["windows"] if w["retained"] >= minimum]


def build(scale_dir, threshold_dir, output, alpha: float = .01) -> dict:
    """Build and save a frozen calibration artifact using only clean stages A/B."""
    a, b = _load_stage(Path(scale_dir)), _load_stage(Path(threshold_dir))
    if a["config"]["seed_namespace"] == b["config"]["seed_namespace"] or a["seeds"] & b["seeds"] or a["hashes"] & b["hashes"]:
        raise ValueError("scale and threshold stages must have disjoint namespaces, seeds and streams")
    if set(a["groups"]) != set(b["groups"]):
        raise ValueError("calibration source stratum/model/settings mismatch")
    strata = {}
    for key in a["groups"]:
        ag, bg = a["groups"][key], b["groups"][key]
        name = ag[0]["result"]["workload"]
        aw = next(w for w in a["config"]["workloads"] if w["name"] == name)
        s = settings(a["config"], aw)
        target = ag[0]["result"]["model_target"]
        exact_target = make_model(model_spec(aw)).target
        if not np.isfinite(target) or target != exact_target or any(p["result"]["model_target"] != target for p in ag + bg):
            raise ValueError("calibration target mismatch")
        if len(ag) < 2:
            raise ValueError("at least two independent scale trials required")
        eligible = [_eligible(p, s["min_measurement"]) for p in ag + bg]
        times = [w["time"] for w in eligible[0]]
        expected_times = [t for t in range(s["window"], s["burn_in"] + s["samples"] + 1, s["window"])
                          if t - s["burn_in"] >= s["min_measurement"]]
        if times != expected_times:
            raise ValueError("calibration application grid differs from declared budget")
        if any([w["time"] for w in rows] != times for rows in eligible) or len(times) < s["confirmation_windows"]:
            raise ValueError("calibration application checkpoint grid mismatch or insufficient looks")
        monitor_times = [m["time"] for m in ag[0]["monitor"]]
        if monitor_times != list(range(s["window"], s["burn_in"] + s["samples"] + 1, s["window"])):
            raise ValueError("calibration monitor grid differs from declared budget")
        if any([m["time"] for m in p["monitor"]] != monitor_times for p in ag + bg):
            raise ValueError("calibration monitor checkpoint grid mismatch")
        floor = max(1e-12, 128 * np.finfo(float).eps * max(1., abs(target)))
        scale, application_cutoffs, order_statistics = {}, {}, {"application": {}, "monitor": {}}
        for criterion in ("cumulative", "recent"):
            values = np.array([[w[criterion + "_mean"] - target for w in _eligible(p, s["min_measurement"])] for p in ag])
            if not np.isfinite(values).all():
                raise ValueError("nonfinite scale estimates")
            sd = np.maximum(values.std(axis=0, ddof=1), floor)
            scale[criterion] = {str(t): float(v) for t, v in zip(times, sd)}
            horizon = []
            for p in bg:
                errors = [w[criterion + "_error"] for w in _eligible(p, s["min_measurement"])]
                if not np.isfinite(errors).all() or any(e < 0 for e in errors):
                    raise ValueError("invalid application calibration errors")
                horizon.append(horizon_score(np.asarray(errors) / sd, s["confirmation_windows"]))
            order = predictive_cutoff(horizon, alpha)
            application_cutoffs[criterion] = order["cutoff"]
            order_statistics["application"][criterion] = order
        monitor_cutoffs = {}
        for rule in RULES:
            horizon = [horizon_score([m[rule] for m in p["monitor"]], s["confirmation_windows"]) for p in bg]
            order = predictive_cutoff(horizon, alpha)
            monitor_cutoffs[rule] = order["cutoff"]
            order_statistics["monitor"][rule] = order
        strata[key] = {"workload": name, "model_spec": model_spec(aw), "settings": s, "target": target,
                       "scale": scale, "scale_floor": floor, "application_cutoffs": application_cutoffs,
                       "monitor_cutoffs": monitor_cutoffs, "alpha": alpha,
                       "trainingcounts": {"scale_trials": len(ag), "threshold_trials": len(bg)},
                       "order_statistics": order_statistics,
                       "provenancehashes": {"scale_manifest": a["manifest_sha256"], "threshold_manifest": b["manifest_sha256"],
                                            "scale_checkpoints": a["checkpoint_payload_sha256"], "threshold_checkpoints": b["checkpoint_payload_sha256"]}}
    created = time.time()
    artifact = {"version": "v03-independent-predictive-calibration-1", "engine_sha256": engine_hash(),
                "calibration_builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "strata": strata,
                "created_unix": created,
                "created_before_fault_runs": {"freeze_time_unix": created, "fault_data_used": False,
                                              "scope": "future experiments referencing this frozen artifact; earlier exploratory/historical studies remain separate"},
                "scale_stage": {"directory": a["directory"], "namespace": a["config"]["seed_namespace"], "manifest_sha256": a["manifest_sha256"]},
                "threshold_stage": {"directory": b["directory"], "namespace": b["config"]["seed_namespace"], "manifest_sha256": b["manifest_sha256"]},
                "decision": "strict score > frozen cutoff for K consecutive eligible looks; no backdating",
                "false_alarm_scope": "alpha separately per detector and separately per application criterion, across its whole planned horizon; no joint all-rule alpha claim",
                "predictive_guarantee": "marginal finite-sample exchangeability bound conditional on independent stage A; realized-cutoff false-alarm probability assessed on fresh stage D",
                "monitor_definitions": SCORE_DEFINITIONS, "fixed_composite_scales": COMPOSITE_SCALES}
    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        raise ValueError("frozen calibration already exists; choose a new artifact path")
    atomic_json(out, artifact)
    return artifact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scale", required=True, type=Path)
    parser.add_argument("--threshold", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--alpha", type=float, default=.01)
    args = parser.parse_args()
    artifact = build(args.scale, args.threshold, args.output, args.alpha)
    print(json.dumps({"strata": len(artifact["strata"]), "alpha": args.alpha,
                      "engine_sha256": artifact["engine_sha256"], "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
