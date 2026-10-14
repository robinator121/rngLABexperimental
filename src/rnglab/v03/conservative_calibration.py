"""Separate, simultaneous tolerance calibration for future confirmation runs.

This never edits the original frozen artifact or core sampler/calibration code.
The stronger realized false-alarm guarantee assumes IID healthy trial horizon
scores from fixed distributions conditional on independent stage-A scales and
settings. Finite exchangeability alone only justifies marginal prediction.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import time
import numpy as np
from scipy.stats import binom

from rnglab.v02.runner import atomic_json, digest
from rnglab.v03.calibration import _load_stage, _eligible, horizon_score
from rnglab.v03.monitor import RULES
from rnglab.v03.runner import engine_hash


def tolerance_rank(n: int, alpha: float = .01, confidence: float = .95,
                   family_size: int = 90) -> dict:
    """Smallest upper order rank with the Bonferroni tolerance coverage bound.

    For cutoff X_(r), the probability that its true strict-exceedance probability
    exceeds alpha is at most P[Binomial(n,1-alpha)>=r]. Choose this <=delta/M.
    Dependence between different metrics does not invalidate the union bound.
    """
    if not isinstance(n, (int, np.integer)) or n < 1 or not isinstance(family_size, (int, np.integer)) or family_size < 1:
        raise ValueError("positive integer sample count and family size required")
    if not 0 < alpha < 1 or not 0 < confidence < 1:
        raise ValueError("alpha and confidence must be in (0,1)")
    allowed = (1 - confidence) / family_size
    ranks = np.arange(1, n + 1)
    tails = binom.sf(ranks - 1, n, 1 - alpha)
    eligible = np.flatnonzero(tails <= allowed)
    if not len(eligible):
        minimum = math.ceil(math.log(allowed) / math.log1p(-alpha))
        raise ValueError(f"{n} trials cannot support requested simultaneous tolerance coverage; need at least {minimum}")
    index = int(eligible[0])
    rank, failure = int(ranks[index]), float(tails[index])
    return {"rank": rank, "training_trials": int(n), "target_realized_false_alarm_upper": alpha,
            "requested_simultaneous_confidence": confidence, "family_size": int(family_size),
            "per_metric_coverage_failure_bound": failure, "per_metric_allowed_failure": allowed,
            "simultaneous_confidence_lower_bound": max(0., 1 - family_size * failure),
            "marginal_predictive_false_alarm_bound": (n + 1 - rank) / (n + 1),
            "decision": "strict score > cutoff"}


def refine(original_file, threshold_dir, output, alpha: float = .01,
           confidence: float = .95) -> dict:
    """Reuse original scales and verified healthy B scores at a fixed higher rank."""
    original_path, output_path = Path(original_file), Path(output)
    original = json.loads(original_path.read_text())
    if original["engine_sha256"] != engine_hash():
        raise ValueError("original calibration engine mismatch")
    source = _load_stage(Path(threshold_dir))
    if source["manifest_sha256"] != original["threshold_stage"]["manifest_sha256"]:
        raise ValueError("threshold manifest differs from original calibration")
    if source["config"]["seed_namespace"] != original["threshold_stage"]["namespace"]:
        raise ValueError("threshold namespace differs from original calibration")
    if set(source["groups"]) != set(original["strata"]):
        raise ValueError("threshold strata differ from original calibration")
    if output_path.resolve() == original_path.resolve() or output_path.exists():
        raise ValueError("preserve original/frozen calibration; output must be a new path")
    family_size = len(original["strata"]) * (len(RULES) + 2)
    artifact = copy.deepcopy(original)
    simultaneous_failure_sum = 0.
    for key, stratum in artifact["strata"].items():
        group = source["groups"][key]
        if stratum["provenancehashes"]["threshold_checkpoints"] != source["checkpoint_payload_sha256"]:
            raise ValueError("threshold checkpoint evidence differs from original calibration")
        if len(group) != stratum["trainingcounts"]["threshold_trials"]:
            raise ValueError("threshold sample count differs from original calibration")
        if any(p["result"]["model_target"] != stratum["target"] or p["result"]["model_spec"] != stratum["model_spec"] for p in group):
            raise ValueError("threshold target/model differs from original calibration")
        selection = tolerance_rank(len(group), alpha, confidence, family_size)
        rank = selection["rank"]
        order_statistics = {"application": {}, "monitor": {}}
        k, minimum = stratum["settings"]["confirmation_windows"], stratum["settings"]["min_measurement"]
        for criterion in ("cumulative", "recent"):
            scores = []
            expected_times = sorted(stratum["scale"][criterion], key=int)
            for p in group:
                rows = _eligible(p, minimum)
                if [str(w["time"]) for w in rows] != expected_times:
                    raise ValueError("threshold application grid differs from original scale")
                standardized = [w[criterion + "_error"] / stratum["scale"][criterion][str(w["time"])] for w in rows]
                scores.append(horizon_score(standardized, k))
            cutoff = float(np.partition(scores, rank - 1)[rank - 1])
            if cutoff < original["strata"][key]["application_cutoffs"][criterion]:
                raise ValueError("refinement must not lower an original cutoff")
            stratum["application_cutoffs"][criterion] = cutoff
            order_statistics["application"][criterion] = {**selection, "cutoff": cutoff,
                "predictive_false_alarm_bound": selection["marginal_predictive_false_alarm_bound"]}
            simultaneous_failure_sum += selection["per_metric_coverage_failure_bound"]
        for rule in RULES:
            scores = [horizon_score([m[rule] for m in p["monitor"]], k) for p in group]
            cutoff = float(np.partition(scores, rank - 1)[rank - 1])
            if cutoff < original["strata"][key]["monitor_cutoffs"][rule]:
                raise ValueError("refinement must not lower an original cutoff")
            stratum["monitor_cutoffs"][rule] = cutoff
            order_statistics["monitor"][rule] = {**selection, "cutoff": cutoff,
                "predictive_false_alarm_bound": selection["marginal_predictive_false_alarm_bound"]}
            simultaneous_failure_sum += selection["per_metric_coverage_failure_bound"]
        stratum["alpha"] = alpha
        stratum["order_statistics"] = order_statistics
        stratum["conservative_tolerance_selection"] = selection
    created = time.time()
    artifact.update({
        "version": "v03-simultaneous-conservative-tolerance-calibration-1",
        "original_calibration_sha256": digest(original),
        "original_calibration_file_sha256": hashlib.sha256(original_path.read_bytes()).hexdigest(),
        "conservative_builder_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "created_unix": created,
        "created_before_fault_runs": {"freeze_time_unix": created, "fault_data_used": False,
            "scope": "future confirmation runs linked to this artifact; existing exploration retains its original calibration"},
        "refinement_trigger": "Healthy holdout validation exposed realized false-alarm variability in the marginal-predictive calibration; this safeguard uses only original healthy B scores and original A scales.",
        "original_artifact_preserved": str(original_path.resolve()),
        "simultaneous_tolerance_guarantee": {
            "target_realized_false_alarm_upper": alpha, "requested_confidence": confidence,
            "family_size": family_size,
            "achieved_confidence_lower_bound": max(0., 1 - simultaneous_failure_sum),
            "scope": "simultaneous confidence that each detector/criterion has healthy whole-horizon confirmed-event probability <=alpha; not an all-rules joint alarm-probability bound",
            "assumptions": "Conditionally IID healthy trial horizon scores from fixed distributions, given independent stage-A scales and frozen settings; finite exchangeability alone is insufficient.",
            "construction": "smallest rank r with binom.sf(r-1,n,1-alpha)<=(1-confidence)/M; strict exceedance, no interpolation",
        },
        "predictive_guarantee": "Each order statistic also has its separately recorded marginal predictive bound; for n999/r998 this is .002, not a claim that realized false-alarm probability equals .002 or .01.",
        "false_alarm_scope": "Simultaneous tolerance confidence for each individual whole-horizon detector/application criterion; a fresh independent D2 stage assesses the realized pipeline.",
    })
    output_path.parent.mkdir(parents=True, exist_ok=True)
    atomic_json(output_path, artifact)
    return artifact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--original", "--calibration", required=True, type=Path)
    parser.add_argument("--threshold", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--alpha", type=float, default=.01)
    parser.add_argument("--confidence", type=float, default=.95)
    args = parser.parse_args()
    threshold = args.threshold
    if threshold is None:
        threshold = Path(json.loads(args.original.read_text())["threshold_stage"]["directory"])
    artifact = refine(args.original, threshold, args.output, args.alpha, args.confidence)
    print(json.dumps({"output": str(args.output), "strata": len(artifact["strata"]),
                      "guarantee": artifact["simultaneous_tolerance_guarantee"]}, indent=2))


if __name__ == "__main__":
    main()
