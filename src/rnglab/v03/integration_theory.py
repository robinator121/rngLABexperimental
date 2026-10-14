"""Bounded independent-trial verification of copy-process integration theory.

Run with ``python -m rnglab.v03.integration_theory --output ... --workers 2``.
The design is frozen in plan.json before trial execution. Cells checkpoint
atomically; completed cells resume without duplicating evidence. This module
does not alter the main physical-workload runner or historical V0.2 results.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import time

import numpy as np
import scipy
from scipy import stats

from rnglab.v03.uniforms import UniformProcess, copy_variance_inflation


MASTER = 202610093
INTEGRANDS = ("x2", "indicator_half", "cos_2pi")
MEANS = np.array([1/3, .5, 0.])
VARIANCES = np.array([4/45, .25, .5])
NS = (128, 4096, 65536)
RHOS = (0., .2, .8, .98)


def _atomic_json(path, data):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def _source_hashes():
    paths = (Path(__file__), Path(__file__).with_name("uniforms.py"))
    return {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def _cell_identity(phase, n, rho):
    return f"{phase}-n{n}-rho{rho:.2f}"


def _cell_trial(task):
    """One complete experimental cell; at most two workers run independently."""
    output, phase, stage, n, rho, trials, hashes = task
    output = Path(output)
    identity = _cell_identity(phase, n, rho)
    path = output / "cells" / f"{identity}.npz"
    fingerprint = hashlib.sha256(json.dumps({"phase": phase, "stage": stage, "n": n,
        "rho": rho, "trials": trials, "hashes": hashes}, sort_keys=True).encode()).hexdigest()
    if path.exists():
        with np.load(path, allow_pickle=False) as old:
            if str(old["fingerprint"]) != fingerprint:
                raise ValueError(f"existing cell {identity} belongs to a different frozen design/source")
        return {"cell": identity, "resumed": True}
    estimates = np.empty((trials, 3))
    seconds = np.empty(trials)
    stream_hashes = np.empty(trials, dtype="U64")
    start = time.perf_counter()
    n_index = NS.index(n)
    for trial in range(trials):
        before = time.perf_counter()
        # Same seed at each rho => paired source substreams within each N.
        # Different N, stage, and trial => distinct independent seed blocks.
        seed = [MASTER, stage, n_index, trial]
        provider = UniformProcess({"kind": "copy", "rho": rho}, seed, 1,
                                  channels=("proposal",))
        values = provider.sample(n)[:, 0]
        estimates[trial] = (np.mean(values*values), np.mean(values < .5),
                            np.mean(np.cos(2*np.pi*values)))
        stream_hashes[trial] = hashlib.sha256(values.astype("<f8").tobytes()).hexdigest()
        seconds[trial] = time.perf_counter() - before
    elapsed = time.perf_counter() - start
    temporary = path.with_suffix(".npz.tmp")
    with temporary.open("wb") as file:
        np.savez_compressed(file, estimates=estimates, elapsed_seconds=seconds,
            stream_sha256=stream_hashes, fingerprint=fingerprint,
            seed_prefix=np.array([MASTER, stage, n_index]), n=n, rho=rho,
            cell_elapsed_seconds=elapsed,
            worker_peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    os.replace(temporary, path)
    return {"cell": identity, "resumed": False, "seconds": elapsed,
            "peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}


def _holm(pvalues):
    pvalues = np.asarray(pvalues, dtype=float)
    order = np.argsort(pvalues)
    adjusted = np.maximum.accumulate((len(pvalues)-np.arange(len(pvalues))) * pvalues[order])
    result = np.empty_like(pvalues)
    result[order] = np.minimum(1., adjusted)
    return result


def _variance_bootstrap(values, replicates, seed):
    """Percentile CI over truly independent trial means, no Gaussian assumption."""
    rng = np.random.default_rng(seed)
    n = len(values)
    estimates = np.empty((replicates, values.shape[1]))
    means = np.empty_like(estimates)
    # Blocks bound temporary allocation to ~3*256*n*8 bytes per process.
    for offset in range(0, replicates, 256):
        count = min(256, replicates-offset)
        indices = rng.integers(0, n, size=(count, n))
        samples = values[indices]
        estimates[offset:offset+count] = samples.var(axis=1, ddof=1)
        means[offset:offset+count] = samples.mean(axis=1)
    return np.quantile(estimates, [.025, .975], axis=0), np.quantile(means, [.025, .975], axis=0)


def analyze(output, plan):
    """Frozen descriptive and inferential analysis, regenerated from cell files."""
    rows, paired = [], []
    trials = plan["trials_per_cell"]
    alpha = plan["family_alpha"]
    family_size = len(plan["cells"]) * 3
    quantile = stats.t.ppf(1-alpha/(2*family_size), trials-1)
    phase = plan["phase"]
    raw_path = output / "trial_results.csv"
    with raw_path.open("w", newline="") as raw:
        writer = csv.writer(raw)
        writer.writerow(["phase", "stage", "cell", "n", "rho", "trial", "seed", "stream_sha256",
            "estimate_x2", "estimate_indicator_half", "estimate_cos_2pi", "runtime_seconds",
            "emitted_uniforms", "base_draws", "decision_draws", "initial_draws",
            "uniform_source_sha256"])
        for cell_index, cell in enumerate(plan["cells"]):
            n, rho = cell["n"], cell["rho"]
            identity = _cell_identity(phase, n, rho)
            with np.load(output / "cells" / f"{identity}.npz", allow_pickle=False) as data:
                estimates = data["estimates"].copy()
                durations = data["elapsed_seconds"].copy()
                hashes = data["stream_sha256"].copy()
                total_elapsed = float(data["cell_elapsed_seconds"])
                rss = float(data["worker_peak_rss_mib"])
            for trial in range(trials):
                seed = [MASTER, plan["stage"], NS.index(n), trial]
                writer.writerow([phase, plan["stage"], identity, n, rho, trial,
                    json.dumps(seed, separators=(",", ":")), hashes[trial],
                    *estimates[trial].tolist(), durations[trial], n, n, n, 9,
                    plan["source_hashes"]["uniforms.py"]])
            inflation = copy_variance_inflation(rho, n)
            predicted = VARIANCES * inflation/n
            measured_var = estimates.var(axis=0, ddof=1)
            empirical_mean = estimates.mean(axis=0)
            se = np.sqrt(measured_var/trials)
            ci_variance, ci_mean_boot = _variance_bootstrap(estimates,
                plan["bootstrap_replicates"], [MASTER, plan["stage"], 9001, cell_index])
            for k, integrand in enumerate(INTEGRANDS):
                statistic = (empirical_mean[k]-MEANS[k])/se[k]
                pvalue = float(2*stats.t.sf(abs(statistic), trials-1))
                single_ci_radius = stats.t.ppf(.975, trials-1)*se[k]
                simultaneous_radius = quantile*se[k]
                mean_ci = [empirical_mean[k]-single_ci_radius, empirical_mean[k]+single_ci_radius]
                simultaneous_ci = [empirical_mean[k]-simultaneous_radius, empirical_mean[k]+simultaneous_radius]
                # Frozen scientific precision goal: absolute systematic bias smaller
                # than 15% of one-trial theoretical SD, simultaneous family-wise CI.
                margin = .15 * np.sqrt(predicted[k])
                row = {"phase": phase, "stage": plan["stage"], "cell": identity,
                    "n": n, "rho": rho, "integrand": integrand, "trials": trials,
                    "analytic_mean": float(MEANS[k]), "sample_mean": float(empirical_mean[k]),
                    "mean_bias": float(empirical_mean[k]-MEANS[k]),
                    "mean_ci95": [float(v) for v in mean_ci],
                    "mean_bootstrap_ci95": ci_mean_boot[:, k].tolist(),
                    "mean_simultaneous_ci": [float(v) for v in simultaneous_ci],
                    "bias_t_statistic": float(statistic), "bias_pvalue": pvalue,
                    "bias_precision_margin": float(margin),
                    "mean_bias_within_precision_margin": bool(
                        simultaneous_ci[0] > MEANS[k]-margin and simultaneous_ci[1] < MEANS[k]+margin),
                    "theoretical_inflation": inflation, "analytic_trial_variance": float(predicted[k]),
                    "sample_trial_variance": float(measured_var[k]),
                    "variance_ratio": float(measured_var[k]/predicted[k]),
                    "variance_ratio_bootstrap_ci95": (ci_variance[:, k]/predicted[k]).tolist(),
                    "variance_ratio_within20pct_precision": bool(
                        ci_variance[0, k]/predicted[k] > .8 and ci_variance[1, k]/predicted[k] < 1.2),
                    "cell_runtime_seconds": total_elapsed, "mean_trial_runtime_seconds": float(durations.mean()),
                    "worker_peak_rss_mib": rss}
                rows.append(row)
            if rho != 0 and any(c["n"] == n and c["rho"] == 0 for c in plan["cells"]):
                with np.load(output / "cells" / f"{_cell_identity(phase, n, 0)}.npz") as control:
                    differences = estimates-control["estimates"]
                for k, integrand in enumerate(INTEGRANDS):
                    delta = float(differences[:, k].mean())
                    stderr = float(differences[:, k].std(ddof=1)/np.sqrt(trials))
                    radius = float(stats.t.ppf(.975, trials-1)*stderr)
                    paired.append({"phase": phase, "n": n, "rho": rho, "integrand": integrand,
                        "paired_trials": trials, "mean_difference_from_iid": delta,
                        "ci95": [delta-radius, delta+radius], "paired_standard_error": stderr,
                        "pairing": "same base/decision/init substreams, same N/trial; independent across trials"})
    adjusted = _holm([row["bias_pvalue"] for row in rows])
    for row, pvalue in zip(rows, adjusted):
        row["bias_pvalue_holm"] = float(pvalue)
        row["reject_zero_bias_holm"] = bool(pvalue < alpha)
    _atomic_json(output / "summary.json", {"plan": plan,
        "executed_cells": len(plan["cells"]), "trial_executions": trials*len(plan["cells"]),
        "independent_seed_blocks": trials*len({cell["n"] for cell in plan["cells"]}),
        "independence_note": "1000 independent trials per cell; rho configurations are paired, integrands share streams",
        "rows": rows, "paired_comparisons": paired,
        "mean_bias_rejections_holm": sum(row["reject_zero_bias_holm"] for row in rows),
        "bias_precision_passes": sum(row["mean_bias_within_precision_margin"] for row in rows),
        "variance20pct_precision_passes": sum(row["variance_ratio_within20pct_precision"] for row in rows),
        "inference_notes": ["Mean Student-t tests rely on the independent-trial ensemble CLT, not IID within-stream samples.",
            "Holm adjustment is over all cell/integrand bias tests in this phase.",
            "Bootstrap variance intervals are percentile intervals over independent trials; 95% intervals are descriptive, not simultaneous.",
            "Failure to reject zero bias is not proof of unbiasedness; the analytic marginal argument establishes the ideal-model mean.",
            "Precision criteria were frozen before execution; no cell was selected or removed based on its result."]})
    table_fields = [key for key in rows[0] if not isinstance(rows[0][key], (list, dict))]
    with (output / "summary.csv").open("w", newline="") as table:
        writer = csv.DictWriter(table, table_fields)
        writer.writeheader()
        writer.writerows({key: row[key] for key in table_fields} for row in rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--phase", choices=("exploration", "confirmation"), default="exploration")
    parser.add_argument("--trials", type=int, default=1000)
    parser.add_argument("--workers", type=int, choices=(1, 2), default=2)
    parser.add_argument("--bootstrap", type=int, default=5000)
    parser.add_argument("--analyze-only", action="store_true")
    args = parser.parse_args()
    if args.trials < 2 or args.bootstrap < 100:
        parser.error("at least two trials and 100 bootstrap replicates required")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "cells").mkdir(exist_ok=True)
    stage = 12001 if args.phase == "exploration" else 12002
    cells = ([{"n": n, "rho": rho} for n in NS for rho in RHOS] if args.phase == "exploration"
             else [{"n": n, "rho": rho} for n in (4096, 65536) for rho in (0., .8, .98)])
    plan = {"version": "v03-integration-theory-1", "phase": args.phase, "master_seed": MASTER,
        "stage": stage, "cells": cells, "trials_per_cell": args.trials,
        "integrands": list(INTEGRANDS), "bootstrap_replicates": args.bootstrap,
        "family_alpha": .05, "mean_bias_precision_sd_fraction": .15,
        "variance_ratio_precision_interval": [.8, 1.2], "seed_scheme": "[master,stage,N-index,trial]; same across rho",
        "source_hashes": _source_hashes(), "marginal_proof": "stationary uniform initialization and independent copy/refresh decisions",
        "precision": "53-bit dyadic PCG base; continuous analytic targets; grid correction <=~1e-16 in means",
        "primary_analysis": "exact finite-N variance, independent-trial bootstrap CI, Holm-adjusted two-sided mean tests",
        "confirmation_rule": "fresh stage12002, rho0/.8/.98 at N4096/65536 fixed before exploratory results",
        "numpy": np.__version__, "scipy": scipy.__version__, "python": platform.python_version(),
        "architecture": platform.machine()}
    plan_path = output / "plan.json"
    if plan_path.exists():
        previous = json.loads(plan_path.read_text())
        if previous != plan:
            raise ValueError("existing frozen plan differs; use a new output directory")
    else:
        _atomic_json(plan_path, plan)
        source_dir = output / "source_snapshot"
        source_dir.mkdir(exist_ok=True)
        for name in plan["source_hashes"]:
            shutil.copyfile(Path(__file__).with_name(name), source_dir / name)
    start = time.perf_counter()
    statuses = []
    if not args.analyze_only:
        tasks = [(str(output), args.phase, stage, cell["n"], cell["rho"], args.trials,
                  plan["source_hashes"]) for cell in cells]
        if args.workers == 1:
            for task in tasks:
                result = _cell_trial(task)
                statuses.append(result)
                print(json.dumps(result), flush=True)
        else:
            with ProcessPoolExecutor(max_workers=args.workers) as executor:
                for result in executor.map(_cell_trial, tasks):
                    statuses.append(result)
                    print(json.dumps(result), flush=True)
    rows = analyze(output, plan)
    _atomic_json(output / "execution.json", {"elapsed_seconds": time.perf_counter()-start,
        "workers": args.workers, "statuses": statuses, "parent_peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        "analyze_only": args.analyze_only, "cpu_policy": "at most two local workers",
        "blas_threads": {key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}})
    print(json.dumps({"phase": args.phase, "rows": len(rows),
        "mean_bias_rejections_holm": sum(row["reject_zero_bias_holm"] for row in rows),
        "bias_precision_passes": sum(row["mean_bias_within_precision_margin"] for row in rows),
        "variance_precision_passes": sum(row["variance_ratio_within20pct_precision"] for row in rows)}), flush=True)


if __name__ == "__main__":
    main()
