"""Stage D pilot design calculations; writes design.json, no sampler data.

    PYTHONPATH=src .venv/bin/python research/v019-d-pilot/design.py
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta, binom

from rnglab.seeds import namespace_seed

HERE = Path(__file__).resolve().parent
SOURCE = "research/v018-h1-detection/results/report_numbers.json"
H, WORKERS, N_CAL = 40960, 8, 999
LENGTHS = (10240, 20480, 40960)
COVERAGES = (0.92, 0.93, 0.94, 0.945, 0.95, 0.955, 0.96, 0.97)
SKILLS = (0.05, 0.10, 0.15, 0.20)
ELIGIBILITY = (0.6, 0.8, 1.0)
DATASETS, BOOTSTRAPS, PROFILES = 1000, 2000, 100000
LEVEL = 1 - 0.05 / 3


def rng(purpose):
    return np.random.default_rng(namespace_seed("rnglab-v019-d-pilot", "design", purpose, "00000"))


def coverage(n):
    k = np.arange(n + 1)
    lo, hi = np.zeros(n + 1), np.ones(n + 1)
    lo[1:] = beta.ppf(0.025, k[1:], n - k[1:] + 1)
    hi[:-1] = beta.ppf(0.975, k[:-1] + 1, n - k[:-1])
    accepted = k[(lo >= 0.92) & (hi <= 0.98)]
    probabilities = []
    for p in COVERAGES:
        passed = float(binom.pmf(accepted, n, p).sum())
        probabilities.append({"true_coverage": p, "p_pass": passed,
                              "all_three_if_independent": passed**3})
    return {"n": n, "cp_level": 0.95, "band": [0.92, 0.98],
            "accepted_counts": accepted.tolist(), "probabilities": probabilities,
            "joint_note": "p_pass^3 is only 'if independent'. Lengths are positively correlated in reality; "
                          "this is not an upper bound or an estimate of their joint pass probability. "
                          "Pairwise coverage correlation alone does not determine joint band acceptance."}


def power(n_eval):
    minimum = 120 * n_eval // 200
    scenarios = []
    for e in ELIGIBILITY:
        n = round(n_eval * e)
        for skill in SKILLS:
            purpose = f"power-eval{n_eval}-e{e:.1f}-skill{skill:.2f}"
            draws = rng(purpose + "-data")
            d = np.sqrt(skill) / 2
            q_model = 0.5 + d * (2 * draws.integers(0, 2, size=(DATASETS, n)) - 1)
            y = draws.binomial(1, q_model)
            correct = (y == (q_model > 0.5)).sum(axis=1)
            # BS_base = .25; BS_model = .25 + d^2 - d*(2*correct/n - 1).
            # Resampling pairs is exactly Binomial(n, correct/n) here: the
            # paired losses take only two values, irrespective of sign or y.
            observed = 4 * d * (2 * correct / n - 1) - 4 * d**2
            boot_correct = rng(purpose + "-bootstrap").binomial(
                n, correct[:, None] / n, size=(DATASETS, BOOTSTRAPS))
            boot_skill = 4 * d * (2 * boot_correct / n - 1) - 4 * d**2
            lower = np.quantile(boot_skill, (1 - LEVEL) / 2, axis=1, method="linear")
            passed = (observed >= 0.10) & (lower > 0) & (n >= minimum)
            p = float(passed.mean())
            scenarios.append({"eligibility_fraction": e, "eligible_count": n,
                              "true_brier_skill": skill, "q_base": 0.5,
                              "q_model_values": [float(0.5 - d), float(0.5 + d)],
                              "passed_datasets": int(passed.sum()), "power": p,
                              "monte_carlo_se": float(np.sqrt(p * (1 - p) / DATASETS)),
                              "data_seed_purpose": purpose + "-data",
                              "bootstrap_seed_purpose": purpose + "-bootstrap"})
    return {"evaluation_per_family": n_eval, "minimum_eligible": minimum,
            "minimum_note": "120 at n_eval=200; scaled proportionally (60 at 100, 240 at 400).",
            "datasets_per_scenario": DATASETS, "bootstrap_resamples_per_dataset": BOOTSTRAPS,
            "interval_level": LEVEL, "lower_percentile": (1 - LEVEL) / 2,
            "interval_method": "Two-sided percentile paired bootstrap, NumPy linear quantile.",
            "pass_rule": "eligible >= minimum AND observed skill >= 0.10 AND lower bound > 0",
            "construction": "Independent equiprobable signs Z; q_model=0.5+Z*d; Y|Z~Bernoulli(q_model); "
                            "q_base=0.5. Both forecasts are calibrated, with marginal outcome rate 0.5. "
                            "E[BS_base]=0.25, E[BS_model]=0.25-d^2, so d=sqrt(skill)/2. "
                            "F0+ is represented by a constant baseline (uninformative alarm time). "
                            "Predictions are fixed, with no fitting uncertainty or refitting.",
            "eligibility_note": "Conditional power at fixed round(e*n_eval) eligible runs, independent of "
                                "outcome; does not model random eligibility counts or selection effects.",
            "scope": "Per-family power, identical assumptions for SHIFT, LAG1 and DUP; "
                     "not the probability that at least one family passes.",
            "bootstrap_note": "Binomial counts resample the two observed loss-pair types exactly, "
                              "equivalent in distribution to resampling eligible trajectories.",
            "scenarios": scenarios}


def ramp_average(length, onset, cap_time):
    # Positive offsets are 1,...,L-1-t0. Sum the uncapped arithmetic sequence
    # and then the capped terms; the t=t0 term is zero, and s_max cancels.
    last = np.maximum(length - 1 - onset, 0)
    uncapped = np.minimum(last, np.floor(cap_time))
    return (uncapped * (uncapped + 1) / (2 * cap_time) + last - uncapped) / length


def ramps():
    draws = rng("ramp-profiles")
    onset = draws.integers(2048, 8193, size=PROFILES)
    cap_time = np.exp(draws.uniform(np.log(2048), np.log(163840), size=PROFILES))
    quantiles = (0.05, 0.25, 0.50, 0.75, 0.95)
    return {"draws": PROFILES, "seed_purpose": "ramp-profiles",
            "onset": "Integer uniform [2048,8192], inclusive",
            "cap_time": "Independent log-uniform [2048,163840], float64",
            "cap_note": "First capped sweep is t0+ceil(T_cap); reaching the cap before L means this is < L. "
                        "Continuous crossing t0+T_cap < L is also reported separately.",
            "quantile_probabilities": list(quantiles),
            "lengths": [{"L": length,
                         "share_reaching_cap_before_L": float(np.mean(onset + np.ceil(cap_time) < length)),
                         "share_continuous_crossing_before_L": float(np.mean(onset + cap_time < length)),
                         "mean_severity_over_cap_quantiles": np.quantile(
                             ramp_average(length, onset, cap_time), quantiles, method="linear").tolist()}
                        for length in LENGTHS]}


def main():
    runtime = json.loads((HERE.parents[1] / SOURCE).read_text())["runtime"]
    stage_c_sweeps = 10240  # Stage C: research/v018-h1-detection/PROTOCOL.md §2.1.
    seconds_per_sweep = runtime["wall_minutes"] * 60 * runtime["workers"] / (
        runtime["runs_stored"] * stage_c_sweeps)
    per_run_bytes = 3 * 8 * H
    options = []
    for healthy, fault, evaluation in ((2000, 500, 200), (1000, 300, 100), (2000, 1000, 400)):
        total = N_CAL + healthy + 3 * fault
        options.append({"healthy": healthy, "fault_per_family": fault,
                        "fitting_per_family": fault - evaluation, "evaluation_per_family": evaluation,
                        "calibration": N_CAL, "total_runs": total,
                        "wall_minutes_extrapolated_8_workers": total * H * seconds_per_sweep / WORKERS / 60,
                        "storage_total_GiB": total * per_run_bytes / 2**30,
                        "criterion_1": coverage(healthy), "criterion_4": power(evaluation)})
    out = {"seed_recipe": 'namespace_seed("rnglab-v019-d-pilot", "design", purpose, "00000")',
           "runtime": {"source": SOURCE, "stage_c_wall_minutes": runtime["wall_minutes"],
                       "stage_c_workers": runtime["workers"], "stage_c_runs_stored": runtime["runs_stored"],
                       "stage_c_sweeps_per_run": stage_c_sweeps, "worker_seconds_per_sweep": seconds_per_sweep,
                       "pilot_sweeps_per_run": H, "pilot_workers": WORKERS,
                       "note": "Linear cost-per-sweep extrapolation from Stage C, not a pilot benchmark; "
                               "assumes unchanged throughput, excludes new engineering and model fitting."},
           "storage": {"bytes_per_run": per_run_bytes, "MiB_per_run": per_run_bytes / 2**20,
                       "note": "Uncompressed CV float64, TD float64, int64 marks: 3*8*H bytes per run; "
                               "excludes other logs, metadata, checkpoint statistics and container overhead."},
           "size_options": options, "ramp_profiles": ramps()}
    (HERE / "design.json").write_text(json.dumps(out, indent=1, allow_nan=False) + "\n", encoding="utf-8")
    print(f"Stage C cost: {seconds_per_sweep:.9f} worker-seconds/sweep; storage: {per_run_bytes / 2**20:.4f} MiB/run.")
    print("Criterion 1 probabilities use p=" + ",".join(f"{p:g}" for p in COVERAGES))
    for option in options:
        c1, c4 = option["criterion_1"], option["criterion_4"]
        counts = c1["accepted_counts"]
        print(f"Healthy={option['healthy']}, fault/family={option['fault_per_family']}, "
              f"eval/family={option['evaluation_per_family']}: {option['total_runs']} total runs; "
              f"{option['wall_minutes_extrapolated_8_workers']:.2f} min on 8 workers (extrapolation); "
              f"{option['storage_total_GiB']:.3f} GiB.")
        print(f"  Criterion 1: accepted k={counts[0]}..{counts[-1]}; "
              "P(pass)=" + ",".join(f"{r['p_pass']:.4f}" for r in c1["probabilities"]))
        print("  All three, if independent: " + ",".join(
            f"{r['all_three_if_independent']:.4f}" for r in c1["probabilities"]))
        print(f"  Criterion 4 per-family power (minimum eligible={c4['minimum_eligible']}; "
              "true skills=0.05,0.10,0.15,0.20):")
        for e in ELIGIBILITY:
            rows = [r for r in c4["scenarios"] if r["eligibility_fraction"] == e]
            print(f"    e={e:.1f}, n={rows[0]['eligible_count']}: " + ",".join(f"{r['power']:.3f}" for r in rows))
    print("Power: 1000 datasets/scenario, 2000 paired bootstraps/dataset; fixed eligible counts; minimum scales with n_eval.")
    print("Independence is a hypothetical comparison, not an upper bound for correlated lengths.")
    for row in out["ramp_profiles"]["lengths"]:
        print(f"Ramp L={row['L']}: cap reached={row['share_reaching_cap_before_L']:.4f}; "
              "Sbar/s_max q05,25,50,75,95=" + ",".join(f"{q:.4f}" for q in row["mean_severity_over_cap_quantiles"]))
    print("Wrote research/v019-d-pilot/design.json (design calculations only).")


if __name__ == "__main__":
    main()
