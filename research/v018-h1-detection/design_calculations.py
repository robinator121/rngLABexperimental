#!/usr/bin/env python3
"""Design numbers for the draft Stage C (H1) protocol (no simulation).

RNGLab v0.18 (written by Claude Code, 2026-11-06). Computes from formulas and
the smoke output only:
- the whole-run false-alarm calibration: order-statistic threshold from N_cal
  healthy runs, its marginal false-alarm rate and the spread of the conditional
  rate given the realized threshold (Beta law of a uniform order statistic);
- Clopper-Pearson precision of a rate at the planned replica counts;
- the Bonferroni level of the H1 existence rule and the probability that its
  three criteria hold at given true rates (exact binomial; the paired difference
  by a normal approximation);
- the healthy calibration gate pass probability at the planned healthy count;
- runtime and storage projections from smoke/smoke_timing.json.

Usage: python research/v018-h1-detection/design_calculations.py > research/v018-h1-detection/design.json
"""
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import beta, binom, norm

HERE = Path(__file__).resolve().parent
ALPHA_RUN = 0.05
N_CAL, R_HEALTHY, R_FAULT = 999, 4096, 256
FAMILIES, SEVERITIES = 4, 6
SWEEPS, WORKERS, PARALLEL_SLOWDOWN = 10240, 8, 1.35  # slowdown observed in Phase 2 (block time vs smoke)


def cp(k, n, level):
    a = 1 - level
    return (0.0 if k == 0 else beta.ppf(a / 2, k, n - k + 1), 1.0 if k == n else beta.ppf(1 - a / 2, k + 1, n - k))


def prob_lower_above(p, n, bound, level):
    ks = [k for k in range(n + 1) if cp(k, n, level)[0] > bound]
    return float(binom.sf(min(ks) - 1, n, p)) if ks else 0.0


def prob_upper_below(p, n, bound, level):
    ks = [k for k in range(n + 1) if cp(k, n, level)[1] < bound]
    return float(binom.cdf(max(ks), n, p)) if ks else 0.0


def main():
    smoke = json.loads((HERE / "smoke" / "smoke_timing.json").read_text())
    runs = smoke["runs"][1:]  # first run includes import and JIT warm-up
    per_run = float(np.mean([r["run_seconds"] + r["battery_seconds"] for r in runs]))
    k = math.ceil((N_CAL + 1) * (1 - ALPHA_RUN))
    conditional = beta(N_CAL + 1 - k, k)
    level = 1 - 0.05 / (FAMILIES * SEVERITIES)
    total_runs = N_CAL + R_HEALTHY + FAMILIES * SEVERITIES * R_FAULT
    rule = {}
    for wrong, conv, gap in ((0.7, 0.05, 0.3), (0.7, 0.1, 0.2), (0.6, 0.1, 0.3), (0.8, 0.15, 0.4)):
        se = math.sqrt((gap + 2 * conv) / R_FAULT)  # paired difference, discordant mass at most gap + 2*conv
        rule[f"wrong{wrong}_conv{conv}_gap{gap}"] = {
            "P_wrong_lower_above_0.5": prob_lower_above(wrong, R_FAULT, 0.5, level),
            "P_conv_upper_below_0.2": prob_upper_below(conv, R_FAULT, 0.2, level),
            "P_gap_lower_above_0_normal_approx": float(norm.sf(norm.ppf(1 - (1 - level) / 2) - gap / se))}
    print(json.dumps({
        "script": "research/v018-h1-detection/design_calculations.py",
        "inputs": {"smoke": "research/v018-h1-detection/smoke/smoke_timing.json"},
        "calibration": {"N_cal": N_CAL, "whole_run_alpha": ALPHA_RUN, "threshold_order_statistic": k,
                        "marginal_false_alarm_rate": (N_CAL + 1 - k) / (N_CAL + 1),
                        "conditional_rate_quantiles": dict(zip(("q025", "q50", "q975"), conditional.ppf([0.025, 0.5, 0.975]).tolist()))},
        "precision": {str(r): {f"p{p}": [float(x) for x in cp(round(p * r), r, 0.95)] for p in (0.05, 0.1, 0.5)}
                      for r in (R_FAULT, R_HEALTHY)},
        "h1_rule": {"existence_intervals": FAMILIES * SEVERITIES, "bonferroni_level": level, "replicas": R_FAULT,
                    "criteria": "wrong-rate lower bound > 0.5; conventional-alarm upper bound < 0.2; paired app-minus-conv alarm lower bound > 0",
                    "probabilities": rule},
        "healthy_gate": {"R_healthy": R_HEALTHY, "pass_probability_at_true_coverage": {
            f"{p:.3f}": float(binom.cdf(max(k2 for k2 in range(R_HEALTHY + 1) if cp(k2, R_HEALTHY, 0.95)[1] <= 0.98), R_HEALTHY, p)
                              - binom.cdf(min(k2 for k2 in range(R_HEALTHY + 1) if cp(k2, R_HEALTHY, 0.95)[0] >= 0.92) - 1, R_HEALTHY, p))
            for p in (0.93, 0.94, 0.95)}},
        "runtime": {"smoke_seconds_per_run": per_run, "total_runs": total_runs,
                    "single_core_hours": total_runs * per_run / 3600,
                    f"wall_minutes_{WORKERS}_workers_with_slowdown_{PARALLEL_SLOWDOWN}": total_runs * per_run * PARALLEL_SLOWDOWN / WORKERS / 60,
                    "pilot_runs": 3 * 10 * 32 + 6 * 32 + 199, "pilot_wall_minutes": (3 * 10 * 32 + 6 * 32 + 199) * per_run * PARALLEL_SLOWDOWN / WORKERS / 60},
        "storage": {"per_run_MiB": (2 * SWEEPS * 8 + (SWEEPS + 1) * 8) / 2**20,
                    "total_GiB": total_runs * (2 * SWEEPS * 8 + (SWEEPS + 1) * 8) / 2**30,
                    "note": "CV and TD float64 plus sweep byte marks; raw bytes regenerable by replay"},
    }, indent=1))


if __name__ == "__main__":
    main()
