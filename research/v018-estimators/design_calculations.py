#!/usr/bin/env python3
"""Design numbers for the v0.18 Phase 1 estimator protocol (no simulation).

RNGLab v0.18 (written by Claude Code, 2026-11-05). Computes, from formulas only:
- AR(1) integrated autocorrelation times and exact finite-n variances of the mean;
- the first-order relative bias of non-overlapping batch means for AR(1);
- the Clopper-Pearson 95% half-width at each planned replica count and the exact
  probability that a coverage gate [0.92, 0.98] fails when true coverage is
  0.95 (chance failure) or 0.93 / 0.90 (sensitivity);
- a runtime projection from research/v018-estimators/smoke/smoke_timing.json.

Usage: python research/v018-estimators/design_calculations.py > research/v018-estimators/design.json
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import beta, binom

HERE = Path(__file__).resolve().parent
PHIS = [0.0, 0.5, 0.9, 0.99, 0.997]
NS = [8192, 65536]
BATCHES = [32, 128, 512, 1024]
GATE = (0.92, 0.98)


def exact_var_mean(phi: float, n: int) -> float:
    """Var of the mean of a stationary unit-variance AR(1) of length n."""
    if phi == 0.0:
        return 1.0 / n
    return ((1 + phi) / (1 - phi) - 2 * phi * (1 - phi ** n) / (n * (1 - phi) ** 2)) / n


def bm_relative_bias(phi: float, b: int) -> float:
    """First-order E[BM]/sigma^2 - 1 for AR(1): -2 sum k rho_k / (b tau)."""
    tau = (1 + phi) / (1 - phi)
    return -2 * phi / (1 - phi) ** 2 / (b * tau)


def cp_interval(k: int, n: int) -> tuple[float, float]:
    lo = 0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1)
    hi = 1.0 if k == n else beta.ppf(0.975, k + 1, n - k)
    return lo, hi


def gate_failure_probability(n: int, p: float) -> float:
    ks = np.arange(n + 1)
    inside = np.array([GATE[0] <= cp_interval(k, n)[0] and cp_interval(k, n)[1] <= GATE[1] for k in ks])
    return float(binom.pmf(ks[~inside], n, p).sum())


def main() -> None:
    ar1 = []
    for phi in PHIS:
        tau = (1 + phi) / (1 - phi)
        for n in NS:
            ar1.append({"phi": phi, "n": n, "tau_int": tau, "n_over_tau": n / tau,
                        "exact_var_mean": exact_var_mean(phi, n),
                        "asymptotic_var_mean": tau / n,
                        "run_adequate_under_rule_n_ge_160_tau": n >= 160 * tau,
                        "bm_first_order_relative_bias": {str(b): bm_relative_bias(phi, b) for b in BATCHES}})
    gates = []
    for units in (1024, 2048, 4096, 32768):
        k = round(0.95 * units)
        lo, hi = cp_interval(k, units)
        gates.append({"units": units, "cp95_halfwidth_at_0.95": (hi - lo) / 2,
                      "P_gate_fail_true_0.95": gate_failure_probability(units, 0.95),
                      "P_gate_fail_true_0.93": gate_failure_probability(units, 0.93),
                      "P_gate_fail_true_0.90": gate_failure_probability(units, 0.90)})
    smoke = json.loads((HERE / "smoke" / "smoke_timing.json").read_text())
    per_chain = {row["L"]: row["seconds_per_chain"] for row in smoke["ising"] if abs(row["T"] - 1.5) < 1e-9}
    per_chain[3] = per_chain[4] * 9 / 16  # not measured; scaled by sites per sweep
    chains = 32768
    ising_seconds = sum(per_chain[L] * chains * 3 for L in (2, 3, 4))
    ar_rows = {(r["n"]): (r["generate_seconds"] + r["estimator_work_seconds"]) / r["reps"] for r in smoke["ar1"]}
    ar_seconds = sum(ar_rows[n] * 4096 * len(PHIS) for n in NS)
    print(json.dumps({
        "script": "research/v018-estimators/design_calculations.py",
        "ar1": ar1,
        "coverage_gate": {"interval": GATE, "method": "Clopper-Pearson 95% must lie inside the interval",
                          "by_units": gates},
        "runtime_projection": {
            "ising_seconds_per_chain_2048_sweeps": per_chain,
            "ising_L3_note": "L3 scaled from L4 by 9/16 sites; not measured",
            "ising_chains_per_cell": chains, "ising_cells": 9,
            "ising_single_core_seconds": ising_seconds,
            "ar1_seconds_per_replica_by_n_excluding_geyer_loops": ar_rows,
            "ar1_single_core_seconds_excluding_geyer_loops": ar_seconds,
            "source": "research/v018-estimators/smoke/smoke_timing.json"},
    }, indent=2))


if __name__ == "__main__":
    main()
