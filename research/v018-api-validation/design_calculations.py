#!/usr/bin/env python3
"""Design numbers for the v0.18 Phase 2 API-validation protocol (no simulation).

RNGLab v0.18 (written by Claude Code, 2026-11-05). Computes, from formulas and
the saved Phase II confirmation outputs only:
- the exact probability that the coverage gate (Clopper-Pearson 95% interval
  inside [0.92, 0.98]) passes, as a function of true coverage, for each
  candidate replica count R;
- an approximate pass probability of the variance-ratio gate (bootstrap 95%
  interval inside [0.75, 1.25]), from a log-normal approximation whose width
  is checked against the Phase II bootstrap intervals;
- the bias-gate interval half-width implied by the Phase II between-replica
  variances;
- the approximate power of the paired refreshed-minus-local comparison
  (lower 95% bound > 0.10), using the Phase II discordant-pair rates;
- the approximate probability that both required conditions pass, under
  perfect calibration and under the Phase II point estimates;
- expected source bytes per replica and role;
- runtime and trace-storage projections from smoke/smoke_timing.json.

Usage: PYTHONPATH=src python research/v018-api-validation/design_calculations.py > research/v018-api-validation/design.json
"""
import json
import math
from pathlib import Path

from scipy.stats import beta, binom, norm, t

from rnglab14.models.pimc_harmonic import exact_finite_bead_energy

HERE = Path(__file__).resolve().parent
PHASE2 = HERE.parents[1] / "research" / "phase2-20261030" / "experiment" / "confirmation-v1" / "primary_decision.json"
RS = [1024, 2048, 4096]
COVERAGE_GATE, RATIO_GATE, BIAS_GATE, MARGIN = (0.92, 0.98), (0.75, 1.25), 0.01, 0.10
BATCHES = 16  # BM512 on n = 8192
WORKERS, N, TRACES_PER_REPLICA = 4, 8192, 4  # CV, TD, centroid, rg2 as float64


def cp(k, n):
    return (0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1), 1.0 if k == n else beta.ppf(0.975, k + 1, n - k))


def coverage_pass(p, r):
    """Exact P(CP95 interval of Binomial(r, p) lies inside the gate)."""
    ks = [k for k in range(r + 1) if COVERAGE_GATE[0] <= cp(k, r)[0] and cp(k, r)[1] <= COVERAGE_GATE[1]]
    return float(binom.cdf(max(ks), r, p) - binom.cdf(min(ks) - 1, r, p)), min(ks), max(ks)


def log_ratio_se(r):
    """SE of log(mean BM variance / empirical variance) for normal replica means:
    empirical variance ~ chi2(r-1)/(r-1), each BM variance ~ chi2(15)/15."""
    return math.sqrt(2 / (r - 1) + 2 / (BATCHES - 1) / r)


def ratio_pass(rho, r):
    se, z = log_ratio_se(r), norm.ppf(0.975)
    lo = (math.log(RATIO_GATE[0]) + z * se - math.log(rho)) / se
    hi = (math.log(RATIO_GATE[1]) - z * se - math.log(rho)) / se
    return float(max(0.0, norm.cdf(hi) - norm.cdf(lo)))


def paired_pass(delta, discordant, r):
    """P(lower 95% bound of the mean paired coverage difference > MARGIN)."""
    se = math.sqrt((discordant - delta ** 2) / r)
    return float(norm.sf((MARGIN + norm.ppf(0.975) * se - delta) / se))


def main():
    p2 = json.loads(PHASE2.read_text())
    rows = {k: v["row"] for k, v in p2["gates"].items()}
    cmp_ = p2["comparison"]
    discordant = (cmp_["improved_only"] + cmp_["local_only"]) / 1024
    smoke = json.loads((HERE / "smoke" / "smoke_timing.json").read_text())
    sec = smoke["mean_seconds_per_replica"]

    coverage = {str(r): {f"{p:.3f}": coverage_pass(p, r)[0] for p in (0.92, 0.93, 0.935, 0.94, 0.945, 0.95, 0.955, 0.96, 0.97)}
                for r in RS}
    accept_k = {str(r): {"min_covered": coverage_pass(0.95, r)[1], "max_covered": coverage_pass(0.95, r)[2]} for r in RS}
    ratio = {str(r): {f"{rho:.3f}": ratio_pass(rho, r) for rho in (0.80, 0.85, 0.90, 0.944, 1.0, 1.1)} for r in RS}
    ratio_check = {c: {"phase2_ratio": rows[c]["estimated_to_empirical_variance_ratio"],
                       "phase2_bootstrap_ci95": rows[c]["variance_ratio_bootstrap_ci95"],
                       "approximation_ci95": [rows[c]["estimated_to_empirical_variance_ratio"] * math.exp(s * norm.ppf(0.975) * log_ratio_se(1024))
                                              for s in (-1, 1)]}
                   for c in ("centroid_refresh", "iid_gaussian")}
    bias = {c: {str(r): float(t.ppf(0.975, r - 1) * math.sqrt(rows[c]["empirical_variance_chain_means"] / r)) for r in RS}
            for c in ("centroid_refresh", "iid_gaussian")}
    paired = {str(r): {f"{d:.3f}": paired_pass(d, discordant, r) for d in (0.12, 0.13, cmp_["paired_coverage_difference"])} for r in RS}
    local_k = rows["local_stationary"]["covered"] / 1024
    local_fail = {str(r): float(binom.cdf(max(k for k in range(r + 1) if cp(k, r)[1] < COVERAGE_GATE[0]), r, local_k)) for r in RS}
    joint = {}
    for r in RS:
        perfect = coverage_pass(0.95, r)[0] * ratio_pass(1.0, r)
        refresh = (coverage_pass(rows["centroid_refresh"]["coverage"], r)[0]
                   * ratio_pass(rows["centroid_refresh"]["estimated_to_empirical_variance_ratio"], r))
        joint[str(r)] = {"both_required_pass_perfect_calibration": perfect ** 2,
                         "refresh_pass_at_phase2_point_estimates": refresh,
                         "both_required_pass_refresh_at_phase2_iid_perfect": refresh * perfect}
    runtime = {str(r): {"single_core_seconds": r * sum(sec.values()), f"wall_minutes_{WORKERS}_workers": r * sum(sec.values()) / WORKERS / 60}
               for r in RS}
    # Bytes per replica from the package mappings: 8 per proposal word, 16 per
    # centroid draw, 16 per Box-Muller pair of standard normals.
    expected_bytes = {"start_two_paths": 16 * ((2 * 16 + 1) // 2), "proposal": 8 * 16 * N,
                      "centroid": 16 * N, "iid_paths": 16 * ((16 * N + 1) // 2), "initialization_in_run": 0}
    storage = {str(r): {"trace_GiB": 3 * r * TRACES_PER_REPLICA * N * 8 / 2**30, "cv_only_GiB": 3 * r * N * 8 / 2**30} for r in RS}

    print(json.dumps({
        "script": "research/v018-api-validation/design_calculations.py",
        "inputs": {"phase2_primary_decision": str(PHASE2.relative_to(HERE.parents[1])),
                   "smoke": "research/v018-api-validation/smoke/smoke_timing.json"},
        "reference_energy": {"finite_bead_E_M": exact_finite_bead_energy(16, 2.0), "continuum": 0.5 / math.tanh(1.0)},
        "phase2_reference": {c: {"covered": rows[c]["covered"], "replicas": rows[c]["replicas"], "coverage": rows[c]["coverage"],
                                 "variance_ratio": rows[c]["estimated_to_empirical_variance_ratio"], "bias": rows[c]["bias"],
                                 "median_estimated_tau": rows[c]["median_estimated_tau"]} for c in rows},
        "phase2_paired": {**cmp_, "discordant_rate": discordant},
        "coverage_gate_pass_probability_exact": coverage,
        "coverage_gate_acceptance_region_k": accept_k,
        "variance_ratio_gate_pass_probability_approx": ratio,
        "variance_ratio_approximation_check_R1024": ratio_check,
        "bias_gate_ci95_half_width_at_phase2_variance": bias,
        "paired_comparison_pass_probability_approx": paired,
        "local_undercoverage_reproduction_probability_at_phase2_rate": local_fail,
        "joint_required_pass_probability_approx": joint,
        "joint_note": "products of component probabilities; bias gate treated as passing with probability 1 (half-width << 0.01)",
        "expected_bytes_per_replica": expected_bytes,
        "runtime_projection": runtime,
        "trace_storage_projection": storage,
    }, indent=1))


if __name__ == "__main__":
    main()
