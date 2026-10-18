"""RNGLab v0.14: Reproducible Integrated Scientific Experiment.

Executes a complete end-to-end multi-domain study:
1. Interchangeable Sources: Reference PCG64, ChaCha20/8, Physical GM, Hash Conditioned, Degraded Faults
2. Decoupled Stochastic Channels across 5 Scientific Applications:
   - MC Pi Estimation
   - 2D Ising Statistical Mechanics (L=8, T=2.269)
   - Variational Monte Carlo (QMC Harmonic Oscillator with Observable Masking check)
   - Path Integral Monte Carlo (PIMC Quantum Ring Polymer)
   - Monte Carlo Decision Support System (DSS Investment Options)
3. Dual Monitoring (Source Health vs Application Integrity in 4-Quadrant Matrix)
4. Corrected RUS Forecasting with Censoring Separation and Uncertainty Abstention
5. Multi-Testing Significance Adjustments (Holm-Bonferroni & Benjamini-Hochberg FDR)
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import datetime
import hashlib
import json
import logging
import math
from pathlib import Path
import sys

from rnglab14.forecast import (
    CUSUMChangePointPredictor,
    LinearExtrapolationPredictor,
    LipschitzStructuralBoundPredictor,
    NoChangeBaselinePredictor,
    TrajectoryFamily,
    TrajectoryGenerator,
    evaluate_forecast_on_trajectory,
    summarize_forecasting_performance,
)
from rnglab14.models import (
    HarmonicOscillatorVMC,
    Ising2D,
    MonteCarloDecisionSupport,
    PathIntegralHarmonicOscillator,
    run_pi_monte_carlo,
)
from rnglab14.monitoring import DualIntegrityMonitor, MonitoringQuadrant
from rnglab14.sources import (
    ChaCha20Source,
    DegradedFaultSource,
    GeigerMullerPulseSource,
    HashBlockConditioner,
    PCG64Source,
    StochasticChannels,
    VonNeumannConditioner,
)
from rnglab14.standards.health_90b import AdaptiveProportionTest, RepetitionCountTest
from rnglab14.statistics import benjamini_hochberg_fdr, holm_bonferroni_correction

logger = logging.getLogger("rnglab14.experiment")


def run_experiment(output_dir: Path, base_seed: int = 101):
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = output_dir / "experiment_manifest.json"
    results_path = output_dir / "experiment_results.json"

    logger.info("Starting RNGLab v0.14 Integrated Scientific Experiment...")
    started_at = datetime.datetime.now(datetime.timezone.utc).isoformat()

    # 1. Source Inventory
    sources = {
        "PCG64_Reference": PCG64Source(seed=base_seed),
        "ChaCha20_FullRounds": ChaCha20Source(rounds=20),
        "ChaCha8_ReducedRounds": ChaCha20Source(rounds=8),
        "GeigerMuller_Pulse": GeigerMullerPulseSource(mean_cps=500.0, seed=base_seed),
        "Hash_Conditioned": HashBlockConditioner(PCG64Source(seed=base_seed + 1), input_block_bytes=64),
        "Degraded_MarkovBias": DegradedFaultSource(
            PCG64Source(seed=base_seed + 2), bias_delta=0.08, autocorrelation_rho=0.40
        ),
    }

    # 2. Application Runs with Decoupled Channels
    app_results = {}

    for src_name, src in sources.items():
        logger.info(f"Evaluating application models on source: {src_name}...")
        channels = StochasticChannels(
            proposal=src,
            acceptance=PCG64Source(seed=base_seed + 10),
            initialization=PCG64Source(seed=base_seed + 20),
        )

        # A. Pi MC
        pi_res = run_pi_monte_carlo(channels, samples=10000)

        # B. 2D Ising
        ising = Ising2D(L=8, temperature=2.269)
        ising_res = ising.run_simulation(channels, sweeps=300, burn_in=50)

        # C. VMC Harmonic Oscillator (Ground state alpha = 0.5)
        vmc = HarmonicOscillatorVMC(alpha=0.5, delta=1.0)
        vmc_res = vmc.run_simulation(channels, steps=3000, burn_in=500)

        # D. PIMC Harmonic Oscillator (Beads M = 16, T = 0.5)
        pimc = PathIntegralHarmonicOscillator(beads=16, temperature=0.5, displacement_delta=0.4)
        pimc_res = pimc.run_simulation(channels, sweeps=1500, burn_in=300)

        # E. Decision Support System
        dss = MonteCarloDecisionSupport()
        dss_res = dss.evaluate_decisions(channels, scenarios=2000)

        app_results[src_name] = {
            "pi_mc": asdict(pi_res),
            "ising_2d": asdict(ising_res),
            "vmc_harmonic": asdict(vmc_res),
            "pimc_harmonic": asdict(pimc_res),
            "decision_support": asdict(dss_res),
        }

    # 3. Dual Integrity Monitoring Simulation (4-Quadrant Study)
    logger.info("Executing Dual Integrity Monitoring analysis...")
    dual_monitor = DualIntegrityMonitor(source_alarm_threshold=0.05, app_failure_threshold=0.10)

    for t in range(1, 101):
        # Simulate ramp fault starting at t = 40
        if t < 40:
            src_defect = 0.01 + 0.005 * (t / 40.0)
            app_error = 0.02 + 0.01 * (t / 40.0)
        elif t < 70:
            # Source defect crosses threshold at t = 50, application crosses at t = 65 -> Early Warning window!
            src_defect = 0.02 + 0.003 * (t - 40)
            app_error = 0.03 + 0.003 * (t - 40)
        else:
            # Full failure
            src_defect = 0.15
            app_error = 0.20
        dual_monitor.record_step(t, src_defect, app_error)

    dual_report = dual_monitor.summarize()

    # 4. Corrected RUS Forecasting Benchmark
    logger.info("Executing corrected RUS Forecasting benchmark...")
    forecasting_reports = {}
    trajectories_to_eval = [
        TrajectoryFamily.LINEAR_DRIFT,
        TrajectoryFamily.STEPWISE_FAILURE,
        TrajectoryFamily.CONSTANT_BIAS,
    ]
    eval_times = [100, 150, 200, 250, 300]
    threshold = 0.40

    predictors = [
        ("NoChangeBaseline", NoChangeBaselinePredictor(threshold=threshold)),
        ("LinearExtrapolation", LinearExtrapolationPredictor(threshold=threshold, window=50)),
        ("CUSUMChangePoint", CUSUMChangePointPredictor(threshold=threshold, baseline_mean=0.0)),
        ("LipschitzStructuralBound", LipschitzStructuralBoundPredictor(threshold=threshold, lipschitz_L=0.001)),
    ]

    for traj_fam in trajectories_to_eval:
        gen = TrajectoryGenerator(traj_fam, seed=base_seed + 50)
        path = gen.generate_path(total_steps=500)
        vals = [p.parameter for p in path]

        for pred_name, pred in predictors:
            fcasts = evaluate_forecast_on_trajectory(
                pred, vals, threshold=threshold, evaluation_times=eval_times, model_name=pred_name
            )
            rep = summarize_forecasting_performance(fcasts)
            forecasting_reports[f"{traj_fam.value}_{pred_name}"] = asdict(rep)

    # 5. Multi-Testing Statistical Significance Control
    logger.info("Executing multi-testing error rate corrections...")
    # Gather nominal p-values from Pi MC Z-scores across all evaluated sources
    p_values = []
    test_labels = []
    for src_name, res in app_results.items():
        z = abs(res["pi_mc"]["z_score"])
        # Two-tailed p-value: 2 * (1 - Phi(z))
        from math import erf
        phi = 0.5 * (1.0 + erf(z / math.sqrt(2.0)))
        p_val = max(1e-12, 2.0 * (1.0 - phi))
        p_values.append(p_val)
        test_labels.append(f"Pi_MC_{src_name}")

    holm_res = holm_bonferroni_correction(p_values, alpha=0.05)
    fdr_res = benjamini_hochberg_fdr(p_values, q_fdr=0.05)

    multi_test_report = {
        "nominal_p_values": {label: round(p, 6) for label, p in zip(test_labels, p_values)},
        "holm_bonferroni_rejections": {
            label: dec for label, dec in zip(test_labels, holm_res.reject_decisions)
        },
        "benjamini_hochberg_fdr_rejections": {
            label: dec for label, dec in zip(test_labels, fdr_res.reject_decisions)
        },
        "total_hypotheses": len(p_values),
        "fwer_alpha": 0.05,
    }

    # Compile Final Structured Result
    completed_at = datetime.datetime.now(datetime.timezone.utc).isoformat()
    final_payload = {
        "experiment_title": "RNGLab v0.14 Integrated Scientific Verification",
        "started_at": started_at,
        "completed_at": completed_at,
        "base_seed": base_seed,
        "sources_evaluated": list(sources.keys()),
        "application_results": app_results,
        "dual_monitoring_report": asdict(dual_report),
        "forecasting_reports": forecasting_reports,
        "multi_testing_report": multi_test_report,
    }

    payload_json = json.dumps(final_payload, indent=2, sort_keys=True)
    payload_hash = hashlib.sha256(payload_json.encode("utf-8")).hexdigest()
    final_payload["provenance_sha256"] = payload_hash

    with open(results_path, "w", encoding="utf-8") as f:
        json.dump(final_payload, f, indent=2)

    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(
            {
                "experiment_id": "v014_integrated_slice",
                "sha256": payload_hash,
                "timestamp": completed_at,
                "output_file": str(results_path),
            },
            f,
            indent=2,
        )

    logger.info(f"Experiment completed successfully. SHA-256: {payload_hash}")
    return final_payload


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(name)s: %(message)s")
    out_dir = Path("runs/v014-integrated-experiment").resolve()
    run_experiment(out_dir)


if __name__ == "__main__":
    main()
