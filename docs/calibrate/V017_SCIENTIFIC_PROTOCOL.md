# CALIBRATE scientific protocol

Local protocol frozen after32-replica runtime/correctness pilots and before fresh calibration/confirmation data. No external registry or blinded custodian. Primary hypotheses and all seed namespaces are specified by the source/config below. Historical v0.16 results inform design only.

H1 uses twelve primary observable/cell endpoints, whole-path5% target, and a deliberately distinct10% upper validation ceiling at simultaneous95% confidence. A passing result does NOT establish conditional5% risk for a frozen threshold. H2 uses exact paired discordance tests with Holm across twelve mild-bias primary endpoints; detection improvement must also exceed20percentagepoints. Excess loss and practical .1/.2 tolerance margins are separate quantities. H3 independently checks physics and measures coverage; it can remain partly unsupported even when H1 passes. H4/H5 specify combined prediction versus fixed development-prevalence baseline; other models are secondary. Forecasting is restricted to the IsingL3/Tc cell only if its adjusted healthy validation passes. Other cell failures cannot be hidden by a pooled rate.

Accuracy events, calibrated persistent excursions, cohort excess loss, coverage underperformance and analytical structural invalidity remain separate. A calibrated anomaly is not itself causal attribution. Formal marginal rank validity requires exchangeable complete healthy trajectories under the unchanged finite execution protocol; does not transfer to changed run lengths, models or shifted faults.

Counts:399 calibration and512 independent healthy confirmation chains/cell;128 paired fault chains/arm. At5% FPR,512 trials give approximate pointwise95% precision±1.9percentagepoints; simultaneous bounds are wider. H2's20pp target has approximate unadjusted80% power at these counts under ordinary paired discordance rates, with reduced power after Holm. Forecast tests with up to768 paths may be inconclusive for a .01 Brier margin, especially with rare events; no sample-size extension after results.

Runtime:two CPU workers,one numerical thread each; estimated peak memory below3GiB,1200second stage budget checked between checkpoints. Cell checkpoints include traces and result JSON. No cloud/GPU/service changes. Pilot data excluded. Finite-cell process limits are not a hard aggregate memory guarantee.

{
  "version": "0.17.0.dev1",
  "root_namespace": "CALIBRATE-20261022-v1",
  "workers": 2,
  "threads_per_worker": 1,
  "Ising": {
    "lattices": [
      2,
      3,
      4
    ],
    "temperatures": [
      1.5,
      2.269185314213022,
      3.5
    ],
    "steps": 2048,
    "burn": 512,
    "window": 512,
    "stride": 128,
    "accuracy_tolerance": [
      0.1,
      0.1,
      0.1,
      0.2,
      0.2
    ]
  },
  "quantum": {
    "models": [
      "vmc",
      "pimc4",
      "pimc16"
    ],
    "steps": 4096,
    "burn": 1024,
    "window": 512,
    "stride": 128
  },
  "counts": {
    "pilot": 32,
    "calibration": 399,
    "healthy_confirmation": 512,
    "fault_confirmation": 128
  },
  "faults": [
    "mild_bias",
    "coupled",
    "correlation",
    "periodic",
    "structural_bias",
    "unused_source"
  ],
  "hypotheses": {
    "H1": "Whole-path healthy false declaration target .05; validation simultaneous upper95<=.10 for twelve primary cells",
    "H2": "Paired harmful mild-bias declarations exceed healthy at Holm alpha.05 AND detection improvement>=.20; excess loss independently reported",
    "H3": "Physical references independently agree; interval coverage validated separately, never assumed from formula validity",
    "H4": "Frozen combined predictor beats train prevalence constant in test Brier by .01 with upper adjusted CI below -.01, and useful continuation>=.20 with risk gate and lower cost than always stop",
    "H5": "Same combined candidate retains H4 score/utility advantage in excluded OOD families"
  },
  "statistics": {
    "alpha": 0.05,
    "primary_cells": 12,
    "calibration_unit": "complete independent trajectory, one max-persistent score per observable",
    "primary_H1_bound": "one-sided exact Bonferroni .05/12 upper",
    "H2_family": "twelve mild-bias primary endpoints, Holm exact paired McNemar",
    "bootstrap": 2000,
    "rank_target": 0.05,
    "tolerance_floor": true,
    "within_path_points_are_not_replicates": true
  },
  "forecast": {
    "L": 3,
    "T": 2.269185314213022,
    "steps": 2048,
    "burn": 512,
    "window": 512,
    "stride": 128,
    "cutoff": 1408,
    "horizon": 512,
    "families": [
      "healthy",
      "abrupt",
      "drift"
    ],
    "OOD": [
      "periodic",
      "correlation"
    ],
    "counts": [
      128,
      128,
      256,
      256
    ],
    "gate": "L3/Tc whole-path H1 upper<=.10; otherwise do not run forecasting",
    "primary_candidate": "combined",
    "primary_baseline": "development prevalence constant",
    "H4_H5_CI": "paired trajectory percentile bootstrap98.75% per endpoint; approximate family adjustment",
    "continued_risk_validation_ceiling": 0.1,
    "min_validation_continued": 32,
    "continue_fraction_target": 0.2,
    "costs": {
      "continue_failure": 10,
      "stop": 1,
      "review_and_stop": 1.3,
      "intervention": 0.5,
      "intervention_failure_loss_no_mitigation": 10,
      "intervention_failure_loss_hypothetical50percent": 5
    },
    "intervention_benefit_measured": false
  },
  "limits": {
    "workers": 2,
    "expected_peak_memory_GiB": 3,
    "wall_budget_seconds": 1200,
    "cloud": false,
    "GPU": false,
    "retry": 0
  },
  "references": [
    "https://proceedings.mlr.press/v25/vovk12/vovk12.pdf",
    "https://arxiv.org/abs/0811.1729"
  ]
}
