# RNGLab post-run scientific audit

23 October 2026 | Evidence-first assessment | Analysis is separate from original experiments

The GPU did useful work. This campaign did not establish a new physical discovery or validate RNGLab as a general integrity-monitoring system. It was a standalone scientific-compute suite: eight-dimensional integration, a heat-equation stencil, and synthetic neural-network regression. The old RNGLab v0.9.1 service is failed and its configured repository is absent on this VM. There are no new Ising, VMC, PIMC, physical-RNG or forecasting experiments in these runs.

## Strongest supported findings

1. Healthy Monte Carlo uncertainty is well calibrated for this smooth IID-style integration problem. CPU coverage is 94.99%, GPU coverage 94.89%; cluster confidence intervals include 95%. Convergence exponents are -0.500 CPU and -0.478 GPU, consistent with the familiar inverse-square-root law.
2. A consumed +0.02 clipped-uniform shift changes the integral from 0.0967713381 to 0.0844191863: a 12.76% reduction. No shifted-batch interval covers the original target. Intervals still mostly cover the deliberately shifted target. This demonstrates the importance of freezing the scientific target; it is not a newly discovered RNG phenomenon or evidence of an implemented operational detector.
3. Heat-solver floating-point error dominates its tiny discretization error, but remains small in relative terms. A new 3.51-second bounded verification reduced 4096-square, 1,024-step RMS residual from 3.82e-6 in FP32 to 1.37e-14 in FP64 using a separately written kernel. No RNG is used by that workload.
4. The L4 provides useful, workload-specific acceleration. Active GPU MC median is 247.4 million samples/s versus 13.58 million per CPU worker and 203.1 million for the measured 15-worker aggregate. The resulting comparisons are 18.22x versus one worker and 1.22x versus the aggregate, with different precision/RNG implementations and timed scopes. BF16 training is 3.86x faster per logged step; equal-step validation error shows no demonstrated improvement.

## Risks and unresolved findings

Restarted campaigns replay the same hybrid stream prefixes. For example, 313,964 of 313,978 matched CPU rows from the preceding stopped run reproduce exactly in the completed run; the remaining partial final records use different sample counts. These restarts are not independent replications. Hashes identify source versions, but earlier harness versions are not fully preserved and one short run has no original run seal. The baseline records a misleading 3,600-second configuration even though its observed execution lasted about 215 seconds.

Shifted GPU intervals over-cover their shifted target: 96.57%, versus nominal 95%, with standardized variance about 0.895. Two related endpoints survive Holm adjustment within an eight-test exploratory family. The effect weakens later in the run. This is an unresolved finite-sample calibration anomaly, not confirmed PRNG failure. A new seed-disjoint, frozen calibration study is the highest-value next experiment.

## Explicit answers

| Question | Evidence-based answer |
| --- | --- |
| What did we discover? | Reproducible engineering behavior: calibrated healthy MC, known bias propagation hidden by retargeted scoring, precision-dominated heat residuals, workload-specific acceleration, and one preliminary shifted-GPU calibration anomaly. |
| What did we fail to demonstrate? | Novel physics, superior RNG quality, entropy guarantees, general detection/forecast reliability, MCMC validity, a new ML precision advantage, or an unexpected research breakthrough. |
| Strongest evidence? | Analytic MC moments with cluster coverage; exact seed replay; analytic heat eigenmode plus paired precision control; repeated measured BF16 throughput. |
| Publishable conclusions? | Not a standalone discovery paper as it stands. A transparent methods/benchmark contribution might become publishable with prespecified matched comparisons, fresh seed replication and a stronger research question. |
| Did L4 accelerate work? | Yes: MC 18.22x one worker / 1.22x measured CPU aggregate, and BF16 about 3.86x FP32 training rate. These are observed workload comparisons, not universal speedups. |
| Principal bottlenecks? | Serial GPU phases, repeated MC allocations/reductions/synchronizations, diagnostic work excluded from stencil rates, fixed ML order and sparse timing, plus incomplete profiling. Exact bottleneck attribution is not measured. |
| Hidden failures? | No malformed JSON, nonfinite rows, nonzero final-run exits, recorded OOM events or seal mismatches. Scientific accounting/provenance limitations and two explicitly stopped runs remain. |
| Contradiction with RNGLab history? | None established. Healthy IID calibration does not contradict previous correlated-MCMC undercoverage or observable masking. Those workloads were not rerun. |
| Thesis relevance? | Useful numerical-validation and target-integrity case studies; limited support for GPU methods. Not new evidence for PIMC/VMC, physical RNG sources or predictive useful horizons. |
| Highest-value next experiment? | Seed-disjoint, matched-input CPU/GPU precision-by-shift calibration with an immutable original target and prespecified variance/coverage endpoints. |

## Preservation and scope

Inventoried 329 readable files, 811,449 JSONL rows across five runs; verified 235 existing per-run checksum entries with no mismatch. All 329 original hashes are unchanged after analysis. One root-only SSH backup was recorded as unreadable and excluded; it is unrelated to experiment science. No heavy campaign was restarted. New bounded verification is clearly separated from historical observations. Implementation checks were performed algorithmically within the validation suite, not by independent external peer review.

See `SCIENTIFIC_FINDINGS.md`, `STATISTICAL_AUDIT.md`, the integrated PDF and machine-readable `evidence/` and `tables/` for numerical details.
