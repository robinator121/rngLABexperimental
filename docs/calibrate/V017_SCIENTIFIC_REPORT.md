# RNGLab v0.17 - Project CALIBRATE

CALIBRATE improves the definition and validation of failure, but is not a generally validated reliability product.

The v0.16 healthy failure counts were independently reconstructed from saved traces:49/128 before cutoff,57/79 later. New definitions separate accuracy excursions,cohort excess risk,coverage failure,structural invalidity and healthy-calibrated persistent excursions. Primary whole-path false alarms met their stated confidence ceiling in10/12 cells,not all12. This is a partial result,not a verified5% conditional guarantee.

Controlled mild-bias discrimination is supported in seven selected cells after both paired multiplicity-adjusted effect tests and healthy-validation gates. Mathematical/physical references agree,but healthy PIMC uncertainty remains wrong. IID reference controls are consistent with intended coverage; longer local chains with short blocks still under-cover.

Restricted IsingL3/Tc incident prediction improves proper scoring:combined test Brier.12078 versus prevalence.24206,with paired adjusted confidence interval for the difference entirely below the required-.01 margin. Theory features provide no incremental benefit. This test mainly separates healthy from a known gradual drift; abrupt paths were already failed. The test policy continues143/460 paths with4 events,but its declared no-mitigation cost1.117 exceeds always-stop1. OOD scoring gain is inconclusive and continuation is zero. Thus full predictive usefulness and generalization remain unsupported.

The next priority for v0.18 is a validated intervention/continuation decision study on a clearly delimited nonstructural defect population,with external independent review and improved MCMC uncertainty. Do not scale forecast sophistication before fixing those acceptance gates.

|Hypothesis|Outcome|
|---|---|
|H1|PARTIAL / GLOBAL CRITERION FAILED:10/12 cells pass|
|H2|SUPPORTED IN SEVEN SELECTED CELLS WITH VALIDATED FALSE-ALARM BOUND; no general claim|
|H3|PARTIAL:references agree; uncertainty coverage remains contradicted for several methods|
|H4|FULL HYPOTHESIS UNSUPPORTED:score improvement survives test; full decision-utility gate fails|
|H5|UNSUPPORTED:OOD gain inconclusive,no useful continuation|

Implemented:rank/persistence endpoint framework,paired causal statistics,frozen seed-separated campaigns,multiple-testing control,four-action decision logic,full traces and explicit audit. Executed:healthy/fault confirmation,conditional forecasting,cluster work controls and PIMC uncertainty probes. Validated only within stated finite protocols. Pending:failed calibration cells,cluster-specific endpoint calibration,nonlinear observable UQ,causal decomposition of PIMC undercoverage,real intervention efficacy and external personnel review.

See all twelve V017 reports,calibration_analysis.json,forecast_analysis.json,MECHANISM_EVIDENCE.json,AUDIT_RESULTS.json and physics-controls/results.json for exact counts,parameters,uncertainty and provenance.

## Healthy susceptibility counterexample

At L4, T=1.5, the exact zero-field susceptibility is 10.403227. Across 512 healthy chains the within-chain plug-in estimate averaged 0.745828, a bias of -9.657399 (between-chain 95% interval [-9.859323, -9.455475]). Subtracting the squared finite-chain magnetization is severely misleading when chains do not sample both magnetized phases. Accurate energy does not validate this response observable. This is healthy finite-run mixing/estimator error, not evidence of defective randomness.

An exploratory known-zero-mean symmetry estimator, N times mean(m²)/T, averaged 10.403680, with bias interval [-0.003270, 0.004176]. It uses exact finite-lattice zero-field spin-flip symmetry; this is not a generic estimator repair or a proof of second-moment equilibration. The post-confirmation calculation reuses saved traces, does not change any registered endpoint, and is not new confirmation. All nine cells are in SUSCEPTIBILITY_PROBE.json. Susceptibility and heat-capacity interval calibration still require separate validation.
