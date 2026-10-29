# v0.16 forensic reproduction

All entries in the existing v0.16 seal verified without changing any old file. Reconstructed the original rolling64/two-crossing endpoint from the saved healthy Ising traces. The reported49/128 pre-cutoff events and57/79 subsequent eligible events are exact. They arise from an oversensitive operational definition applied to correlated finite windows, not evidence that healthy PRNGs broke.

The saved test Brier contrast is recomputed from raw probabilities; theory adds no improvement. OOD deterioration and zero permitted continuations are directly preserved in analysis.json. These are failed predictive-usefulness/utility tests, not a software arithmetic defect.

Strong acceptance bias can make two magnetized states absorbing: loss of ergodicity is genuine structural damage caused by that intervention. Interpreting it through a unique stationary equilibrium is a model-assumption failure. Corrected harmonic PIMC formulas remain valid, but26/32 healthy CV coverage in the last v0.16 sample is below the nominal95% target; this is an uncertainty/mixing concern, not an RNG defect. VMC constant local energy at the exact state is known mathematics and cannot certify position sampling.

The two-state perturbation identity and augmented-driver controls passed prior checks within their stated domain. v0.15 remains a proposal in the inspected archive, not an accepted experimental release. New CALIBRATE evidence uses separate sources and seeds; no old result is reclassified as independent confirmation.

{
  "checks": [
    "all765_existing_sealed_entries_verified"
  ],
  "healthy_pre_cutoff_failures": 49,
  "healthy_eligible": 79,
  "healthy_future_failures": 57,
  "total_healthy_trajectories": 128,
  "theory_minus_combined_brier": 2.183129554352934e-05,
  "source_checksums_unchanged": true,
  "causal_interpretation": "v016 label records repeated finite-window deviation, not a causal defect finding"
}
