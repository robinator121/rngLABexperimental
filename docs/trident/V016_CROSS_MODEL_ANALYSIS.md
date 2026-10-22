# Cross-model comparison

VMC uses alpha=.5,3000 steps,500 burn-in,32 chains per arm; energy is exactly.5 under every sampled position. x² and x⁴ target.5 and.75. PIMC uses beta2,M4/16,3000 sweeps,500 burn-in; the corrected inherited TD/CV estimator and finite-bead reference are imported directly. Position x² targets finite-M energy; centroid² targets1/beta. Full sequential traces and last500 versus retained-window means are saved.

Proposal shift +.08 and proposal/acceptance coupling are applied to the actual consumed uniforms. Ising acceptance bias uses+.15; these numerical severities are not physically matched. Across models the common taxonomy is the affected channel, not a universal defect-strength scale.

Every metric bias, between-chain CI, finite-run interval coverage and binomial CI is in `CROSS_MODEL_TABLE.json`. Within-chain IAT/ESS use a truncated positive ACF-pair estimate and remain exploratory, especially with deterministic bead scan, short chains, nonstationarity or severe defects. Constant VMC energy has undefined IAT/ESS, stored null. An exactly correct constant energy does not establish correct sampling.

![Coverage](figures/cross_model_coverage.png)

The old recovery already established VMC masking and PIMC undercoverage. New seeds extend numerical replication, not mathematical novelty. Severe asymmetric proposals may lack a useful stationary regime; cumulative/reset differences are descriptive, not a recovery certificate.
