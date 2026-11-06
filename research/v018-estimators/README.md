# v0.18 Phase 1: estimator repair

Completed 2026-11-05. Pre-registered protocol `rnglab-v018-estimators-v1` ([PROTOCOL.md](PROTOCOL.md), frozen at `44676eb`), pre-generation manifest at `f3be97a`, [Amendment 01](AMENDMENT_01.md) committed before analysis. Full report: **[REPORT.md](REPORT.md)**.

## What was tested

Before RNGLab can attribute an anomaly to a random source, its own uncertainty estimates have to be trustworthy. This campaign tested that on two systems with exact answers:

- **Zero-field Ising susceptibility** on 2×2, 3×3 and 4×4 lattices at T = 1.5, T_c and 3.5. 32,768 healthy chains per cell, three estimators: the plug-in **P**, the symmetry-aware **S**, and a generic 16-chain jackknife **G**.
- **Long-run variance of an AR(1) process** with known finite-n variance. Batch means, overlapping batch means, two spectral estimators and Geyer's estimator, plus a batch-length diagnostic that abstains when the run is too short. 4,096 replicas per setting.

## Results

| Finding | Evidence |
|---|---|
| The v0.17 plug-in failure replicates: relative bias −0.917 at L4 T1.5 | gate S-F1 |
| **S is calibrated in all 9 cells** (coverage 0.925–0.949; bias within ±0.0007) | gates S-S |
| **G is unbiased but its interval fails in 5 of 9 cells.** At T = 1.5 coverage is 0.836–0.916 despite an inflated jackknife variance; at T_c the variance is inflated beyond tolerance | gates S-G (failed) |
| Short batches fail at long τ: BM(32) covers 0.547 at τ ≈ 199 and 0.332 at τ ≈ 662 | gates L-F1, L-F2 |
| BM, OBM, SV-B, SV-TH (b = 512) and Geyer are calibrated for τ ≤ 19 | 15 gates L-A, L-G |
| The diagnostic either abstains or gives calibrated intervals (21 of 21 evaluable gates pass). It is conservative: it abstains at τ ≈ 199, n = 8192, where BM(512) would cover 0.926 | gates L-D |
| Split-R̂ of m flags the sign non-mixing behind P's failure (2,048 of 2,048 ensembles at L4 T1.5) | gate D-F1 |
| **Control gate "controls valid" is FALSE:** 1 of 20 control intervals missed. The pre-specified adjusted view flags nothing | gates L-C, Amendment 01 |
| The independent recount agrees on 28 of 28 primary counts | `results/recount.json` |

Verdicts: H-S **not supported** (because of G), H-L **supported**, H-D **supported**.

## Scope

These are instrument checks on small lattices and an idealised AR(1), at fixed run lengths. They are not RNG-fault results, and they do not transfer to PIMC without Phase 2. S depends on an exact symmetry of the Ising model; no generic susceptibility repair is established.

## Files and reproduction

[REPRODUCE.md](REPRODUCE.md) gives the steps. The Ising integer traces (`traces/`, about 2 GB) are not in Git; their SHA-256 values are recorded in `results/ising/*.npz`. Code, protocol drafting, analysis and documents were produced with Claude Code (Anthropic) under the author's direction.
