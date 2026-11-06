# v0.18 Phase 2: the Phase II core result through the package API

Completed 2026-11-06. The pre-registered protocol is `rnglab-v018-api-validation-v1` ([PROTOCOL.md](PROTOCOL.md), frozen at `0cbefb4`). The pre-generation manifest was pushed at `ff948e2` before the first draw. Full report: **[REPORT.md](REPORT.md)**.

## What was tested

Phase II found three things with a standalone NumPy/Numba implementation of harmonic path-integral Monte Carlo (M = 16, β = 2, δ = 0.4, exact starts, n = 8192):
- local Metropolis under-covers;
- an exact centroid refresh restores calibrated uncertainty intervals;
- an IID Gaussian-path control is calibrated.

This campaign asked whether the repository's own package API reproduces that. It used `PathIntegralHarmonicOscillator.run_simulation` with the reviewed v0.18 run options, the package's PCG64 sources and the `rnglab12` uncertainty code, with fresh seeds and 4,096 replicas per condition. The primary interval was CV with batch means of length 512 (16 batches, t₁₅).

## Results

| Finding | Evidence |
|---|---|
| **The refreshed sampler is calibrated:** 3846/4096 covered (0.9390, CP95 [0.9312, 0.9461]), variance ratio 0.971 | gate G-R, pass |
| **The IID control is calibrated:** 3884/4096 (0.9482), variance ratio 0.993 | gate G-I, pass |
| Local Metropolis under-covers: 3323/4096 (0.8113) | gate F-L, pass |
| The refresh improves paired coverage by 0.1277 [0.1138, 0.1414] | gate P-RL, pass |
| All three conditions agree with Phase II within the Wald 95% interval | descriptive |
| Byte ledger, start pairing, trace consistency and replay are all exact | checks I-1 to I-4 |
| The independent recount agrees on 21 of 21 quantities | `results/recount.json` |
| **The batch-length diagnostic accepts every refreshed run but covers only 0.8828:** the "10·τ̂" rule misses a slow shape mode | secondary, negative |
| **The package's default batch sizes under-cover with the refresh:** 0.8611 at b = 32, 0.9260 at b = 128 | secondary, negative |
| **The refresh condition has a small bias whose interval excludes zero:** CV −0.000321 [−0.000567, −0.000075]; unexplained, inside the gate | secondary, open |

Verdicts: H-API **supported**, H-REP **supported**. The Phase II core result is reproduced through the package API.

## Scope

This is a healthy-baseline result at one parameter point, with exact starts and one run length. It is not an RNG-fault result, and it tests none of H1, H2 or H3. Agreement with Phase II is distributional, not bit-for-bit.

## Files and reproduction

[REPRODUCE.md](REPRODUCE.md) gives the steps.
- **Traces** (about 3 GiB, `design.json`) are not in Git; their SHA-256 values are in `results/run_manifest.json`.
- **AI disclosure:**
  - The PIMC run options and their 23 tests were written by OpenAI Codex and ported unchanged after review.
  - The golden test, the start-path provenance, the protocol, the scripts, the analysis and these documents were produced with Claude Code (Anthropic).
  - All of this was done under the author's direction.
