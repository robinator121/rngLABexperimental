# v0.18 Stage C: application-aware detection at matched false alarms (H1)

Completed 2026-11-08. Pre-registered protocol `rnglab-v018-h1-detection-v1` ([PROTOCOL.md](PROTOCOL.md), frozen at `015bed6`). The pre-generation manifest was pushed at `0aca45a` before the first draw. One pre-analysis amendment, approved by the author: [AMENDMENT_01.md](AMENDMENT_01.md). Full report: **[REPORT.md](REPORT.md)**.

## What was tested

**H1-C, per fault family:** at a whole-run false-alarm rate of 0.05 for both arms, is there a severity at which:
- the application's uncertainty interval is wrong in most runs;
- a conventional battery, run online on the same consumed byte stream, stays quiet;
- reference-free application monitors alarm more often than the battery?

The application was harmonic path-integral Monte Carlo (M = 16, β = 2, δ = 0.4) through the package API. The conventional arm was nine tests: seven from SP 800-22 plus lag-1 word correlation and adjacent identical words. The application arm combined two monitors:
- the TD − CV estimator invariant;
- the uphill-proposal fraction.

There were four fault families: SHIFT, LAG1 and DUP (synthetic, word-level) and GFSR (six real lag pairs, R250 included). Each family had 6 severities with 256 runs each, against 999 calibration runs and 4,096 healthy runs.

## Results

| Finding | Evidence |
|---|---|
| **H1-C is NOT SUPPORTED in any family**: no H1 cell and no undetected-harm cell among the 24 | `results/gates.json` |
| **The battery alarmed after onset in every run of all 18 synthetic cells**, including 7 "conventional first" cells where the application's wrong-rate stayed below the class bound of 0.2 (0.027 to 0.113; healthy 0.056) | T2, `results/report_numbers.json` |
| The application arm never alarmed alone in a synthetic cell. When both alarmed, it was later by a median of 512 to 5,376 sweeps | T4 |
| LAG1 and DUP were caught first by the two word-level tests aimed at them. SHIFT was caught first by SP 800-22 tests | secondary |
| The SP 800-90B health tests alone stayed at about their healthy rate for LAG1, and for DUP up to severity 0.19 | secondary |
| **GFSR, all six lag pairs: harmless and quiet** (wrong-rate 0.027 to 0.066; healthy 0.056) | T2 |
| V1 3867/4096 (0.944092); V2 false alarms C 0.048584, A 0.0424805; I-1 to I-4 pass | gates |
| The independent recount agrees 319/319. The superseded pre-amendment analysis gives identical output except the 24 added descriptive intervals | `results/recount.json`, `results/amendment01_comparison.json` |

Verdicts, all interpretable: H1-C **not supported** for SHIFT, LAG1, DUP and GFSR. This negative result is the main finding of Stage C.

## Scope

This covers one application at one parameter point, three synthetic fault families and one real generator family. The battery includes tests that target LAG1 and DUP directly. GFSR showing no harm here says nothing about cluster algorithms. Nothing here tests H2 or H3. See REPORT §8.

## Files and reproduction

[REPRODUCE.md](REPRODUCE.md) gives the steps.
- **Traces** (about 2.6 GiB, `design.json`) are not in Git. Their SHA-256 values are in `results/run_manifest.json`.
- **AI disclosure:**
  - Produced with Claude Code (Anthropic) under the author's direction: protocol, library, campaign, analysis, recount, report.
  - Written by OpenAI Codex (`gpt-6-astra`) and reviewed by Claude Code:
    - `report_tables.py`, `plot.py` and `report_numbers.py`;
    - the seven Amendment 01 fixes, from Codex's own review.
  - Contributions are listed one per line in `docs/v018/AI_LOG.md`.
