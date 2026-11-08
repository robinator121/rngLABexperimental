# v0.18 results summary: positive and negative findings side by side

Running summary, updated after each phase. Every entry names its source. Produced by Claude Code (Opus 5.5); the decisions are the author's. Last update: Stage C, 2026-11-08.

## Phase 0: audit and consolidation

| Positive | Negative or limiting |
|---|---|
| Phase II coverage and fault counts, the test counts, crypto BLUE results and the detector inventory are backed by repository scripts and outputs; the detector figures were re-derived from the raw export (`AUDIT.md` §1) | The `SOURCE_COMPARISON.md` benchmark matrix is unbacked: throughputs 29–156× above what the code achieves for 4 of 5 sources, and an "all 15" SP 800-22 claim for a 3-test implementation (`AUDIT.md` §2.1) |
| The v0.17 healthy susceptibility counterexample matches its saved output exactly, and its deterministic replay is now a test | The previous README VMC numbers do not reproduce (re-run: 1.656 / −0.248 / 0.296) |
| One v0.18 namespace re-exports the audited components unchanged; 9 identity and anchor tests | v0.16/v0.17 outputs exist only off-repo; GitHub CI has never executed in this repository (billing lock) |

## Phase 1: estimator repair (`research/v018-estimators/REPORT.md`)

| Positive | Negative or limiting |
|---|---|
| Symmetry-aware susceptibility S is calibrated in all 9 Ising cells (coverage 0.925–0.949) | The generic 16-chain jackknife G fails its interval gates in 5 of 9 cells (coverage 0.836–0.916 at T = 1.5); no generic repair is established, so **H-S is not supported** |
| The plug-in failure (−0.917 relative bias at L4 T1.5) and the short-batch failures (BM(32): 0.547, 0.332) reproduce under fresh seeds | S relies on an exact symmetry most models lack; S coverage at T = 1.5 sits near the gate's lower edge (0.925–0.930) for an unexplained reason |
| BM, OBM, SV-B, SV-TH and Geyer are calibrated for τ ≤ 19 (15 of 15 gates); **H-L supported** | Spectral estimators run 6–8% low at b/n = 1/16; at long τ only Geyer stays near nominal, in an AR(1) setting that favours it |
| The batch-length diagnostic never accepts a run and then mis-covers (21 of 21 evaluable gates); split-R̂ flags the Ising sign non-mixing; **H-D supported** | The diagnostic is conservative: it abstains at τ ≈ 199, n = 8192, where BM(512) would cover 0.926 |
| The independent recount agrees on 28 of 28 counts; Amendment 01 changed R̂ by at most 4.9 × 10⁻⁵ | The frozen control verdict is FALSE (1 of 20 intervals missed), exposing a protocol-design flaw that Amendment 01 documents |

## Phase 2: the Phase II core result through the package API (`research/v018-api-validation/REPORT.md`)

| Positive | Negative or limiting |
|---|---|
| Through the package API, the refreshed sampler (3846/4096, 0.9390) and the IID control (3884/4096, 0.9482) pass the Phase II healthy gates; **H-API supported** | The batch-length diagnostic accepts every refreshed run but its intervals cover only 0.8828: the AR(1)-derived 10·τ̂ rule misses the slow shape mode |
| Local Metropolis under-covers (3323/4096, 0.8113) and the refresh improves paired coverage by 0.1277 [0.1138, 0.1414]; **H-REP supported**; all three conditions agree with Phase II | The package's default batch sizes under-cover even with the refresh (0.8611 at b = 32, 0.9260 at b = 128) |
| Byte ledger, start pairing, trace consistency and replay are exact (I-1 to I-4); the independent recount agrees 21/21 | Refresh-condition bias intervals exclude zero (CV −0.000321, TD +0.001713), far inside the gate; unexplained (post-hoc check in `posthoc_bias_check.json`) |
| The PIMC run options were ported from Codex unchanged; a golden test captured before the port shows the defaults are bit-identical | A healthy baseline at one parameter point only; tests none of H1–H3 |

## Stage C: application-aware detection at matched false alarms, H1 (`research/v018-h1-detection/REPORT.md`)

| Positive | Negative or limiting |
|---|---|
| All gates pass: I-1 to I-4; V1 3867/4096 (0.944092); V2 false alarms C 0.048584, A 0.0424805. The independent recount agrees 319/319 | **H1-C NOT SUPPORTED in all four families** (SHIFT, LAG1, DUP, GFSR): no H1 cell and no undetected-harm cell among 24 |
| Amendment 01 was approved before analysis. The superseded analysis gives identical output except the 24 added descriptive intervals | **The battery alarmed in every run of all 18 synthetic cells**, including 7 "conventional first" cells where the application's wrong-rate stayed below the class bound of 0.2 (0.027 to 0.113; healthy 0.056) |
| The application-arm alarm rate rises with severity together with harm, while the battery is saturated (descriptive; a hypothesis for Stage D) | The application arm never alarmed alone. When both alarmed it was later by a median of 512 to 5,376 sweeps |
| GFSR, all six lag pairs including R250: harmless and quiet in this application (wrong-rate 0.027 to 0.066; healthy 0.056) | LAG1 and DUP were caught by the two word-level tests aimed at them. Whether SP 800-22 alone would catch them first was not measured |
| | The SP 800-90B health tests alone stayed near their healthy rate for LAG1, and for DUP up to severity 0.19 |

## Stage C sensitivity 01: the battery without the targeted word tests (`research/v018-h1-detection/SENSITIVITY_01_RESULTS.md`)

Post hoc and exploratory; plan frozen at `c7231fa` before the analysis; stored Stage C data only.

| Positive | Negative or limiting |
|---|---|
| With only the seven in-repository SP 800-22 tests (C′, recalibrated; healthy 182/4096 whole-run alarms), no cell is H1 or undetected harm; the Stage C verdicts stand in all four families | **For LAG1, "the battery detects first" depended on the lag-1 test aimed at it.** C′ raised 10 to 41 post-onset alarms per 256 runs, and up to 186 of 256 runs were wrong with no C′ alarm (none with the full battery) |
| DUP: the SP 800-22 subset still detected before harm (DUP-3: 164 alarms, 17/256 wrong) | LAG1 came close to the H1 class: the application arm alarmed in up to 256 of 256 runs, the A − C′ interval excluded zero in all six cells, and LAG1-4 missed H1 only by its wrong-rate lower bound (0.496 against 0.5). Descriptive, not a confirmed H1 result |
| The self-check reproduced `results/gates.json` exactly | Post hoc, on the same data, with the Stage C classes; C′ is a 7-of-15 subset of SP 800-22 |

## Central hypothesis (H1–H3, `ROADMAP.md`)

Evidence status:
- **H1, as narrowed** (calibrated, matched false alarms, same stream): not supported in Stage C for SHIFT, LAG1, DUP and GFSR with the nine-test battery. Phase II's family-dependent evidence included channel-specific faults that a single-stream battery cannot see by construction; those were excluded from Stage C (PROTOCOL §9) and remain untested under matched false alarms. A post-hoc sensitivity analysis shows that for LAG1 the negative result relied on the targeted lag-1 test: against the SP 800-22 subset alone, LAG1 comes close to the H1 pattern (exploratory).
- **H2 and H3** are untested.
