# v0.18 results summary: positive and negative findings side by side

Running summary, updated after each phase. Every entry names its source. Produced by Claude Code (Opus 5.5); the decisions are the author's. Last update: Phase 1, 2026-11-05.

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

## Central hypothesis (H1–H3, `ROADMAP.md`)

No v0.18 phase so far tests H1, H2 or H3. Evidence status: H1 is mixed and depends on the fault family (Phase II); H2 and H3 are untested.
