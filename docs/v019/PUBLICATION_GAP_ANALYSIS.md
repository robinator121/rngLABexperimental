# Publication gap analysis: what is established, what is missing, in what order

**Status: 2026-11-11, for the author's planning.** Companion to `SCIENTIFIC_ASSESSMENT.md`. "Established" means: produced by a committed script from a frozen protocol, with an independent recount where one exists.

## 1. Claims and their evidence today

| # | Claim | Evidence in the repository | Status | Gap before a TDK paper | Gap before a journal |
|---|---|---|---|---|---|
| 1 | The refreshed sampler's fixed-batch intervals cover at the nominal rate | Phase 2: 3846/4096; Stage C V1: 3867/4096 | established at batches of 480 and 512 | criterion 1 of the pilot extends it to 640, 1,280 and 2,560 | a second parameter point |
| 2 | An AR(1)-calibrated batch-length rule fails on that sampler | Phase 2: 3616/4096 = 0.8828 | established as an observation | fit of the two-scale reading; baselines of `SCIENTIFIC_ASSESSMENT.md` §1 | a second system; the failure region as a proposition; prospective confirmation |
| 3 | H1-C: application monitors do not detect before a calibrated battery | Stage C, 24 cells, recount 319/319 | established **negative** for four synthetic families | none; state the dependence on the targeted tests (Sensitivity 01) | an external battery; a real weak generator with documented application failures |
| 4 | Forecasts need a declared continuation; the loss otherwise is bounded below | Phase II counterexample; bound and tests in v0.19 | established (elementary) | none | on its own it is not a paper result |
| 5 | A declared-model forecaster is calibrated and skilful under its family, and not outside it | surrogate study, R1 to R7 | established **on a surrogate only** | say "surrogate" wherever it is quoted | the same in the real application |
| 6 | The validity of a batch-means interval under drift follows the doubly noncentral form | quadrature against simulation in the Gaussian model | established for the model | check batch-mean normality and independence in the sampler | that check, published |
| 7 | H2-R: validity can be forecast after the alarm in the application | none | **not tested** | the pilot's criteria | a campaign with held-out and unseen conditions |
| 8 | H3: an intervention extends the useful run | none | **not tested** | none planned for the TDK | the whole of Stage E |

## 2. Before the first pilot draw

From `IMPLEMENTATION_AUDIT.md`. The pilot must not start until the author says so.

| Step | Owner | State |
|---|---|---|
| Commit the measurement layer | pilot | uncommitted in the main checkout |
| Campaign driver, calibration, fitting with the evaluation barrier, scores, criteria | pilot | not written |
| Decide F-01: accept that the oracle is unavailable where the frozen call returns NaN, or amend before analysis | **author** | open |
| Decide where the run records live (19.4 GB) and whether writes are atomic (F-05, F-14) | **author** | open |
| Record the TDK deadline in the pre-generation manifest | **author** | missing: criterion 5 is otherwise INSUFFICIENT_INFORMATION |
| Smoke run in the smoke namespace, then `verify.py` on it | pilot, audit | not run |
| Pre-generation manifest with code hashes and the libc version (F-08) | pilot | missing |
| `preflight.py` shows no FAIL | audit | currently NOT READY: code, manifest |
| Plan about nine hours of computing and 23 GiB (F-05, F-06) | author | measured once, at six workers |

## 3. Order of work proposed

1. **Close Section 2 and run the frozen pilot.** It is already designed, it costs hours, and its stop rule decides whether Stage D continues. Its independent recount exists (`research/v019-audit/verify.py`).
2. **In parallel, a pre-registered Direction A study on stored healthy traces**, then a fresh confirmation set. No new protocol exists yet; it needs one before any analysis that counts.
3. **Write the TDK around what is established:** claims 1 to 4 as results, claim 5 as methods on a surrogate, claims 7 and 8 as the declared programme with the pilot's outcome, whichever it is.
4. **Only if the pilot passes:** the campaign protocol of `FORECASTING_METHODS.md` §6.

## 4. Inventory of unverified or weakly verified statements

| Statement | Where | Why it is not verified |
|---|---|---|
| The slow "shape mode" explains the 0.8828 coverage | Phase 2 report, `SCIENTIFIC_ASSESSMENT.md` §1 | a reading; no model was fitted |
| F-B-lite's SHIFT features are near noise at issuance | audit F-03 | one full-length run |
| The pilot takes about nine hours | audit F-06 | one benchmark at six workers; eight not measured |
| Bit-for-bit replay fails on another platform | audit F-08 | inferred from the libm dependence; no second platform was tried |
| The surrogate's scales resemble LAG1 in the sampler | `FORECASTING_METHODS.md` §3 | chosen from four probe runs whose outputs are not in the repository |
| No prior work forecasts Monte Carlo validity under a degrading generator | `SCIENTIFIC_ASSESSMENT.md` §2 | web searches only; no database or citation-graph search |
| What the cited papers contain | both assessments | abstracts and secondary descriptions; full texts not read |
| The second measurement-layer snapshot is correct | audit §6 | compared through outputs, not re-read |
| Batch means in the sampler are normal and independent enough for the validity model | `FORECASTING_METHODS.md` §7 | not tested |
| Anything about real generators | everywhere | all faults are synthetic adapters |

## 5. What must not be claimed

- That application monitors detect what conventional tests miss: Stage C found the opposite for these faults.
- That remaining useful life can be forecast: only a conditional, surrogate-level statement exists.
- That the batch-length failure is a new phenomenon: slow modes that defeat window rules are documented ([SSV11], [Wol04]).
- That any result is certified, proven in general, or established as novel.
