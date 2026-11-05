# RNGLab central hypothesis: how v0.18 supports it, what is missing, staged roadmap

Status: **adopted by the author, 2026-11-05.** Adopted: the staged roadmap (Section 6), the Phase 3 redesign as the H1 campaign (decision 1), and the widened Phase 3.5 review (decision 2). Still open: decisions 3–5 (external test battery, README hypothesis table, v0.19/v0.20 versioning). Drafted by Claude Code (Opus 5.5) from the author's specification of the central hypothesis. It changes no frozen protocol: the Phase 1 campaign (`research/v018-estimators/`) was running when this was written and is unaffected.

## 1. The central hypothesis (author's specification)

RNGLab asks whether RNG degradation can be **detected through application-level behaviour while conventional RNG tests still pass**, whether the **remaining useful lifetime** of the affected application can then be **predicted**, and whether **targeted interventions** extend it by a measurable amount. The contribution is the full sequence:

degradation → conventional tests still acceptable → application anomaly → early detection → lifetime forecast → intervention → reassessment → measured lifetime extension.

- **H1, application-aware detection:** application observables (bias, coverage, variance, convergence) reveal degradation before conventional RNG tests detect it, distinguished from sampler failure, numerical error and application defects.
- **H2, remaining-useful-lifetime prediction:** after detection, predict the remaining budget (bytes, requests, sweeps, time) before a predefined, justified failure threshold, with uncertainty, and verify it against observed crossings.
- **H3, intervention and lifetime extension:** interventions (independent true-random source, reseed from independently assessed entropy, alternative generation strategy, state reset) delay failure. Extension = failure budget with intervention − without, in matched experiments. Report failures and counterproductive cases.

Supporting work (estimator calibration, convergence diagnostics, Ising/PIMC/AR(1) campaigns) establishes that application metrics can be trusted. It does not replace the hypothesis.

## 2. What the existing evidence already says

| Question | Evidence | Source | Reading |
|---|---|---|---|
| H1, strong acceptance-uniform shift (PIMC, refreshed sampler) | The marginal source-mean rule raised an alarm before every one of the 470 observed post-onset failures (median 2,304 sweeps earlier). The application-shape rule raised an alarm for 0 of 470. | `research/phase2-20261030/reports/Controlled_RNG_fault_report.md`, `experiment/fault-study-v1/alarm_failure_times.csv` | **Against H1 for this fault:** the source test won. |
| H1, proposal/acceptance coupling | The marginal source rule missed it; a joint-channel source rule detected it. The application rule crossed in 431 of 484 failures: 96 earlier, 115 at the same endpoint, 220 later; 53 never. | same | **Weak for H1:** single-stream tests are blind to coupling *by construction*, and the application rule mostly alarmed at or after failure. |
| H1, v0.14 "0.0% hazardous gap" | Superseded: one uncalibrated run. | `docs/v018/AUDIT.md` §2.2 | No evidence. |
| H2, forecasting (v0.16, v0.17) | Theory features gave no benefit. Abrupt faults had already failed at issuance. A Brier gain only for healthy vs gradual drift. Decision utility (H4) and unseen families (H5) failed. | `docs/trident/V016_NEGATIVE_FINDINGS.md`, `docs/calibrate/V017_NEGATIVE_FINDINGS.md` (outputs off-repo) | **Unsupported so far.** |
| H2, identifiability | Two degradation processes with identical observed history and different futures. | `research/phase2-20261030/rus/` | Forecasting needs a **declared continuation model**. Unrestricted H2 is provably impossible. |
| H3 | Never measured. v0.17 used only a hypothetical 50% mitigation in a cost sensitivity. | `docs/calibrate/V017_SCIENTIFIC_PROTOCOL.md` (`intervention_benefit_measured: false`) | **No evidence.** |
| Prior work | Application-level Ising tests catching generators that pass standard tests are established (Ferrenberg, Landau & Wong 1992; Coddington 1996). | Phase II next-milestone report | H1 is new only as a *calibrated, latency-matched* comparison, not as an observation. |

Current status: H1 has mixed, fault-dependent evidence that leans against a general claim. H2 and H3 are untested, and H2 has a known impossibility boundary. That should be stated wherever the hypothesis is presented.

## 3. How the v0.18 plan maps onto H1–H3

| v0.18 phase | Role for the central hypothesis | Covers | Gap |
|---|---|---|---|
| 0 audit, namespace, README | Credibility of every number later used as evidence | prerequisite | The README leads with Phase II, as the brief required; it does not yet state H1–H3 or their status |
| 1 estimator repair (frozen, running) | Calibrated uncertainty and convergence diagnostics; separates healthy sampler failure from RNG faults | prerequisite for the H1 "not sampler failure" control | none for its own purpose |
| 2 API validation | Healthy baseline through the package API; byte-accounted sources | prerequisite; budget accounting in bytes exists (`RandomSource.consumed_bytes`) | no mid-run source switch, checkpoint or branching (needed for H3) |
| 3 severity response and monitors (not yet frozen) | The closest current test of H1: where the simulation breaks versus where monitors notice, at matched false alarms | H1, partially | no conventional RNG battery arm; synthetic adapters only; constant-severity faults only; monitors not required to be reference-free |
| 3.5 novelty check | Decides whether the compute is justified | gate | scope must expand to the whole sequence and to prognostics/RUL and intervention literature |
| 4 decision pilot (exploratory) | Offers an application-defined failure threshold (decision loss) | H2 threshold, optional | exploratory only |

## 4. Constraints the roadmap must respect

1. **"Conventional tests pass" needs an operational definition.** Name the battery, the stream it sees (what a practitioner would test: the generator's output as consumed, channels interleaved), the sample size and the matched whole-run false-alarm rate. Some faults are invisible to single-stream tests by construction: coupling uses two individually uniform streams. For those, "the application detects while tests pass" is true by design and carries little information. The informative regime is faults a battery *can* eventually detect, where the question is which signal comes first.
2. **Separate evaluation truth from monitor information.** Exact references (finite-bead energy, Ising enumeration, AR(1) variance) score correctness. Monitors may use only what a practitioner has: healthy calibration runs and invariants such as E[m] = 0 at zero field, agreement of the thermodynamic and virial energy estimators, and acceptance-rate stationarity. Otherwise H1 is not realistic.
3. **A cumulative Monte Carlo estimate under a persistent fault fails predictably, and that gives H2 an identifiable form.** If each post-onset sample carries bias β and the long-run standard deviation is σ, the running mean after onset at n₀ has bias of roughly β(n − n₀)/n, while its interval half-width shrinks as z·σ/√n. Coverage is therefore lost at an n* that follows from (β, σ, n₀), so failure is eventually certain and its timing is computable once β is estimable. Drifting β(t) integrates the same way. Unannounced abrupt onsets remain unpredictable (RUS counterexample), so forecasts can only be issued after detection. This approximation must be derived properly and validated in the H2 protocol before use.
4. **Measure budgets in the units the hypothesis names.** `RandomSource` already counts consumed bytes per source, so megabytes and requests come directly from the package API path. The vectorized kernels count uniforms instead.
5. **In Monte Carlo, an RNG intervention does not undo damage already done.** Switching to a healthy source stops new bias, but contaminated samples stay in the running average, and the chain state may be off target (for example an absorbed Ising state). The lifetime gained therefore depends on whether accumulators are reset or the contaminated window discarded, and whether the state re-equilibrates. Discarding data and re-equilibrating cost budget and can be counterproductive. Reseeding from the compromised source is the required negative control.
6. **Failure thresholds need a declared justification.** A coverage-based threshold (the nominal interval stops covering at the calibrated rate) is justified statistically without a user. Tolerance-based thresholds need a user or an application loss (V018_PRIORITIES item 2; the Phase 4 decision loss).
7. **Attribution to the RNG** needs paired common-random-number designs, a source-swap counterfactual, and the Phase 1/2 healthy calibration. Phase II showed that a healthy sampler can fail on its own (811 / 1,024 coverage).

## 5. Missing experimental capabilities

| ID | Capability | Needed for | Exists today |
|---|---|---|---|
| C1 | Conventional RNG battery on the consumed stream, run in streaming windows and calibrated to the same whole-run false-alarm rate as the application monitors | H1 | Partly: SP 800-22 monobit, block frequency and runs only; SP 800-90B RCT/APT. No external battery (e.g. TestU01, PractRand) and no streaming calibration |
| C2 | Package-level fault adapters with severity **and time profile** (constant, step, linear drift, random-walk ageing, intermittent), byte-accounted | H1, H2 | Constant-severity adapters only in frozen Phase II research code; Phase 3 starts this |
| C3 | Real generators with tunable known weaknesses beside the synthetic adapters (candidates to verify in Phase 3.5: small-modulus or low-bit LCGs, short-lag shift-register/lagged-Fibonacci types, reduced-round ciphers) | H1 external validity | `ChaCha20Source(rounds=…)` only; untested as a weak generator |
| C4 | Reference-free application monitors (invariants, estimator-pair consistency, acceptance-rate and autocorrelation drift, the Phase II shape monitor) | H1 | Phase II shape monitor (research code); no package monitors |
| C5 | Online alarm times in sweeps and bytes; sequential rules (e.g. CUSUM) with calibrated thresholds | H1 latency | Whole-run maxima only |
| C6 | Pre-declared failure thresholds and online failure-time measurement | H1, H2, H3 | Phase II toy rule (rolling-window error > 0.05) |
| C7 | Lifetime forecaster under declared degradation models, with prediction intervals; evaluation by interval coverage, calibration and error against observed crossings; baselines (no change, the analytic n*, empirical hazard); a held-out and an unseen degradation family | H2 | No; v0.16/v0.17 forecasters target a different, confounded endpoint |
| C8 | Intervention operators in the API: mid-run source switch, reseed from an independent source, reset or discard of accumulators, state reset with re-equilibration; each with a cost in samples and bytes | H3 | No |
| C9 | Checkpoint and branch: copy the chain and source state at the intervention time and run matched with/without branches from an identical prefix | H3 matched design | No; `RandomSource` hashes consumption but has no state copy |
| C10 | A second application for generality (Ising with Metropolis and Wolff; the Wolff kernel has no unit test) | H1–H3 | Kernels exist; `wolff_vectorized` untested |

Out of scope until a separate SP 800-90B assessment exists: physical sources (the detector data carries no entropy claim).

## 6. Staged roadmap

Each stage is its own versioned, frozen protocol with a fresh seed namespace. Each stage's negative outcome is a reportable result.

| Stage | Content | Hypothesis | Main deliverable | Depends on |
|---|---|---|---|---|
| **A (now)** | Finish Phase 1 (running); Phase 2 API validation exactly as specified | foundation | calibrated estimators; validated healthy PIMC baseline through the API | — |
| **B** | Phase 3.5 novelty review, widened to the full sequence: application-based RNG testing, online RNG health tests, MCMC diagnostics, prognostics/RUL estimation, restart/intervention in Monte Carlo | gate | `RELATED_WORK.md`; a stop/go on novelty **before** compute | — |
| **C** | Phase 3 redesigned as the **H1 campaign**: severity curves (acceptance shift, partial coupling, lag-1) plus C1 battery arm, C4 reference-free monitors, C5 latencies in sweeps and bytes, one C3 real weak generator, C6 thresholds. Primary endpoint per fault family: the severity region where the application is wrong while the battery passes at matched false alarms, and the latency difference with confidence intervals | H1 | severity × detection maps; "pass-but-broken" regions; RNG-tests-first regions reported equally | A, B |
| **D** | **H2 campaign**: C2 drift and ageing profiles; forecasts issued at the H1 alarm time; checks against observed crossings on held-out trajectories and one unseen profile family; the analytic n* as benchmark; explicit non-identifiable cases | H2 | forecast calibration and error, in bytes and sweeps | C |
| **E** | **H3 campaign**: C8 and C9 branched matched design (identical prefix, intervention at the alarm, common random numbers after the branch where possible); interventions: independent source switch, independent reseed, accumulator reset, state reset; negative control: reseed from the compromised source | H3 | lifetime extension with confidence intervals; failure and counterproductive conditions | C, D |
| **F** | Integration: the full sequence on two models (PIMC, Ising), independent replication, archived data, paper | all | submission package | E |

**Focus recommendation.** The brief asks for one focused, defensible contribution. For the **TDK**, stages A–C give a complete story: a validated baseline, plus a calibrated test of H1 whose answer may be "it depends on the fault family", which is still a result. H2 and H3 then appear as the declared research programme, with the identifiability argument (constraint 3) as a concrete plan. For a **paper**, H1 plus the scoped, declared-model H2 (stage D) is the smallest contribution likely to be new; H3 strengthens it but is the most engineering-heavy. All stages fit CPU-only local compute. The costly part is C8/C9 engineering, not simulation.

## 7. Decisions this asks of the author

1. **Phase 3 design:** add C1, C4, C5, C6 and one C3 generator, and make "application wrong while the battery passes, at matched false alarms" the primary endpoint. Phase 3 is not frozen, so this violates no protocol.
2. **Phase 3.5 scope:** widen it as in stage B. Its stop/go then applies to the whole sequence, not only to severity maps.
3. **Conventional battery:** keep the in-repo SP 800-22 subset plus 90B health tests, or add an external battery (TestU01 or PractRand compiled locally, with the version pinned). The external battery is stronger against a referee but adds a dependency.
4. **README and later deliverables:** add the central hypothesis with a status table (H1 mixed, H2 and H3 untested) above the Phase II result. `TDK_OUTLINE.md` and `PUBLICATION_PLAN.md` should be framed around H1–H3 with the stage boundaries above.
5. **Versioning:** stages D and E could become v0.19 and v0.20 rather than being added to v0.18.

## 8. Guardrails

- Phase 1 and Phase 2 results are not to be relabelled as H1–H3 evidence; they establish trust in the instruments.
- The v0.16/v0.17 forecasting results stay negative. Stage D is a new design with a declared model, not a reanalysis.
- Synthetic adapter faults are labelled as such; any statement about real generators rests on C3 only.
- "Passes conventional tests" always names the battery, stream, sample size and false-alarm rate. It is never presented as evidence that a source is random or secure.
