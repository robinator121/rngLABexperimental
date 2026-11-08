# Stages D and E revisited after Stage C

**Status: DECIDED by the author, 2026-11-08.** The first option of every decision in Section 5 was adopted. Order: with about six months to the TDK submission, the Stage D feasibility pilot (v0.19) runs alongside the TDK write-up of stages A to C, so that its result can enter the TDK. The draft was revised after an adversarial review whose findings were checked against the Stage C files. No data were generated for this document, and it changes no frozen protocol. It revisits the "conditional GO" for Stages D (H2) and E (H3) in `RELATED_WORK.md` §1, as `ROADMAP.md` asks after Stage C. Every number is copied from the Stage C file named beside it.

## 0. Recommendation in brief

- **Stage D: GO for a bounded feasibility pilot of a reframed H2 (H2-R, Section 2), not yet for a campaign.**
  - The forecast is issued at the conventional battery's operational alarm, not at an application alarm.
  - The forecast target: the probability that the run's result is still valid at each of a few declared run lengths. This is the run's survival curve at fixed horizons.
  - A battery-based forecaster, and the same forecaster with the application monitors added, are compared against baselines that ignore the run's fault signals.
  - The pilot has stop criteria (2.7). A D protocol is written only if it passes them.
- **Stage E: stays conditional; build nothing yet (Section 3).** Switching at the alarm is the obvious comparator, but its benefit is untested, and E should be gated on a cost trade-off and on D's pilot showing usable lead time.
- **Alternative: stop after C.** Write up stages A to C for TDK, with D and E as the declared programme. `ROADMAP.md` §6 already calls A to C a complete TDK story. This is the right choice if the TDK deadline comes before a D pilot and campaign could finish.

## 1. What Stage C changes for D and E

**Findings** (`research/v018-h1-detection/REPORT.md`; tables T1 to T5 in `results/report_tables.md`; numbers also in `results/report_numbers.json`):

1. **The battery alarmed in every synthetic run.** In SHIFT, LAG1 and DUP it alarmed after onset in every run of all 18 cells. That includes the 7 "conventional first" cells, where the wrong-rate stayed below the class bound of 0.2 (T2; REPORT §1).
   - This depends on the battery. A post-hoc sensitivity analysis removed the two word-level tests, leaving the SP 800-22 subset. That subset raised only 10 to 41 post-onset alarms per 256 LAG1 runs, while up to 186 of 256 runs were wrong (`research/v018-h1-detection/SENSITIVITY_01_RESULTS.md`).
2. **In severity, detection set in below harm.** In the pilot, the battery alarmed in at least 31 of 32 runs at lower severities than those where the wrong-rate first reached 0.5, in all three families (REPORT §6.4, `pilot/summary.json`).
   - This compares constant-severity cells. Stage C did not measure when harm sets in within a run: harm was assessed over a fixed post-onset window (`PROTOCOL.md` §2.4).
3. **The first-firing test is associated with the family.** In the campaign, lag-1 correlation fired first in all LAG1 runs, adjacent identical words in all DUP runs, and only SP 800-22 tests in SHIFT runs (T5; REPORT §2).
   - This is descriptive, not a validated classifier. The recorded name is the lowest-p test at the first alarming checkpoint (`analyze.py`), and healthy and GFSR alarms also produce names.
4. **The application arm's alarm rate rose with severity, along with the wrong-rate** (REPORT §6.1, descriptive). Not every application signal tracks harm, though:
   - the shape monitor and the CUSUM largely missed SHIFT (REPORT §6.2);
   - in LAG1, A1 stayed near its healthy level, even where the CV bias reached 0.0294 and 224 of 256 runs were wrong (LAG1-5; T2, T5).
5. **GFSR showed no detectable excess harm** in this setting: wrong-rates 0.027 to 0.066, against 0.056 for healthy runs (T2). Its individual runs were wrong at about the healthy rate.

**Consequences:**

- **The roadmap's sequence does not occur for these faults.** The roadmap's sequence is "conventional tests still pass → application anomaly → early detection → forecast". Here the battery itself alarmed in every run. So a D forecast should be issued at the **battery alarm**, not at an application alarm as `ROADMAP.md` §6 assumed.
- **A binary alarm says *that* the stream is faulty, not *whether* the result will be valid.** The battery alarmed in every run of DUP-2 (wrong-rate 0.027) and of DUP-5 (wrong-rate 1) (T2). That gap is what D can fill.
  - The practical question after an alarm: when restarting is costly (a long run, an embedded source), will the current run's result still be valid at its planned end?
- **Battery statistics respond to the fault parameters, but are not yet measurements of them.**
  - The lag-1 test measures the correlation of transformed uniforms, while the LAG1 parameter is a latent Gaussian correlation (`h1lib.py::_transform`).
  - The duplicate count needs normalising by the number of adjacent pairs (`h1lib.py::word_tests`).
  - SHIFT shows only in the frequency-type statistics.
  - Drift, and windows that straddle the onset, change all of these.

  A battery-based forecaster therefore needs a calibrated observation model built from the raw, signed statistics and their sample sizes. Saturated or floored p-values are not severity measurements. With that model it is a strong reference-free competitor, and D should include it.
- **Stage C does not show that the application histories carry no extra information.** It shows no incremental *detection* benefit under Stage C's endpoint, where the battery's binary alarm was saturated. Whether the application signal adds *forecasting* information is open; it is D's continuation of the H1 question.
- **D's faults are synthetic.** GFSR is left out of the degradation campaign on the limited basis of finding 5; healthy runs remain the control. Every D result is about synthetic adapters and declared degradation models (`ROADMAP.md` §8).

## 2. Stage D reframed: H2-R

**H2-R.** After an operational alarm of the conventional battery on a degrading source, can one forecast whether a running Monte Carlo estimate will still be valid at declared run lengths? Valid means that the estimate's uncertainty interval covers the exact value. The forecasts must be calibrated under a declared continuation model and beat baselines that use the same non-fault information.

**Sub-question, continuing H1:** do the application monitors add forecasting skill to the battery (F-AB against F-B, Section 2.4)?

### 2.1 Endpoint (evaluation truth only)

- **Recommended: validity at declared run lengths.**
  - For each run and each declared length L₁ < … < L_k up to the horizon, record whether the interval of the estimate at L covers E_M.
  - A forecast gives a probability for each L, so it is the survival curve at fixed horizons.
  - Every length is observed in every run, so nothing is censored inside the horizon. A proper score for binary outcomes (Brier or logarithmic [GR07]) applies directly.
- **Not recommended: last exit** (the first checkpoint after which the interval excludes E_M at every checkpoint up to the horizon).
  - It is a retrospective, horizon-dependent label. A recovery beyond the horizon would undo the event, and an exclusion at the final checkpoint alone counts as a failure.
  - If a time-to-event endpoint is wanted as a secondary, it must be a prospectively ascertainable sustained-exclusion event, administratively censored at the horizon and reported with a horizon-sensitivity check. The Brier score has documented problems under administrative censoring [KB23].
- **The measurement procedure is new and must be validated.**
  - Stage C validated one fixed post-onset window: sweeps 2,560 to 10,239, batch means of 16 batches of 480 (`PROTOCOL.md` §2.4).
  - For each L, the D protocol must fix the measurement start, the batching, incomplete batches, and the handling of statuses other than `ESTIMATED`.
  - Healthy coverage must be validated at each L. Interval behaviour under drift must be checked in the pilot: differences between batches under drift widen the interval, so exclusion need not grow monotonically with severity.

### 2.2 Degradation profiles and the declared continuation model

- **The generative model is declared in full:**
  - onset distribution, initial severity, ramp-rate distribution, parameter bounds and saturation;
  - the clock that drives degradation: per sweep, or per consumed byte. If per byte, the forecaster must also forecast future byte consumption.

  The protocol states which of these the forecaster is given and which it must infer.
- **Primary profile:** a linear ramp in the fault parameter, with a per-run rate drawn from a declared range.
  - The range keeps slow rates under which runs stay valid to the horizon. Choosing only rates that fail inside the horizon would remove the cases that make the forecast hard.
- **Continuation model:** the declared generative model itself, with its parameter uncertainty propagated into the forecast.
  - At the first alarm, the post-onset history may be too short to identify a slope. A broad forecast is then the correct one. Abstention is allowed and is reported.
- **Held-out conditions, reported separately:**
  - ramp rates or onsets outside the training range (extrapolation, not interpolation);
  - an unseen family, a step followed by a hold or a random walk. Under it the declared model is wrong, so it tests robustness, not calibration.
- **Identifiability:** the RUS counterexample (`research/phase2-20261030/rus/`) shows that without a declared continuation, the same history allows different futures. Forecasts are claimed only under the declared model.

### 2.3 Issuance, inputs and training

- **Issuance at an operational alarm.** Stage C counted only alarms after the known onset (`PROTOCOL.md` §2.7), but a deployed monitor cannot discard earlier false alarms. D therefore needs:
  - a declared battery that issues the alarm. With the SP 800-22 subset alone, most LAG1 runs would trigger no forecast (finding 1);
  - a rule that uses no onset knowledge: the first alarm, with a declared re-arming policy;
  - thresholds recalibrated over D's longer horizon. Stage C's healthy whole-run alarm rate for the battery, 0.048584 (V2, `results/gates.json`), holds only for Stage C's monitoring schedule;
  - reports of pre-onset alarms, runs that never alarm, and runs already invalid at issuance.
- **Deployment inputs, with no new reference query:**
  - the battery's windowed statistics, raw and signed, with sample sizes;
  - the application monitors A1 and A2;
  - elapsed sweeps and bytes.

  A2 is computed from the byte ledger (`PROTOCOL.md` §2.6). So any dynamics derived from the ledger count as application information, and the battery side gets only the elapsed counts.
- **Training with the exact reference** makes the claim a supervised laboratory characterisation of this one application. It needs the author's explicit decision (Section 5, decision 5). If allowed:
  - every learned transformation is frozen before the test data are opened;
  - whole trajectories go to one split, together with the paired runs that share their pre-onset prefix (the Stage C pairing, `PROTOCOL.md` §3.2);
  - pilot, fitting, calibration and final evaluation use separate sets with fresh seed namespaces.

### 2.4 Forecasters

| ID | Information | Role |
|---|---|---|
| F0 | the training distribution of validity at each L given an alarm, per diagnosed family | baseline that ignores the run's signals |
| F0+ | F0 plus the elapsed sweeps at the alarm and the remaining length to each L | baseline with the same non-signal information as the others |
| F-B | battery statistics, through the calibrated observation model, forward-simulated under the declared model | reference-free, battery side |
| F-AB | F-B plus A1, A2 and any ledger dynamics | combined |
| F-A | the application monitors alone | exploratory only |
| Oracle | a predictive distribution given the true severity path, which keeps the future Monte Carlo randomness | benchmark; approximate, not a proven bound |

- **What F-B and F-AB must model.** A mean-bias crossing is not enough. They must model the cumulative estimator: dilution by healthy samples, chain relaxation, autocorrelation and the interval width (constraint 3 of `ROADMAP.md` §4). Stage C's fixed-window mean biases do not establish this map.
- **Contrasts:**
  - primary: F-B against F0+, asking whether there is forecasting skill at all;
  - key secondary: F-AB against F-B, asking whether the application monitors add information. This is the H1-continuation question; F-A against F-B does not test it.
  - F-B and F-AB get equal fitting and tuning budgets.
- **The multiplicity rule must be fixed in the D protocol.** It must name the primary forecaster, score, issuance time, aggregation and smallest meaningful improvement, and a testing hierarchy. Everything outside the hierarchy is exploratory.

**The main technical risk is quantification.** Reference-free signals may detect a fault without measuring the bias it causes (finding 4). F-AB could then fail to beat F-B. That would be a reportable negative, not a design failure.

### 2.5 Evaluation

- **Unit and timing:** one trajectory. Forecasts are issued at the alarm (primary) and at declared later checkpoints (secondary). Comparisons are aggregated per trajectory, so long-lived runs are not over-weighted.
- **Metrics:**
  - a proper score for the validity probabilities (Brier or logarithmic [GR07]);
  - calibration (reliability) and sharpness [GBR07];
  - the prognostics metrics of [SCSSG10] where they apply. Their exact names and definitions must be checked against the full text before freezing; only its abstract was checked in Stage B.
- **Also reported:** the fractions of runs with an alarm, with no alarm, invalid at issuance, and abstained.
- **Negative outcomes, reported with equal weight:**
  - no forecaster beats F0+;
  - the forecasts are miscalibrated;
  - F-AB is no better than F-B;
  - the skill vanishes on the held-out conditions.

### 2.6 Guardrails carried over

- **Not a rescue of v0.16/v0.17:** those forecasting results stay negative (`ROADMAP.md` §8). D is a new design with a declared model.
- **Synthetic only:** nothing about real generators.
- **Narrow scope:** one application (PIMC at one parameter point).

### 2.7 Feasibility pilot and its stop criteria

A bounded pilot comes before any D protocol, in its own seed namespace. Its thresholds are fixed before it runs. It asks:

1. **Healthy validity:** does the interval at each declared L cover E_M at close to the nominal rate in healthy runs?
2. **Something to forecast:** at the operational alarm, is validity at the planned end uncertain for a non-negligible share of runs? If it is almost always already lost, or almost always kept, the forecast is trivial.
3. **Lead time:** does the alarm come early enough before invalidity for an action to matter?
4. **Signal:** does any reference-free statistic move with the realised validity beyond what F0+ predicts (descriptive)?
5. **Cost:** do the engineering and the campaign compute fit the schedule?

**Stop rule:** stop D if criterion 1 or 2 fails. If criteria 3 to 5 fail, reconsider the scope.

**Engineering the pilot needs:**
- time-profile fault adapters (C2) on top of the Stage C adapters, keeping the byte ledger and the sweep clock;
- endpoint measurement at the declared lengths (C6);
- forecaster code (C7).

**Compute, for scale:** Stage C stored 11,239 runs in 110.4 wall-clock minutes on 8 workers (`results/report_numbers.json`). Each run had 10,240 sweeps (`PROTOCOL.md` §2.1; sweep-clock check I-1 in `results/gates.json`). D runs are longer; the sizes come from a design script after the pilot.

## 3. Stage E revisited

- **Switching at the alarm is the obvious comparator, and its benefit is not established.**
  - The battery alarmed in every synthetic run (finding 1), so switching at the alarm would act early.
  - But Stage C did not measure when harm sets in within a run (finding 2).
  - Switching also leaves the contaminated samples in the accumulator and may leave the chain off target (constraint 5 of `ROADMAP.md` §4).

  Prevention, and recovery after switching, are untested.
- **The H3 questions that remain informative:**
  - **Contaminated accumulators:** after a switch, keep the running estimate, reset it, or discard a window. Discarding costs budget and can be counterproductive.
  - **Policies at matched cost:** compare never switching, switching at the alarm, switching on a simple battery-severity rule, and switching when a forecast predicts invalidity at the planned end.
    - The branches start from the same state.
    - The costs include discarded work, re-equilibration and the replacement source.
    - The costs also include acting on healthy alarms, at a rate recalibrated over D's horizon (Stage C: 0.048584, V2, `results/gates.json`).
  - **Negative control:** reseeding from the compromised source.
- **What gates E.** A forecast with a better score does not by itself give a better intervention. Conversely, a simple battery rule may improve decisions even if the forecasters fail. So E is gated on:
  - a declared cost and loss trade-off;
  - D's pilot showing a usable lead time (2.7, criterion 3).

  It is not gated on a forecaster beating F0.
- **Engineering:** mid-run source switch (C8) and checkpoint-and-branch (C9) are not built, and they are the most expensive part of the roadmap.
- **Recommendation:** keep E conditional, and decide after D's pilot.

## 4. Literature added for D and E

**Search.** A model's own knowledge (no web), 20 candidates. Each DOI was checked against Crossref with `scripts/v018/verify_dois.py` (`docs/v018/evidence/agy_literature_de_check.json`):
- 10 match title, year and first author;
- 7 DOIs point to other works or carry other metadata;
- 3 DOIs do not exist.

Only matches are used. One more reference, [KB23], came from the review and was checked on the publisher's page; it has no DOI. As in Stage B, all were verified bibliographically, and their full texts were not read.

**New for D:**
- **Degradation models with first-passage times:** a nonlinear diffusion model [SWPZ12] and gamma processes [vN09]. Either could serve as the declared generative model in place of a deterministic ramp.
- **Forecast evaluation:**
  - strictly proper scoring rules [GR07];
  - calibration and sharpness of probabilistic forecasts [GBR07];
  - scoring rules for continuous distributions [MW76];
  - problems of the Brier score under administrative censoring [KB23], relevant only to a time-to-event secondary endpoint (2.1).

**New for E:** the optimum checkpoint interval for long-running computations [Dal06], for costing checkpoint-and-branch (C9).

**Not found by this search** (recall is weak, since the search used the model's own knowledge):
- a model of time-varying drift in a true random number generator;
- an analytic bound on bias accumulation in Monte Carlo estimates caused by a defective generator, beyond the established failures [FLW92].

**Existing anchors** (`RELATED_WORK.md`):
- RUL methods [SWHZ11];
- prognostics metrics [SCSSG10];
- entropy-source ageing [WSMG20];
- rejuvenation [HKKF95];
- reseeding recovery [DPRVW13];
- coupled chains [JOA20].

**References** (new; volume, issue and pages from Crossref, or from the publisher's page for [KB23]):
- [Dal06] J. T. Daly. A higher order estimate of the optimum checkpoint interval for restart dumps. *Future Gener. Comput. Syst.* 22(3), 303–312, 2006. doi:10.1016/j.future.2004.11.016
- [GBR07] T. Gneiting, F. Balabdaoui, A. E. Raftery. Probabilistic forecasts, calibration and sharpness. *J. R. Stat. Soc. B* 69(2), 243–268, 2007. doi:10.1111/j.1467-9868.2007.00587.x
- [GR07] T. Gneiting, A. E. Raftery. Strictly proper scoring rules, prediction, and estimation. *J. Am. Stat. Assoc.* 102(477), 359–378, 2007. doi:10.1198/016214506000001437
- [KB23] H. Kvamme, Ø. Borgan. The Brier score under administrative censoring: problems and a solution. *J. Mach. Learn. Res.* 24(2), 1–26, 2023. https://jmlr.org/papers/v24/19-1030.html
- [MW76] J. E. Matheson, R. L. Winkler. Scoring rules for continuous probability distributions. *Manage. Sci.* 22(10), 1087–1096, 1976. doi:10.1287/mnsc.22.10.1087
- [SWPZ12] X.-S. Si, W. Wang, C.-H. Hu, D.-H. Zhou, M. G. Pecht. Remaining useful life estimation based on a nonlinear diffusion degradation process. *IEEE Trans. Reliab.* 61(1), 50–67, 2012. doi:10.1109/TR.2011.2182221
- [vN09] J. M. van Noortwijk. A survey of the application of gamma processes in maintenance. *Reliab. Eng. Syst. Saf.* 94(1), 2–21, 2009. doi:10.1016/j.ress.2007.03.019

## 5. Decisions for the author (first option recommended)

**Adopted 2026-11-08: the first option of each decision.** Order: with about six months to the TDK submission, the D pilot runs as v0.19 alongside the TDK write-up of stages A to C.

1. **Stage D:** feasibility pilot for H2-R with the stop criteria of 2.7 / full H2-R campaign without a pilot / D as originally written (forecast at an application alarm) / STOP after C and write up A to C.
2. **Endpoint:** validity at declared run lengths / a prospective sustained-exclusion time with administrative censoring.
3. **Profiles:** a linear ramp as primary with slow rates kept in, held-out rates and onsets, and step-and-hold or random walk as the unseen family / another set.
4. **Fault families:** SHIFT, LAG1 and DUP (GFSR left out: no detectable excess harm in Stage C) / SHIFT only, which is cheaper.
5. **Exact reference in training:** allowed, with the claim stated as a supervised laboratory characterisation of this application / not allowed, so F-B and F-AB need an analytic observation model.
6. **Forecasters:** F0, F0+, F-B, F-AB and the oracle, with F-A exploratory / F0+, F-B and F-AB only.
7. **Stage E:** decide after D's pilot / plan it now in parallel.
8. **Versioning:** D as v0.19 (open decision 5 in `ROADMAP.md`) / within v0.18.

## 6. Limits of this review

- **No new data.** Statements about D's feasibility rest on Stage C outputs and on the design reasoning above. The pilot in 2.7 is the first test of the main risk.
- **Stage C's endpoint differs from D's.** Stage C validated one fixed window, not intervals at several run lengths, and not intervals under drift (2.1).
- **Literature:** recall is limited (Section 4). Before a D paper, a database search is still needed on two topics: drift models for random sources, and RUL forecasts for computations.
