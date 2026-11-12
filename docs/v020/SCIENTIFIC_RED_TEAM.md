# Stage D pilot: scientific red-team review of the frozen design

**Status: review for the author, 2026-11-12.** It challenges the design frozen on 2026-11-10 (`research/v019-d-pilot/PROTOCOL.md`, `protocol.json`, SHA-256 `5329464210019e968116e8c3c3a13a09f33b12e0368ecfaa69752a1eb43c5e97`). It changes no protocol, no threshold and no result. No pilot, smoke or campaign run exists, and none was made for this review.

**What the review rests on.** The frozen documents; the pipeline and its tests; the F-01 study (`F01_NUMERICAL_RESOLUTION.md`); the v0.19 audit and surrogate study (`docs/v019/IMPLEMENTATION_AUDIT.md`, `FORECASTING_METHODS.md`); and exact calculations in the design's own simplified models (`research/v019-audit/redteam_numbers.py`, output `docs/v020/evidence/redteam_numbers.json`). Nothing here is evidence about the real sampler under a ramp: that evidence does not exist yet.

**Limits.** The review was written by the same hand that built the pipeline and the independent recount. It is an internal challenge, not an external review. Every proposal below is a proposal: nothing is applied, and anything that touches the frozen protocol would be a versioned amendment that only the author can approve.

## 0. Findings by severity

| ID | Severity | Finding | Question |
|---|---|---|---|
| RT-01 | **Blocker until decided** | Under the frozen oracle rule, which runs have an oracle value depends on the last bits of δ, and the missing runs are the ones the oracle would call "almost surely invalid". The oracle score is then unavailable or, as an available-case figure, biased. The choice between amendment options A, B and C must be made before fitting, and before any outcome is seen | 6 |
| RT-02 | **Blocker until decided** | Any analysis the author wants beyond the frozen ones (RT-04, RT-05) has to be declared before the first pilot draw. Added afterwards it is exploratory and must be labelled so | 2, 4, 7 |
| RT-03 | Major | The target "valid at L" is interval coverage of a known exact value. It is not observable in use, it is noisy even for a healthy run, and drift can make a biased run "valid" by widening its interval | 1 |
| RT-04 | Major | Criterion 4 cannot tell prediction from recognition of a state that is already decided. No baseline uses the state at issuance, and the score is not stratified by it | 7 |
| RT-05 | Major | The baseline F0+ is a straight line in the alarm sweep. Features that correlate with the alarm sweep can earn "skill" by repairing that line rather than by adding battery information | 4 |
| RT-06 | Major | Fitting and evaluation runs come from the same declared ramp distribution, and the forecasters are told the family. A pass shows skill inside the declared family only | 3, 7 |
| RT-07 | Major | A family that passes criterion 4 at the edge reports an inflated skill: at a true skill of 0.10 the mean observed skill of passing families is 0.124 to 0.129 | 5 |
| RT-08 | Major | The oracle leaves out the widening of the interval by drift (F-02) and assumes a bias that is linear in the mean severity. It is a reference, not an upper bound; "the forecaster beats the oracle" would mean nothing | 6 |
| RT-09 | Major | One failed or interrupted calibration run makes the threshold unavailable (protocol, `failed_runs.calibration`), and criteria 2 to 4 then count as failed. Under the policy `record_failure` a power cut can therefore decide the pilot | 5, 8 |
| RT-10 | Major | A negative criterion 4 is weak evidence: the probability of passing is about 0.5 at a true skill of 0.10. It supports "not demonstrated", not "no predictive information" | 8 |
| RT-11 | Minor | Short intervals at issuance under-cover in healthy runs if the v0.18 variance ratios carry over (0.881 at an alarm at sweep 512, 0.935 at 2,048). "Already invalid at issuance" then partly measures the estimator, not the fault | 2, 7 |
| RT-12 | Minor | The statistic at the alarm is selected by the alarm: it sits at the threshold when its own test fires. The slope uses four windows that overlap by half | 3 |
| RT-13 | Minor | Summaries of fault runs contain evaluation outcomes. Run before fitting, they show them to the analyst. The models have no free choice left, so the channel is the analyst's later decisions | 2 |
| RT-14 | Minor | Criterion 1 can end as "insufficient precision" with noticeable probability even when coverage is as in Stage C (0.82 to 0.94 probability that all three lengths pass at a true coverage of 0.94) | 8 |
| RT-15 | Minor | For SHIFT the target statistic is not what fires the alarm (F-03, one run). Criterion 4 is then unlikely to pass for SHIFT, which leaves two families for "at least one" | 4 |
| RT-16 | Future work | Held-out rates, onsets and families; a second application; a tolerance-based target; field evidence that a source degrades along a ramp | 1, 3, 7 |

No finding shows that the pilot, read as what it says it is (a feasibility screen with descriptive results), is invalid. RT-03 to RT-10 limit what a pass or a fail may be said to mean.

## 1. Is the forecasting target operationally meaningful?

**The target.** y = 1 if the 95% batch-means interval of the CV mean over sweeps [0, L) contains the exact finite-bead energy E_M.

- **It needs the exact answer.** A user never knows whether an interval covers the truth. The target is a laboratory construct, and the protocol says so ("a known-family laboratory experiment").
- **It is noisy for a healthy run.** Healthy coverage in this sampler was 3846/4096 and 3867/4096 in earlier stages. So about 6% of healthy runs are "invalid", and no forecaster can predict which. This sets a floor on the Brier score that is unrelated to the generator.
- **Drift can buy validity.** A ramp makes the batch means trend. The spread between batches grows, the interval widens, and a biased estimate can be covered. The protocol states that validity is not monotone in L and "does not imply an unbiased or accurate estimate". A family could therefore pass criterion 2 (valid fraction in [0.2, 0.8]) because its intervals inflate, not because its outcomes are undecided in any useful sense. The signed error and half-width are stored per run and length, so this can be seen in the results; the criteria do not look at it.
- **What an operator needs** is closer to "will the error exceed a tolerance" or "how much run is left before it does". Criterion 3 measures run left after the alarm, not time before invalidity, as the protocol notes.

**Assessment (RT-03, major).** The target is well defined and exactly scorable, which is its merit. Its operational meaning is limited, and any report should say "forecast of interval coverage against a known reference", not "forecast of result reliability". A tolerance-based target is campaign material (RT-16).

## 2. Are calibration and evaluation genuinely separated?

**Separated by construction, and checked:**
- The five conditions have distinct seed labels; fitting and evaluation are distinct replicas of the same condition. The test `test_analysis_seeds_are_distinct_from_each_other_and_from_run_seeds` checks distinctness on the frozen plan.
- Healthy calibration (threshold, application constants, σ_eff) and healthy validation (criterion 1, gate V2-D) are different runs.
- The fitting stage opens calibration and fitting records only. The test `test_fit_opens_fitting_runs_only_and_evaluation_is_opened_once` records every run the stage loads. The model file is hashed before the evaluation stage starts, and the evaluation stage verifies the hash and refuses a second run.
- The bootstrap does not refit. Its intervals are conditional on the fitted models, which the protocol states.

**Where separation is procedural only:**
- **RT-13.** Criteria 2 and 3 use all 1,000 fault runs per family, and the fault summary reports validity at L for them. Whoever reads that summary before fitting has seen the evaluation outcomes in aggregate. The models cannot use this: every feature, option and fallback is frozen. The remaining channel is a human decision taken afterwards, above all an amendment. The order tested here is generate, calibrate, fit, issue, evaluate, summaries, criteria; keeping to it closes the channel.
- **RT-01 and RT-02** follow from the same point: the oracle rule and any added analysis must be fixed before outcomes exist.
- **RT-11.** "Already invalid at issuance" uses 16 batches of n/16 sweeps, where n is the alarm sweep. For early alarms the batches are short. If the batch-means variance ratios measured in v0.18 carry over (0.601 at 32 sweeps, 0.875 at 128, 0.971 at 512), healthy coverage of these intervals would be 0.881 for an alarm at sweep 512, 0.935 at 2,048 and 0.947 at 8,192 (`redteam_numbers.json`, normal batch means). The protocol checks these intervals on healthy runs, so the size will be measured. Criterion 2's bound (fewer than half already invalid) has room for it.

**Assessment.** Separation of the data is sound. The open risk is the order of human decisions, not the code.

## 3. Does the design create hidden information leakage?

**Checked and not found:**
- Features use battery checkpoints up to and including the alarm, and the last application checkpoint at or before it. The test `test_features_use_nothing_after_the_alarm` changes the target statistic and the application checkpoints after the alarm and requires identical features.
- The profile draws (t₀, T_cap) and the severity path are stored in the record for I-2, I-3 and the oracle. The feature code does not read them. The oracle does, and is labelled as using future information.
- Standardisation uses fitting-run means and deviations only.

**Present by design, and to be stated in any report:**
- **RT-06. In-distribution skill.** Fitting and evaluation profiles are drawn from the same law (onset uniform on [2,048, 8,192], T_cap log-uniform on [2,048, 163,840]). A model learns that law. The alarm sweep alone then says a great deal: an early alarm means a steep ramp. This is why F0+ is the baseline, and it is why skill here is skill under a known prior. The surrogate study found that forecasters fitted on the declared family were miscalibrated and worse than the alarm-time baseline on two other families (skill −0.256 and −0.218). Held-out rates and families are deferred to the campaign.
- **Known family and aimed statistics.** One model per family, with a target statistic chosen for that family's mechanism. No classifier is tested.
- **RT-12. Selection by the alarm.** The alarm is the first strict exceedance of S, the largest −log₁₀ p over ten p-values. When the family's own test fires it, the target statistic is at its threshold value by construction and carries little beyond the alarm sweep. It varies mainly when another test fires. The slope is computed from at most four windows of 512 sweeps at a stride of 256, so neighbouring values share half their data. Both effects are the same in fitting and evaluation: they weaken the features; they do not leak.
- **Pre-onset alarms** are kept and scored. If false alarms were spread evenly over the 159 checkpoints, about 0.6% of fault runs would alarm before onset (`redteam_numbers.json`; an assumption, not a measurement). Their features hold no fault information. The effect on scores is small.

**Assessment.** No leakage from evaluation into fitting was found. The design's information advantage is the declared ramp law, which is disclosed, and RT-06 limits the claim.

## 4. Are the baseline comparisons appropriate?

- **F0 (base rate) and F0+ (alarm sweep)** are the right kind of baseline: they use what is known at issuance without the battery statistic.
- **RT-05. F0+ is linear on the logit scale.** The relation between alarm sweep and validity under a capped ramp with random onset need not be. F-B-lite keeps the same linear term and adds two features that are correlated with the alarm sweep. Part of a measured skill could be a repair of the baseline's shape. A flexible function of the alarm sweep alone (binned, or with a quadratic term) would separate the two. It could be added as a descriptive baseline without touching criterion 4, but only by an amendment written before generation (RT-02).
- **RT-04. No state baseline.** See Section 7.
- **RT-15. SHIFT.** In one full-length run the frequency statistic was not what fired the alarm (F-03). If that holds, F-B-lite's target feature is weak for SHIFT. Criterion 4 needs one family of three, at a level corrected for three, so the criterion is still usable; the family-wise chance of passing is lower than the design's power table suggests if one family cannot pass.
- **The oracle** is a reference built with future information, with two known defects (RT-01, RT-08). It should be presented as "the declared analytic model", not as the best achievable score.
- **F-AB-lite against F-B-lite** is descriptive. The surrogate has no application monitor, so nothing is known in advance about this contrast.

## 5. Can the stopping rule introduce selection bias?

**No optional stopping.** Sample sizes are fixed, evaluation happens once, thresholds were fixed before generation, and "insufficient information" counts as failed. That part is sound.

**RT-07. Selection of the passing family.** Criterion 4 passes a family only if its observed skill is at least 0.10. In the simplified model that `design.py` uses for the power (exact enumeration, `redteam_numbers.json`):

| Eligible runs | True skill | P(family passes) | Mean observed skill | Mean observed skill of passing families |
|---|---|---|---|---|
| 240 | 0.05 | 0.034 | 0.050 | 0.112 |
| 240 | 0.10 | 0.527 | 0.100 | 0.129 |
| 240 | 0.15 | 0.758 | 0.150 | 0.169 |
| 400 | 0.05 | 0.009 | 0.050 | 0.108 |
| 400 | 0.10 | 0.492 | 0.100 | 0.124 |
| 400 | 0.15 | 0.919 | 0.150 | 0.156 |

A pilot that passes must not hand its skill estimate to the campaign as the expected effect. The campaign needs its own seeds, must not pool pilot runs, and should size itself for a skill below the pilot's estimate. The model has no fitting uncertainty and fixed eligibility, as in `design.py`.

**Selection by eligibility.** Only runs with an alarm before L are scored. Steep ramps alarm early and are over-represented among the eligible. All models are scored on the same runs, so the comparison between models is fair; the population it describes is "alarmed runs under the declared law", and the counts of runs without an alarm are reported per L.

**RT-09. Failures and interruptions.**
- A run with a failure record is never run again and keeps its place in every denominator. That prevents selection by re-drawing.
- The threshold needs exactly 999 usable calibration maxima. One missing calibration run makes it unavailable, issuance is then unknown for every run, and criteria 2 to 4 are "insufficient information", which counts as failed.
- A run that was started and left no record (a killed process, a power cut) is not covered by the protocol. The pipeline refuses to continue until the author has chosen a policy. Under `regenerate_same_id` the same seed labels give the same bytes, so nothing is selected and nothing is lost. Under `record_failure` the run is lost; if it is a calibration run, the consequence above follows. The cause of an interruption is unrelated to the outcome of the run, so neither policy biases an estimate. They differ in what an accident costs. The choice is the author's.

**What "reconsider the scope" allows** after a failed criterion 3 or 4 is not defined. That leaves room for a redesign informed by the pilot's data. A campaign protocol should be frozen on fresh seeds, and it should say which pilot results it used.

## 6. Could oracle failures bias the conclusions?

**The stop criteria do not use the oracle.** No oracle defect can change the result of criteria 1 to 5. What is at stake is the reference line in the report, and the agreement between the pipeline and the independent recount.

**RT-01. Missing values are not random, and not reproducible.**
- `scipy.stats.nct.cdf(−c, 15, δ)` returns NaN for a set of noncentralities between 8.53 and 43.61. Under the frozen rule such a value is "unavailable with reason", and one unavailable run makes the cell's primary oracle score unavailable.
- The runs concerned have a large bias relative to the standard error. They are the runs for which the oracle's answer would be "invalid with probability above 1 − 10⁻⁸". Leaving them out removes the oracle's easiest cases. An available-case score would make the oracle look worse than it is.
- Whether a value exists flips between neighbouring doubles: on 20,000 points in [8.6, 43.5], 8,655 were not finite, and 521 changed state at the next representable δ (`research/v019-audit/results/f01_oracle.json`). Two correct implementations that sum in a different order will disagree on which runs are unavailable. On the synthetic campaign the pipeline and the recount agreed on 1,972 of 1,974 fields; the two that differed were unavailable counts. With option A forced in memory they agreed on all 1,974 (`evidence/synthetic_agreement.json`).
- On the synthetic campaign the frozen rule left 6 of 9 oracle cells without a primary score. This says how the rule behaves when large δ occurs. It does not say how often large δ will occur in the pilot.

**RT-08. Model defects that remain under any option.**
- The formula ignores the widening of the interval by drift (F-02). In the surrogate study the formula without widening was miscalibrated at L = 20,480. The direction in the real sampler is not known.
- β is one slope through the origin, pooled over three lengths. Linearity of the bias in the mean severity is assumed, not shown, for this sampler.

**Assessment.** The oracle cannot bias the stop decision. It can bias the comparison "how far are the forecasters from the reference", and under the frozen rule it makes the recount gate depend on rounding. The decision between A, B and C is needed before fitting (`F01_NUMERICAL_RESOLUTION.md`, Section 8). Under C the report needs a declared exception for oracle availability in the recount, which is also the author's decision.

## 7. Can Stage D tell predictive information from recognising an already degraded state?

**What the protocol provides.** The share already invalid at issuance (criterion 2 requires less than half), the run left after the alarm (criterion 3), and the table of transitions from validity at issuance to validity at L. These describe the population. They do not enter the score.

**RT-04. What is missing.**
- Criterion 4's skill is one number over all eligible evaluation runs. It can come entirely from runs whose outcome was settled at the alarm: already invalid and staying so, or already at the severity cap. A model that recognises "this run is far gone" earns the same skill as one that foresees a change.
- No baseline carries the state at issuance. A logistic model on the alarm sweep and the indicator "valid at issuance" is the natural control. It uses the exact reference, so it is not an operational forecaster; it is the right yardstick for this question.
- Two descriptive additions would answer the question without changing a criterion: the Brier scores of every model split by validity at issuance and by whether the severity had reached the cap at the alarm; and the state baseline above. Both need an amendment before generation (RT-02). If the author prefers not to amend, they can still be computed afterwards and labelled exploratory.

**A limit no analysis removes (RT-06).** Under a ramp, the future is a deterministic function of the present state (onset, rate, cap) and the family's law. Estimating the state and applying the declared continuation is the whole forecast. Two continuations with the same history cannot be told apart by any forecaster (`FORECASTING_METHODS.md`, Section 1; the excess Brier score is at least E[(p₁ − p₂)²]/4 in one of them). So the pilot can show that the signals at the alarm carry information about the outcome beyond the alarm time, under the declared law. It cannot show prediction in a sense that would survive a different continuation. That needs held-out families.

**RT-11** matters here too: with short issuance intervals, "already invalid" includes healthy under-coverage, so a split by that label is slightly blurred. The healthy short-interval check gives the size.

## 8. Would a negative outcome still support the pre-registered claims?

**What is pre-registered.** The pilot decides whether a campaign is worth running. It is declared "not a test of H2-R", and every result is descriptive. So no outcome of the pilot, positive or negative, supports or refutes H2-R. What each negative supports:

| Outcome | What it supports | What it does not support |
|---|---|---|
| Criterion 1: outside the band | The 16-batch interval mis-covers in healthy runs at these lengths. The endpoint itself is unreliable. Relevant to the uncertainty-estimation line of the programme | Anything about forecasting |
| Criterion 1: straddles an edge | Insufficient precision (RT-14). Not miscalibration, and not a stop | — |
| Criterion 2 fails | Under the declared ramp law the outcomes at L = 20,480 are mostly decided, or mostly decided at the alarm | That forecasting is trivial or impossible in general (the protocol says this) |
| Criterion 3 fails | The battery alarm comes too late to leave run to act on | Anything about validity |
| Criterion 4 fails | These three-feature models, at this sample size, did not show a skill of 0.10 over the alarm sweep | Absence of predictive information (RT-10). The interval for the skill is reported; only an upper end well below 0.10 would support "little skill for these models" |
| Gate V2-D fails | The alarm is not calibrated; criteria 2 to 4 are uninterpretable | — |
| A failed implementation check | Nothing. Interpretation is blocked | — |
| A missing calibration run (RT-09) | Nothing scientific. The pilot fails for an engineering reason | — |

**What stands whatever the pilot gives:** the Stage C result on detection; the identifiability bound; the surrogate results, stated as surrogate results; and, new here, the numerical findings on the oracle formula.

**The asymmetry to avoid.** A pass leads to a campaign; a fail leads to "reconsider the scope". If every fail leads to a new design, H2-R is never exposed to refutation. Before generation the author could state in one sentence what a failed criterion 2 or 4 will mean for the thesis, for example "reported as not shown feasible under the declared ramp with these models". This is advice on framing, not a protocol change.

## 9. Decisions this review puts to the author

None is taken here.

1. **Oracle rule (RT-01):** option A, B or C of proposed Amendment 01, before fitting. Under C, also the wording of the recount exception.
2. **Supplementary descriptive analyses (RT-04, RT-05):** amend before generation to add a state baseline, a flexible alarm-time baseline and score splits by state at issuance; or leave the protocol as it is and treat any such analysis as exploratory.
3. **Interrupted-run policy (RT-09):** `regenerate_same_id` or `record_failure`, recorded in the pre-generation manifest.
4. **Order of stages (RT-13):** confirm that summaries of fault runs are produced after fitting.
5. **Wording of the outcome (Sections 1, 7, 8):** "coverage against a known reference", "information beyond the alarm time under the declared law", and what a failed criterion will be reported as.

## Reproduction

    python research/v019-audit/redteam_numbers.py          # docs/v020/evidence/redteam_numbers.json
    PYTHONPATH=src python research/v019-audit/f01_oracle.py  # research/v019-audit/results/f01_oracle.json
    PYTHONPATH=src python tests/v019/synthetic_agreement.py  # docs/v020/evidence/synthetic_agreement.json (needs the measurement layer)
