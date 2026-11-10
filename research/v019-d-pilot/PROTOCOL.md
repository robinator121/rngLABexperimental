# v0.19 Stage D feasibility pilot (H2-R): protocol

**Status: FROZEN DESIGN, 2026-11-10.** The author approved the design ("Approve as recommended": the first option of every decision in Section 9), then the implementation specification `protocol.json` ("Accept and freeze"). Its SHA-256 is in `protocol.sha256`. Before the freeze, two wordings were corrected to the Stage C conventions (the LAG1 first word in Section 2.2, gate V2-D in Section 2.4), and criterion 3 received an eligibility minimum (`protocol.json` open questions Q01, Q08 and Q09). Revised before approval after an adversarial review (24 findings, checked and addressed). Nothing has been generated. It implements Section 2.7 of `docs/v018/STAGE_DE_REVIEW.md` (decided 2026-11-08).

**Freeze procedure:**
1. The author approves this design. Done 2026-11-10.
2. The implementation specification `protocol.json` is written (Section 7). Done 2026-11-10: drafted, then checked against the code and in a second, independent review. It fixes every field, formula and convention that this document leaves to it, and it must agree with this document. The design numbers of Sections 2.2, 3 and 5 are already computed by `design.py` (output `design.json`), and are recomputed if a decision changes them.
3. The author approves the design numbers, and `protocol.json` is frozen with its SHA-256 recorded. Done 2026-11-10.
4. The engineering is built and tested (Section 7).
5. A smoke run is made in its own namespace and is excluded from every result.
6. A pre-generation manifest, with the hashes of all code including the imported Stage C modules, is pushed before the first pilot draw.
7. Generation, then fitting, then evaluation (Section 4.4).

Any change after step 3 is a versioned amendment, written before the affected analysis. Failed runs are reported, never silently replaced.

## 1. Purpose and status of the results

- **What the pilot decides:** whether a Stage D campaign is worth running.
- **What it is not:** a test of H2-R. Every result is descriptive, and is compared only with the stop criteria of Section 5. Passing the pilot does not validate the campaign features deferred in Section 6.
- **A known-family laboratory experiment:** the forecasters are told the injected fault family. No operational family classifier is tested.
- **Negative outcomes** are reported with the same prominence as positive ones.

## 2. Definitions

### 2.1 Application

Unchanged from Stage C (`research/v018-h1-detection/PROTOCOL.md` §2.1):
- `PathIntegralHarmonicOscillator(beads=16, temperature=0.5, displacement_delta=0.4)`;
- the centroid refresh and an exact start;
- no burn-in.

Every run lasts **H = 40,960 sweeps**, numbered t = 0, …, H − 1.

### 2.2 Degradation profile (C2)

- **Severity per sweep:** s(t) = min(s_max, r · (t − t₀)) for t ≥ t₀, and 0 before. It applies to every word consumed during sweep t. So s(t₀) = 0, and the first positive severity is at sweep t₀ + 1.
- **The transformation starts at sweep t₀.** It starts at the first word consumed in sweep t₀, which defines the onset byte offset B₀. The Stage C state rules then apply:
  - LAG1 initialises z = Φ⁻¹(u + 2⁻⁵⁴) at that word, the midpoint rule Stage C uses for every word (its `protocol.json`, `implementation_notes_added_at_freeze.lag1_first_word`), keeps its AR(1) state across severity changes, and uses ρ_t in both coefficients;
  - DUP draws one fault-generator value for every word from B₀ on, including while s(t) = 0, and "the previous word" means the previously emitted word;
  - SHIFT is unchanged apart from s(t).
- **Caps:** each family's highest Stage C severity, at full precision from `research/v018-h1-detection/results/gates.json`. These are SHIFT-5, LAG1-5 and DUP-5. Held constant, these severities left 245/256, 224/256 and 256/256 runs wrong over the Stage C window (T2). Under a ramp, the outcomes are not predicted by those counts.
- **Per-run draws**, from the `profile` generator (NumPy `default_rng`), in this order:
  - t₀, an integer uniform on [2,048, 8,192];
  - T_cap, log-uniform on [2,048, 163,840] sweeps, as a float64, with r = s_max / T_cap.

  The two are independent. The range of T_cap includes slow ramps that never come near the cap within H. In `design.json`, 0.19, 0.46 and 0.65 of the profiles reach the cap before L = 10,240, 20,480 and 40,960.
- **Recorded per run:** t₀, T_cap, r, s_max, B₀, and the SHA-256 of the prefix before B₀.

### 2.3 Endpoint: validity at declared run lengths (C6)

- **Lengths:** L ∈ {10,240, 20,480, 40,960}. The estimate is the CV mean over sweeps [0, L).
- **Why sweep 0:** the exact start removes initialisation bias. Whether the intervals cover at the nominal rate is a separate question, checked by criterion 1. The repository's healthy results for this sampler show no reason to expect failure: Phase 2 (3846/4096) and Stage C V1 (3867/4096).
- **Interval:** batch means with exactly 16 batches of L/16 sweeps, and t₁₅. The batch result is `summarize_series(cv[0:L], reference=E_M, batch_sizes=(L//16,))["batch_means"][str(L//16)]`, with E_M = `exact_finite_bead_energy(16, 2.0)`.
- **Valid at L:** that batch result has status `ESTIMATED` and `covered` is true. Anything else is invalid, and the reason is recorded. An exception or non-finite input counts as invalid and is reported separately. Every declared length is divisible by 16.
- **Validity is not monotone in L.** Drift between batches can widen the interval, so a run invalid at one L can be valid at a later one. All three-length patterns are reported. So, per run and per L, are:
  - the signed error (estimate − E_M) and the interval half-width;
  - the 16 batch means;
  - the severity path.

  Validity here means interval coverage. It does not imply an unbiased or accurate estimate.
- **Stored per run:** the per-sweep CV and TD traces and byte marks to H. The design script gives the storage.

### 2.4 Operational alarm

- **Issuing battery:** the full Stage C battery of nine tests, ten p-values (decision 1). The SP 800-22 subset of Sensitivity 01 is recalibrated on the same calibration runs and reported alongside, descriptively.
- **Windows:**
  - battery: 512-sweep windows every 256 sweeps, so 159 checkpoints with ends 512, 768, …, 40,960;
  - application: 1,024-sweep windows on the same stride, so 157 checkpoints from 1,024.
- **Window cap:** if a battery window is longer than 1,000,000 bits, the SP 800-22 bit tests use its first 1,000,000 bits, the limit of the in-repository implementation. The word-level tests use the whole window. The original and used byte counts are recorded. A 512-sweep window can reach 1,048,576 bits, at 256 bytes per sweep.
- **Statistic:**
  - at each checkpoint, S = −log₁₀(max(min p over the available tests, `P_FLOOR`)), as in Stage C;
  - unavailable p-values are omitted;
  - if none is available, S = 0.
- **Threshold:** the 950th order statistic of the whole-horizon maxima of 999 healthy calibration runs. An alarm is a strict exceedance. This gives a nominal marginal whole-horizon false-alarm rate of 0.05.
- **Gate V2-D:** the upper end of the two-sided 95% Clopper–Pearson interval of the healthy-validation whole-horizon alarm rate must be at most 0.10, as for Stage C's V2. Otherwise every issuance-based result is reported as uninterpretable.
- **Issuance:** the first exceedance, with no onset knowledge and no re-arming. Pre-onset alarms are kept, because they are operational alarms, and are reported separately.
- **Eligibility for a forecast at L:** the alarm sweep is below L. Runs with no alarm, or with an alarm at or after L, are reported per L and are not scored.
- **Already invalid at issuance:** the same validity rule as Section 2.3, applied to cv[0:n], where n is the alarm sweep rounded down to a multiple of 16. These much shorter intervals are checked separately on healthy runs.

### 2.5 Signals stored for forecasters (C7, logging part)

**Per battery checkpoint**, in addition to the p-values and their statuses:
- the signed lag-1 correlation r, with the word count and pair count;
- the adjacent-duplicate count, with the pair count;
- the signed standardised frequency imbalance and the signed runs residual;
- the remaining test statistics, each with its formula, sign, sample count and unavailable-value convention fixed in `protocol.json`;
- the bit count, and whether the window was cut.

**Per application checkpoint:**
- A1, the signed window mean of TD − CV over 1,024 sweeps, with its count;
- A2, the uphill-proposal count over the 16 × 1,024 proposals, computed from the complete local-byte ledger.

## 3. Data

| Set | Condition key | Replicas | Use |
|---|---|---|---|
| Calibration, healthy | `calibration` | 0–998 | Alarm thresholds, application constants, σ_eff |
| Healthy validation | `healthy` | 0–1,999 | Criterion 1, gate V2-D, short-interval check |
| SHIFT, LAG1, DUP | `fault/SHIFT` etc. | 0–999 each | Fitting: replicas 0–599. Evaluation: replicas 600–999 |

- **Seeds:** `namespace_seed("rnglab-v019-d-pilot", condition, role, f"{i:05d}")` for the roles `local`, `centroid`, `start`, `fault` and `profile`. The smoke uses `rnglab-v019-d-pilot-smoke`. Neither namespace appears in the repository or in local working folders (searched 2026-11-10).
- **No pairing:** families do not share seeds, so the fitting and evaluation sets are disjoint trajectories.
- **Size:** 5,999 runs. `design.py` extrapolates from Stage C's cost per sweep (`research/v018-h1-detection/results/report_numbers.json`) to about 236 minutes on 8 workers and 5.5 GiB of uncompressed traces (`design.json`).

## 4. Analysis (descriptive)

### 4.1 Healthy runs

- coverage at each L, with 95% Clopper–Pearson intervals;
- the alarm rate (gate V2-D);
- the coverage of the short issuance-time intervals.

### 4.2 Fault runs, per family

- issuance outcomes: pre-onset alarm, post-onset alarm, no alarm;
- the alarm sweep minus t₀;
- validity at each L, and the transition from validity at issuance to validity at L;
- the quantities listed in Section 2.3, against the severity path.

### 4.3 Oracle reference

1. **Cumulative severity:** S̄(L) = L⁻¹ Σ_{t<L} s(t), the exact discrete sum including the cap.
2. **Bias slope:** a per-family slope β, fitted by least squares through the origin of the signed error (estimate − E_M) on S̄(L), pooled over the three lengths, on the fitting runs only.
3. **Noise scale:** σ_eff, the long-run standard deviation of CV. It is estimated per L as √L times the mean batch-means standard error over the healthy calibration runs, so that the standard error at L is σ_eff/√L.
4. **Probability:** P(valid at L) = F_{t₁₅,δ}(c) − F_{t₁₅,δ}(−c), the noncentral t distribution, where δ = β·S̄(L)/(σ_eff/√L) and c is the 0.975 quantile of t₁₅.

This ignores the widening of the interval by drift. It is a benchmark, not a bound. Its own calibration on the evaluation runs is reported separately from the forecasters.

### 4.4 Feasibility forecasters

**Scope.** One model per family and per L. Fitted on fitting runs only, and scored once on eligible evaluation runs.

**Fitting procedure.** L2-regularised logistic regression (C = 1, lbfgs) on inputs standardised with fitting-run means and standard deviations. If a fitting set has only one class, or the fit fails to converge, the model falls back to F0 and this is reported.

**Order.** Fitting writes the model file, and its hash is recorded, before the evaluation outcomes are read.

| Model | Inputs |
|---|---|
| F0 | none: the fitting-run fraction valid at L among eligible runs |
| F0+ | the alarm sweep |
| F-B-lite | F0+ plus the family's target battery statistic at the alarm checkpoint, and its least-squares slope over the last four checkpoints. SHIFT: the signed frequency imbalance. LAG1: the signed lag-1 r · √pairs. DUP: duplicates per pair |
| F-AB-lite | F-B-lite plus A1/σ_D and the standardised A2 deviation, at the last application checkpoint at or before the alarm. If none exists (an alarm at 512 or 768), F-AB-lite returns the F-B-lite prediction, and this is reported |

**Scores, per L:**
- **Brier score**, and Brier skill 1 − BS_model/BS_F0+. If BS_F0+ = 0, skill is undefined: insufficient information.
- **Calibration** in five bins of width 0.2, with the count in each bin.
- **Paired-bootstrap percentile intervals** over evaluation trajectories, with 10,000 resamples. The seed comes from the `analysis` role of the pilot namespace, and models are not refitted.
- **Contrasts:**
  - F-B-lite against F0+ (criterion 4);
  - F-AB-lite against F-B-lite, descriptive only.

## 5. Stop criteria (thresholds fixed before generation)

| # | Criterion | Passes if | If it fails |
|---|---|---|---|
| 1 | Healthy validity | For every L, the 95% CP interval of healthy coverage (2,000 runs) lies inside [0.92, 0.98], the Stage C V1 band. Per L, this accepts 1,864 to 1,946 covered runs. The probability of passing is 0.384 at a true coverage of 0.93, 0.938 at 0.94 and 0.9998 at 0.95 (`design.json`) | Report whether the interval straddles a band edge (insufficient precision) or lies outside it (miscalibration). Only miscalibration **stops** D as designed |
| 2 | Outcome-balance screen | For at least two families, among all fault runs eligible at L = 20,480 (at least 300 per family): the fraction valid at 20,480 lies in [0.2, 0.8], and fewer than half were already invalid at issuance | **Stop.** The declared ramp model gives mostly decided outcomes. This does not show that forecasting is intrinsically trivial |
| 3 | Remaining run after issuance | For at least two families, the median of 20,480 − alarm sweep among eligible runs (at least 300 per family, as for criterion 2) is at least 5,120 sweeps | Reconsider the scope. This criterion measures run left to act on, not time before invalidity |
| 4 | Signal | For at least one family, F-B-lite has a Brier skill over F0+ at L = 20,480 of at least 0.10, with the lower bound of its paired-bootstrap interval at level 1 − 0.05/3 above 0, on at least 240 eligible evaluation runs. The power is about 0.5 at a true skill of 0.10, because the observed skill must reach 0.10; 0.78 to 0.91 at 0.15; and 0.92 to 0.99 at 0.20, over eligible fractions 0.6 to 1 (`design.json`, a simplified calibrated model without fitting uncertainty) | Reconsider the scope: skill was not demonstrated by these models at this sample size |
| 5 | Schedule | The engineering of Section 7 is complete within six weeks of approval, and pilot results exist at least three months before the TDK deadline. The author records the deadline in the pre-generation manifest; without it, the result is "insufficient information" | Reconsider the scope, or keep D as the declared programme in the TDK |

**Criterion 2 partly reflects the declared ramp distribution.** It checks that the declared model yields a problem worth forecasting; it does not show that real degradation would.

**Eligibility minimums.** If an eligibility minimum is not met, the criterion's result is "insufficient information", and it counts as failed for the stop rule.

## 6. Deliberately deferred to the campaign

The pilot leaves out several features the campaign would need:
- held-out ramp rates and onsets;
- the unseen profile family;
- the Bayesian forecaster with forward simulation (Section 8);
- an oracle predictive distribution;
- sharpness beyond the calibration bins;
- abstention;
- forecasts at checkpoints after issuance;
- the multiplicity hierarchy of `STAGE_DE_REVIEW.md` §2.4.

## 7. Engineering and checks before the smoke

All new code lives in `research/v019-d-pilot/`. It imports or copies the Stage C library and leaves the frozen Stage C files unchanged. The effort ratings come from an engineering review of the Stage C code (2026-11-10).

| Item | Content | Effort |
|---|---|---|
| `protocol.json` | Every field, formula, sign, unit, dtype and unavailable-value convention left open above | Small |
| `design.py` | Runtime, storage, the criterion 1 operating characteristics, the criterion 4 power scenarios, and the share of profiles reaching the cap by each L | Done (`design.json`) |
| C2 | The profile inside the Stage C adapter transform, driven by a sweep counter advanced at the end-of-sweep callback; the replay path receives the sweep byte marks | Medium |
| C6 | The interval at each L and at issuance, from the stored trace | Small, plus validation |
| Window cap, C7 logging, onset-blind alarm | Sections 2.4 and 2.5 | Small |
| Forecasters and scores | Section 4.4 | Small to medium |

**Implementation checks.** These carry forward Stage C's I-1 to I-4, and they block interpretation if any fails:
- **I-1, sweep clock:** H marks per run, monotone, with every sweep consuming 128 + 8k bytes.
- **I-2, profile:** an independent recomputation of s(t) and B₀ from the recorded draws and marks.
- **I-3, fault realisation:** an independent implementation re-applies the profile to the recorded base stream. **A step profile** (s(t) = s_max from t₀) reproduces the Stage C transform byte for byte, as a separate test from the ramp tests.
- **I-4, replay:** a sample of runs regenerated bit for bit.

**Independent recount.** A separate script recomputes the endpoint, the issuance, the eligibility and the scores.

## 8. Sources

- **Design:** `docs/v018/STAGE_DE_REVIEW.md` (decided 2026-11-08); `research/v018-h1-detection/SENSITIVITY_01_RESULTS.md`.
- **Engineering map and adversarial review:** read-only reviews of the Stage C code and of the draft, 2026-11-10.
- **Method suggestions** (2026-11-10). These were the Bayesian update of the ramp rate with forward simulation (deferred to the campaign), the analytic validity approximation (corrected here for the cap and the t distribution), and draft thresholds (adjusted).
- **References:** checked against Crossref (`docs/v018/evidence/dpilot_refs_check.json`: 5 match, 1 DOI does not exist and was dropped). As before, they were verified bibliographically; their full texts were not read.
  - [GLLR05] N. Z. Gebraeel, M. A. Lawley, R. Li, J. K. Ryan. Residual-life distributions from component degradation signals: a Bayesian approach. *IIE Trans.* 37(6), 543–557, 2005. doi:10.1080/07408170590929018. This is the Bayesian rate updating that the campaign forecaster would follow.
  - [Bri50] G. W. Brier. Verification of forecasts expressed in terms of probability. *Mon. Weather Rev.* 78(1), 1–3, 1950. doi:10.1175/1520-0493(1950)078<0001:VOFEIT>2.0.CO;2
  - [Sch82] B. Schmeiser. Batch size effects in the analysis of simulation output. *Oper. Res.* 30(3), 556–568, 1982. doi:10.1287/opre.30.3.556
  - Already in the programme: [GR07], [GBR07], [SWHZ11], [KB23].

## 9. Decisions for the author (first option recommended)

**Decided 2026-11-10:** the first option of all seven.


1. **Issuing battery:** the full Stage C battery, with the SP 800-22 subset reported alongside / the SP 800-22 subset / both co-primary.
2. **Horizon and lengths:** H = 40,960 with L ∈ {10,240, 20,480, 40,960} / H = 20,480 with L ∈ {10,240, 20,480}.
3. **Ramp parametrisation:**
   - onset uniform on [2,048, 8,192];
   - T_cap log-uniform on [2,048, 163,840];
   - caps at the Stage C maxima;

   or another range.
4. **Stop criteria and thresholds:** as in Section 5 / changed.
5. **Clock:** per sweep / per consumed word (the smaller code change, but future byte consumption must then be forecast).
6. **Sizes**, all from `design.json`:
   - 2,000 healthy validation runs and 1,000 runs per family (600 fitting, 400 evaluation): about 236 minutes;
   - 2,000 and 500 (300/200): about 177 minutes, with criterion 4 power 0.47 to 0.70 at a true skill of 0.15;
   - 1,000 and 300 (200/100): about 114 minutes, with criterion 1 passing with probability 0.685 at a true coverage of 0.94.
7. **Branch:** v0.19 on a new branch `feature/v019` from `feature/v018` / continue on `feature/v018`.
