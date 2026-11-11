# Forecasting the remaining useful run under a degrading random input: what is possible, methods, and a validation study

**Status: 2026-11-11.** The code is in `src/rnglab/prognosis/` (tests: `tests/v019/test_prognosis.py`). The validation study is in `research/v019-forecast-methods/`: protocol frozen at `459d144` (SHA-256 `4b3b6a24…0dba`), run once, results at `a568973`. Every study number below is copied from `results/tables.md` or `results/verdicts.json` there.

**Scope.** The study uses a Gaussian surrogate workload. It says how the *methods* behave. It is not a test of H2 or H2-R, and it says nothing about the path-integral sampler or about any generator. The frozen Stage D pilot is unchanged: its forecasters and endpoints are not replaced by anything here.

## 1. What cannot be forecast without an assumption about the future

Let the history be everything observed up to issuance. Take two laws of the run ("worlds") under which the history has the same distribution but the future does not.

**Statement.** A forecaster that sees only the history gives the same forecast q in both worlds. Let p₁ and p₂ be the probabilities of the outcome given the full severity path in each world. The forecaster's Brier score exceeds that of the path-informed forecast by E[(q − p_w)²] in world w. Since (q − p₁)² + (q − p₂)² ≥ (p₁ − p₂)²/2 for every q,

max over the two worlds of the excess ≥ E[(p₁ − p₂)²] / 4.

No method avoids this; only information about which world holds does.

**Sizes found here:**
- **Phase II Bernoulli counterexample** (`research/phase2-20261030/rus/`): the failure probabilities inside the horizon are 0.002472 and about 1, so the bound is 0.2488. A constant forecast of 1/2 scores 0.25: the history is worth nothing there (`test_minimax_excess_brier_for_the_bernoulli_counterexample`).
- **Ramp that continues or freezes at the alarm** (study scenarios E): the bound is 0.0079, 0.0701 and 0.1268 at L = 10,240, 20,480 and 40,960 (Section 4.4).

**What follows.** These forecasts are not available from the history alone:
1. whether a fault seen at the alarm will grow, persist or stop;
2. a failure time or a validity probability that is calibrated across continuation laws that share the history;
3. any warning before the onset of an unannounced abrupt fault;
4. a failure time beyond the horizon: it is censored, and nothing here extrapolates it.

A forecast is therefore always "under declared family F". Inside a declared family, onset and rate are still only partly identified at issuance. The gap between the posterior forecaster and the oracle in Section 4.1 measures that.

## 2. Methods

| Module | Content |
|---|---|
| `families` | severity paths with exact cumulative sums: the Stage D ramp, step-and-hold, ramp-then-hold, sampled random walks; the declared prior and its grid |
| `validity` | probability that a 16-batch batch-means interval covers the reference under a drifting mean |
| `surrogate` | Gaussian workload with the Stage D geometry: alarm calibration, validity at declared lengths, failure time with censoring |
| `forecasters` | base rate, alarm-time logistic, lite logistic, persistence, declared-model posterior (`BayesRamp`), oracle, empirical failure-time baseline |
| `scoring` | Brier and logarithmic scores, reliability, Spiegelhalter's Z [Spi86], paired bootstrap, ranked probability score [Eps69], interval coverage and median error under censoring, Kaplan–Meier [KM58] |

**Validity model.** With independent normal batch means, the interval at length L covers the reference with probability

P(δ, λ) = E_V[Φ(c·√(V/15) − δ) − Φ(−c·√(V/15) − δ)], V ~ noncentral χ²(15, λ),

where δ = β·S̄·√L/σ_eff carries the bias and λ = (b/σ_eff²)·β²·Σ(s̄_j − S̄)² carries the spread of the batch means caused by drift. The statistic is doubly noncentral t ([Kri68], [BA68]). With λ = 0 it is the Stage D oracle formula. The code checks the quadrature against adaptive integration and against simulation.

**Posterior forecaster.** A grid over the declared onset and rate carries the prior. The likelihood is the joint Gaussian density of every statistic observed up to the alarm, with correlation 1/2 between adjacent overlapping windows. The forecast mixes the validity model over the posterior. Because the likelihood uses everything the stopping rule saw, issuing at the first threshold crossing needs no selection correction.

**Information barrier.** A forecaster receives a `History`: a read-only copy of the statistics up to the alarm. Tests rewrite the whole future of each run and require identical forecasts.

**Censoring.** A failure-time forecast has an explicit last cell, "not failed inside the horizon". Scores use only thresholds inside the horizon, where the outcome is observed for every run.

## 3. Validation study: design

Surrogate with β = 0.06, σ_eff = 0.613, κ = 103, horizon 40,960; alarm at the 950th of 999 healthy maxima. 999 calibration runs, 4,000 training runs, 2,000 evaluation runs per scenario, 10,000 bootstrap resamples. Scenarios: A declared ramp; B1 slower rates and B2 later onsets than declared; C step-and-hold; D random walk; E the declared ramp continued or frozen at each run's own alarm, with the same noise. Reading rules R1 to R7 were fixed before generation (`PROTOCOL.md` §7).

Estimated from calibration and training: σ_eff 0.6125, β 0.05984, κ 102.998, threshold 3.4038.

## 4. Results

| Rule | Outcome |
|---|---|
| R1 calibration when the declared family is true | **MET** |
| R2 skill over the alarm-time baseline at L = 20,480 | **MET**: 0.228 [0.189, 0.268] |
| R3 better than the lite features at L = 20,480 | **MET**: Brier difference −0.0132 [−0.0181, −0.0084] |
| R4 drift-widening term | matters at 20,480 and 40,960; not at 10,240 |
| R5 expectation: miscalibrated under unseen families | **MET** in C and in D, at every length |
| R6 identical histories | forecasts identical; bound 0.0079 / 0.0701 / 0.1268 |
| R7 failure-time skill over the empirical baseline | **MET**: 0.111 [0.092, 0.129] |

### 4.1 Declared family true (A)

| L | eligible | valid | alarm-time | lite | posterior | oracle | lite skill | posterior skill | oracle skill |
|---|---|---|---|---|---|---|---|---|---|
| 10,240 | 1,734 | 0.850 | 0.1120 | 0.1068 | 0.1034 | 0.0985 | 0.047 [0.026, 0.067] | 0.076 [0.044, 0.109] | 0.120 [0.084, 0.158] |
| 20,480 | 2,000 | 0.503 | 0.1812 | 0.1531 | 0.1399 | 0.1153 | 0.155 [0.123, 0.186] | 0.228 [0.189, 0.268] | 0.364 [0.317, 0.409] |
| 40,960 | 2,000 | 0.183 | 0.0950 | 0.0927 | 0.0799 | 0.0577 | 0.025 [0.001, 0.049] | 0.160 [0.121, 0.196] | 0.392 [0.337, 0.445] |

Columns 4 to 7 are Brier scores; skill is relative to alarm-time.

- The posterior forecaster and the oracle pass the calibration rule at every length (|Z| at most 1.27 and 2.14).
- **The lite features (the form of the pilot's F-B-lite) capture part of the available skill:** 0.155 against 0.228 at 20,480, and 0.025 against 0.160 at 40,960.
- **Persistence fails under a ramp:** skill −0.21, −1.36 and −5.85.
- **Without the widening term** the posterior forecaster is miscalibrated at 20,480 (Z = 3.55) and scores worse at 20,480 and 40,960 (differences −0.0033 [−0.0056, −0.0009] and −0.0014 [−0.0026, −0.0003]).
- **Failure time:** 271 of 2,000 runs censored, 241 already failed at issuance. Ranked probability score 0.1232 against 0.1385. The 80% interval covers 0.841 [0.824, 0.857], above nominal. Median absolute error 4,352 sweeps against 5,376.

### 4.2 Rates and onsets outside the declared range (B1, B2)

| Scenario | L | posterior skill | calibrated |
|---|---|---|---|
| B1 slow rates | 20,480 | 0.093 [0.050, 0.137] | yes |
| B1 slow rates | 40,960 | 0.111 [0.087, 0.133] | **no** (Z = −4.25) |
| B2 late onsets | 20,480 | 0.010 [−0.007, 0.027] | **no** (Z = 19.8) |
| B2 late onsets | 40,960 | 0.060 [0.050, 0.070] | **no** (Z = 56.7) |

- In B1 persistence comes close to the oracle (skill 0.338 against 0.346 at 40,960): a very slow ramp is nearly constant.
- In B2 the failure-time skill is **negative**: −0.113 [−0.126, −0.099].

### 4.3 Unseen families (C, D)

| Scenario | L | posterior skill | lite skill | persistence skill | posterior Z |
|---|---|---|---|---|---|
| C step-and-hold | 20,480 | **−0.256** [−0.299, −0.214] | −0.203 | 0.367 | 70.8 |
| C step-and-hold | 40,960 | −0.017 [−0.030, −0.003] | −0.016 | 0.511 | 116.7 |
| D random walk | 20,480 | **−0.218** [−0.252, −0.185] | −0.162 | 0.049 | 49.2 |
| D random walk | 40,960 | −0.069 [−0.096, −0.042] | −0.015 | −3.58 | 33.7 |

- Under a wrong family the posterior forecaster and the lite logistic are **worse than the alarm-time baseline**.
- Failure-time skill is −0.098 (C) and −0.096 (D). The nominal 80% interval covers 0.450 [0.427, 0.473] and 0.522 [0.500, 0.544].
- Persistence is the right family for C and does best there, but it is not calibrated either (Z between 9 and 20).

### 4.4 Identical histories, two futures (E)

| L | valid if the ramp continues | valid if it freezes | bound | posterior: excess if continues | excess if freezes | persistence: excess if continues | excess if freezes |
|---|---|---|---|---|---|---|---|
| 10,240 | 0.841 | 0.939 | 0.0079 | 0.0034 | 0.0279 | 0.0350 | 0.0001 |
| 20,480 | 0.521 | 0.919 | 0.0701 | 0.0245 | 0.2575 | 0.3092 | 0.0014 |
| 40,960 | 0.189 | 0.862 | 0.1268 | 0.0200 | 0.4860 | 0.5767 | 0.0097 |

Every forecast was identical in the two worlds, as it must be. Each forecaster is good in the world it assumes and poor in the other. None is good in both.

### 4.5 An unusual calibration draw

The healthy calibration set covered 930 of 999 runs at L = 10,240 (0.9309). The surrogate's law is exactly 0.95; a count this low has probability 0.005. It entered only the threshold and σ_eff (0.6125 against 0.613). It is reported because it happened.

## 5. Meaning for Stage D

- **Nothing in the frozen pilot changes.** These are observations for reading it and for a later campaign.
- A posterior forecaster over the whole statistic history is the campaign's planned "F-B". On the surrogate it adds skill over the lite features, most at the longest length.
- The oracle should use the doubly noncentral form. The frozen formula stays in the pilot as written (`IMPLEMENTATION_AUDIT.md`, F-01 and F-02).
- **Skill under the declared ramp is not evidence of skill under real degradation.** Under two other families the same forecaster was worse than a baseline that ignores the run's signals.

## 6. Proposed elements of a campaign forecasting protocol (not run, not frozen)

1. **Declared family and prior** as in the pilot. State in the claim that every score is conditional on it.
2. **Observation model per family**, fitted on training runs: mean and variance of the target statistic against window-mean severity, and the correlation of adjacent windows. For SHIFT use the saturated-word count (audit F-03).
3. **Bias and noise model checked first:** the batch-mean normality and independence that the validity model assumes, on healthy and fault training runs.
4. **Forecasters:** alarm-time baseline, lite logistic, posterior forecaster, persistence, oracle with the doubly noncentral form.
5. **Primary contrast:** posterior against alarm-time at L = 20,480 by paired bootstrap; calibration by Z and reliability.
6. **Robustness arms scored separately:** held-out rates and onsets, step-and-hold, random walk, and the continue-or-freeze pair with its bound.
7. **Failure time as a secondary endpoint** with administrative censoring, scored as in Section 2.
8. Separate seed namespaces for training, calibration and evaluation; forecasters hashed before evaluation data exist.

## 7. Limits

- The surrogate is Gaussian with independent blocks. The sampler's batch means may be neither; that is untested.
- The scales (β, σ_eff, κ) are a choice of the order seen in four audit probe runs, not a fit.
- The surrogate has no application monitor, so it says nothing about F-AB against F-B.
- One study, one seed set per scenario. The R1 thresholds were set on development seeds.
- References were checked against Crossref for author, title, venue and year (`docs/v019/evidence/refs_check.json`, `crossref_lookup.json`). Their full texts were not read.

## References

- [BA68] W. G. Bulgren, D. E. Amos. A note on representations of the doubly non-central t distribution. *J. Am. Stat. Assoc.* 63(323), 1013–1019, 1968. doi:10.1080/01621459.1968.11009328
- [Eps69] E. S. Epstein. A scoring system for probability forecasts of ranked categories. *J. Appl. Meteorol.* 8(6), 985–987, 1969. doi:10.1175/1520-0450(1969)008<0985:ASSFPF>2.0.CO;2
- [KM58] E. L. Kaplan, P. Meier. Nonparametric estimation from incomplete observations. *J. Am. Stat. Assoc.* 53(282), 457–481, 1958. doi:10.1080/01621459.1958.10501452
- [Kri68] M. Krishnan. Series representations of the doubly noncentral t-distribution. *J. Am. Stat. Assoc.* 63(323), 1004–1012, 1968. doi:10.1080/01621459.1968.11009327
- [Spi86] D. J. Spiegelhalter. Probabilistic prediction in patient management and clinical trials. *Stat. Med.* 5(5), 421–433, 1986. doi:10.1002/sim.4780050506
- Already in the programme: [GR07], [GBR07], [KB23], [GLLR05], [SWHZ11] (`docs/v018/STAGE_DE_REVIEW.md`, `research/v019-d-pilot/PROTOCOL.md`).
