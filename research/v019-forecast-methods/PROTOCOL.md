# v0.19 forecasting methods validation on a surrogate workload: protocol

**Status: FROZEN before generation, 2026-11-11.** `protocol.json` holds every parameter and reading rule; its SHA-256 is in `protocol.sha256` and `run_study.py` refuses to draw if they differ. Nothing has been drawn from the study namespace. A smoke run in its own namespace (`smoke/`, reduced sizes) exercised the code and is excluded from every result.

## 1. What this study is, and what it is not

- **It is** a check of forecasting *methods* (`src/rnglab/prognosis/`) against a workload whose laws are known exactly, so that calibration and skill can be judged against exact answers.
- **It is not** a test of H2 or H2-R, not evidence about the path-integral sampler, and not evidence about any random number generator. It is separate from the Stage D pilot (`research/v019-d-pilot/`), whose frozen forecasters and endpoints it does not touch.
- **Negative outcomes** are reported with the same prominence as positive ones. The study runs once.

## 2. Workload

A Gaussian surrogate with the Stage D pilot's geometry (`rnglab.prognosis.surrogate`):

- horizon 40,960 sweeps; checkpoints every 256 sweeps; monitoring windows of 512 sweeps; declared lengths 10,240, 20,480 and 40,960;
- the estimate at length L is a mean with a 16-batch batch-means interval; valid means the interval covers the reference;
- the observable has bias β·s(t) and independent Gaussian noise with long-run scale σ_eff;
- the monitoring statistic of a window is Gaussian with unit variance and mean κ times the window-mean severity; adjacent windows overlap by half and have correlation 1/2;
- β = 0.06, σ_eff = 0.613, κ = 103, s_max = 0.5. These are of the order seen for the LAG1 family in four audit probe runs and in Stage C; they are a choice, not a fit.

**Declared degradation family:** the Stage D ramp, s(t) = min(s_max, s_max·(t − t₀)/T_cap) from t₀, with t₀ integer-uniform on [2,048, 8,192] and T_cap log-uniform on [2,048, 163,840].

**Alarm:** first strict exceedance of the 950th order statistic of 999 healthy whole-horizon maxima, with no onset knowledge.

**Failure time:** the first checkpoint from which the interval excludes the reference at three consecutive checkpoints. It is right-censored at sweep 40,448 (administrative censoring).

## 3. Sets and seeds

| Set | Size | Use |
|---|---|---|
| calibration (healthy) | 999 | alarm threshold, σ_eff estimate |
| training (declared ramp) | 4,000 | β and κ estimates, logistic fits, empirical failure-time baseline |
| forward simulation | 100 per hypothesis | failure-time distribution of each ramp hypothesis |
| evaluation, per scenario | 2,000 | scoring only |

Seeds are `namespace_seed("rnglab-v019-forecast-methods", set, role)`. No set shares a seed with another. Forecasters are fixed before any evaluation set is drawn; every evaluation set is opened once.

## 4. Scenarios

| Scenario | Evaluation law | Question |
|---|---|---|
| A_declared | the declared ramp | calibration and skill when the declared family is true |
| B1_slow_rates | ramp, T_cap log-uniform on [163,840, 655,360] | rates outside the declared range |
| B2_late_onsets | ramp, t₀ uniform on [8,193, 16,384] | onsets outside the declared range |
| C_step_hold | step to a level log-uniform on [0.005, 0.5], then constant | an unseen family |
| D_random_walk | reflected random walk with drift from a uniform onset | an unseen family |
| E_continue, E_freeze | the declared ramp; the same paths and the same noise with the ramp frozen at each run's own alarm | identical histories, different futures |

## 5. Forecasters

Base rate; logistic regression on the alarm sweep (alarm-time); logistic regression on the alarm sweep, the statistic at the alarm and its slope over the last four checkpoints (lite-logistic, the form of the Stage D pilot's F-B-lite); persistence of the severity seen at the alarm; the declared-model posterior forecaster (bayes-ramp), with and without the drift-widening term of the validity model; and the validity model at the true path (oracle, a reference that uses evaluation truth). Failure time: bayes-ramp by forward simulation, against the empirical distribution by alarm-time quintile.

All parameters (β, σ_eff, κ, logistic coefficients) are estimated from the calibration and training sets. The exact reference is used in training only.

## 6. Scores

- Validity at each L, on runs whose alarm is below L: Brier score, logarithmic score, Brier skill over alarm-time with a paired-bootstrap 95% interval (10,000 resamples of runs), Spiegelhalter's Z, reliability in ten bins.
- Failure time, on runs with an alarm: ranked probability score over the in-horizon thresholds, coverage of the 80% central interval (decidable for censored runs), median absolute error where forecast median and outcome are both inside the horizon, the three other cases counted, and Kaplan–Meier medians.

## 7. Reading rules (fixed before generation)

| Rule | Statement |
|---|---|
| R1 | A: bayes-ramp and oracle are calibrated at every L (\|Z\| < 3 and weighted reliability gap < 0.05 over bins with at least 50 runs). |
| R2 | A, L = 20,480: the lower 95% bound of the skill of bayes-ramp over alarm-time is above 0. |
| R3 | A, L = 20,480: the 95% interval of Brier(bayes-ramp) − Brier(lite-logistic) lies below 0. |
| R4 | A, each L: whether the 95% interval of Brier(bayes-ramp) − Brier(no-widening) lies below 0. Reported per L. |
| R5 | B1, B2, C, D: R1's calibration rule and skill for every forecaster. Expectation stated in advance: in C and in D, bayes-ramp fails R1's rule at one L at least. |
| R6 | E: forecasts are identical in both worlds (exact equality, otherwise an implementation error). Reported: every forecaster's excess over the path oracle in each world, and the bound mean((p_continue − p_freeze)²)/4. |
| R7 | A, failure time: the lower 95% bound of the RPS skill of bayes-ramp over the empirical baseline is above 0. |

The R1 thresholds were set on development seeds before the freeze (`protocol.json`, R1). R6's bound is a theorem (`docs/v019/FORECASTING_METHODS.md`, section 1); the study reports its size.

## 8. Reproduction

```bash
PYTHONPATH=src .venv/bin/python research/v019-forecast-methods/run_study.py
.venv/bin/python research/v019-forecast-methods/report_tables.py
```

`run_study.py` writes `results/results.json` and refuses to overwrite it. `report_tables.py` writes `results/tables.md` and `results/verdicts.json` from that file alone.
