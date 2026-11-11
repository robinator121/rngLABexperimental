# v0.19 Stage D pilot: independent implementation audit

**Status: audit of the state on 2026-11-11, branch `feature/v019-audit` (from `33f701e`).** No pilot, smoke or campaign data were generated. Engineering runs used the audit-only namespace `rnglab-v019-audit-smoke`. No frozen file was changed. Numbers come from `research/v019-audit/results/audit_numbers.json` (script `audit_numbers.py`) and `results/benchmark.json` (script `bench.py`) unless another file is named.

## 1. State found

| Item | Found |
|---|---|
| Protocol | `research/v019-d-pilot/protocol.json` matches `protocol.sha256` (`53294642…5c97`), frozen at `b33846f` |
| Profile and adapter (C2) | committed at `33f701e` (`dlib.py`). Its 26 tests pass; with `tests/v018` that is 106, the figure reported at the checkpoint |
| Measurement layer | `measure.py` and its tests exist only as uncommitted files in the main checkout. Two snapshots were examined: SHA-256 `b7c25b67…` (2026-11-10 23:33) and `f8d5ac52…` (2026-11-11 00:05) |
| Not in the tree | campaign driver, calibration, issuance file, forecaster fitting with the evaluation barrier, scores, criteria, the pilot's own recount, smoke run, pre-generation manifest |
| TDK deadline | not recorded: stop criterion 5 is INSUFFICIENT_INFORMATION until it is |
| Baseline tests | 451 passed at `33f701e` before any change |

## 2. Findings

Type: **P** needs the author's decision under the amendment procedure (nothing was changed); **R** implementation repair that leaves frozen semantics unchanged; **I** information.

| ID | Type | Finding | Evidence |
|---|---|---|---|
| F-01 | P | **The frozen oracle call can return NaN.** `scipy.stats.nct.cdf` (SciPy 1.18.1) is not finite at 312 of 901 noncentralities on a 0.05 grid, from δ = 8.55 to 43.6. Where it is finite its error is below 10⁻¹³. The true probability at the first NaN is 1.2 × 10⁻⁹. Under `oracle.unavailable` a NaN makes that run's oracle probability unavailable, and one unavailable value makes the oracle score of that family and length unavailable. Stop criteria 1 to 4 do not use the oracle | `oracle_call`; test `test_frozen_oracle_call_is_nan_or_correct` |
| F-02 | I | **The oracle formula leaves out the widening of the interval by drift** (the protocol says so). Under normal batch means the exact probability is doubly noncentral (`src/rnglab/prognosis/validity.py`). In the surrogate study the formula without widening was miscalibrated at L = 20,480 (Z = 3.55) and scored worse at 20,480 and 40,960 | `research/v019-forecast-methods/results/tables.md`, R4 |
| F-03 | I | **The SHIFT target statistic is not what fires the alarm.** In the five windows of one full-length SHIFT run with mean severity below 0.002 the frequency z averaged 0.80, while those windows held 10.2 saturated words on average against none before onset. The saturated fraction equals the mean severity (0.996 per unit). F-B-lite's SHIFT features are therefore near noise at issuance, so a failed criterion 4 for SHIFT would not show that no signal exists | `full_length_runs.SHIFT.shift_saturation` (one run) |
| F-04 | I | **The 1,000,000-bit cap is correct but unreachable.** Windows held 849,344 to 876,928 bits in four full-length runs; a cut needs an uphill fraction above 0.907, against 0.635 to 0.660. The smoke will not exercise the cap. The cut path is covered by synthetic windows only. The longest-run test always uses 10,000-bit blocks | `full_length_runs`; `test_bit_cap_boundary_is_exactly_one_million_bits` |
| F-05 | I | **Storage is about 22.8 GiB**, not the 5.49 GiB baseline (which the protocol says covers traces only): 3.23 MB of run record and 0.84 MB of trace per run. The run records alone are 19.4 GB, too large to commit | `campaign_cost` |
| F-06 | I | **Wall time is about 8.8 hours at six workers** on this machine (7.3 with F-07), against the 3.93-hour extrapolation. Each job ran 2.38 times slower in parallel than alone. Eight workers were not measured | `campaign_cost`, `benchmark.json` |
| F-07 | R | The I-3 transformer takes 2.6 to 5.3 s per run; the independent vectorised transform takes 0.03 to 0.64 s and gave identical bytes in every test. `dlib.severity` is called once per word (10% of generation time). Neither was changed | `benchmark.json` |
| F-08 | I | **Bit-for-bit replay holds on one platform only.** The sampler squares scalars with libm `pow`, which differs from the rounded product in 912 of 10⁶ cases; `exp`, `log1p` and `cos` also come from libm. I-4 must run in the generating environment, and the manifest should record the libc version. `preflight.py` compares golden digests | `power_versus_product` |
| F-09 | I | The LAG1 "midpoint" u + 2⁻⁵⁴ is exact only below 1/2; above it the sum rounds. One word value (mantissa 2⁵³ − 1) is outside the inverse-CDF domain and gives a recorded run failure. Inherited from Stage C | `test_lag1_midpoint_rounding_and_domain_failure_set` |
| F-10 | R | A window of identical words gives r = 1 in the Stage C word test (rounding in `np.corrcoef`), where the signal schema says unavailable. The first measurement snapshot logged it as available; the second handles it. A regression test now guards it. Unreachable for the pilot's families | `test_measure_lag1_is_unavailable_for_a_constant_window` |
| F-11 | I | "Already invalid at issuance" uses 16 batches of n/16 sweeps, 32 sweeps at n = 512. For this sampler Phase 2 measured coverage 0.8611 at batch 32 and 0.9260 at 128, so its healthy base rate is above 5% for early alarms | `research/v018-api-validation/REPORT.md` §4 |
| F-12 | I | σ_eff averages standard errors, which is 0.9835 of the long-run value for 15 degrees of freedom; δ is 1.7% too large | `gates` |
| F-13 | I | Gate V2-D passes if at most 173 of 2,000 healthy runs alarm. Criterion 1 accepts 1,864 to 1,946, with pass probabilities 0.384, 0.938 and 0.9998, as `design.json` states. Exact cap shares are 0.1936, 0.4583 and 0.6529, within two standard errors of the simulated values the protocol quotes | `gates` |
| F-14 | R | `write_run` writes the trace and then the record without temporary files. A process killed between them leaves a trace without a record. `integrity.py` classifies that as PARTIAL and offers an atomic commit | `test_killed_writer_never_leaves_a_run_that_passes_as_complete` |
| F-15 | I | L-BFGS-B in SciPy 1.18.1 ends with "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL" or "…RELATIVE REDUCTION OF F…"; the iteration limit gives `success=False`, status 1. The frozen call and a Newton optimum differed by less than 2 × 10⁻⁴ in probability, so the recount tolerance of 10⁻¹² applies only with saved parameters | `test_newton_optimum_and_frozen_lbfgsb_agree` |

## 3. What was verified

- **Whole-run regeneration.** `replay.py` rebuilds a run from its seed labels with its own sweep loop and NumPy's PCG64. On eight full-length runs (two per family) marks, CV, TD, severity and all stream hashes were bit-identical to the pilot driver's output.
- **Stored outputs.** `verify.py` made 6,818 comparisons per set of four stored runs (both snapshots) with no discrepancy. Every tampering in the tests was reported.
- **Stage C regression.** The independent functions reproduce the stored threshold 3.8177723774232715, 199/4096 healthy alarms and 3867/4096 coverage.
- **Interruption.** A writer killed with SIGKILL at arbitrary moments never left a run classified as complete with wrong content.

## 4. Traceability

Status: **V** verified independently; **S** verified on the uncommitted snapshot only; **N** not implemented in the pilot (an independent version exists in `research/v019-audit/`).

| Requirement (`protocol.json`) | Implementation | Independent check | Status |
|---|---|---|---|
| application, sweep_clock | `dlib.generate_run`, `SweepClock` | `replay.py`, bitwise | V |
| profile | `dlib.draw_profile`, `severity` | `recount.draw_profile`, `severity_path` | V |
| fault_transforms, step equivalence | `dlib.FaultyWords`, `reapply` | `recount.apply_fault`; existing step tests | V |
| data.seeds, no pairing | `dlib.seeds_for` | 30,445 seeds disjoint | V |
| endpoint, issuance length | `measure.endpoint` | `recount.endpoint`, exact sums | S |
| battery, bit cap, statistic | `measure.battery_checkpoint` | `recount.battery_checkpoint`, NIST examples | S |
| first alarm, onset class, eligibility | `measure.issuance` and others | `recount.first_alarm` and others | S |
| application signals | `measure.application_checkpoint` | `verify.verify_run` | S |
| storage per run | `measure.write_run`, `load_run` | `verify.load_stored`, `integrity.py` | S, F-14 |
| thresholds, V2-D, constants, σ_eff | none | `verify.recount_campaign` | N |
| oracle | none | `recount`, F-01, F-02, F-12 | N |
| forecasters, scores, criteria 1 to 4 | none | `verify.fit_and_predict`, `recount` | N |
| run manifest, failed-run rules | none | `integrity.build_manifest` | N |
| I-1 to I-3 | `dlib.check_i1/2/3` | `verify.verify_run` | V |
| I-4, recount | none | `verify.deep_sample`, this directory | N |
| smoke, pre-generation manifest | none | `preflight.py` reports NOT READY | N |

## 5. Reproduction

Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1, glibc 2.44.

```bash
PYTHONPATH=src .venv/bin/python -m pytest tests/v019 -q
PYTHONPATH=src .venv/bin/python research/v019-audit/audit_numbers.py
PYTHONPATH=src .venv/bin/python research/v019-audit/preflight.py
PYTHONPATH=src .venv/bin/python research/v019-audit/verify.py --root <output root> --namespace <namespace>
```

Tests that need the measurement layer or the Stage C traces skip unless `RNGLAB_DPILOT_MEASURE` or `RNGLAB_STAGEC_TRACES` points to them.

## 6. Not verified

- The second measurement snapshot was not re-read line by line; it was compared through its outputs.
- F-03 rests on one SHIFT run. F-06 was measured once, at six workers.
- The whole suite passes on this branch: 561 passed and 5 skipped as committed; with the measurement layer and the Stage C traces supplied, the skipped tests run (`docs/v019/evidence/pytest_v019_audit.txt`). The measurement layer itself is not on this branch.
- Nothing downstream of the per-run records exists in the pilot, so none of it could be audited.
