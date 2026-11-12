# Stage D pilot: implementation status

**Status: 2026-11-12, branch `feature/v020-stage-d-readiness`.** The execution path from the pre-generation manifest to `criteria.json` is implemented and tested on synthetic fixtures. **Nothing has been generated:** no pilot run, no smoke run, no campaign result. The frozen protocol (`research/v019-d-pilot/protocol.json`, SHA-256 `5329464210019e968116e8c3c3a13a09f33b12e0368ecfaa69752a1eb43c5e97`) is unchanged, and no amendment has been applied.

The pilot cannot be started from this branch as it is. Section 3 lists what is missing and which decisions only the author can take.

## 1. What is implemented

All primary code is in `research/v019-d-pilot/`. The independent verification is in `research/v019-audit/` and never imports the primary code.

| Protocol component | Where | Tests | State |
|---|---|---|---|
| Sampler adapter, ramp profile, fault transforms, I-2, I-3 (`dlib.py`) | already in v0.19 | `test_d_pilot_lib.py`, `test_audit_d_pilot.py` | Done before this work; unchanged |
| Measurement layer: checkpoints, endpoints, issuance rule, run format (`measure.py`) | **another developer's uncommitted file; not in this branch** | its own tests, and the tests here when the file is supplied | Used read-only. See B-1 |
| Run plan, canonical JSON, atomic write, context (`pilot_common.py`) | new | `test_pilot_driver.py` | Done |
| Generation driver: one record per planned run, failure records, staging and rename with the record last, journal, shards, resume (`pilot_generate.py`) | new | `test_pilot_driver.py` | Done |
| Pre-generation manifest: protocol and design hashes, every executed file's hash, versions and libc, planned runs and seed labels, I-4 selection, commit, workers, policy, deadline, approved amendments | `pilot_generate.build_pre_generation_manifest` | `test_pre_generation_manifest_contents` | Done. Pushing it is the author's step |
| Authorisation gate: no draw in the pilot or smoke namespace without the manifest, unchanged code hashes, a named policy, the author's `authorization.json`, and thread variables at 1 | `pilot_generate.require_generation_allowed` | `test_generation_is_refused_without_the_authors_authorization` | Done. The authorisation file is an engineering addition (B-9) |
| Step check of I-3 for replica 0 of each fault family | `pilot_generate.step_check` | `test_step_check_agrees_with_stage_c_and_detects_a_wrong_stream` | Done |
| I-4 replay of the 25 selected runs | `pilot_generate.replay_check` | `test_replay_detects_a_changed_value` | Done |
| Run manifest and failure list | `pilot_generate.finalize` | `test_pilot_driver.py` | Done |
| Calibration: 999 ordered maxima, full battery and subset, thresholds, ties, application constants, σ_eff | `pilot_analysis.calibrate` | `test_threshold_is_the_prescribed_order_statistic`, `test_incomplete_calibration_…` | Done |
| Issuance, eligibility, onset class, issuance endpoint | `pilot_analysis.issue` | `test_issuance_rows_respect_the_boundaries` | Done |
| Forecasters F0, F0+, F-B-lite, F-AB-lite; the frozen optimiser call and convergence test; fallbacks; standardisation | `pilot_forecast.fit` | `test_the_fit_is_the_frozen_call_…`, `test_fallbacks_follow_the_protocol`, `test_prediction_rules` | Done |
| Evaluation barrier: fitting opens calibration and fitting runs only; model hash before evaluation; evaluation once | `pilot_forecast.fit`, `evaluate` | `test_fit_opens_fitting_runs_only_…`, `test_fit_refuses_a_run_that_is_not_a_fitting_run` | Done |
| Oracle: S̄, β, σ_eff, probability by the frozen call; unavailable values kept unavailable | `pilot_forecast.oracle_probability` | `test_oracle_rule_is_frozen_…`, `test_f01_oracle.py` | Done as frozen. See B-3 |
| Scores: Brier, skill, calibration bins, paired bootstrap with the saved seeds, contrasts | `pilot_forecast.evaluate` | `test_scores_equal_a_direct_recomputation`, `test_a_missing_forecast_…` | Done |
| Healthy and fault summaries, V2-D, short-interval check | `pilot_analysis.healthy_summary`, `fault_summary` | `test_pilot_pipeline.py` | Done |
| Implementation checks file (I-1 to I-4) | `pilot_analysis.implementation_checks` | `test_pilot_pipeline.py` | Done |
| Stop criteria 1 to 5 with precedence | `pilot_analysis.criteria` | `test_criterion_1_…`, `test_criteria_2_and_3_…`, `test_criterion_4_and_5`, `test_precedence_rules_…` | Done |
| Command line | `run_pilot.py` | `test_command_line_stops_without_the_measurement_layer` | Done |
| Independent recount and comparison with every result file | `research/v019-audit/verify.py --compare-primary`, `recount.py` | `test_audit_verify.py`, `test_audit_recount.py`, `test_independent_recount_agrees_…` | Done for criteria 1 to 4. See Section 2 |
| Store integrity scan and preflight | `research/v019-audit/integrity.py`, `preflight.py` | `test_audit_integrity.py` | Done |

## 2. What is not implemented or not done

- **`results/report_numbers.json` and `RESULTS.md`.** They hold the author's report. The pipeline stops at `criteria.json` (B-8).
- **`results/recount.json` is not written by the pipeline.** It is the output of the independent script (`verify.py --out`). The script records the `--root` argument as given: pass a relative path, so that no local path enters a committed file.
- **Criterion 5 is not recounted.** The protocol asks that its dates be "independently checked against manifest dates and the declared deadline". The recount reports criterion 5 as INSUFFICIENT_INFORMATION and compares criteria 1 to 4 only.
- **The smoke has not been run,** so `smoke/smoke_timing.json` does not exist. The command exists (`run_pilot.py smoke-timing --smoke`).
- **The pipeline has never run on real pilot data.** The analysis stages were exercised on 207 synthetic runs in the audit-only namespace; the generation path was timed on full-length runs in that namespace and the runs were deleted (`STAGE_D_READINESS_REPORT.md`, Section 3).
- **No amendment is applied.** `amendments.json` does not exist, so the oracle follows the frozen rule.

## 3. Blocking items and decisions for the author

The pipeline has no default for any of these. Where a value is needed and missing, the command stops and says which.

| ID | Item | Options or what is needed | What the pipeline does until then |
|---|---|---|---|
| B-1 | **The measurement layer is not in the repository.** `research/v019-d-pilot/measure.py` is an uncommitted file of another developer | Its author commits it. Its hash then enters the pre-generation manifest | `run_pilot.py` stops: "measure.py is not in this tree". Preflight reports NOT READY |
| B-2 | **Policy for a run that was started and left no record** (killed process, power cut). The protocol is silent | `regenerate_same_id`: run the same identity again; same seed labels, same bytes. `record_failure`: write a failure record and never run it. Consequence of the second: an interrupted calibration run makes the threshold unavailable (protocol, `failed_runs.calibration`) | `manifest` refuses to write without `--interrupted-run-policy`; `generate` refuses to pass an interrupted run |
| B-3 | **Oracle rule (F-01).** The frozen call returns NaN for part of the noncentralities between 8.53 and 43.61, and which ones depends on the last bits of δ | Proposed Amendment 01, option A (replace only a non-finite lower term at δ > 0 by 0), option B (a stable formula everywhere), or option C (no amendment, with a declared exception in the recount for oracle availability). To be decided before `fit` | Frozen rule. Option A exists behind a switch that is off unless `amendments.json` records an approved amendment `A01` with the hash of its document. The rule in force is written into `models.json` |
| B-4 | **Criterion 5 dates** | TDK deadline (`manifest --tdk-deadline`), date the engineering was complete, date the results were available (`criteria --engineering-complete-on --results-available-on`) | Criterion 5 is INSUFFICIENT_INFORMATION |
| B-5 | **Number of worker processes.** The protocol says 8. The v0.19 timing used 6 | 8, or another number recorded as a deviation. Throughput at 6, 8 and 12 workers is in the readiness report | Default 8. Another value sets `workers_deviate_from_protocol` in the manifest |
| B-6 | **Where the run records and traces live.** About 23 GiB, too large for Git | A directory with enough free space, and a decision on what is archived and where | The driver writes under `research/v019-d-pilot/`; `results/runs/` and `traces/` must stay outside Git |
| B-7 | **Smoke directory layout.** The protocol names `smoke/traces/` and `smoke/runs/`; the measurement layer writes `smoke/traces/` and `smoke/results/runs/` | Accept the layout of the measurement layer, or change that file | The driver follows the measurement layer |
| B-8 | **Report files** (`report_numbers.json`, `RESULTS.md`) | Written by the author after the criteria | Not produced |
| B-9 | **`authorization.json`.** An interlock added here: generation needs a file in which the author names the SHA-256 of the manifest. Not in the protocol | Keep it, or remove the check | Generation is refused without it |
| B-10 | **Supplementary descriptive analyses** raised by the red-team review (state baseline, flexible alarm-time baseline, scores split by state at issuance) | Amend before generation, or treat as exploratory afterwards | Not implemented |
| B-11 | **A storage error during generation** (disk full, input/output error). The protocol lists serialization exceptions among the failures that get a record, and a run with a record is never run again | Keep this reading, or declare that an error of the storage medium stops the driver and is handled as an interruption under B-2 | The run gets a failure record with phase `serialization`. Preflight requires about 31 GiB free to keep this from happening |

## 4. What the pipeline refuses

Invalid input is reported and stops the stage. Nothing is repaired, clipped or filled in.

- A stored run that does not verify (changed record, changed trace, record without trace, wrong identity, wrong namespace): generation and every analysis stage stop; the files are left as found.
- A trace without a record: moved to `interrupted/`, never deleted, and only under the author's policy.
- An existing record or trace is never replaced. A run with a record, including a failure record, is never run again.
- A non-finite trace: refused by the storage layer and recorded as a failure.
- A sampler or measurement exception: a failure record with the phase, the exception and the progress made.
- Incomplete calibration (fewer than 999 usable maxima): no threshold; issuance UNKNOWN; criteria 2 to 4 INSUFFICIENT_INFORMATION or UNINTERPRETABLE.
- A non-finite or out-of-range oracle probability: unavailable with a reason; the cell's primary score is unavailable.
- A missing forecast: the model's primary score is unavailable; criterion 4 is INSUFFICIENT_INFORMATION for that family.
- A model file whose hash differs from the model manifest, or evaluation outputs that already exist: evaluation stops.
- A smoke plan: no threshold, model or criterion is computed from it.
- A pilot or smoke namespace in a test plan, or a test engine with a pilot plan: refused.
- An amendment without a recorded approval, or whose document is missing or changed after approval: the pipeline stops.

## 5. Separation of implementation and verification

- `research/v019-audit/recount.py`, `replay.py` and `verify.py` use the standard library, NumPy and SciPy. They do not import `research/v019-d-pilot/` or `src/rnglab`. They read the stored files and the saved model parameters and recompute endpoints, thresholds, issuance, features, predictions, scores, bootstrap intervals and criteria 1 to 4. They do not refit.
- The comparison requires exact equality of flags, counts, labels and null patterns, and `np.isclose(rtol=1e-12, atol=1e-12)` for finite floats, as the protocol prescribes.
- **Limit.** Both code paths were written by the same hand in the same period. They share no code, but they can share a misreading of the protocol. Agreement between them is evidence of consistent arithmetic, not of a correct reading. A reading of the protocol by a second person remains worthwhile.

## 6. Tests

Evidence: `docs/v020/evidence/pytest_v020.txt` (commit `425dc61`, Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1).

| Run | Passed | Failed | Skipped |
|---|---|---|---|
| Whole suite, tree as committed | 614 | 0 | 31 |
| `tests/v019` with the measurement layer and the Stage C traces supplied | 220 | 0 | 0 |

- The 31 skips are tests that need the measurement layer (30) or the Stage C traces (1). Neither is in this branch. With both supplied, every test of `tests/v019` runs and passes.
- Before this work the same two runs gave 561 passed with 5 skipped, and 141 passed (`docs/v019/evidence/pytest_v019_audit.txt`).
- New test files: `test_f01_oracle.py` (13 tests), `test_pilot_driver.py` (29), `test_pilot_pipeline.py` (36); extended: `test_audit_verify.py`, `test_audit_integrity.py`.
- The measurement layer used had SHA-256 `f8d5ac52…1b9be`. It is an uncommitted file and can change; the counts hold for that version.

## 7. Commands

Order of a pilot, once the items of Section 3 are settled (nothing below has been run for the pilot or smoke namespace):

    PYTHONPATH=src python research/v019-audit/preflight.py --workers 8
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py status
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py manifest --interrupted-run-policy <policy> --tdk-deadline <YYYY-MM-DD> --workers <n>
    # the author reviews and pushes pre_generation_manifest.json, then writes authorization.json
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py generate
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py replay
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py calibrate
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py fit
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py issue
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py evaluate
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py summaries
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py checks
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py criteria --engineering-complete-on <YYYY-MM-DD> --results-available-on <YYYY-MM-DD>
    python research/v019-audit/verify.py --root research/v019-d-pilot --namespace rnglab-v019-d-pilot --deep sample --compare-primary --out research/v019-d-pilot/results/recount.json

The smoke uses the same `status`, `manifest`, `generate` and `replay` with `--smoke`, then `smoke-timing --smoke`.

Verification without any pilot data:

    python -m pytest tests/v019 -q -p no:cacheprovider
    RNGLAB_DPILOT_MEASURE=<path to measure.py> python -m pytest tests/v019 -q -p no:cacheprovider
    PYTHONPATH=src RNGLAB_DPILOT_MEASURE=<path to measure.py> python tests/v019/synthetic_agreement.py
    PYTHONPATH=src python research/v019-audit/f01_oracle.py
