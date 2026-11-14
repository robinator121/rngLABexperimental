# Stage D pilot: implementation status

**Status: 2026-11-14, branch `feature/v020-stage-d-readiness`.** The execution path from the pre-generation manifest to `criteria.json` is implemented and tested on synthetic fixtures. **Nothing has been generated:** no pilot run, no smoke run, no campaign result. The frozen protocol (`research/v019-d-pilot/protocol.json`, SHA-256 `5329464210019e968116e8c3c3a13a09f33b12e0368ecfaa69752a1eb43c5e97`) is unchanged, and no amendment has been applied.

Changes of 2026-11-13: the author's decision on interrupted runs (`regenerate_same_id`, Section 3a) is implemented on both sides; the measurement layer is in the branch (B-1); preflight also checks the smoke namespace.

Changes of 2026-11-14: the author confirmed the readings R-1 to R-5 with two clarifications (Section 3a) and decided B-12; no test samples in the pilot or smoke namespace any more; the whole test suite runs without deselection.

Neither the smoke nor the pilot can be started yet. Section 3 lists what is missing and which decisions only the author can take.

## 1. What is implemented

All primary code is in `research/v019-d-pilot/`. The independent verification is in `research/v019-audit/` and never imports the primary code.

| Protocol component | Where | Tests | State |
|---|---|---|---|
| Sampler adapter, ramp profile, fault transforms, I-2, I-3 (`dlib.py`) | already in v0.19 | `test_d_pilot_lib.py`, `test_audit_d_pilot.py` | Done before this work; unchanged |
| Measurement layer: checkpoints, endpoints, issuance rule, run format (`measure.py`) | its developer's file, copied into this branch byte for byte (SHA-256 `f8d5ac52…1b9be`) | `test_d_pilot_measure.py` (its own 23 tests), and the pipeline tests | In the branch; `measure.py` unchanged. Its tests use the audit-only namespace since 2026-11-14 (B-12) |
| Run plan, canonical JSON, atomic write, context (`pilot_common.py`) | new | `test_pilot_driver.py` | Done |
| Generation driver: one record per planned run, failure records, staging and rename with the record last, shards, resume (`pilot_generate.py`) | new | `test_pilot_driver.py` | Done |
| Interrupted runs: same identity run again, at most three retries, then an explicit unresolved status; append-only chained audit log; one `generate` per output tree; stop on request without losing an attempt | `pilot_generate.prepare`, `run_one`, `journal`, `generate` | Section 3a | Done |
| Pre-generation manifest: protocol and design hashes, every executed file's hash, versions and libc, planned runs with seed labels, seeds and the hash of their random input, I-4 selection, commit, workers, the interrupted-run rules, deadline, approved amendments | `pilot_generate.build_pre_generation_manifest`; `run_pilot.py manifest [--preview]` | `test_pre_generation_manifest_contents`, `test_manifest_preview_…` | Done. Pushing it is the author's step |
| Authorisation gate: no draw in the pilot or smoke namespace without the manifest, unchanged code hashes, the decided interrupted-run rules, the author's `authorization.json`, and thread variables at 1 | `pilot_generate.require_generation_allowed` | `test_generation_is_refused_without_the_authors_authorization` | Done. The authorisation file is an engineering addition (B-9) |
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
| Command line | `run_pilot.py` | `test_manifest_preview_shows_the_decision_and_writes_nothing` | Done |
| Independent recount and comparison with every result file | `research/v019-audit/verify.py --compare-primary`, `recount.py` | `test_audit_verify.py`, `test_audit_recount.py`, `test_independent_recount_agrees_…`, `test_criterion_5_…` | Done. Criterion 5: dates only, see Section 2 |
| Independent check of the audit log: run identity, completion status, retry history | `research/v019-audit/verify.py` (`check_attempt_log`, part of every run of the script), `integrity.read_journal`, `recount.rng_input_sha256` | `test_the_independent_check_of_the_attempt_log_finds_what_was_changed`, `test_the_attempt_log_of_the_campaign_is_checked_independently` | Done |
| Store integrity scan and preflight, for the pilot and for the smoke | `research/v019-audit/integrity.py`, `preflight.py [--namespace rnglab-v019-d-pilot-smoke]` | `test_audit_integrity.py` | Done |

## 2. What is not implemented or not done

- **`results/report_numbers.json` and `RESULTS.md`.** They hold the author's report. The pipeline stops at `criteria.json` (B-8).
- **`results/recount.json` is not written by the pipeline.** It is the output of the independent script (`verify.py --out`). The script records the `--root` argument as given: pass a relative path, so that no local path enters a committed file.
- **Criterion 5 is checked for consistency, not for truth.** The protocol asks that its dates be "independently checked against manifest dates and the declared deadline". With `--compare-primary`, `verify.py` now recomputes the label and both deadlines from the frozen protocol and the manifest, and compares them with `criteria.json`. It also reports a date that is not a calendar date, an engineering date later than the manifest, and a results date earlier than the manifest; the last two follow from the order of the freeze steps and do not change the label. The two evidence dates themselves are the author's statement. Nothing in the repository can confirm them.
- **The smoke has not been run,** so `smoke/smoke_timing.json` does not exist. The command exists (`run_pilot.py smoke-timing --smoke`).
- **The pipeline has never run on real pilot data.** The analysis stages were exercised on 207 synthetic runs in the audit-only namespace; the generation path was timed on full-length runs in that namespace and the runs were deleted (`STAGE_D_READINESS_REPORT.md`, Section 3).
- **No amendment is applied.** `amendments.json` does not exist, so the oracle follows the frozen rule.
- **Generation with several worker processes has not been run under test.** The tests run the same functions (`prepare`, `work`, `run_one`, `finalize`) in one process, and kill and resume that process. The command that starts eight processes exists and was used for the timing runs before the audit log was added; the smoke will be its first run with the log.

## 3. Blocking items and decisions for the author

Where a value is needed and missing, the command stops and says which.

| ID | Item | Options or what is needed | What the pipeline does until then |
|---|---|---|---|
| B-1 | **The measurement layer.** Settled in this branch on 2026-11-13: `measure.py` and `tests/v019/test_d_pilot_measure.py` were copied byte for byte from the main checkout, where they are still untracked. The test file was changed on 2026-11-14 (B-12); `measure.py` is as copied | The author confirms that this is the version to freeze. The untracked copies in the main checkout will stand in the way of a later merge there until they are removed or committed | The file's hash enters the pre-generation manifest. A change of the file after the manifest stops `generate` |
| B-2 | **Policy for a run that was started and left no record.** Decided by the author on 2026-11-13: `regenerate_same_id`, at most three retries | Section 3a lists the implementation and the five readings, which the author confirmed on 2026-11-14 | Implemented. No other policy exists in the code |
| B-3 | **Oracle rule (F-01).** The frozen call returns NaN for part of the noncentralities between 8.53 and 43.61, and which ones depends on the last bits of δ | Proposed Amendment 01, option A (replace only a non-finite lower term at δ > 0 by 0), option B (a stable formula everywhere), or option C (no amendment, with a declared exception in the recount for oracle availability). To be decided before `fit` | Frozen rule. Option A exists behind a switch that is off unless `amendments.json` records an approved amendment `A01` with the hash of its document. The rule in force is written into `models.json` |
| B-4 | **Criterion 5 dates** | TDK deadline (`manifest --tdk-deadline`), date the engineering was complete, date the results were available (`criteria --engineering-complete-on --results-available-on`) | Criterion 5 is INSUFFICIENT_INFORMATION |
| B-5 | **Number of worker processes.** The protocol says 8. The v0.19 timing used 6 | 8, or another number recorded as a deviation. Throughput at 6, 8 and 12 workers is in the readiness report | Default 8. Another value sets `workers_deviate_from_protocol` in the manifest |
| B-6 | **Where the run records and traces live.** About 23 GiB, too large for Git | A directory with enough free space, and a decision on what is archived and where | The driver writes under `research/v019-d-pilot/`; `results/runs/` and `traces/` must stay outside Git |
| B-7 | **Smoke directory layout.** The protocol names `smoke/traces/` and `smoke/runs/`; the measurement layer writes `smoke/traces/` and `smoke/results/runs/` | Accept the layout of the measurement layer, or change that file | The driver follows the measurement layer |
| B-8 | **Report files** (`report_numbers.json`, `RESULTS.md`) | Written by the author after the criteria | Not produced |
| B-9 | **`authorization.json`.** An interlock added here: generation needs a file in which the author names the SHA-256 of the manifest. Not in the protocol | Keep it, or remove the check | Generation is refused without it |
| B-10 | **Supplementary descriptive analyses** raised by the red-team review (state baseline, flexible alarm-time baseline, scores split by state at issuance) | Amend before generation, or treat as exploratory afterwards | Not implemented |
| B-11 | **A storage error during generation** (disk full, input/output error). Settled by the decision of 2026-11-13, requirement 5, as read here: an error of the operating system is an interruption of the infrastructure, not a failure of the run | Confirmed by the author on 2026-11-14 (Section 3a, reading R-1). Before the decision such a run got a failure record with phase `serialization` | No record is written; the attempt is logged as interrupted and the worker stops. An encoding failure that is not an operating-system error (a non-finite array) still gets its failure record |
| B-12 | **Tests that made sampler runs with seeds of the smoke namespace.** Decided by the author on 2026-11-14: tests may not sample in the pilot or smoke namespace before the manifest. Two test files did. `test_d_pilot_measure.py` made one full-length run of `fault/SHIFT` replica 0 with a supplied profile. `test_d_pilot_lib.py`, in the repository since v0.19, made 64-sweep runs of `healthy` replica 0, of replica 0 of each fault family and of `fault/DUP` replica 1 with hand-built profiles, and one run of `fault/SHIFT` replica 0 with its seeded profile up to 16 sweeps past the onset | Both files now use the audit-only namespace `rnglab-v019-audit-smoke`; one check in `test_audit_verify.py` that drew a smoke profile uses another label. No test samples in the pilot or smoke namespace. What remains is seed derivation, the prescribed I-4 selection, and the generator fingerprint of the manifest (R-5) | **The short runs of `test_d_pilot_lib.py` were made by every run of the test suite up to 2026-11-13,** including those reported in the earlier evidence files; the statement of `pytest_v020_policy.txt` that no smoke-seeded run was made is wrong in this respect. Nothing was stored, and the smoke is excluded from every result |

## 3a. Interrupted runs: the decision of 2026-11-13

The frozen protocol says what happens to a run that failed (`failed_runs`), and is silent on a run that was started and left no record because its process ended. The author decided `regenerate_same_id` with eight requirements. The decision changes no seed, sample size, threshold, criterion or label of the protocol. It is frozen into the pre-generation manifest as `interrupted_run_policy` and `interrupted_run_rules`; `generate` refuses a manifest that carries anything else.

| # | Requirement | Implementation | Tests (`tests/v019`) |
|---|---|---|---|
| 1 | An interrupted run without a committed result is run again with its run ID, seed labels, generator configuration and protocol parameters | `prepare` puts the same identity back on the work list with the next attempt number. `run_one` derives the seeds from the same labels and refuses to draw if the seeds or the hash of the random input differ from the manifest. Code and protocol hashes are compared with the manifest at every `generate` | `test_a_run_that_was_started_and_left_no_record_is_run_again_with_the_same_identity`, `test_seeds_and_random_input_are_checked_against_the_manifest_before_anything_is_drawn`, `test_killed_generation_resumes_to_the_same_store` |
| 2 | Every attempt and every interruption reason in an append-only audit log | `results/journal.jsonl`: opened for appending, locked for each write, flushed to disk. Each entry has a sequence number and the SHA-256 of the entry before it. Entries: `generate_begin`, `start`, `commit`, `interrupted` (with the reason and who detected it), `restart`, `unresolved`, `trace_set_aside`, `stop_requested`, `worker_exit`, `generate_end`. A log that does not verify stops the driver | `test_the_audit_log_is_append_only_and_chained`, `test_a_write_cut_short_is_left_in_place_and_declared` |
| 3 | Atomic result commits; a partial output is never a result | Unchanged: files are staged, flushed and renamed, the record last. A trace without a record is moved to `interrupted/` and logged. A trace without a record and without a logged attempt stops the driver, because the driver cannot have left it | `test_killed_generation_resumes_to_the_same_store`, `test_a_trace_without_a_record_is_set_aside_not_deleted`, `test_a_record_stored_before_the_log_said_so_is_reconciled_not_rerun` |
| 4 | At most three infrastructure retries; exhaustion gives an explicit unresolved status and keeps the affected criteria from passing | Four attempts at most. After the fourth interruption the run gets a failure record with phase `infrastructure` and `infrastructure_status: UNRESOLVED_INFRASTRUCTURE`, listing its attempts, and is never run again. `criteria.json` names the runs and the criteria they touch; a `PASS` of such a criterion becomes `INSUFFICIENT_INFORMATION` | `test_three_retries_then_the_run_is_unresolved_and_never_run_again`, `test_an_unresolved_run_keeps_the_criteria_that_use_it_from_passing` |
| 5 | Infrastructure interruptions are kept apart from numerical, statistical and application failures; the latter stay as evidence | Infrastructure: the process ended without a record; an operating-system error or a memory error anywhere in the chain of the exception; a memory error that the measurement layer recorded inside a checkpoint. Everything else is the run's failure: it gets a failure record and is never retried | `test_an_error_of_the_system_leaves_no_record_and_is_retried`, `test_the_runs_own_failure_is_evidence_and_is_never_retried`, `test_a_memory_error_recorded_inside_the_measurement_is_an_interruption` |
| 6 | No seed, parameter, threshold or criterion is changed to obtain output | No code path takes a seed or parameter from anywhere but the frozen labels and the protocol. The manifest comparison of requirement 1 and the independent recomputation of requirement 8 would report it | as 1 and 8 |
| 7 | Deterministic generator input is verified with hashes; no bit-identical floating-point output is assumed across environments | `rng_input_sha256` per planned run: per generator, the SHA-256 of its first 128 64-bit outputs; integer arithmetic only. It is in the manifest, is checked before every attempt, and is in every `start` entry. Every attempt also logs the hash of its software versions. The I-4 replay stays a same-platform comparison | `test_the_random_input_hash_is_reproduced_by_the_independent_generators` |
| 8 | The independent recount verifies run identity, completion status and retry history | `verify.py` (`check_attempt_log`): the chain of the log; seeds and random-input hash recomputed with its own generators against the manifest, every attempt and the stored record; attempt numbers without gaps, four at most; every attempt but the last logged as interrupted for an infrastructure reason before the next starts; no attempt after a record exists; exactly one commit entry naming the hash of the record file; the unresolved runs and blocked criteria of `criteria.json`; the same history in `run_manifest.json`. Preflight recomputes every planned run of the manifest | `test_the_independent_check_of_the_attempt_log_finds_what_was_changed`, `test_the_attempt_log_of_the_campaign_is_checked_independently`, `test_preflight_manifest_and_store_rules` |

**Stopping on purpose.** The first Ctrl-C or SIGTERM lets every process finish the run in progress and then stops; no attempt is consumed. Only one `generate` can work on an output tree (`.generate.lock`).

**Readings applied here.** Each follows from the eight requirements as read, and each could have been read otherwise. **The author confirmed R-1 to R-5 as implemented on 2026-11-14,** with two clarifications, which are part of the rules that the manifest freezes (`unresolved_status`, `rng_input_use`):

- *R-4.* `UNRESOLVED_INFRASTRUCTURE` describes an unresolved execution status, not a proven infrastructure cause. The full attempt history and the original error evidence are preserved. Implemented: an `interrupted` entry written by the attempt holds every exception of the chain with its message and the frames in which the system error was raised; the unresolved run's record repeats, for each of its four attempts, the start (time, process, software hash, log entry) and the interruption as logged, and states what the status means. In error text the output root and the repository root are written as `<root>` and `<repo>`; nothing else is altered.
- *R-5.* The author explicitly authorises the deterministic fingerprint of the generators at manifest time. It is identity verification only, not experimental sampling. No generated value may influence parameters, selection, measurements or scientific decisions. Implemented: the hash is compared for equality and written down; nothing else reads it.

| ID | Reading | Alternative |
|---|---|---|
| R-1 | An operating-system error (full disk, quota, input/output error, lost mount) and a memory error are infrastructure interruptions. This reverses the earlier behaviour (B-11), under which a failed write gave a failure record | Keep storage errors as recorded failures |
| R-2 | An unresolved run changes only a `PASS`. A criterion that is already `FAIL`, `INSUFFICIENT_INFORMATION` or `UNINTERPRETABLE` keeps its label, and the file names the unresolved runs beside it. The existing label `INSUFFICIENT_INFORMATION` is used; no label was added. Criterion 5 (dates) is not touched. Which criteria a condition touches: calibration 2 to 4, healthy 1 to 4, fault 2 to 4 | Replace every label of a touched criterion |
| R-3 | An unresolved run still gets a record, because the protocol requires one for every planned ID, and it counts as a failed run in every denominator, as the protocol prescribes for failed runs. One unresolved calibration run therefore leaves the thresholds unavailable | None seen that respects the protocol |
| R-4 | Every start is an attempt, whatever ended it. A run whose process ends every time for a reason of its own (a crash in native code) cannot be told apart from an infrastructure fault at this level: it ends as `UNRESOLVED_INFRASTRUCTURE` after four attempts, with its history, and blocks its criteria | Treat a repeated identical crash as an experimental failure |
| R-5 | The manifest step hashes the first 128 outputs of each planned run's five generators. This reads generator output of the pilot namespace when the manifest is written, before the first draw. No sampler is run, no statistic is computed and no outcome is seen | Hash the seeds only, which would not test that the generators reproduce |

A retry happens only at the next `generate`: an interruption ends the work of its process, since its cause would meet the next run too. With several worker processes the others go on.

## 4. What the pipeline refuses

Invalid input is reported and stops the stage. Nothing is repaired, clipped or filled in.

- A stored run that does not verify (changed record, changed trace, record without trace, wrong identity, wrong namespace): generation and every analysis stage stop; the files are left as found.
- A trace without a record after a logged attempt: moved to `interrupted/`, never deleted, and the run is run again (Section 3a). Without a logged attempt: the driver stops and leaves it.
- An audit log whose chain does not verify, or a second `generate` on the same output tree: the driver stops.
- Seeds or a random-input hash that differ from the pre-generation manifest: nothing is drawn.
- An existing record or trace is never replaced. A run with a record, including a failure record, is never run again.
- A non-finite trace: refused by the storage layer and recorded as a failure.
- A sampler or measurement exception: a failure record with the phase, the exception and the progress made. It is never retried.
- A fifth attempt at a run: never made. The run is recorded as unresolved.
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

Evidence: `docs/v020/evidence/pytest_v020_confirmed.txt` (Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1). Earlier counts: `pytest_v020_policy.txt` (commit `8725fc4`, three tests deselected): 715 passed with 2 skipped, and 292 passed; `pytest_v020_criterion5.txt` (commit `0b1addc`, before the interrupted-run policy and without the measurement layer in the tree): 638 passed with 42 skipped, and 255 passed; `pytest_v020.txt` (commit `425dc61`): 614 passed with 31 skipped, and 220 passed.

| Run | Passed | Failed | Skipped | Deselected |
|---|---|---|---|---|
| Whole suite, tree as committed | 718 | 0 | 2 | 0 |
| `tests/v019` with the Stage C traces supplied and the full-length storage test enabled | 295 | 0 | 0 | 0 |

- The two skips need the Stage C traces, which are not in this branch, and the switch for a 30-second storage test. With both supplied, every test of `tests/v019` runs and passes.
- **No test is deselected, and no test samples in the pilot or smoke namespace** (B-12). The earlier evidence files were produced while `test_d_pilot_lib.py` still made short runs with smoke-namespace seeds.
- Before this work the same two runs gave 561 passed with 5 skipped, and 141 passed (`docs/v019/evidence/pytest_v019_audit.txt`).
- New test files: `test_f01_oracle.py` (13 tests), `test_pilot_driver.py` (43, of which 16 on interrupted runs), `test_pilot_pipeline.py` (51); extended: `test_audit_verify.py`, `test_audit_integrity.py`, `test_audit_recount.py` (23 tests of the criterion 5 dates). From its developer: `test_d_pilot_measure.py` (23).
- The measurement layer has SHA-256 `f8d5ac52…1b9be`, the version that was used read-only before. It is now in the branch; its hash enters the pre-generation manifest.

## 7. Commands

Order of a pilot, once the items of Section 3 are settled (nothing below has been run for the pilot or smoke namespace):

    PYTHONPATH=src python research/v019-audit/preflight.py --workers 8
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py status
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py manifest --preview --tdk-deadline <YYYY-MM-DD> --workers <n>   # writes nothing
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py manifest --tdk-deadline <YYYY-MM-DD> --workers <n>
    # the author reviews and pushes pre_generation_manifest.json, then writes authorization.json
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py generate    # again after an interruption; Ctrl-C once to stop cleanly
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py replay
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py calibrate
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py fit
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py issue
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py evaluate
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py summaries
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py checks
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py criteria --engineering-complete-on <YYYY-MM-DD> --results-available-on <YYYY-MM-DD>
    python research/v019-audit/verify.py --root research/v019-d-pilot --namespace rnglab-v019-d-pilot --deep sample --compare-primary --out research/v019-d-pilot/results/recount.json

The smoke comes first (freeze step 5). It uses the same `status`, `manifest`, `generate` and `replay` with `--smoke`, then `smoke-timing --smoke`; its preflight is `preflight.py --namespace rnglab-v019-d-pilot-smoke`, and its verification is `verify.py --root research/v019-d-pilot/smoke --namespace rnglab-v019-d-pilot-smoke --replicas 2 --deep all --no-campaign`.

Verification without any pilot data:

    PYTHONPATH=src python -m pytest tests -q -p no:cacheprovider --import-mode=importlib
    PYTHONPATH=src python tests/v019/synthetic_agreement.py
    PYTHONPATH=src python research/v019-audit/f01_oracle.py
