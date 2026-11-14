# Stage D pilot: readiness report

**Status: for the author's review, 2026-11-14.** Branch `feature/v020-stage-d-readiness`, not pushed. The frozen protocol is unchanged (`research/v019-d-pilot/protocol.json`, SHA-256 `5329464210019e968116e8c3c3a13a09f33b12e0368ecfaa69752a1eb43c5e97`). **No pilot, smoke or campaign run was generated.** No pilot execution is permitted until the author has reviewed this report and approved it.

Companion documents: `F01_NUMERICAL_RESOLUTION.md`, `STAGE_D_IMPLEMENTATION_STATUS.md`, `SCIENTIFIC_RED_TEAM.md`.

**Update of 2026-11-13.** The author decided the interrupted-run policy: `regenerate_same_id`, at most three retries. It is implemented and tested on both sides (`STAGE_D_IMPLEMENTATION_STATUS.md`, Section 3a), the measurement layer is in the branch, and preflight can check the smoke. Conditions C-1 and C-2 below are met. The assessment is unchanged.

**Update of 2026-11-14.** The author confirmed the readings R-1 to R-5 as implemented, with two clarifications: the status `UNRESOLVED_INFRASTRUCTURE` describes an unresolved execution, not a proven cause, and keeps its whole attempt history and error evidence; the generator fingerprint computed for the manifest is authorised as identity verification only. The author also decided that no test may sample in the pilot or smoke namespace before the manifest. Two test files did; both now use the audit-only namespace (status document, B-12). The assessment is unchanged.

## 1. Assessment: CONDITIONAL GO, for the smoke first

- **Today the pilot cannot start, and should not.** Several decisions that only the author can take are open (Section 2.5), the smoke has not been run, and no manifest is authorised. The driver refuses to run without them.
- **No finding requires a redesign.** The execution path is complete from the pre-generation manifest to `criteria.json`, it rejects invalid input instead of repairing it, and it agrees with an independent recount on a synthetic campaign. The one numerical defect found in a frozen formula (F-01) concerns the oracle reference only, cannot change a stop criterion, and has a documented resolution that awaits a decision.
- **The go is conditional and staged.** Passing tests do not make the pilot ready: the pipeline has never read a real pilot run. The smoke (ten runs, its own namespace, excluded from every result) is the first real test of generation, storage, replay and the independent verification together. The pilot proper should start only after the smoke has passed these checks.

**Conditions for the smoke** (all are actions or decisions of people, none is new engineering):

| # | Condition | Who |
|---|---|---|
| C-1 | `research/v019-d-pilot/measure.py` and its tests are committed; the full test suite passes with the file in the tree. **Met in this branch.** Since 2026-11-14 the file's tests use the audit-only namespace (B-12) and the whole suite runs without deselection | its developer |
| C-2 | Interrupted-run policy chosen. **Decided on 2026-11-13: `regenerate_same_id`, at most three retries** | author |
| C-3 | Storage location with about 31 GiB free for the pilot (the preflight requirement); `results/runs/` and `traces/` kept outside Git | author |
| C-4 | Worker count chosen; a value other than 8 is recorded as a deviation from the protocol | author |
| C-5 | Preflight reports READY in the environment that will generate | operator |
| C-6 | Smoke manifest written, reviewed, pushed; `authorization.json` written by the author | author |

**Further conditions for the pilot:**

| # | Condition | Who |
|---|---|---|
| C-7 | Smoke: 10 of 10 runs stored with status SUCCESS, I-1 to I-3 passed, all ten replayed bit for bit, the independent verification reports no discrepancy | operator, then author |
| C-8 | Oracle rule decided before fitting: Amendment 01 option A or B recorded in `amendments.json`, or option C with a declared recount exception | author |
| C-9 | Decision on supplementary descriptive analyses (red-team RT-04, RT-05): amend before generation, or declare them exploratory | author |
| C-10 | TDK deadline recorded in the pilot's manifest, or criterion 5 accepted as INSUFFICIENT_INFORMATION | author |
| C-11 | Pilot manifest written from a clean working tree, reviewed, pushed; `authorization.json` names its hash | author |
| C-12 | The machine is left to the pilot for its duration (Section 4), on mains power, with sleep disabled | operator |

If any of C-1 to C-6 is not met, the assessment for the smoke is NO-GO. If C-7 fails, the pilot is NO-GO until the cause is found and reported; a failed smoke is not repaired by editing its output.

## 2. What is known, and how

### 2.1 Mathematically established

These follow from definitions and do not depend on a computation.

- The oracle quantity q(δ) = P(|T| ≤ c), T noncentral t with 15 degrees of freedom, is even in δ, equals 0.95 at δ = 0 and decreases in |δ|.
- For δ ≥ 0 the lower term F(−c; δ) is at most Φ(−δ), because T ≤ −c < 0 requires Z + δ < 0. Where the library fails (δ ≥ 8.53) this is at most 7.3 × 10⁻¹⁸. Replacing a failed lower term by 0 therefore changes q by less than that.
- A forecast of the outcome needs a declared continuation: two continuations with the same history give any history-based forecaster an excess Brier score of at least E[(p₁ − p₂)²]/4 in one of them (`docs/v019/FORECASTING_METHODS.md`, Section 1).
- In the simplified model that `design.py` uses for criterion 4, a family that passes at a true skill of 0.10 reports a mean skill of 0.124 to 0.129 (exact enumeration, `evidence/redteam_numbers.json`). The model is an idealisation; the arithmetic is exact.

### 2.2 Verified computationally

On this machine (Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1, glibc 2.44), by scripts whose outputs are committed.

- **F-01.** `scipy.special.nctdtr(15, δ, −c)` is not finite at 1,554 of 12,001 grid points, all between δ = 8.53 and 43.61, and only for δ > 0. Three stable formulations agree with a 70-digit reference to better than 3 × 10⁻¹⁵. On 20,000 points in [8.6, 43.5], 521 change between NaN and a number at the next double (`research/v019-audit/results/f01_oracle.json`).
- **The pipeline's rules,** on deterministic synthetic fixtures: atomic storage, resume after a killed process, failure records, refusal of damaged stores, the authorisation gate, the threshold order statistic, issuance boundaries, the evaluation barrier, the frozen optimiser call, fallbacks, scores, criteria 1 to 5 and their precedence (Section 6).
- **Agreement with the independent recount** on a synthetic campaign of 207 runs: 1,972 of 1,974 compared fields under the frozen oracle rule; the two that differ are counts of unavailable oracle values. 1,974 of 1,974 with option A forced in memory (`evidence/synthetic_agreement.json`).
- **The generation path on real full-length runs** in the audit-only namespace: every timed run was stored with status SUCCESS and passed I-1 to I-3 (`evidence/pipeline_benchmark.json`). These runs were deleted; they are not results.
- **Carried over from v0.19:** bit-for-bit replay of golden runs on this platform; the independent implementations of the transforms and of the endpoint.
- **Determinism.** By construction a run is a function of its seed labels, so the same identity gives the same bytes and the number of worker processes cannot change a result. The replay of golden runs confirms this on this platform; a killed and resumed synthetic generation ends in the same store (`test_killed_generation_resumes_to_the_same_store`).

### 2.3 Supported only by surrogate or synthetic experiments

- Forecast skill of 0.228 [0.189, 0.268] under the declared ramp family, and negative skill under two other families: a Gaussian surrogate workload, not the sampler (`docs/v019/FORECASTING_METHODS.md`).
- The oracle formula without the widening term was miscalibrated at L = 20,480 in that surrogate (F-02).
- On the synthetic campaign the frozen oracle rule left 6 of 9 cells without a primary score. The fixture's bias is large enough to reach large δ; this says how the rule behaves there, not how often the pilot will be there.
- The behaviour of every analysis stage: exercised on synthetic records in the measurement layer's schema, never on stored pilot runs.

### 2.4 Hypothetical

Nothing in the repository supports or contradicts these yet. They are what the pilot is for, or what it assumes.

- That the real sampler under the declared ramp gives undecided outcomes at L = 20,480 (criterion 2), an alarm early enough (criterion 3), and skill over the alarm sweep (criterion 4).
- That the estimate's bias is linear in the mean severity, as the oracle assumes.
- That δ reaches the failing range in the pilot. With the surrogate's scales it would, just, at L = 40,960; that is an illustration.
- That short issuance intervals under-cover as the v0.18 variance ratios suggest (0.881 for an alarm at sweep 512).
- That throughput measured over one to three minutes holds for seven hours on this machine (Section 3).
- That the faults and the ramp resemble anything a generator does in use.

### 2.5 Decisions that need the author's approval

| Decision | Options | Needed before |
|---|---|---|
| Interrupted-run policy | decided 2026-11-13: `regenerate_same_id`, three retries. The readings R-1 to R-5 were confirmed on 2026-11-14 (status document, Section 3a) | the smoke manifest |
| Tests that sampled under smoke-namespace seeds (B-12) | decided 2026-11-14: not authorised. Both test files now use the audit-only namespace | done |
| Worker count | 8 (protocol) / another value, recorded as a deviation | the smoke manifest |
| Storage location and archive | — | the smoke |
| Smoke directory layout | as the measurement layer writes it / as the protocol names it | the smoke |
| `authorization.json` interlock | keep / remove | the smoke |
| A storage error during generation (disk full) | settled by the decision of 2026-11-13 as read here (R-1): an interruption, retried under the same identity | the smoke manifest |
| Oracle rule (F-01) | Amendment 01 option A / option B / option C with a recount exception | `fit` |
| Supplementary descriptive analyses | amend before generation / exploratory afterwards / none | the pilot manifest |
| Criterion 5 dates | TDK deadline; engineering-complete date; results-available date | the pilot manifest; `criteria` |
| Report files (`report_numbers.json`, `RESULTS.md`) | the author's text | after `criteria` |
| What a failed criterion 2 or 4 will be reported as | wording | the pilot manifest (advised) |

## 3. Resources, measured

**What was measured.** 195 full-length runs went through the pipeline's own `run_one` (sampler, I-2, I-3, measurement, atomic storage) in the audit-only namespace, with 1 to 12 worker processes, the pilot's mix of conditions and two repetitions. Then twelve more were stored and read back. All were deleted. Scripts and outputs: `tests/v019/pipeline_benchmark.py` → `evidence/pipeline_benchmark.json`; `tests/v019/pipeline_read_timing.py` → `evidence/pipeline_read_timing.json`.

**The machine.** 8 cores, 16 hardware threads, 15.7 GiB of memory, 92 GiB free on the file system at the start. It was in desktop use during the timing, with about 5 GiB of memory available. The timing of each configuration lasted one to three minutes.

| Workers | Runs timed | Runs per hour | Pilot generation, projected | Seconds per run | Slower than alone | Workers on CPU | Sum of peak memory |
|---|---|---|---|---|---|---|---|
| 1 | 12 | 272–275 | 21.5–21.8 h | 13.0 | 1.00 | 0.99 | 0.30 GiB |
| 4 | 24 | 659 | 8.4–8.5 h | 20.3 | 1.57 | 0.93 | 1.14 GiB |
| 6 | 36 | 804–819 | **6.6–6.8 h** | 24.2 | 1.86 | 0.90 | 1.70 GiB |
| 8 (protocol) | 48 | 747–774 | **6.8–7.0 h** | 33.0 | 2.54 | 0.86 | 2.14 GiB |
| 12 | 72 | 752–777 | 7.0–7.2 h | 51.1 | 3.93 | 0.86 | 3.40 GiB |
| 16 | — | not measured | — | — | — | — | — |

The ranges are the two repetitions. The projection is the pilot's count of each condition times the mean time of that condition's runs, divided by the number of workers. Sixteen workers were skipped both times: the machine had 5.1 to 5.2 GiB available against the 5.9 GiB the script requires for them.

**Runtime.**
- Six, eight and twelve workers give the same generation time within the spread of the measurement: 6.6 to 7.2 hours. Six gave the highest throughput. More than six processes add nothing on this machine.
- The v0.19 estimate of 8.8 hours at six workers came from one measurement of a different code path. On the pipeline's path the figure is 6.6 to 6.8 hours. The design's extrapolation (3.9 hours at eight workers) is too low by a factor of 1.7 to 1.8.
- **The protocol's eight workers are as good as any other choice.** Nothing measured here argues for deviating from the protocol.

**Storage.** 22.5 GiB for 5,999 runs: 3.9 to 4.1 MB per run, of which the JSON run record is 3.2 to 3.4 MB and the compressed trace 0.7 to 0.9 MB. This confirms the v0.19 figure (22.8 GiB) and not the design's 5.5 GiB, which covers three uncompressed arrays only. The preflight check requires 1.25 × 5,999 × 4.4 MB, about 31 GiB, free.

**Memory.** Peak resident memory was 293 to 303 MiB per worker at every worker count. Eight workers need about 2.2 GiB. Available memory did not fall below 3.7 GiB during the timing. The preflight check requires 300 MiB of available memory per worker, which is the measured peak without a margin.

**Scheduling overhead.** At eight workers a run takes 2.5 times as long as alone, yet the workers were on a processor for 86% of the wall time, and system time was 2 to 5% of processor time. So most of the slowdown is slower execution while running, not waiting: the processes share caches and memory bandwidth, and beyond eight they share cores. Involuntary context switches per run rose from 425 alone to 596 at six workers, 1,222 at eight and 7,495 at twelve.

**Where the time goes** (one process, mean over the pilot's mix):

| Stage | Seconds | Share |
|---|---|---|
| Sampler and fault transform | 7.3 | 56% |
| Measurement (159 battery and 157 application checkpoints, endpoints) | 3.2 | 25% |
| Checks I-2 and I-3 | 2.4 | 18% |
| Atomic storage | 0.13 | 1% |

- A healthy run takes 9.4 to 9.6 s, a fault run 16.1 to 17.2 s. The difference is the fault transform in the sampler (2.4 to 3.0 s) and the independent re-application of the profile in I-3 (4.6 to 5.1 s).
- The step check for replica 0 of each fault family adds 6.6 to 9.2 s, three times in the pilot.
- **No optimisation was made.** v0.19 found that a vectorised transform reproduces the I-3 bytes in a fraction of the time (F-07). Adopting it would change verification code whose hash enters the manifest, not a scientific output. It still needs the author's approval, and it would save roughly an hour.

**After generation** (one process, reading only; the arithmetic of the stages is not included and was not measured at full size):

| Step | Per run | For the pilot |
|---|---|---|
| Store scan with trace hashes (resume, run manifest, completeness gate) | 24 ms | 2.4 min each |
| Analysis stages, about 22,000 run loads in all | 48 ms | 17.6 min |
| Independent verification: read, shallow check of 5,974 runs, regeneration of 25 | 47 ms, 78 ms, 1.7 s | 13.2 min |
| I-4 replay by the production engine | 25 runs, each about as long as a run alone | 5 to 6 min |

## 4. Execution plan

Nothing in this section has been run for the pilot or smoke namespace. Commands are in `STAGE_D_IMPLEMENTATION_STATUS.md`, Section 7.

**Before anything.** Conditions C-1 to C-6. Preflight READY. A clean working tree at the commit that the manifest will name.

**Smoke** (ten runs, namespace `rnglab-v019-d-pilot-smoke`, excluded from every result).
1. `preflight.py --namespace rnglab-v019-d-pilot-smoke`; `manifest --smoke` (`--preview` shows it without writing), review, push, `authorization.json`.
2. `generate --smoke`: a few minutes. About 40 MB.
3. `replay --smoke`: all ten runs regenerated and compared.
4. `smoke-timing --smoke`; then `verify.py --root research/v019-d-pilot/smoke --namespace rnglab-v019-d-pilot-smoke --replicas 2 --no-campaign --deep all`.
5. Stop and report if any run failed, any check failed, or the verification lists a discrepancy. No threshold, model or criterion is computed from the smoke; the pipeline refuses to.

**Pilot** (5,999 runs).

| Resource | Requirement | Basis |
|---|---|---|
| Disk | 31 GiB free where `results/runs/` and `traces/` live; 22.5 GiB expected use | Section 3; preflight |
| Memory | 2.3 GiB available for eight workers; more than 4 GiB advised | 300 MiB per worker |
| Generation time | 6.6 to 7.2 hours projected; **plan for up to 11 hours** | Section 3; allowance of 1.5 on the largest projection |
| After generation | about one hour: replay, six analysis stages, independent verification | Section 3 |

The allowance is not a measurement. The projection rests on minutes of load on a machine in desktop use. Hours of sustained load, a lower clock frequency when the processor is hot, or other programs can lengthen it. The driver runs the conditions in the protocol's order, so the 2,999 healthy-type runs come first and the 3,000 fault runs, each of which takes almost twice as long, come last: the count of stored runs grows faster in the first part than in the second. `status` shows the count at any time.

**Order.** `manifest`, review, push, `authorization.json`, `generate`, `replay`, `calibrate`, `fit`, `issue`, `evaluate`, `summaries`, `checks`, `criteria`, then the independent verification with `--compare-primary`. This is the order the tests use. It keeps the summaries of fault runs, which contain evaluation outcomes, after fitting.

**Recovery.**

| Event | What is on disk | What to do | What the pipeline will not do |
|---|---|---|---|
| A worker or the machine dies during generation | Every stored run is whole: a run is written to a staging directory and renamed into place, the record last. Runs in progress have no record; the audit log holds their `start` entry | Run `generate` again. It scans the store (2.4 min), verifies every trace hash and the chain of the audit log, logs each interrupted attempt, and runs those identities again with the same seeds | Run a stored run again; make a fifth attempt at a run |
| The operator wants to stop | — | Ctrl-C once (or SIGTERM): every process finishes its run in progress and stops. No attempt is consumed. `generate` continues later | — |
| A run raises in the sampler, a check or the measurement | A failure record with phase, exception and progress | Nothing. It is reported in `failures.jsonl` and keeps its place in every denominator | Run it again or replace it |
| **The disk fills up, a write fails, or memory runs out** | No record for the run in progress; an `interrupted` entry with the reason. The worker that met the error stops | Remove the cause (free space, memory), then run `generate` again. Each run has three retries | Write a failure record for it; go on past the cause |
| **A run is interrupted four times** | A failure record with phase `infrastructure` and the status `UNRESOLVED_INFRASTRUCTURE`, listing the four attempts | Nothing in the pipeline. The run counts as failed, `criteria.json` names it, and the criteria that use it cannot be PASS. One such calibration run leaves the thresholds unavailable. The author decides what to report | Run it a fifth time; drop it from a denominator |
| A stored file is changed or damaged | The file, as found | The driver and every stage stop and name the run. The author decides; nothing is repaired | Delete, re-create or skip it |
| An analysis stage dies | Either its output file, complete, or nothing: each result file is written through a temporary file and renamed | Run the stage again if its output does not exist | Overwrite an existing result file |
| `evaluate` dies between its two files | `predictions.jsonl` without `scores.json` | The stage refuses to run again, because the evaluation runs are opened once. The author decides whether to remove the file and repeat, and records it | Repeat the evaluation by itself |
| A failed calibration run | A failure record | Nothing can be done: the threshold needs 999 usable maxima. Criteria 2 to 4 are then INSUFFICIENT_INFORMATION or UNINTERPRETABLE | Lower the order statistic or use another run |

**Reproducibility.** Bit-for-bit replay holds on one platform (v0.19, F-08). `replay` and the deep part of the independent verification must run in the environment that generated, before any library is updated. The manifest records Python, NumPy, SciPy, the C library and the machine type.

## 5. Risks that remain after the conditions are met

| Risk | Consequence | Mitigation in place |
|---|---|---|
| The analysis stages meet a real record they cannot handle | A stage stops with an error; no result is repaired | The smoke; the synthetic fixtures use the measurement layer's own schema and functions |
| A calibration run is interrupted four times | Threshold unavailable; criteria 2 to 4 count as failed | Three retries, each after the operator has removed the cause; stopping on purpose costs no attempt |
| The driver is killed again and again while the same runs are in progress (a script that restarts it, an unstable machine) | Those runs use up their three retries and end unresolved | Every start is an attempt (R-4). Stop with Ctrl-C, not with a kill; find the cause before the next `generate` |
| Disk fills during generation | The runs in progress are interrupted and retried later; nothing is recorded as failed | Preflight requires 31 GiB; reading R-1 |
| Generation with several worker processes meets the audit log for the first time in the smoke | A fault of the driver would show there | The same functions are tested in one process, with kills; the smoke is ten runs |
| Generation takes longer than planned | Delay only: generation resumes | Audit log, resume |
| A library update between generation and replay | I-4 fails for a reason unrelated to the pilot | Manifest records versions; replay directly after generation |
| The pipeline and the recount share a misreading of the protocol | Agreement between them would not show it | A second person reads the protocol against `STAGE_D_IMPLEMENTATION_STATUS.md`, Section 1 |
| The oracle under the frozen rule | Unavailable oracle scores, and a recount that differs in oracle availability | The decision C-8 |

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

## 7. Reproduction

Everything below runs without pilot data, and nothing below samples in the pilot or smoke namespace.

    # numerical study of the oracle call (F-01)
    PYTHONPATH=src python research/v019-audit/f01_oracle.py
    python -m pytest tests/v019/test_f01_oracle.py -q

    # pipeline and independent recount on synthetic fixtures
    PYTHONPATH=src python -m pytest tests/v019 -q -p no:cacheprovider --import-mode=importlib
    PYTHONPATH=src python tests/v019/synthetic_agreement.py

    # timing (generates and deletes full-length runs in the audit-only namespace; about 25 minutes).
    # The figures of Section 3 were measured before the audit log of attempts was added; they were not measured again.
    PYTHONPATH=src python tests/v019/pipeline_benchmark.py --workers 1,4,6,8,12,16 --repeat 2
    PYTHONPATH=src python tests/v019/pipeline_read_timing.py

    # numbers of the red-team review
    python research/v019-audit/redteam_numbers.py

    # whole suite
    PYTHONPATH=src python -m pytest tests -q -rs -p no:cacheprovider --import-mode=importlib

    # preflight for the smoke, and a preview of its manifest; neither writes anything
    PYTHONPATH=src python research/v019-audit/preflight.py --namespace rnglab-v019-d-pilot-smoke
    PYTHONPATH=src python research/v019-d-pilot/run_pilot.py manifest --smoke --preview
