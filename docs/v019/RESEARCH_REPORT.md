# v0.19 research report: audit, independent verification, forecasting foundation, assessment

**Date: 2026-11-11. Branch `feature/v019-audit`, from `33f701e`. Not pushed.** No pilot, smoke or campaign data of Stage D were generated. No frozen file and no historical result was changed.

## 1. State of Stage D when the work began

- The protocol is frozen and its hash matches (`b33846f`). The profile and fault adapter are committed (`33f701e`); their 26 tests pass, 106 with `tests/v018`.
- The measurement layer was being written by someone else: `measure.py` existed only as uncommitted files in the main checkout. This work therefore ran in its own worktree and never touched those files.
- Nothing after the per-run measurements exists in the pilot: no driver, calibration, fitting, scores, criteria, recount, smoke or pre-generation manifest. The TDK deadline is not recorded.

## 2. What was implemented

| Area | Artefact | Tests |
|---|---|---|
| Independent recount | `research/v019-audit/recount.py` | `test_audit_recount.py` (46) |
| Independent regeneration of a run | `replay.py` | `test_audit_d_pilot.py` (24) |
| Stored-output verifier and campaign recount | `verify.py` | `test_audit_verify.py` (7) |
| Campaign integrity and preflight | `integrity.py`, `preflight.py`, `golden.json` | `test_audit_integrity.py` (16) |
| Benchmark and audit numbers | `bench.py`, `audit_numbers.py`, `results/` | — |
| Forecasting components | `src/rnglab/prognosis/` | `test_prognosis.py` (22) |
| Methods validation study | `research/v019-forecast-methods/` (frozen, run once) | — |

Documents: `IMPLEMENTATION_AUDIT.md`, `FORECASTING_METHODS.md`, `SCIENTIFIC_ASSESSMENT.md`, `PUBLICATION_GAP_ANALYSIS.md`, this report. Test evidence: `docs/v019/evidence/pytest_v019_audit.txt`.

## 3. Questions investigated

1. Does the committed Stage D code do what the frozen protocol says, byte for byte?
2. Can every stored quantity be reproduced by code that shares nothing with the pilot?
3. Which forecasts are impossible without an assumption about the future, and by how much?
4. Under a declared family, does a forecaster built on that family beat simple baselines, and what happens when the family is wrong?
5. Is the system ready for a campaign of hours on this machine?
6. Where is RNGLab's strongest original contribution?

## 4. What was demonstrated

- **The profile, clock and fault adapter match the protocol.** An independent regeneration reproduced marks, traces, severity and stream hashes bit for bit on eight full-length runs.
- **The measurement layer's per-run outputs match independent recomputation**: 6,818 comparisons per set of four stored runs, no discrepancy, for both snapshots examined.
- **The independent functions reproduce Stage C**: threshold 3.8177723774232715, 199/4096 alarms, 3867/4096 coverage.
- **A bound on forecasting without a declared future**: any history-based forecaster loses at least E[(p₁ − p₂)²]/4 in one of two worlds with the same history. It is 0.2488 in the Phase II counterexample.
- **On a surrogate, under the declared ramp**, the declared-model forecaster is calibrated and has Brier skill 0.228 [0.189, 0.268] over an alarm-time baseline at L = 20,480; features of the pilot's "lite" form reach 0.155.
- **The interval-validity probability under drift is doubly noncentral.** Leaving the widening term out, as the pilot's oracle formula does by design, gave a miscalibrated forecaster at L = 20,480 on the surrogate.
- **A killed writer never leaves a run that passes as complete** under the atomic commit provided.

## 5. What failed, and what was found wrong

- **The frozen oracle call returns NaN** for 312 of 901 noncentralities from 8.55 (audit F-01). It needs the author's decision before analysis.
- **Outside the declared family the forecasters fail**: skill −0.256 (step-and-hold) and −0.218 (random walk) at L = 20,480; nominal 80% failure-time intervals cover 0.45 and 0.52; with late onsets the failure-time skill is −0.113.
- **Persistence**, the simplest model-based baseline, is far worse than the base rate under a ramp.
- **Readiness**: storage about 22.8 GiB against a 5.49 GiB baseline; about 8.8 hours at six workers against a 3.93-hour extrapolation; the bit cap cannot be reached and the smoke will not exercise it.
- **Faults in this work, found and fixed before use**: the verifier crashed on a truncated trace (found by its own test); a draft calibration rule had its threshold at the noise floor (changed before the freeze, on development seeds); three recalled bibliographic details were wrong (corrected against Crossref).
- **Not done**: no static type checker is installed, so the type hints are unchecked; no eight-worker benchmark; no bit-exact fast sampler (the reference path costs 6 to 9 s per run and is fixed by the protocol).

## 6. What remains unknown

- Whether validity can be forecast in the real application at all: the pilot has not run.
- Whether the sampler's batch means are normal and independent enough for the validity model.
- Whether any real generator degrades in a way the declared ramp describes.
- Whether the batch-length failure of Phase 2 is the two-time-scale effect it is read as, and whether established remedies already fix it.
- Everything about H3.

## 7. The next experiment

1. **Run the frozen Stage D pilot, once Section 2 of `PUBLICATION_GAP_ANALYSIS.md` is closed and the author authorises it.** It is the cheapest decisive experiment in hand: its stop rule says whether Stage D continues, and both outcomes are reportable. Verify it with `verify.py`.
2. **Then, or alongside, a pre-registered Direction A study on the stored healthy traces**: fit the two-scale model to the variance-ratio curve, compare the established interval methods against the exact energy, and confirm on fresh seeds. It is the shortest path to a result that could be new.

## 8. Is there a defensible original contribution today?

- **For a TDK: yes**, as a rigorous negative and methodological study: H1-C not supported at matched false alarms with an independent recount, and a documented transfer failure of an uncertainty rule against an exact reference.
- **For a journal: not yet.** No positive result is established; the forecasting results are on a surrogate; the identifiability bound is elementary.
- Details and the comparison of the two directions: `SCIENTIFIC_ASSESSMENT.md`.

## 9. Reproduction

Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1, glibc 2.44.

```bash
PYTHONPATH=src .venv/bin/python -m pytest tests -q
```

```bash
PYTHONPATH=src .venv/bin/python research/v019-audit/audit_numbers.py
```

```bash
.venv/bin/python research/v019-forecast-methods/report_tables.py
```

Source commits on `feature/v019-audit`: `4b55378` recount and replay; `326f795` verifier; `3cacdfa` integrity and preflight; `4e38355` forecasting components and benchmark; `459d144` methods protocol, frozen; `a568973` its results; `783c1c2` audit numbers; `6afaf69` audit; `23e08ab` forecasting document; then the assessment documents.
