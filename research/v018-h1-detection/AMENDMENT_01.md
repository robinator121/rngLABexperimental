# Stage C Amendment 01: analysis fixes before any analysis

**Status: APPROVED by the author, 2026-11-08** (proposed 2026-11-07, while generation was running under the unchanged manifest `manifest_pregeneration.json`, commit `0aca45a`). Generation finished before the approval. No analysis, replay or recount had been run on campaign data, and no campaign result had been looked at, when the approval was given and recorded.

## Origin

A read-only review by OpenAI Codex (`gpt-6-astra`, high reasoning) compared `analyze.py` and `recount.py` with the frozen protocol (`PROTOCOL.md` §§2.4–5.4), and found seven gaps. Codex then implemented the fixes in a separate worktree. Claude Code reviewed the diff line by line and accepted it unchanged (`docs/v018/AI_LOG.md`).

## Changes

| # | File | Protocol | Before | After |
|---|---|---|---|---|
| 1 | `analyze.py` | §5.1 (a failure means stop) | `interpretable` depended on V1 and V2 only | it also requires I-1 to I-4 to pass |
| 2 | `analyze.py` | §5.1 | I-3 and I-4 passed if `replay.json` simply lacked entries | missing entries fail: I-4 needs every cell, I-3 every fault cell |
| 3 | `analyze.py`, `recount.py` | §§4, 5.3 (Bonferroni over 24) | divisor = number of fault cells found | frozen constant `BONFERRONI_CELLS = 24` |
| 4 | `recount.py` | §2.4 (non-`ESTIMATED` counts as wrong) | a zero or non-finite batch variance was not counted as wrong | it is counted as wrong |
| 5 | `analyze.py`, `recount.py` | §2.4 (onset exclusions) | a cell with every run excluded would divide by zero | null rate and interval, class "indeterminate" |
| 6 | `recount.py` | §5.4 (recount every cell rate) | counts compared, rates and denominators not | denominators and the wrong, C-alarm and A-alarm rates also compared |
| 7 | `analyze.py`, `recount.py` | §4 (paired intervals at 95% and Bonferroni) | Bonferroni interval only | adds the descriptive 95% interval `A_minus_C.ci95` (seed `gap95/<cell>`), recomputed by the recount |

## Expected effect on this campaign

**None on any verdict, class or count**, given facts fixed before generation:
- `grid.json` has exactly 24 fault cells, so change 3 gives the same level.
- Change 4 concerns an event of probability zero for continuous data.
- The latest onset sweep seen in the pilot was 2,485 against the limit of 2,560 (`pilot/summary.json`), so change 5 is not expected to apply.
- Changes 1 and 2 matter only if a check fails or the replay is incomplete.
- Changes 6 and 7 add checks and a descriptive interval.

The amendment closes gaps that could have hidden a failure; it does not tune the analysis.

## Hashes

| File | Manifest (before) | Amended |
|---|---|---|
| `research/v018-h1-detection/analyze.py` | `d172f146…af5e42` | `7d36f099…1930e5` |
| `research/v018-h1-detection/recount.py` | `fd44ed3a…03f1e0` | `56c127cc…23755c` |

Tests after the change: 25 passed (`tests/v018/test_h1_campaign.py`, `tests/v018/test_h1_lib.py`). That run was Codex's; the full suite is re-run before commit.

## Approval

The author's approval, verbatim:

> **Stage C Amendment 01: APPROVED, 2026-11-08**
>
> I approve the seven proposed analysis and recount corrections as protocol-conformance amendments, subject to the documented test suite passing and verification that the generation process and campaign data remain unchanged.
>
> Approval is granted before inspection or analysis of campaign results, based on the recorded blind review.
>
> The amendment must not alter the prespecified hypothesis, statistical decision thresholds, experimental grid, generation parameters, or interpretation criteria beyond correcting the identified protocol violations.
>
> The additional 95% confidence intervals are descriptive and must not replace the prespecified multiplicity-adjusted inference.
>
> Claims about the amendment having no effect on campaign outcomes remain expectations, not established findings. Any observed changes must be reported transparently.
>
> Record the approved source hashes and test results before executing campaign analysis.

## Pre-analysis record

Written by `amendment01_record.py` to `results/amendment01_preanalysis.json` (2026-11-08T05:09:16Z), after the amended scripts were installed and before any analysis.

**Approved source hashes**

| File | Superseded (`0aca45a`, = manifest) | Approved |
|---|---|---|
| `analyze.py` | `d172f14659f3c61bbf51f98c0159d5e77b36e681d41a3b4e9117f065e3af5e42` | `7d36f0999f0ef8b23c299f0440f4f3b589b9458aad7ad0e1770b631b251930e5` |
| `recount.py` | `fd44ed3ad8473ab9781f1bb7ad105acf01d6d7c7c67cfc2949b72c3a6803f1e0` | `56c127cccec034fd663cc54e0ccdddd800232d12fa89f81e8b050d371d23755c` |

**Conditions of the approval**

| Condition | Check (record field) | Result |
|---|---|---|
| Documented test suite passes | full suite, `docs/v018/evidence/pytest_stageC_amendment01.txt` (`tests`) | 425 passed |
| Generation process unchanged | of the 14 manifest-listed files, only the two amended ones differ (`only_amended_files_changed`); protocol and grid hashes match the manifest | yes; yes |
| Campaign data unchanged | every cell `.npz` and `.jsonl` matches `results/run_manifest.json` (`campaign_file_mismatches`); 11,239 runs stored of 11,239 planned | 26 cells, no mismatch |
| No change to thresholds or decision criteria | top-level constants of both scripts compared (`constants_changed`, `constants_removed`); the only addition is `BONFERRONI_CELLS = 24` | none changed or removed |
| 95% intervals descriptive only | the cell class is computed from the Bonferroni intervals alone (`analyze.py`, the `classify(...)` call on `ci_bonferroni`); `ci95` is only reported | by code inspection |
| Before analysis | none of `replay.json`, `gates.json`, `secondary.json`, `recount.json` exists (`analysis_outputs_present`) | none present |

**How a change in outcome will be detected.** That the amendment changes no verdict, class or count is an expectation, not a finding. To test it, the superseded `analyze.py` (`0aca45a`) will also be run on the same data, with its output kept apart in `results/amendment01_superseded/`. The approved version's output, `results/gates.json`, is the analysis of record. Every difference between the two in verdicts, classes, counts, rates and intervals, including a failure of the superseded version to run, will be reported in the Stage C report.
