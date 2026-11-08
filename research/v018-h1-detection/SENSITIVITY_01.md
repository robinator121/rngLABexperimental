# Stage C sensitivity analysis 01: the battery without the targeted word tests

**Status: analysis plan, frozen by the commit that adds it together with `sensitivity01.py`, before any result for arm C′ is computed.** Planned on 2026-11-08, *after* the Stage C results were known, so this analysis is post hoc and exploratory. It changes no Stage C file, class or verdict. It uses only stored Stage C data; nothing is generated.

## 1. Question

In Stage C, the first alarm in every LAG1 and DUP run came from one of the two word-level tests aimed at those faults: lag-1 correlation and adjacent identical words (REPORT §2; T5 in `results/report_tables.md`). So for those two families, "the battery detects first" may be true by construction.

**Question:** if the conventional arm is restricted to the seven in-repository SP 800-22 tests (arm C′) and recalibrated to the same whole-run false-alarm rate, is there any cell where the application is harmed while C′ stays quiet?

## 2. What was known when this plan was written

- **Seen:** every Stage C output (`results/gates.json`, `results/secondary.json`, `results/report_tables.md`, `REPORT.md`). That includes the first-firing counts of T5:
  - in LAG1 runs, lag-1 correlation fired first in every run;
  - in DUP runs, adjacent identical words fired first in every run;
  - in SHIFT runs, only SP 800-22 tests fired first;
  - the healthy-run first-firing counts.
- **Not seen:** per-test p-values or alarm rates of the SP 800-22 tests in any LAG1 or DUP cell, and any threshold or alarm count for C′.

## 3. Arm C′ and everything else

- **Inputs:** the stored p-values `conv_p` in `results/cells/*.npz`. C′ keeps columns 0 to 7 of `h1lib.CONV_P_KEYS`: frequency, block_frequency, runs, serial_1, serial_2, approximate_entropy, cumulative_sums, longest_run. These are the seven SP 800-22 tests of `PROTOCOL.md` §2.5; serial contributes two p-values. Dropped: lag1_correlation, adjacent_duplicates.
- **Statistic:** at each checkpoint, −log₁₀ of the smallest p-value over the eight columns, floored at `h1lib.P_FLOOR`, with missing values treated as no evidence. This is exactly `analyze.conv_series` restricted to those columns. The whole-run statistic is the maximum over the 39 checkpoints.
- **Threshold τ_C′:** the 950th order statistic of the 999 healthy calibration whole-run maxima, the same rule as `PROTOCOL.md` §2.5.
- **Unchanged from Stage C, computed by `analyze.gates` with only the conventional series replaced:**
  - detection and latency (`PROTOCOL.md` §2.7);
  - arm A, its constants and τ_A;
  - the harm endpoint and V1;
  - the cell classes and the Bonferroni level 1 − 0.05/24.
- **Resampling:** the bootstrap intervals use the namespace `rnglab-v018-h1-sensitivity01` (unused: `git grep`, 2026-11-08), so no Stage C resampling stream is reused.

## 4. Self-check, before any C′ result

`sensitivity01.py` first runs `analyze.gates` with all ten columns and the Stage C namespace. Its output must equal `results/gates.json` exactly, after a JSON round trip. If it does not, the script stops and computes no C′ result. The check can be run on its own (`--self-check`) before this plan is frozen, since it reproduces only what Stage C already reported.

## 5. Outputs

`results/sensitivity01.json`, holding:
- the self-check result;
- the git commit and the SHA-256 of this plan and of the script at run time, with a check that both match the freezing commit;
- τ_C′, and the healthy whole-run false-alarm rate of C′ with its 95% Clopper–Pearson interval. The Stage C V2 rule applies: if the upper bound exceeds 0.10, the C′ results are reported as uninterpretable;
- **per cell:**
  - C′ post-onset alarms k/n with intervals;
  - the A − C′ gap with intervals;
  - pass_but_broken, the class under C′, and the Stage C class under C from `results/gates.json`;
  - the counts only C′, only A and neither, and the median A − C′ latency;
  - the first-firing SP 800-22 test counts under C′ (descriptive);
- the per-family verdict under C′: SUPPORTED if any cell is in class H1.

## 6. Interpretation, fixed now

- **No cell of a family in class H1 or undetected_harm under C′:** for that family, the Stage C negative result does not depend on the two targeted word tests.
- **A cell in class H1 or undetected_harm under C′:** the in-repository SP 800-22 subset alone missed that fault where the application was harmed. This is an exploratory finding. It marks a candidate region for a future pre-registered test, not a confirmed H1 result.
- **Every class change between C and C′ is reported**, with the pass_but_broken counts, for every cell.
- **The Stage C verdicts stand.** The word-level tests belong to the frozen Stage C battery, and C′ does not replace it.

## 7. Limits

- **Post hoc:** the same data, thresholds rule and classes as Stage C, and a plan written after its results were known.
- **Not the full SP 800-22 suite:** C′ is the in-repository subset, seven of the suite's fifteen tests. Whether the full suite, or a battery such as TestU01, would catch LAG1 and DUP first is not answered.
- **Scope:** one application at one parameter point. The severity grids came from the pilot's wrong-rate (`PROTOCOL.md` §3.4), so they span the harm transition, which is where the classes in Section 6 would appear.
