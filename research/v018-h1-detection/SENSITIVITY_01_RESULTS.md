# Stage C sensitivity analysis 01: results

**Post hoc and exploratory.** The plan is [SENSITIVITY_01.md](SENSITIVITY_01.md), frozen at `c7231fa` together with `sensitivity01.py`, before any result for arm C′ was computed. Output: `results/sensitivity01.json` (run log `results/sensitivity01.log`). Every number below is from `results/sensitivity01_table.md`, generated from those outputs by `sensitivity01_table.py`. The Stage C verdicts, files and classes are unchanged.

**Question.** Arm C′ is the conventional battery without its two word-level tests (lag-1 correlation, adjacent identical words), leaving the seven in-repository SP 800-22 tests. It is recalibrated to the same whole-run false-alarm rate. Is there then a cell where the application is harmed while C′ stays quiet?

## 1. Checks

- **Self-check:** with all ten columns, the script reproduced `results/gates.json` exactly.
- **Freeze record:** both files were last changed at `c7231fa`, the HEAD at run time.
- **Thresholds:** τ_C′ = 3.80732, against Stage C's τ_C = 3.81777; τ_A unchanged.
- **Healthy whole-run C′ alarms:** 182/4096 (0.0444336), CP95 [0.0383, 0.0512]. The V2 rule passes; V1 and the implementation checks pass. The results are interpretable.

## 2. The pre-specified result

- **No cell is in class H1 or undetected_harm under C′.** The verdict is NOT SUPPORTED in all four families.
- **Reading, by the rule fixed in plan §6:** the Stage C negative verdict does not depend on the two targeted word tests, in any family.
- **Five cells change class:**
  - DUP-0, DUP-1, DUP-2: conventional_first → harmless_quiet;
  - LAG1-0: conventional_first → indeterminate;
  - LAG1-5: both_detect_harmed → indeterminate.

## 3. What the class rule does not show: LAG1 comes close to H1

The rule in plan §6 is coarse, and for LAG1 it hides a strong pattern. Descriptively, under C′:

| LAG1 cell | Wrong k/n | C′ alarms (Stage C: C) | A alarms | Wrong with no C′ alarm (Stage C) | A − C′ mean [Bonferroni CI] |
|---|---|---|---|---|---|
| LAG1-0 | 26/256 | 10 (256) | 40 | 23 (0) | 0.117 [0.045, 0.195] |
| LAG1-1 | 35/256 | 11 (256) | 69 | 34 (0) | 0.227 [0.138, 0.320] |
| LAG1-2 | 53/256 | 16 (256) | 148 | 51 (0) | 0.516 [0.402, 0.625] |
| LAG1-3 | 98/256 | 21 (256) | 236 | 92 (0) | 0.840 [0.766, 0.906] |
| LAG1-4 | 152/256 | 27 (256) | 255 | 137 (0) | 0.891 [0.828, 0.948] |
| LAG1-5 | 224/256 | 41 (256) | 256 | 186 (0) | 0.840 [0.763, 0.905] |

- **The SP 800-22 subset barely sees LAG1** at these severities: 10 to 41 post-onset alarms in 256 runs, against 182 whole-run alarms in 4,096 healthy runs.
- **The application arm sees it:** its alarm count rises from 40 to 256, and the A − C′ interval excludes zero in all six cells.
- **The application is often wrong while C′ is quiet:** up to 186 of 256 runs at LAG1-5, against none with the full battery. For comparison, in the GFSR cells, where wrong-rates are near healthy, this count is 7 to 17.
- **Why no cell is H1:**
  - LAG1-4 misses only on the wrong-rate criterion: its lower bound is 0.496 against the required 0.5, while the C′ bound (upper 0.177 < 0.2) and the gap bound (lower 0.828 > 0) both pass.
  - LAG1-5 misses only on the C′ criterion: its upper bound is 0.241 against the required 0.2.

**Exploratory reading.** Against the in-repository SP 800-22 subset alone, LAG1 shows the pattern H1 describes: the application is harmed, the battery mostly stays quiet, and the application monitors alarm far more often. It falls short of the frozen class thresholds by narrow margins. In Stage C, the battery won on LAG1 because of the lag-1 test aimed at it. This marks a candidate region for a future pre-registered test. It is not a confirmed H1 result.

## 4. DUP, SHIFT and GFSR

- **DUP: the SP 800-22 subset still detects before harm.**
  - At the three lowest severities, C′ is quiet (9 to 24 alarms) and the wrong-rate is near healthy, so these cells become harmless_quiet. The Stage C alarms there, where the application was not harmed, came from the adjacent-identical-words test.
  - DUP-3 stays conventional_first under C′: 164 alarms with 17/256 wrong.
  - At DUP-4, C′ alarms in all 256 runs with 75/256 wrong.
  - Approximate entropy fires first most often (`results/sensitivity01_table.md`, first-firing list).
- **SHIFT: unchanged.** C′ alarmed in all 256 runs of every cell, as expected, since only SP 800-22 tests fired first in Stage C.
- **GFSR: unchanged,** harmless_quiet in all six cells.

## 5. What this changes

- **The Stage C verdicts stand,** as frozen. The word-level tests belong to the Stage C battery.
- **"The battery detects first" depends on the battery, for LAG1.** It holds with the targeted lag-1 test, and fails with the SP 800-22 subset alone. Statements of the Stage C result should name the battery.
- **For a future H1 test:** LAG1 against SP 800-22 alone, and against a stronger external battery (`ROADMAP.md` §7, decision 3), with the H1 criteria pre-registered, is now a motivated design.
- **For Stage D:** the forecast is issued at the battery alarm (`docs/v018/STAGE_DE_REVIEW.md` §2.3), so which battery issues it must be declared. With the SP 800-22 subset, most LAG1 runs would never trigger a forecast.

## 6. Limits

- **Post hoc:** a plan written after the Stage C results were known, on the same data, with the Stage C classes and thresholds. The near-H1 reading in Section 3 is descriptive.
- **Not the full SP 800-22 suite:** C′ is the in-repository subset, seven of fifteen tests. A fuller suite might catch LAG1.
- **Scope:** one application at one parameter point.
