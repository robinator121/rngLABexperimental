"""Numbers quoted in docs/v020/SCIENTIFIC_RED_TEAM.md. Exact calculations; nothing is simulated or drawn.

1. criterion_4_selection: in the simplified model that design.py uses for the power of criterion 4
   (two equiprobable forecast values 0.5 +/- d, baseline 0.5, true Brier skill 4 d^2), the number of
   eligible evaluation runs whose outcome falls on the forecast side is binomial, and the percentile
   bootstrap of the skill is a binomial quantile. So the probability that a family passes and the
   mean observed skill of the passing families follow by enumeration. No fitting uncertainty and no
   random eligibility, as in design.py.
2. short_intervals: coverage of a 16-batch t interval when the batch-means variance is a fraction
   `ratio` of the long-run variance (normal batch means): P(|t_15| < c sqrt(ratio)). The ratios are the
   ones measured in v0.18 on the refreshed path-integral sampler (research/v018-api-validation/REPORT.md,
   section 4.1); whether they carry over to the pilot's runs is what its short-interval check measures.
3. pre_onset_alarms: share of fault runs with an alarm before onset if the nominal whole-horizon
   false-alarm probability 0.05 were spread evenly over the 159 battery checkpoints.
4. criterion_1: the pass probabilities of design.json, with the two bounds for all three lengths.

    python research/v019-audit/redteam_numbers.py      writes docs/v020/evidence/redteam_numbers.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy import stats

REPO = Path(__file__).resolve().parents[2]
OUTPUT = REPO / "docs/v020/evidence/redteam_numbers.json"
DESIGN = REPO / "research/v019-d-pilot/design.json"
LEVEL = 1 - 0.05 / 3
MINIMUM_SKILL = 0.10
VARIANCE_RATIOS = {32: 0.601, 128: 0.875, 512: 0.971}   # batch length in sweeps: BM variance ratio (v0.18 report)


def criterion_4_selection(n: int, true_skill: float) -> dict:
    d = np.sqrt(true_skill) / 2
    hit, miss = (0.5 - d) ** 2, (0.5 + d) ** 2            # loss when the outcome is on the forecast side, or not
    k = np.arange(n + 1)
    weight = stats.binom.pmf(k, n, 0.5 + d)
    skill = 1 - (k * hit + (n - k) * miss) / (n * 0.25)   # increasing in k
    # lower end of the two-sided percentile interval: the (1-level)/2 quantile of the resampled count
    k_low = stats.binom.ppf((1 - LEVEL) / 2, n, np.clip(k / n, 0, 1))
    lower = 1 - (k_low * hit + (n - k_low) * miss) / (n * 0.25)
    passed = (skill >= MINIMUM_SKILL) & (lower > 0)
    p = float(weight[passed].sum())
    return {"eligible": n, "true_skill": true_skill, "p_pass": p,
            "p_observed_at_least_0.10": float(weight[skill >= MINIMUM_SKILL].sum()),
            "mean_observed_skill": float((weight * skill).sum()),
            "mean_observed_skill_given_pass": float((weight[passed] * skill[passed]).sum() / p) if p > 0 else None,
            "p_at_least_one_of_three_families": 1 - (1 - p) ** 3}


def short_intervals() -> list[dict]:
    c = stats.t.ppf(0.975, 15)
    return [{"alarm_sweep": 16 * batch, "batch_sweeps": batch, "variance_ratio": ratio,
             "implied_healthy_coverage": float(2 * stats.t.cdf(c * np.sqrt(ratio), 15) - 1)}
            for batch, ratio in VARIANCE_RATIOS.items()]


def pre_onset_alarms() -> dict:
    onsets = np.arange(2048, 8193)
    ends = np.arange(512, 40961, 256)
    before = np.array([(ends <= t0).sum() for t0 in onsets])
    return {"battery_checkpoints": int(len(ends)), "mean_checkpoints_ending_at_or_before_onset": float(before.mean()),
            "share_of_fault_runs_if_false_alarms_were_uniform": float(0.05 * before.mean() / len(ends)),
            "assumption": "nominal 0.05 whole-horizon rate, spread evenly over checkpoints; not a measurement"}


def criterion_1() -> list[dict]:
    rows = json.loads(DESIGN.read_text(encoding="utf-8"))["size_options"][2]["criterion_1"]["probabilities"]
    return [{"true_coverage": r["true_coverage"], "p_pass_one_length": r["p_pass"],
             "all_three_if_independent": r["p_pass"] ** 3, "all_three_if_fully_dependent": r["p_pass"]} for r in rows]


def main() -> int:
    out = {"schema": "rnglab-v020-redteam-numbers-v1",
           "criterion_4_selection": [criterion_4_selection(n, s) for n in (240, 400) for s in (0.0, 0.05, 0.10, 0.15, 0.20)],
           "criterion_4_note": "design.py's simplified model; bootstrap replaced by its exact binomial quantile",
           "short_intervals": short_intervals(), "pre_onset_alarms": pre_onset_alarms(), "criterion_1": criterion_1()}
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for row in out["criterion_4_selection"]:
        print(f"n={row['eligible']} true skill {row['true_skill']:.2f}: P(pass) {row['p_pass']:.3f}, mean observed "
              f"{row['mean_observed_skill']:.3f}, given pass {row['mean_observed_skill_given_pass']}")
    for row in out["short_intervals"]:
        print(f"alarm at {row['alarm_sweep']}: batch {row['batch_sweeps']}, implied healthy coverage {row['implied_healthy_coverage']:.4f}")
    print(out["pre_onset_alarms"])
    for row in out["criterion_1"]:
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
