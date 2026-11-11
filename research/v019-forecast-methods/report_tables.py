"""Tables and reading-rule outcomes of the methods validation study, from results.json only.

    python research/v019-forecast-methods/report_tables.py [--smoke]

Writes tables.md and verdicts.json next to results.json. Every number in
docs/v019/FORECASTING_METHODS.md that concerns the study is copied from these two files.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ORDER = ("base-rate", "alarm-time", "lite-logistic", "persistence", "bayes-ramp-no-widening", "bayes-ramp", "oracle")


def f(value, digits=4):
    return "n/a" if value is None else f"{value:.{digits}f}"


def interval(pair, digits=3):
    return "n/a" if pair is None else f"[{pair[0]:.{digits}f}, {pair[1]:.{digits}f}]"


def verdicts(results: dict) -> dict:
    a = results["scenarios"]["A_declared"]
    out = {}
    cells = {L: a["lengths"][L]["forecasters"] for L in a["lengths"]}
    r1 = {L: {k: cells[L][k]["calibration"]["calibrated"] for k in ("bayes-ramp", "oracle")} for L in cells}
    out["R1"] = {"per_length": r1, "result": "MET" if all(all(v.values()) for v in r1.values()) else "NOT MET"}
    skill = cells["20480"]["bayes-ramp"]["skill_over_alarm_time"]
    out["R2"] = {"skill": skill["skill"], "interval": skill["interval"],
                 "result": "MET" if skill["interval"] and skill["interval"][0] > 0 else "NOT MET"}
    difference = a["lengths"]["20480"]["brier_difference:bayes-ramp-minus-lite-logistic"]
    out["R3"] = {**difference, "result": "MET" if difference["interval"][1] < 0 else "NOT MET"}
    out["R4"] = {L: {**a["lengths"][L]["brier_difference:bayes-ramp-minus-bayes-ramp-no-widening"],
                     "widening_matters": a["lengths"][L]["brier_difference:bayes-ramp-minus-bayes-ramp-no-widening"]["interval"][1] < 0}
                 for L in a["lengths"]}
    r5 = {}
    for name in ("C_step_hold", "D_random_walk"):
        rows = results["scenarios"][name]["lengths"]
        flags = {L: rows[L]["forecasters"]["bayes-ramp"]["calibration"]["calibrated"] for L in rows if rows[L]["forecasters"]}
        r5[name] = {"bayes_ramp_calibrated": flags,
                    "result": "EXPECTATION MET" if not all(flags.values()) else "EXPECTATION NOT MET"}
    out["R5"] = r5
    pair = results["identifiability_pair"]
    out["R6"] = {"forecasts_identical": pair["alarms_identical"] and all(v["forecasts_identical"] for v in pair["lengths"].values()),
                 "bound": {L: v["bound_worst_world_excess"] for L, v in pair["lengths"].items()}}
    ft = a["failure_time"]
    out["R7"] = {"rps_skill": ft["rps_skill"], "interval": ft["rps_skill_interval"],
                 "result": "MET" if ft["rps_skill_interval"][0] > 0 else "NOT MET"}
    return out


def tables(results: dict) -> str:
    lines = [f"# Methods validation study: tables ({'SMOKE, excluded from results' if results['smoke'] else 'frozen study'})", "",
             f"Source: results.json, namespace `{results['namespace']}`, protocol SHA-256 `{results['protocol_sha256']}`.", "",
             f"- alarm threshold {results['calibration']['threshold']:.4f}; estimated sigma_eff {results['calibration']['sigma_eff_hat']:.4f}, "
             f"beta {results['training']['beta_hat']:.5f}, kappa {results['training']['kappa_hat']:.3f}",
             "- healthy valid fraction per length: " + ", ".join(f"{L}: {v:.4f}" for L, v in results["calibration"]["healthy_valid_fraction"].items()), ""]
    for name, row in results["scenarios"].items():
        lines += [f"## {name}", "", f"Runs {row['runs']}, alarms {row['alarms']}, no alarm {row['no_alarm']}, "
                  f"pre-onset alarms {row['pre_onset_alarms']}.", ""]
        for L, cell in row["lengths"].items():
            lines += [f"**L = {L}**: eligible {cell['eligible']}, late alarm {cell['late_alarm']}, valid fraction {f(cell['valid_fraction'])}", ""]
            if not cell["forecasters"]:
                lines += ["(fewer than 50 eligible runs: not scored)", ""]
                continue
            lines += ["| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |",
                      "|---|---|---|---|---|---|---|"]
            for k in ORDER:
                v = cell["forecasters"][k]
                s = v["skill_over_alarm_time"]
                lines.append(f"| {k} | {f(v['brier'])} | {f(v['log_score'])} | {f(s['skill'], 3)} {interval(s['interval'])} | "
                             f"{f(v['calibration']['spiegelhalter_z'], 2)} | {f(v['calibration']['weighted_gap'], 3)} | "
                             f"{'yes' if v['calibration']['calibrated'] else 'no'} |")
            for key in sorted(k for k in cell if k.startswith("brier_difference:")):
                d = cell[key]
                lines.append(f"\n{key.split(':')[1]}: {d['difference']:+.5f} {interval(d['interval'], 5)}")
            lines.append("")
        ft = row["failure_time"]
        lines += [f"**Failure time**: {ft['runs']} runs with an alarm, {ft['censored']} censored at the horizon, "
                  f"{ft['already_failed_at_issuance']} already failed at issuance; Kaplan-Meier median failure sweep "
                  f"{ft['km_median_failure_sweep']}, median remaining sweeps {ft['km_median_remaining_sweeps']}.", "",
                  "| forecaster | RPS | 80% interval coverage [CP95] | upper end censored | median abs. error (sweeps) | both inside | "
                  "forecast inside, outcome censored | forecast censored, outcome inside | both censored |", "|---|---|---|---|---|---|---|---|---|"]
        for k in ("empirical", "bayes-ramp"):
            c, e = ft["interval_coverage"][k], ft["median_error"][k]
            lines.append(f"| {k} | {f(ft['rps'][k])} | {f(c['coverage'], 3)} {interval(c['cp95'])} | {c['upper_end_censored']} | "
                         f"{e['median_absolute_error_sweeps']} | {e['both_inside']} | {e['forecast_inside_outcome_censored']} | "
                         f"{e['forecast_censored_outcome_inside']} | {e['both_censored']} |")
        lines += ["", f"RPS skill of bayes-ramp over empirical: {f(ft['rps_skill'], 3)} {interval(ft['rps_skill_interval'])}", ""]
    pair = results["identifiability_pair"]
    lines += ["## Identical histories, two continuations (E_continue, E_freeze)", "",
              f"Alarms identical: {pair['alarms_identical']}.", "",
              "| L | forecasts identical | mean p (continue) | mean p (freeze) | bound on worst-world excess | bound holds for every forecaster |",
              "|---|---|---|---|---|---|"]
    for L, v in pair["lengths"].items():
        lines.append(f"| {L} | {v['forecasts_identical']} | {f(v['mean_p_continue'])} | {f(v['mean_p_freeze'])} | "
                     f"{f(v['bound_worst_world_excess'], 5)} | {v['bound_holds_for_every_forecaster']} |")
    lines += ["", "Excess Brier score over the path oracle (mean of (q - p)^2), per world:", "",
              "| L | forecaster | continue | freeze | worse of the two |", "|---|---|---|---|---|"]
    for L, v in pair["lengths"].items():
        for k, e in (v["excess_over_path_oracle"] or {}).items():
            lines.append(f"| {L} | {k} | {f(e['continue'], 5)} | {f(e['freeze'], 5)} | {f(max(e.values()), 5)} |")
    return "\n".join(lines) + "\n"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args(argv)
    folder = HERE / ("smoke" if args.smoke else "results")
    results = json.loads((folder / "results.json").read_text())
    outcome = verdicts(results)
    (folder / "verdicts.json").write_text(json.dumps(outcome, indent=1, sort_keys=True) + "\n")
    (folder / "tables.md").write_text(tables(results))
    for rule, value in outcome.items():
        print(rule, value.get("result") if isinstance(value, dict) and "result" in value else json.dumps(value)[:300])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
