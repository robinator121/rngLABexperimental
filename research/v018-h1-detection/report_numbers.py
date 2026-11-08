"""Summary numbers for the Stage C (H1) report, from stored outputs only.

Reads results/{gates,secondary,recount,replay,run_manifest,amendment01_comparison,
amendment01_preanalysis}.json, results/run.log and design.json. Writes and prints
results/report_numbers.json; only selections, counts, ranges, sums and ratios
are used. Latency ranges describe stored cell medians, not pooled-run medians.

    .venv/bin/python research/v018-h1-detection/report_numbers.py
"""
import json
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FAMILIES = ("SHIFT", "LAG1", "DUP", "GFSR")


def read(name):
    return json.loads((RESULTS / (name + ".json")).read_text())


def get(value, path):
    for key in path.split("."):
        if isinstance(value, dict):
            value = value.get(key)
        elif isinstance(value, list) and key.isdigit() and int(key) < len(value):
            value = value[int(key)]
        else:
            return None
    return value


def values(cells, path):
    return [v for c in cells.values() if (v := get(c, path)) is not None]


def bounds(items):
    return {"min": min(items, default=None), "max": max(items, default=None)}


def ratio(k, n):
    return k / n if k is not None and n else None


def summarize(cells):
    classes = values(cells, "class")
    lat = {key + "_sum": sum(values(cells, "latency." + key))
           for key in ("only_A", "only_C", "neither")}
    lat.update({"median_" + unit: bounds(values(cells, "latency." + unit + ".median_A_minus_C"))
                for unit in ("sweeps", "bytes")})
    out = {"cells": list(cells), "class_counts": {k: classes.count(k) for k in classes},
           "wrong_rate": bounds(values(cells, "wrong.rate")),
           "wrong_cp95_upper_max": max(values(cells, "wrong.cp95.1"), default=None),
           "cells_with_C_alarm_in_every_run": sum(c.get("runs") is not None and
               get(c, "C_alarm.k") == c["runs"] for c in cells.values()),
           "pass_but_broken": {**bounds(values(cells, "pass_but_broken")),
                               "sum": sum(values(cells, "pass_but_broken"))},
           "onset_flagged_sum": sum(values(cells, "onset_flagged")),
           "all_runs_assessed": all(c.get("runs") is not None and
               get(c, "wrong.n") == c["runs"] for c in cells.values()), "latency": lat}
    out.update({arm + "_alarm_rate": bounds(values(cells, arm + "_alarm.rate")) for arm in ("C", "A")})
    return out


def secondary(cells):
    field = "health_test_post_onset_alarm_runs"
    firing = {}
    for c in cells.values():
        for test, count in c.get("first_firing_conventional_test", {}).items():
            firing[test] = firing.get(test, 0) + count
    return {field: bounds(values(cells, field)), field + "_per_cell":
            {name: c.get(field) for name, c in cells.items()},
            "first_firing_conventional_test": firing,
            **{key: bounds(values(cells, key)) for key in ("mean_cv_bias", "td_wrong")}}


def runtime():
    manifest, pre = read("run_manifest"), read("amendment01_preanalysis")
    plan = json.loads((HERE / "design.json").read_text()).get("runtime", {})
    lines = (RESULTS / "run.log").read_text().splitlines()
    start = manifest.get("first_start_utc")
    end = lines[-1][4:] if lines and lines[-1].startswith("end ") else None
    minutes = (datetime.fromisoformat(end) - datetime.fromisoformat(start)).total_seconds() / 60 if start and end else None
    planned = plan.get("wall_minutes_8_workers_with_slowdown_1.35")
    return {"start": start, "end": end, "wall_minutes": minutes, "workers": manifest.get("workers"),
            "blocks": len(manifest.get("blocks", [])),
            "planned_wall_minutes": planned, "ratio": ratio(minutes, planned),
            "runs_planned": plan.get("total_runs"), "runs_stored": pre.get("campaign_runs")}


def main():
    g, s, rc, replay, amendment = [read(n) for n in
        ("gates", "secondary", "recount", "replay", "amendment01_comparison")]
    cells = g.get("cells", {})
    groups = {f: {n: c for n, c in cells.items() if c.get("family") == f} for f in FAMILIES}
    v1, v2 = g["validity"]["V1"], g["validity"]["V2"]
    coverage = ratio(v1.get("covered"), v1.get("runs"))
    healthy = {k: v1.get(k) for k in ("runs", "covered", "cp95")}
    healthy.update(coverage=coverage, wrong_rate=1 - coverage if coverage is not None else None,
                   paired_A_minus_C=v2.get("paired_A_minus_C"), paired_ci95=v2.get("paired_ci95"))
    for arm in ("C", "A"):
        q = v2.get("rate_" + arm, {})
        healthy["rate_" + arm] = {**{k: q.get(k) for k in ("alarms", "runs", "cp95")},
                                  "rate": ratio(q.get("alarms"), q.get("runs"))}
    pooled = {n: c for n, c in cells.items() if c.get("family") in FAMILIES[:3]}
    summary = summarize(pooled)
    synthetic = {k: summary[k] for k in ("cells", "cells_with_C_alarm_in_every_run")}
    synthetic.update(only_A_sum=summary["latency"]["only_A_sum"],
                     pass_but_broken_sum=summary["pass_but_broken"]["sum"],
                     latency={"median_sweeps": summary["latency"]["median_sweeps"]})
    for label in ("conventional_first", "indeterminate", "both_detect_harmed"):
        subset = summarize({n: c for n, c in pooled.items() if c.get("class") == label})
        synthetic[label] = {k: subset[k] for k in ("cells", "wrong_rate")}
        if label == "indeterminate":
            synthetic[label].update(C_alarm_rate_min=subset["C_alarm_rate"]["min"],
                                    A_alarm_rate=subset["A_alarm_rate"])
    sc = s.get("cells", {})
    sh = sc.get("healthy", {})
    secondaries = {"healthy": {k: sh.get(k) for k in
                  ("health_test_post_onset_alarm_runs", "mean_cv_bias", "td_wrong")},
                  "families": {f: secondary({n: sc.get(n, {}) for n in group}) for f, group in groups.items()}}
    secondaries["healthy"]["health_test_post_onset_alarm_rate"] = ratio(
        sh.get("health_test_post_onset_alarm_runs"), v1.get("runs"))
    comparisons = rc.get("comparisons", [])
    checks = {"recount": {"agreeing": sum(c.get("agree") is True for c in comparisons),
                          "total": len(comparisons), "all_agree": rc.get("all_agree")},
              "replay": {"cells": len(replay), "I-4_true": sum(c.get("I-4") is True for c in replay.values()),
                         "I-3_cells": sum("I-3" in c for c in replay.values()),
                         "I-3_true": sum(c.get("I-3") is True for c in replay.values())},
              "amendment_comparison": {"superseded_ran": amendment.get("superseded_ran"),
                  **{k: len(amendment[k]) if amendment.get(k) is not None else None
                     for k in ("different", "only_approved", "only_superseded")}},
              "thresholds": {k: get(g, "thresholds." + k) for k in ("C", "A")},
              "bonferroni_level": g.get("bonferroni_level")}
    out = {"healthy": healthy, "families": {f: summarize(c) for f, c in groups.items()},
           "synthetic": synthetic, "secondary": secondaries, "checks": checks, "runtime": runtime()}
    text = json.dumps(out, indent=1) + "\n"
    (RESULTS / "report_numbers.json").write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
