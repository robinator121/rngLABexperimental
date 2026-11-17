"""Tables of the B2 verification, from its result file only.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/report_tables.py [--smoke]

Reads results/verification.json (smoke/verification.json with --smoke) and
writes tables.md next to it. Computes nothing new: it formats recorded values.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def fmt(value) -> str:
    if value is None:
        return "-"
    value = float(value)
    if value == 0:
        return "0"
    return f"{value:.3g}" if 1e-3 <= abs(value) < 1e4 else f"{value:.2e}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    folder = HERE / ("smoke" if parser.parse_args().smoke else "results")
    r = json.loads((folder / "verification.json").read_text())
    v = r["verdicts"]
    lines = [f"# B2 verification: tables from `{folder.name}/verification.json`", "",
             f"Namespace `{r['namespace']}`; code `{r['code']['git_commit']}` (dirty: {r['code']['git_dirty']}); "
             f"started {r['started_utc']}; {r['elapsed_seconds']:.1f} s in total.", "", "## Verdicts", "",
             f"- R1: **{v['R1']}**" + (f" (failed: {', '.join(v['R1_failed_items'])})" if v["R1_failed_items"] else ""),
             f"- R2: **{v['R2']}**", "", "| Deterministic check | Holds |", "|---|---|"]
    lines += [f"| {k} | {ok} |" for k, ok in v["deterministic"].items()]
    laws = r["laws"]
    lines += ["", f"D8: {laws['laws_checked']} interval laws; largest |decimal - path-tree expm| "
                  f"{fmt(laws['max_abs_decimal_minus_tree'])}; largest |end states - Layer A expm| "
                  f"{fmt(laws['max_abs_end_state_minus_layer_a'])}; largest |sum - 1| {fmt(laws['max_abs_sum_minus_one'])}; "
                  f"brackets hold: {laws['brackets_hold']}; most digits used: {laws['max_digits']}; non-finite "
                  f"references: {', '.join(laws.get('nonfinite_references', [])) or 'none'}.", ""]
    lines += ["| XREF case | Parts | Digests equal | Arrays equal | Streams equal | Reference engine seconds |",
              "|---|---:|---|---|---|---:|"]
    for x in r["xref"]:
        lines.append(f"| {x['scheme']}, part_max {x['part_max']} | {x['parts']} | {x['digests_equal']} | "
                     f"{x['arrays_equal']} | {x['stream_digests_equal']} | {x['reference_seconds']:.2f} |")
    lines += ["", "## Families (Holm, family-wise alpha 0.01)", "",
              "| Family | Tests | Smallest p | Smallest Holm-adjusted p | Rejections |", "|---|---:|---:|---:|---|"]
    for name, fam in r["families"].items():
        lines.append(f"| {name} | {len(fam['tests'])} | {fmt(min(fam['p_values']))} | "
                     f"{fmt(min(fam['holm']['adjusted_p_values']))} | {', '.join(fam['rejected']) or 'none'} |")
    members = {}
    for family in ("NC1", "NC2"):
        for entry in r["runs"][family].values():
            members.update(entry["tests"])
    members.update(r["ensemble"]["tests"])
    members.update(r["ensemble"]["comparison_tests"])
    members.update(r["big_tests"])
    for name, fam in r["families"].items():
        lines += ["", f"### {name}: every test", "", "| Test | Statistic (z for binomial and moment tests) | p | Holm-adjusted p |",
                  "|---|---:|---:|---:|"]
        for test, p, adjusted in zip(fam["tests"], fam["p_values"], fam["holm"]["adjusted_p_values"]):
            d = members[test]
            lines.append(f"| {test} | {fmt(d.get('statistic', d.get('z')))} | {fmt(p)} | {fmt(adjusted)} |")
    not_run = [name for name, d in members.items() if "not_applicable" in d]
    errored = [f"{name} ({d['error']})" for name, d in members.items() if "error" in d]
    lines += ["", f"Tests not applicable (fewer than two cells after pooling): {', '.join(not_run) or 'none'}.",
              "", f"Planned tests that could not be computed (failed items): {', '.join(errored) or 'none'}."]
    calibration = folder / "calibration.json"
    if calibration.exists():
        c = json.loads(calibration.read_text())
        lines += ["", "## R0: calibration (complete null batteries)", "",
                  f"`{folder.name}/calibration.json`: {c['batteries']} batteries, namespace `{c['namespace']}`, "
                  f"code `{c['code']['git_commit']}` (dirty: {c['code']['git_dirty']}). C1 holds: **{c['C1']}**.", "",
                  "| Family | Tests | Batteries that FAIL under the null | Rate | 95% interval | C1 |",
                  "|---|---:|---:|---:|---|---|"]
        for name, f in c["families"].items():
            lines.append(f"| {name} | {f['tests']} | {f['failed_batteries']} | {f['false_fail_rate']:.4f} | "
                         f"{f['false_fail_ci95'][0]:.4f} to {f['false_fail_ci95'][1]:.4f} | {f['C1_holds']} |")
        for name, f in c["families"].items():
            lines += ["", f"### {name}: test types at the Holm level ({fmt(f['levels'][-1])})", "",
                      "| Type | Tests | Batteries with a rejection | Share of the union bound | 95% interval |",
                      "|---|---:|---:|---:|---|"]
            for kind, row in sorted(f["types"].items()):
                level = row["levels"][-1]
                low, high = level["battery_rate_ci95_over_union_bound"]
                lines.append(f"| {kind} | {row['instances']} | {level['batteries_with_rejection']} | "
                             f"{level['battery_rate_over_union_bound']:.2f} | {low:.2f} to {high:.2f} |")
    pc = r["positive_controls"]
    lines += ["", "## Positive controls", "", f"Level: p < {fmt(pc['level'])}.", "",
              "| Control | Target | p | Rejected |", "|---|---|---:|---|"]
    lines += [f"| {k} | {pc[k]['target']} | {fmt(pc[k]['p_value'])} | {pc[k]['rejected']} |"
              for k in ("PC-B1", "PC-B2 single run", "PC-B2 comparison")]
    for key in ("PC-B1", "PC-B2 single run"):
        lines += ["", f"{key}, every test (descriptive, outside the families):", "", "| Test | p |", "|---|---:|"]
        lines += [f"| {t} | {fmt(d.get('p_value'))} |" for t, d in pc[key]["all_tests"].items()]
    lines += ["", "## Huge populations", "", "| Initial nuclei | Parts | Seconds | Final populations |", "|---:|---:|---:|---|"]
    for size, entry in r["big"].items():
        lines.append(f"| {size} | {entry['parts']} | {entry['seconds']:.2f} | {entry['populations'][-1]} |")
    lines += ["", "## Single-machine timings (seconds, wall clock; not HPC measurements)", "",
              "| Run | Seconds |", "|---|---:|"]
    lines += [f"| {k} | {t:.3f} |" for k, t in r["timings_seconds"].items()]
    (folder / "tables.md").write_text("\n".join(lines) + "\n")
    print(f"wrote {folder / 'tables.md'}")


if __name__ == "__main__":
    main()
