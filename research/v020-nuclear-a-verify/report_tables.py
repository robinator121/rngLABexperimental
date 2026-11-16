"""Tables of the Layer A verification, from its result file only.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-a-verify/report_tables.py [--smoke]

Reads results/verification.json (smoke/verification.json with --smoke) and
writes tables.md next to it. Computes nothing new: it formats recorded values.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def fmt(value: float) -> str:
    if value == 0:
        return "0"
    return f"{value:.3g}" if 1e-3 <= abs(value) < 1e4 else f"{value:.2e}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--smoke", action="store_true")
    folder = HERE / ("smoke" if parser.parse_args().smoke else "results")
    r = json.loads((folder / "verification.json").read_text())
    v = r["verdicts"]
    lines = [f"# Layer A verification: tables from `{folder.name}/verification.json`", "",
             f"Namespace `{r['namespace']}`; code `{r['code']['git_commit']}` (dirty: {r['code']['git_dirty']}); "
             f"started {r['started_utc']}; {r['elapsed_seconds']:.1f} s in total.", "",
             "## Verdicts", "", f"- R1: **{v['R1']}**" + (f" (failed: {', '.join(v['R1_failed_items'])})"
                                                          if v["R1_failed_items"] else ""),
             f"- R2: **{v['R2']}**", ""]
    lines += ["| Deterministic check | Holds |", "|---|---|"]
    lines += [f"| {k} | {ok} |" for k, ok in v["deterministic"].items()]
    a1, x, rp = r["analytic_crosscheck"], r["xref"], r["replay"]
    lines += ["", f"A1: largest |expm − Bateman path sum| = {fmt(a1['max_abs_expm_minus_paths'])}; "
              f"largest |generator column sum| = {fmt(a1['max_abs_generator_column_sum'])}.",
              f"XREF: {x['nuclei']} nuclei, events {x['events'][0]} and {x['events'][1]}, discrete columns equal: "
              f"{x['discrete_columns_equal']}, final states equal: {x['final_states_equal']}, stream digests equal: "
              f"{x['stream_digests_equal']}, largest relative time difference {fmt(x['max_relative_time_difference'])}.",
              f"REPLAY (chunk {rp['chunk']}): events digest equal: {rp['events_sha256'] == rp['original_events_sha256']}, "
              f"stream digest equal: {rp['stream_sha256'] == rp['original_stream_sha256']}.", ""]
    lines += ["## Families (Holm, family-wise alpha 0.01)", "", "| Family | Tests | Smallest p | Holm rejections |",
              "|---|---:|---:|---|"]
    for name, fam in r["families"].items():
        lines.append(f"| {name} | {len(fam['tests'])} | {fmt(min(fam['p_values']))} | "
                     f"{', '.join(fam['rejected']) or 'none'} |")
    for name, fam in r["families"].items():
        lines += ["", f"### {name}: every test", "",
                  "| Test | Statistic (z for binomial tests) | p | Holm-adjusted p |",
                  "|---|---:|---:|---:|"]
        members = {**{f"NC1 {k} {t}": d for k in ("S1", "S2", "S3") for t, d in r["runs"]["NC1"][k]["tests"].items()},
                   **{f"NC2 {k} {t}": d for k in ("S1", "S2", "S3") for t, d in r["runs"]["NC2"][k]["tests"].items()},
                   **{f"ENS {t}": d for t, d in r["ensemble"]["tests"].items()}}
        for test, p, adjusted in zip(fam["tests"], fam["p_values"], fam["holm"]["adjusted_p_values"]):
            d = members[test]
            stat = d.get("statistic", d.get("z"))
            lines.append(f"| {test} | {fmt(stat)} | {fmt(p)} | {fmt(adjusted)} |")
    pc = r["positive_controls"]
    lines += ["", "## Positive controls", "", f"Level: p < {fmt(pc['level'])}.", "",
              "| Control | Target | p | Rejected |", "|---|---|---:|---|"]
    lines += [f"| {k} | {pc[k]['target']} | {fmt(pc[k]['p_value'])} | {pc[k]['rejected']} |" for k in ("PC1", "PC2")]
    for k, key in (("PC1", "S1"), ("PC2", "S2")):
        tests = r["runs"][k][key]["tests"]
        lines += ["", f"{k}, every test (descriptive, outside the families):", "",
                  "| Test | p |", "|---|---:|"]
        lines += [f"| {t} | {fmt(d['p_value'])} |" for t, d in tests.items()]
    r3 = r["runs"]["NC1"]["S3"]["R3"]
    lines += ["", "## R3 (descriptive): transient equilibrium, NC1 S3", "",
              f"Decays in (20, 40]: Ya-410 {r3['decays_Ya_20_40']}, Yb-410 {r3['decays_Yb_20_40']}; ratio "
              f"{r3['ratio']:.4f}, approximate 95% interval [{r3['interval95'][0]:.4f}, {r3['interval95'][1]:.4f}] "
              f"(standard error of the log ratio {r3['log_ratio_se']:.5f}; correlation of the per-nucleus window "
              f"indicators {r3['indicator_correlation']:.3f}); "
              f"exact expectation {r3['exact_expected_ratio']:.6f}; equilibrium limit "
              f"{r3['transient_equilibrium_limit']:.6f}.", "",
              "## Single-machine timings (seconds, wall clock; not HPC measurements)", "",
              "| Run | Seconds |", "|---|---:|"]
    lines += [f"| {k} | {t:.3f} |" for k, t in r["timings_seconds"].items()]
    lines += [f"| XREF reference engine | {x['reference_engine_seconds']:.3f} |", ""]
    (folder / "tables.md").write_text("\n".join(lines))
    print(f"wrote {folder / 'tables.md'}")


if __name__ == "__main__":
    main()
