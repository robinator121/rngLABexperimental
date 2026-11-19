"""Tables of the Layer A / B2 benchmark, from its result file only.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-bench/report_tables.py

Reads results/benchmark.json and writes results/tables.md. Measured values and projections are kept in
separate sections; nothing is computed here beyond unit conversions and ratios of recorded values.
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
MB = 1024.0 * 1024.0


def spread(row: dict, key: str = "wall_seconds") -> str:
    """The smallest and largest of a configuration's three timed samples."""
    values = row["samples"][key]
    return f"{min(values):.4g}-{max(values):.4g}"


def main() -> None:
    r = json.loads((HERE / "results" / "benchmark.json").read_text())
    m = r["machine"]
    out = ["# Layer A and B2 benchmarks: tables from `results/benchmark.json`", "",
           f"Code `{r['code']['git_commit']}` (dirty: {r['code']['git_dirty']}); started {r['started_utc']}; "
           f"{r['elapsed_seconds']:.0f} s in total.", "",
           f"Machine: {m.get('cpu_model', 'unknown CPU')}, {m['logical_cpus']} logical CPUs, "
           f"{m.get('memory_total_bytes', 0) / MB / 1024:.1f} GiB RAM; {m['platform']}; Python {m['python']}, "
           f"NumPy {m['numpy']}, SciPy {m['scipy']}. One workstation; none of this is an HPC measurement.", "",
           "Every configuration ran in fresh processes: times are medians of three untraced processes; "
           "\"Python peak\" is tracemalloc's peak from a fourth process; \"peak RSS\" includes the interpreter and "
           f"the imported libraries (baseline: {r['baseline']['max_rss_bytes'] / MB:.0f} MiB, import only). "
           "\"samples\" gives the smallest and largest of the three wall times. Fields that differ by seed "
           "(Layer A's events and output size) are medians of the three processes.", "",
           "## 1. Measured: Layer A (S2, horizon 12, events observed on the six-interval grid)", "",
           "| nuclei | wall s | samples | CPU s | events | events/s | words/s | output MiB | Python peak MiB | peak RSS MiB |",
           "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in r["layer_a"]:
        out.append(f"| {row['config']['n']:.0e} | {row['wall_seconds']:.4f} | {spread(row)} | {row['cpu_seconds']:.4f} | "
                   f"{row['events']} | {row['events'] / row['wall_seconds']:.3g} | {row['words_per_second']:.3g} | "
                   f"{row['output_bytes'] / MB:.1f} | {row['traced_peak_bytes'] / MB:.1f} | "
                   f"{row['max_rss_bytes'] / MB:.0f} |")
    out += ["", "## 2. Measured: B2 against the population (S2, six intervals, plan construction included)", "",
            "| nuclei | parts | words/interval | wall s | samples | CPU s | plan s | s per interval | Python peak MiB | peak RSS MiB |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in r["b2_population"]:
        out.append(f"| {row['config']['n']:.4g} | {row['parts']} | {row['words_per_interval']} | "
                   f"{row['wall_seconds']:.4f} | {spread(row)} | {row['cpu_seconds']:.4f} | {row['plan_seconds']:.4f} | "
                   f"{(row['wall_seconds'] - row['plan_seconds']) / 6:.5f} | {row['traced_peak_bytes'] / MB:.2f} | "
                   f"{row['max_rss_bytes'] / MB:.0f} |")
    out += ["", "## 3. Measured: B2 against the number of intervals (S2, 1e6 nuclei, uniform grid on [0, 12])", "",
            "| intervals | outputs kept | wall s | samples | CPU s | plan s | ms per interval | output MiB | Python peak MiB | peak RSS MiB |",
            "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in r["b2_intervals"]:
        k = row["config"]["intervals"]
        out.append(f"| {k} | {row['config']['keep']} | {row['wall_seconds']:.3f} | {spread(row)} | {row['cpu_seconds']:.3f} | "
                   f"{row['plan_seconds']:.3f} | {1e3 * (row['wall_seconds'] - row['plan_seconds']) / k:.3f} | "
                   f"{row['output_bytes'] / MB:.2f} | {row['traced_peak_bytes'] / MB:.2f} | "
                   f"{row['max_rss_bytes'] / MB:.0f} |")
    out += ["", "## 4. Measured: checkpoints (B2, 1e6 nuclei, 2,000 uniform intervals, outputs not kept)", "",
            "| checkpoint every | wall s | samples | CPU s | ms per interval | Python peak MiB | peak RSS MiB |",
            "|---:|---:|---:|---:|---:|---:|---:|"]
    for row in r["b2_checkpoints"]:
        every = row["config"]["checkpoint_every"]
        out.append(f"| {every if every else 'never'} | {row['wall_seconds']:.3f} | {spread(row)} | "
                   f"{row['cpu_seconds']:.3f} | {1e3 * (row['wall_seconds'] - row['plan_seconds']) / 2000:.3f} | "
                   f"{row['traced_peak_bytes'] / MB:.2f} | {row['max_rss_bytes'] / MB:.0f} |")
    out += ["", "## 5. Measured: plan construction against distinct interval lengths (B2, 1e6 nuclei, 200 intervals)", "",
            "One set of interval laws serves every interval of a uniform grid; a grid whose 200 lengths all differ "
            "needs 200 sets.", "",
            "| grid | plan s | plan samples | wall s | Python peak MiB | peak RSS MiB |", "|---|---:|---:|---:|---:|---:|"]
    for row in r["b2_plan_lengths"]:
        grid = "uniform" if row["config"].get("uniform") else "all lengths distinct"
        out.append(f"| {grid} | {row['plan_seconds']:.4f} | {spread(row, 'plan_seconds')} | {row['wall_seconds']:.4f} | "
                   f"{row['traced_peak_bytes'] / MB:.2f} | {row['max_rss_bytes'] / MB:.0f} |")
    out += ["", "## 6. Measured: thread settings (B2, 1e6 nuclei, 1,000 intervals)", "",
            "| setting | wall s | CPU s | CPU/wall |", "|---|---:|---:|---:|"]
    for row in r["threads"]:
        out.append(f"| {row['setting']} | {row['wall_seconds']:.3f} | {row['cpu_seconds']:.3f} | "
                   f"{row['cpu_seconds'] / row['wall_seconds']:.2f} |")
    out += ["", "## 7. Measured: independent processes started together (single-thread settings)", "",
            "Wall time runs from the start of the first process to the last result, so it includes interpreter "
            "start-up and imports; the median task time is each worker's own timing, which excludes them.", "",
            "| workload | processes | wall s | runs/s | median task wall s |", "|---|---:|---:|---:|---:|"]
    for name, rows in r["parallel"].items():
        for row in rows:
            out.append(f"| {name} | {row['processes']} | {row['wall_seconds']:.3f} | {row['runs_per_second']:.2f} | "
                       f"{row['median_task_wall_seconds']:.3f} |")
    p = r["projections"]
    out += ["", "## 8. Projections (NOT measured)", "", p["label"] + ".", "",
            f"- Layer A: {p['layer_a_seconds_per_nucleus_S2_horizon_12']:.3e} s per nucleus, hence about "
            f"{p['layer_a_seconds_for_1e9_nuclei']:.0f} s for 1e9 nuclei (memory would not fit in one piece; "
            "Layer A is chunked over nuclei but keeps every event).",
            f"- B2: {p['b2_seconds_per_interval_S2_n_1e6']:.3e} s per interval at 1e6 nuclei, hence about "
            f"{p['b2_seconds_for_1e6_intervals']:.0f} s for 1e6 intervals.", ""]
    (HERE / "results" / "tables.md").write_text("\n".join(out))
    print(f"wrote {HERE / 'results' / 'tables.md'}")


if __name__ == "__main__":
    main()
