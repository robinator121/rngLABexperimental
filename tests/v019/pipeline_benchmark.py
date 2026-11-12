"""Timing of the Stage D generation path: sampler, checks I-2 and I-3, measurement, atomic storage.

Runs full-length runs through the pipeline's own run_one (research/v019-d-pilot/pilot_generate.py)
with 1 to N worker processes, in the audit-only namespace. The runs are timing material: their
seeds are outside the pilot and smoke namespaces, no analysis stage reads them, and the store is
deleted when the timing ends. Nothing here is a pilot, smoke or campaign result.

For each worker count the job list has the pilot's mix of conditions (1 calibration : 2 healthy :
1 per fault family) and starts at replica 1, so that the additional step check of replica 0 of the
fault families (3 of 5,999 pilot runs) is not over-represented; that check is timed separately.

Needs the measurement layer: RNGLAB_DPILOT_MEASURE=<path to measure.py>, or the file in the tree.

    PYTHONPATH=src python tests/v019/pipeline_benchmark.py --workers 1,4,6,8,12,16 --repeat 2

Writes docs/v020/evidence/pipeline_benchmark.json.
"""
from __future__ import annotations

import os

THREAD_VARIABLES = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
for _name in THREAD_VARIABLES:  # as the driver does, before NumPy is imported
    os.environ[_name] = "1"

import argparse  # noqa: E402
import importlib.util  # noqa: E402
import json  # noqa: E402
import platform  # noqa: E402
import resource  # noqa: E402
import shutil  # noqa: E402
import statistics  # noqa: E402
import subprocess  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import scipy  # noqa: E402

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("_pilot_fixture", HERE / "_pilot_fixture.py")
fx = importlib.util.module_from_spec(_spec)
sys.modules["_pilot_fixture"] = fx
_spec.loader.exec_module(fx)
au, pc, pg = fx.au, fx.pc, fx.pg

OUTPUT = au.REPO / "docs/v020/evidence/pipeline_benchmark.json"
DEFAULT_ROOT = au.REPO / ".devtmp/pipeline_benchmark"
CONDITIONS = (("calibration", "healthy"), ("healthy", "healthy"), ("fault/SHIFT", "SHIFT"), ("fault/LAG1", "LAG1"),
              ("fault/DUP", "DUP"))
PLAN = pc.Plan(fx.NAMESPACE, tuple((key, family, 64) for key, family in CONDITIONS), range(0, 1), range(1, 2), 64, 64, "test")
PILOT_COUNTS = {key: count for key, _, count in pc.FROZEN_PLAN.conditions}
STAGES = ("generation_seconds", "check_seconds", "measurement_seconds")
RSS_MARGIN = 1.25  # a worker count is skipped if workers x largest measured peak x margin exceeds MemAvailable


def jobs_for(count: int, first_replica: int = 1) -> list[tuple[str, str, int]]:
    """`count` jobs (a multiple of 6) in the pilot's mix, interleaved so that every shard gets every condition."""
    units = count // 6
    jobs = []
    for unit in range(units):
        block = [("calibration", "healthy", first_replica + unit), ("healthy", "healthy", first_replica + 2 * unit),
                 ("fault/SHIFT", "SHIFT", first_replica + unit), ("healthy", "healthy", first_replica + 2 * unit + 1),
                 ("fault/LAG1", "LAG1", first_replica + unit), ("fault/DUP", "DUP", first_replica + unit)]
        jobs += block[unit % 6:] + block[:unit % 6]  # rotated, so that a shard of 6 or 12 workers is not one condition
    return jobs


def memory_available_mib() -> float:
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemAvailable:"):
            return int(line.split()[1]) / 1024
    return float("nan")


def worker(root: Path, shard: str) -> int:
    started = time.perf_counter()
    context = fx.context(PLAN, root, need_measure=False)
    jobs = pg.shard([tuple(j) for j in pc.read_json(root / ".pending.json")], shard)
    loaded = time.perf_counter()
    rows = []
    engine = pg.load_engine("production", context)
    for condition, family, replica in jobs:
        t0 = time.perf_counter()
        entry = pg.run_one(context, engine, condition, family, replica)
        rows.append({"condition": condition, "replica": replica, "status": entry.get("status"),
                     "job_seconds": time.perf_counter() - t0})
    usage = resource.getrusage(resource.RUSAGE_SELF)
    report = {"shard": shard, "jobs": rows, "startup_seconds": loaded - started, "seconds": time.perf_counter() - started,
              "user_seconds": usage.ru_utime, "system_seconds": usage.ru_stime, "peak_rss_mib": usage.ru_maxrss / 1024,
              "voluntary_switches": usage.ru_nvcsw, "involuntary_switches": usage.ru_nivcsw}
    target = root / ".bench"
    target.mkdir(exist_ok=True)
    (target / f"worker_{shard.replace('/', '_of_')}.json").write_text(json.dumps(report), encoding="utf-8")
    return 0


def mean(values: list[float]) -> float | None:
    return float(statistics.fmean(values)) if values else None


def one_configuration(root: Path, workers: int, jobs: list[tuple[str, str, int]], *, time_scan: bool) -> dict:
    shutil.rmtree(root, ignore_errors=True)
    root.mkdir(parents=True)
    context = fx.context(PLAN, root, need_measure=False)
    pc.atomic_write(root / ".pending.json", pc.canonical_json([list(j) for j in jobs]), overwrite=True)
    before = {"load_average_1min": os.getloadavg()[0], "memory_available_mib": memory_available_mib()}
    environment = {**os.environ, **{name: "1" for name in THREAD_VARIABLES}}
    started = time.perf_counter()
    processes = [subprocess.Popen([sys.executable, str(Path(__file__).resolve()), "--worker", f"{k}/{workers}", "--root", str(root)],
                                  env=environment) for k in range(workers)]
    lowest_available = before["memory_available_mib"]
    while any(p.poll() is None for p in processes):
        time.sleep(0.5)
        lowest_available = min(lowest_available, memory_available_mib())
    wall = time.perf_counter() - started
    codes = [p.returncode for p in processes]
    reports = [json.loads(p.read_text(encoding="utf-8")) for p in sorted((root / ".bench").glob("worker_*.json"))]
    rows = []
    for report in reports:
        for job in report["jobs"]:
            record = pc.read_json(context.record_path(job["condition"], job["replica"]))
            trace = context.trace_path(job["condition"], job["replica"])
            stages = {name: float(record[name]) for name in STAGES}
            rows.append({**job, **stages, "storage_seconds": job["job_seconds"] - sum(stages.values()),
                         "record_bytes": context.record_path(job["condition"], job["replica"]).stat().st_size,
                         "trace_bytes": trace.stat().st_size if trace.exists() else 0,
                         "checks_passed": all(record[name]["passed"] for name in ("I-1", "I-2", "I-3", "I-3-step") if name in record)})
    per_condition = {}
    for key, _ in CONDITIONS:
        selected = [r for r in rows if r["condition"] == key]
        per_condition[key] = {"jobs": len(selected), **{name: mean([r[name] for r in selected]) for name in
                                                        ("job_seconds", *STAGES, "storage_seconds", "record_bytes", "trace_bytes")}}
    cpu = sum(r["user_seconds"] + r["system_seconds"] for r in reports)
    out = {"workers": workers, "jobs": len(rows), "worker_exit_codes": codes, "wall_seconds": wall,
           "runs_per_hour": 3600 * len(rows) / wall, "all_runs_succeeded": all(r["status"] == "SUCCESS" for r in rows) and len(rows) == len(jobs),
           "all_implementation_checks_passed": all(r["checks_passed"] for r in rows),
           "mean_job_seconds": mean([r["job_seconds"] for r in rows]),
           "mean_stage_seconds": {name: mean([r[name] for r in rows]) for name in (*STAGES, "storage_seconds")},
           "per_condition": per_condition,
           "cpu_seconds": cpu, "system_share_of_cpu": sum(r["system_seconds"] for r in reports) / cpu if cpu else None,
           "cpu_utilisation_of_workers": cpu / (wall * workers),
           "startup_seconds_mean": mean([r["startup_seconds"] for r in reports]),
           "idle_tail_share": 1 - sum(r["seconds"] for r in reports) / (wall * workers),
           "peak_rss_mib_max": max(r["peak_rss_mib"] for r in reports), "peak_rss_mib_sum": sum(r["peak_rss_mib"] for r in reports),
           "involuntary_switches_per_job": sum(r["involuntary_switches"] for r in reports) / max(len(rows), 1),
           "before": before, "lowest_memory_available_mib": lowest_available,
           "pilot_projection_hours": sum(PILOT_COUNTS[key] * per_condition[key]["job_seconds"] for key, _ in CONDITIONS
                                         if per_condition[key]["jobs"]) / workers / 3600
           if all(per_condition[key]["jobs"] for key, _ in CONDITIONS) else None,
           "pilot_projection_gib": sum(PILOT_COUNTS[key] * (per_condition[key]["record_bytes"] + per_condition[key]["trace_bytes"])
                                       for key, _ in CONDITIONS if per_condition[key]["jobs"]) / 2**30
           if all(per_condition[key]["jobs"] for key, _ in CONDITIONS) else None}
    if time_scan:  # what resume, the run manifest and the completeness gate cost per stored run
        t0 = time.perf_counter()
        states = [pg.classify(context, c, r, verify_hash=True)["state"] for c, _, r in jobs]
        seconds = time.perf_counter() - t0
        out["scan"] = {"runs": len(jobs), "all_complete": all(s == "COMPLETE" for s in states), "seconds_per_run": seconds / len(jobs),
                       "pilot_projection_seconds": 5999 * seconds / len(jobs)}
    shutil.rmtree(root, ignore_errors=True)
    return out


ALLOWANCE = 1.5  # planning allowance on the largest projection: sustained load, other use of the machine; not a measurement


def summarise(out: dict) -> dict:
    """Per worker count over the repetitions, and the figures the readiness report quotes."""
    by_workers: dict[int, list[dict]] = {}
    for row in out["configurations"]:
        by_workers.setdefault(row["workers"], []).append(row)
    alone = mean([r["mean_job_seconds"] for r in by_workers.get(1, [])])
    table = []
    for workers, rows in sorted(by_workers.items()):
        table.append({"workers": workers, "repetitions": len(rows), "runs_timed": sum(r["jobs"] for r in rows),
                      "runs_per_hour": [min(r["runs_per_hour"] for r in rows), max(r["runs_per_hour"] for r in rows)],
                      "pilot_hours": [min(r["pilot_projection_hours"] for r in rows), max(r["pilot_projection_hours"] for r in rows)],
                      "mean_job_seconds": mean([r["mean_job_seconds"] for r in rows]),
                      "slowdown_per_job": mean([r["mean_job_seconds"] for r in rows]) / alone if alone else None,
                      "cpu_utilisation_of_workers": mean([r["cpu_utilisation_of_workers"] for r in rows]),
                      "system_share_of_cpu": mean([r["system_share_of_cpu"] for r in rows]),
                      "involuntary_switches_per_job": mean([r["involuntary_switches_per_job"] for r in rows]),
                      "peak_rss_mib_max": max(r["peak_rss_mib_max"] for r in rows),
                      "peak_rss_mib_sum": max(r["peak_rss_mib_sum"] for r in rows),
                      "lowest_memory_available_mib": min(r["lowest_memory_available_mib"] for r in rows),
                      "all_runs_succeeded": all(r["all_runs_succeeded"] for r in rows),
                      "all_implementation_checks_passed": all(r["all_implementation_checks_passed"] for r in rows)})
    single = by_workers.get(1, [])
    stages = {name: mean([r["mean_stage_seconds"][name] for r in single]) for name in (*STAGES, "storage_seconds")} if single else {}
    total = sum(stages.values()) if stages else None
    per_condition = {key: {name: mean([r["per_condition"][key][name] for r in single])
                           for name in ("job_seconds", *STAGES, "storage_seconds", "record_bytes", "trace_bytes")}
                     for key, _ in CONDITIONS} if single else {}
    step = out.get("replica_zero_step_check", {}).get("per_condition", {})
    extra = {key: step[key]["check_seconds"] - per_condition[key]["check_seconds"]
             for key in step if step[key]["jobs"] and key in per_condition} if per_condition else {}
    candidates = [row for row in table if row["workers"] >= 6]
    hours = [h for row in candidates for h in row["pilot_hours"]]
    scans = [r["scan"]["seconds_per_run"] for r in out["configurations"] if "scan" in r]
    return {"per_worker_count": table,
            "single_process_stage_seconds": stages,
            "single_process_stage_share": {name: value / total for name, value in stages.items()} if total else {},
            "single_process_per_condition": per_condition,
            "step_check_extra_seconds_replica_zero": extra,
            "bytes_per_run": {key: per_condition[key]["record_bytes"] + per_condition[key]["trace_bytes"] for key in per_condition},
            "pilot_storage_gib": mean([r["pilot_projection_gib"] for r in out["configurations"]]),
            "pilot_hours_at_6_or_more_workers": [min(hours), max(hours)] if hours else None,
            "planning_allowance_factor": ALLOWANCE,
            "pilot_hours_planning_allowance": ALLOWANCE * max(hours) if hours else None,
            "scan_seconds_per_run": mean(scans), "scan_pilot_minutes": 5999 * mean(scans) / 60 if scans else None,
            "runs_timed_in_total": sum(r["jobs"] for r in out["configurations"]) + out.get("replica_zero_step_check", {}).get("jobs", 0)}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", default="1,4,6,8,12,16")
    parser.add_argument("--repeat", type=int, default=2)
    parser.add_argument("--jobs-per-worker", type=int, default=3)
    parser.add_argument("--root", default=str(DEFAULT_ROOT))
    parser.add_argument("--worker", default=None, help=argparse.SUPPRESS)
    parser.add_argument("--summarise", action="store_true", help="recompute the summary of the stored output; times nothing")
    args = parser.parse_args()
    root = Path(args.root)
    if au.measure_path() is None:
        print("the measurement layer is not available: set RNGLAB_DPILOT_MEASURE")
        return 2
    if args.worker:
        return worker(root, args.worker)
    if args.summarise:
        out = json.loads(OUTPUT.read_text(encoding="utf-8"))
        out["summary"] = summarise(out)
        OUTPUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(json.dumps(out["summary"], indent=1))
        return 0
    counts = [int(w) for w in args.workers.split(",")]
    usage = shutil.disk_usage(root.parent if root.parent.exists() else au.REPO)
    out = {"schema": "rnglab-v020-pipeline-benchmark-v1",
           "what": "timing of the generation path in the audit-only namespace; the runs are deleted and are not pilot, "
                   "smoke or campaign results",
           "namespace": fx.NAMESPACE, "sweeps_per_run": pc.H, "pilot_runs": pc.FROZEN_PLAN.total(), "pilot_counts": PILOT_COUNTS,
           "job_mix": "1 calibration : 2 healthy : 1 SHIFT : 1 LAG1 : 1 DUP, replicas from 1",
           "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                           "machine": platform.machine(), "logical_cpus": os.cpu_count(),
                           "memory_total_mib": os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / 2**20,
                           "disk_free_gib_at_start": usage.free / 2**30, "thread_variables": 1},
           "configurations": [], "skipped": []}
    largest_rss = 0.0
    for repetition in range(args.repeat):
        for workers in counts:
            available = memory_available_mib()
            if largest_rss and workers * largest_rss * RSS_MARGIN > available:
                out["skipped"].append({"workers": workers, "repetition": repetition, "memory_available_mib": available,
                                       "needed_mib": workers * largest_rss * RSS_MARGIN,
                                       "reason": "not enough available memory on the machine at that moment"})
                print(f"workers {workers}: skipped, {available:.0f} MiB available")
                continue
            count = 6 if workers == 1 else 6 * -(-args.jobs_per_worker * workers // 6)
            row = one_configuration(root, workers, jobs_for(count), time_scan=True)
            row["repetition"] = repetition
            largest_rss = max(largest_rss, row["peak_rss_mib_max"])
            out["configurations"].append(row)
            print(f"workers {workers:2d} rep {repetition}: {row['jobs']} runs in {row['wall_seconds']:.0f} s, "
                  f"{row['mean_job_seconds']:.1f} s per run, {row['runs_per_hour']:.0f} runs/h, pilot {row['pilot_projection_hours']:.2f} h, "
                  f"{row['pilot_projection_gib']:.1f} GiB, peak RSS {row['peak_rss_mib_max']:.0f} MiB", flush=True)
            OUTPUT.parent.mkdir(parents=True, exist_ok=True)
            OUTPUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    # replica 0 of each fault family: the step check of I-3 is added to check_seconds
    step = one_configuration(root, 1, [(key, family, 0) for key, family in CONDITIONS if family != "healthy"], time_scan=False)
    out["replica_zero_step_check"] = {"jobs": step["jobs"], "per_condition": step["per_condition"],
                                      "all_implementation_checks_passed": step["all_implementation_checks_passed"]}
    out["summary"] = summarise(out)
    OUTPUT.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote", OUTPUT.relative_to(au.REPO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
