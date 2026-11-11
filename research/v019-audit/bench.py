"""Benchmark of full-length Stage D runs in the audit-only namespace. Timings only; no pilot data.

    PYTHONPATH=src python research/v019-audit/bench.py [--workers 6] [--measure path/to/measure.py] --out benchmark.json

Per family, one full-length run (40,960 sweeps) through each stage, timed alone:

    generate          pilot driver (dlib.generate_run): sampler, fault adapter, sweep clock
    check_i2, check_i3  the pilot's implementation checks
    vector_transform  the independent vectorised transform that I-3 could use instead
    measure, write    the measurement layer, when present
    replay            independent regeneration of the run (research/v019-audit/replay.py)
    verify            stored-output verification with regeneration, every eighth checkpoint

Then `--workers` identical jobs run at once to measure parallel efficiency and
peak memory per process, and a profile of a 4,096-sweep generation names the
functions where the time goes. Every result is reported; nothing is tuned here.
"""
from __future__ import annotations

import argparse
import cProfile
import importlib.util
import json
import os
import platform
import pstats
import resource
import sys
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

_HERE = Path(__file__).resolve().parent
REPO = _HERE.parents[1]
NAMESPACE = "rnglab-v019-audit-smoke"
FAMILIES = ("healthy", "SHIFT", "LAG1", "DUP")


def _load(name: str, path: Path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def _modules(measure_path: str | None):
    verify = _load("v019_audit_verify", _HERE / "verify.py")
    if measure_path:
        measure = _load("v019_bench_measure", Path(measure_path))
        return verify, measure.dlib, measure
    return verify, _load("v019_bench_dlib", REPO / "research/v019-d-pilot/dlib.py"), None


def one_run(job: tuple[str, int, str | None, bool]) -> dict:
    family, replica, measure_path, stages = job
    verify, dlib, measure = _modules(measure_path)
    rc, rp = verify.rc, verify.rp
    condition = "healthy" if family == "healthy" else "fault/" + family
    out: dict = {"family": family, "replica": replica, "pid": os.getpid()}

    def timed(name, call):
        start = time.perf_counter()
        value = call()
        out[name + "_seconds"] = time.perf_counter() - start
        return value

    run = timed("generate", lambda: dlib.generate_run(NAMESPACE, condition, family, replica))
    out["rss_after_generate_kb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    out["local_bytes"] = int(run["local_bytes"])
    if measure is not None:
        record, arrays = timed("measure", lambda: measure.measure_run(run))
        out["rss_after_measure_kb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        if stages:
            with tempfile.TemporaryDirectory(dir=os.environ.get("RNGLAB_BENCH_TMP")) as tmp:
                entry = timed("write", lambda: measure.write_run(record, arrays, Path(tmp)))
                sizes = {name: (Path(tmp) / name).stat().st_size for name in entry["files"]}
                out["trace_file_bytes"] = next(v for k, v in sizes.items() if k.endswith(".npz"))
                out["record_file_bytes"] = next(v for k, v in sizes.items() if k.endswith(".json"))
                stored_record, stored_arrays = verify.load_stored(Path(tmp), condition, replica)
                result = timed("verify", lambda: verify.verify_run(stored_record, stored_arrays, NAMESPACE,
                                                                  deep=True, checkpoint_stride=8))
                out["verify_comparisons"], out["verify_discrepancies"] = result["compared"], len(result["issues"])
    if stages:
        out["I-2"] = timed("check_i2", lambda: dlib.check_i2(run))["passed"]
        out["I-3"] = timed("check_i3", lambda: dlib.check_i3(run))["passed"]
        profile = run["profile"]

        def vector():
            base = rc.pcg64_words(run["seeds"]["local"], run["local_bytes"] // 8)
            emitted = rc.apply_fault(family, base, rc.sweep_of_word(run["marks"]),
                                     rc.severity_path(profile["t0"], profile["T_cap"], profile["s_max"]),
                                     profile["t0"], run["seeds"]["fault"])
            return rc.words_to_bytes(base) == run["raw_base"] and rc.words_to_bytes(emitted) == run["emitted"]

        out["vector_transform_equal"] = timed("vector_transform", vector)
        mine = timed("replay", lambda: rp.replay(NAMESPACE, condition, family, replica))
        out["replay_bit_identical"] = bool(
            all(np.array_equal(mine[k], run[k]) for k in ("marks", "CV", "TD", "severity"))
            and all(mine[k] == run[k] for k in ("emitted_sha256", "base_sha256", "prefix_sha256", "start_sha256")))
    out["rss_peak_kb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return out


def profile_generation(sweeps: int = 4096) -> list[dict]:
    dlib = _load("v019_bench_dlib_profile", REPO / "research/v019-d-pilot/dlib.py")
    profiler = cProfile.Profile()
    profiler.enable()
    dlib.generate_run(NAMESPACE, "fault/LAG1", "LAG1", 1, sweeps=sweeps,
                      profile={"t0": 64, "T_cap": np.float64(1024.0), "r": np.float64(0.5 / 1024.0), "s_max": np.float64(0.5)})
    profiler.disable()
    stats = pstats.Stats(profiler)
    total = sum(row[2] for row in stats.stats.values())
    rows = sorted(stats.stats.items(), key=lambda item: -item[1][2])[:12]
    return [{"function": f"{Path(key[0]).name}:{key[2]}", "calls": value[1], "own_seconds": round(value[2], 4),
             "share_of_profiled_time": round(value[2] / total, 4)} for key, value in rows]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--workers", type=int, default=6)
    parser.add_argument("--measure", default=None, help="path of the measurement layer; default: in-tree if present")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    in_tree = REPO / "research/v019-d-pilot/measure.py"
    measure_path = args.measure or (str(in_tree) if in_tree.is_file() else None)
    import hashlib
    import scipy
    report = {"schema": "rnglab-v019-audit-benchmark-v1", "namespace": NAMESPACE,
              "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                              "platform": platform.platform(), "cpu_count": os.cpu_count(),
                              "load_average_before": os.getloadavg(),
                              "thread_variables": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS")}},
              # The file is identified by its hash; a local path would not help a reader and is not recorded.
              "measure_layer": None if measure_path is None else {
                  "file": Path(measure_path).name, "in_tree": Path(measure_path).resolve() == in_tree.resolve(),
                  "sha256": hashlib.sha256(Path(measure_path).read_bytes()).hexdigest()},
              "dlib_sha256": hashlib.sha256((REPO / "research/v019-d-pilot/dlib.py").read_bytes()).hexdigest()}
    with ProcessPoolExecutor(1) as pool:  # one fresh process per family, alone on the machine
        report["single"] = [pool.submit(one_run, (family, 20, measure_path, True)).result() for family in FAMILIES]
    start = time.perf_counter()
    jobs = [(FAMILIES[i % 4], 30 + i, measure_path, False) for i in range(args.workers)]
    with ProcessPoolExecutor(args.workers) as pool:
        parallel = list(pool.map(one_run, jobs, chunksize=1))
    wall = time.perf_counter() - start
    stage = lambda row: row["generate_seconds"] + row.get("measure_seconds", 0.0)
    alone = {row["family"]: stage(row) for row in report["single"]}
    report["parallel"] = {"workers": args.workers, "wall_seconds": wall, "jobs": parallel,
                          "worker_seconds": sum(stage(row) for row in parallel),
                          "slowdown_per_job": [stage(row) / alone[row["family"]] for row in parallel],
                          "max_rss_kb": max(row["rss_peak_kb"] for row in parallel)}
    report["profile_generation_4096_sweeps"] = profile_generation()
    report["environment"]["load_average_after"] = os.getloadavg()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for row in report["single"]:
        print(row["family"], {k: round(v, 2) for k, v in row.items() if k.endswith("_seconds")},
              "rss MiB", row["rss_peak_kb"] // 1024)
    print("parallel", args.workers, "workers: wall", round(wall, 1), "s; slowdown", [round(s, 2) for s in report["parallel"]["slowdown_per_job"]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
