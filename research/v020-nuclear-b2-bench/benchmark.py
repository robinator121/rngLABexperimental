"""Single-machine benchmarks of Layer A and B2: measured, not projected (projections are labelled).

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-bench/benchmark.py

Every configuration runs in fresh Python processes, so peak resident memory
(ru_maxrss), CPU time and wall time belong to that configuration alone. Times
come from three untraced processes (median reported); Python-level allocation
peaks come from one further process run under tracemalloc, whose times are
discarded because tracing slows allocation-heavy code. Random numbers come from the
engineering namespace rnglab-v020-nuclear-bench. Writes results/benchmark.json.
Nothing here is an HPC measurement.

Also measured: an import-only baseline (the resident memory of the interpreter
and libraries, to read peak RSS against), the B2 time under the default BLAS /
OpenMP thread settings and with one thread, and independent processes started
together with one thread each (as array jobs would run them).
"""
from __future__ import annotations

import datetime
import json
import os
import platform
import statistics
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
NAMESPACE = "rnglab-v020-nuclear-bench"
S2_TIMES = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]
THREAD_VARIABLES = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")


def child_env(threads: str | None = None) -> dict:
    env = dict(os.environ, PYTHONPATH=str(REPO / "src"))
    for name in THREAD_VARIABLES:
        env.pop(name, None)
        if threads is not None:
            env[name] = threads
    return env


def worker(config: dict) -> dict:
    """One measured configuration (runs in its own process)."""
    import resource
    import tracemalloc

    import numpy as np

    from rnglab.nuclear import DecayScheme, WordStream, labelled_source, simulate_decays
    from rnglab.nuclear.aggregate import AggregatePlan, simulate_aggregate

    protocol = json.loads((REPO / "research/v020-nuclear-a-verify/protocol.json").read_text())
    s2 = DecayScheme.from_dict(protocol["schemes"]["S2"])
    stream = WordStream(labelled_source(NAMESPACE, config["engine"], config["n"], config.get("intervals", 6),
                                        config.get("repeat", 0)))
    traced = bool(config.get("trace"))
    if traced:
        tracemalloc.start()
    cpu0, wall0 = time.process_time(), time.perf_counter()
    out: dict = {}
    if config["engine"] == "baseline":
        out.update(parts=None, imported=["numpy", "scipy.special", "rnglab.nuclear", "rnglab.nuclear.aggregate"])
    elif config["engine"] == "layer_a":
        run = simulate_decays(s2, {"Xa-400": config["n"]}, 12.0, stream)
        pops = run.populations(S2_TIMES)
        out.update(events=run.n_events, output_bytes=int(sum(v.nbytes for v in run.events.values())),
                   populations_sum=int(pops[-1].sum()))
    else:
        intervals = config.get("intervals")
        if config.get("distinct"):          # every interval length differs: one set of laws per interval
            steps = np.arange(1, intervals + 1, dtype=np.float64)
            times = (np.concatenate([[0.0], np.cumsum(steps)]) * (12.0 / steps.sum())).tolist()
        else:
            times = [12.0 * k / intervals for k in range(intervals + 1)] if config.get("uniform") else S2_TIMES
        plan_start = time.perf_counter()
        plan = AggregatePlan(s2, times, config["n"])
        out["plan_seconds"] = time.perf_counter() - plan_start
        every = config.get("checkpoint_every", 0)
        run = simulate_aggregate(s2, {"Xa-400": config["n"]}, times, stream, plan=plan,
                                 keep_outputs=config.get("keep", True), checkpoint_every=every,
                                 on_checkpoint=(lambda checkpoint: None) if every else None)
        out.update(parts=plan.parts, words_per_interval=plan.words_per_interval,
                   output_bytes=int(sum(getattr(run, k).nbytes for k in ("populations", "path_counts", "decays",
                                                                       "photons", "neutrons"))))
    wall, cpu = time.perf_counter() - wall0, time.process_time() - cpu0
    traced_peak = None
    if traced:
        _, traced_peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
    out.update(wall_seconds=wall, cpu_seconds=cpu, words=stream.words_read, bytes=stream.bytes_read,
               traced_peak_bytes=traced_peak,
               max_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024)
    return out


DETERMINISTIC = ("words", "bytes", "parts", "words_per_interval", "imported")  # fixed by the configuration
DETERMINISTIC_BY_ENGINE = {"layer_a": ("populations_sum",),   # conserved; its output_bytes follow the events
                           "b2": ("output_bytes",)}           # array sizes fixed by the grid and the scheme


def fixed_fields(config: dict) -> tuple[str, ...]:
    """Fields that every repeat of a configuration must reproduce exactly."""
    return DETERMINISTIC + DETERMINISTIC_BY_ENGINE.get(config.get("engine"), ())


def measure(config: dict, repeats: int = 3, threads: str | None = None) -> dict:
    samples = []
    for repeat in range(repeats):
        cfg = dict(config, repeat=repeat)
        done = subprocess.run([sys.executable, __file__, "worker", json.dumps(cfg)], capture_output=True, text=True,
                              env=child_env(threads), check=True)
        samples.append(json.loads(done.stdout.strip().splitlines()[-1]))
    timed = [k for k in ("wall_seconds", "cpu_seconds", "plan_seconds", "max_rss_bytes") if k in samples[0]]
    summary = {k: statistics.median(s[k] for s in samples) for k in timed}  # every timing is a median
    summary["samples"] = {k: [s[k] for s in samples] for k in timed}
    for key in samples[0]:
        if key in timed or key in ("repeat", "traced_peak_bytes"):
            continue
        values = [s.get(key) for s in samples]
        if all(value == values[0] for value in values):
            summary[key] = values[0]
        elif key in fixed_fields(config):   # fixed by the configuration: a difference means another computation
            raise RuntimeError(f"repeats disagree on {key} ({values}); they must describe the same computation")
        else:                           # stochastic (each repeat has its own seed): median, samples kept
            summary[key] = statistics.median(values)
            summary["samples"][key] = values
    traced = subprocess.run([sys.executable, __file__, "worker", json.dumps(dict(config, repeat=50, trace=True))],
                            capture_output=True, text=True, env=child_env(threads), check=True)
    summary["traced_peak_bytes"] = json.loads(traced.stdout.strip().splitlines()[-1])["traced_peak_bytes"]
    summary["cpu_utilisation"] = summary["cpu_seconds"] / summary["wall_seconds"]
    summary["words_per_second"] = summary["words"] / summary["wall_seconds"]
    summary["thread_setting"] = "default (unset)" if threads is None else f"{threads} thread(s)"
    return {"config": config, **summary}


def parallel(config: dict, processes: int) -> dict:
    """Wall time for `processes` independent copies started together, one BLAS/OpenMP thread each."""
    start = time.perf_counter()
    jobs = [subprocess.Popen([sys.executable, __file__, "worker", json.dumps(dict(config, repeat=100 + k))],
                             stdout=subprocess.PIPE, text=True, env=child_env("1"))
            for k in range(processes)]
    results = [json.loads(job.communicate()[0].strip().splitlines()[-1]) for job in jobs]
    wall = time.perf_counter() - start
    return {"processes": processes, "wall_seconds": wall, "runs_per_second": processes / wall,
            "median_task_wall_seconds": statistics.median(r["wall_seconds"] for r in results),
            "timing_note": "wall_seconds spans process start to the last result, so it includes interpreter "
                           "start-up and imports; median_task_wall_seconds is each worker's own timing, "
                           "which excludes them"}


def machine() -> dict:
    info = {"platform": platform.platform(), "python": sys.version.split()[0], "logical_cpus": os.cpu_count()}
    try:
        text = Path("/proc/cpuinfo").read_text()
        info["cpu_model"] = next(line.split(":", 1)[1].strip() for line in text.splitlines()
                                 if line.startswith("model name"))
        info["memory_total_bytes"] = int(next(line.split()[1] for line in Path("/proc/meminfo").read_text().splitlines()
                                              if line.startswith("MemTotal"))) * 1024
    except (OSError, StopIteration):
        pass
    import numpy
    import scipy
    info.update(numpy=numpy.__version__, scipy=scipy.__version__)
    return info


def main() -> int:
    out = HERE / "results" / "benchmark.json"
    if out.exists():
        raise SystemExit("results/benchmark.json exists; results are never overwritten")
    sys.path.insert(0, str(REPO / "src"))
    from rnglab.nuclear.manifest import code_identity
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    clock = time.perf_counter()
    result = {"generated_by": "research/v020-nuclear-b2-bench/benchmark.py", "started_utc": started,
              "code": code_identity(REPO), "machine": machine(), "namespace": NAMESPACE,
              "note": "single workstation, measured; projections are in the 'projections' block only"}
    result["baseline"] = measure({"engine": "baseline", "n": 0})
    result["layer_a"] = [measure({"engine": "layer_a", "n": n}) for n in (10 ** 3, 10 ** 4, 10 ** 5, 10 ** 6, 2 * 10 ** 6)]
    print("layer A done", flush=True)
    result["b2_population"] = [measure({"engine": "b2", "n": n})
                               for n in (10 ** 3, 10 ** 5, 10 ** 6, 10 ** 9, 10 ** 12, 2 ** 50, 2 ** 53, 2 ** 56, 2 ** 60)]
    print("B2 population sweep done", flush=True)
    result["b2_intervals"] = []
    for intervals in (10, 100, 1000, 10000):
        for keep in (True, False):
            result["b2_intervals"].append(measure({"engine": "b2", "n": 10 ** 6, "intervals": intervals,
                                                   "uniform": True, "keep": keep}))
    print("B2 interval sweep done", flush=True)
    result["b2_checkpoints"] = [measure({"engine": "b2", "n": 10 ** 6, "intervals": 2000, "uniform": True,
                                         "keep": False, "checkpoint_every": every}) for every in (0, 100, 1)]
    result["b2_plan_lengths"] = [measure({"engine": "b2", "n": 10 ** 6, "intervals": 200, "keep": False, **grid})
                                 for grid in ({"uniform": True}, {"distinct": True})]
    print("checkpoint and plan-length configurations done", flush=True)
    threads_config = {"engine": "b2", "n": 10 ** 6, "intervals": 1000, "uniform": True, "keep": False}
    result["threads"] = [dict(measure(threads_config), setting="default (unset)"),
                         dict(measure(threads_config, threads="1"), setting="1 thread")]
    print("thread settings done", flush=True)
    result["parallel"] = {"layer_a_n_1e5": [parallel({"engine": "layer_a", "n": 10 ** 5}, k) for k in (1, 2, 4, 8)],
                          "b2_n_1e6_1000_intervals": [parallel({"engine": "b2", "n": 10 ** 6, "intervals": 1000,
                                                                "uniform": True, "keep": False}, k)
                                                       for k in (1, 2, 4, 8)]}
    print("parallel done", flush=True)
    per_nucleus = statistics.median(r["wall_seconds"] / r["config"]["n"] for r in result["layer_a"][2:])
    per_interval = statistics.median(r["wall_seconds"] / r["config"]["intervals"]
                                     for r in result["b2_intervals"] if not r["config"]["keep"])
    result["projections"] = {
        "label": "PROJECTION, not measured: linear extrapolation of the medians above, this machine, one process",
        "layer_a_seconds_per_nucleus_S2_horizon_12": per_nucleus,
        "layer_a_seconds_for_1e9_nuclei": per_nucleus * 1e9,
        "b2_seconds_per_interval_S2_n_1e6": per_interval,
        "b2_seconds_for_1e6_intervals": per_interval * 1e6,
    }
    result["elapsed_seconds"] = time.perf_counter() - clock
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(result, indent=1) + "\n")
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "worker":
        print(json.dumps(worker(json.loads(sys.argv[2]))))
    else:
        sys.exit(main())
