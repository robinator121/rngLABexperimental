#!/usr/bin/env python3
"""Stage C pilot: severity range finding under the frozen grid rule (PROTOCOL.md section 3.4).

RNGLab v0.18 (written by Claude Code, 2026-11-06). Namespace
rnglab-v018-h1-detection-pilot; nothing here counts toward any endpoint.

    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python research/v018-h1-detection/run_pilot.py run --workers 8
    PYTHONPATH=src python research/v018-h1-detection/run_pilot.py summarize

`run` refuses to start unless protocol.json matches protocol.sha256 and is
frozen. Per-run records go to pilot/records/ (gitignored); completed cells are
skipped on a restart and never overwritten. `summarize` writes
pilot/summary.json: per cell the wrong-rate (the grid input) and, descriptively,
post-onset alarm rates of both arms at pilot-calibrated thresholds (190th of
199 healthy whole-run maxima).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.stats import beta

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import h1lib  # noqa: E402

PILOT = HERE / "pilot"
RECORDS = PILOT / "records"
CAL_CHUNK = 25
CODE_FILES = ["research/v018-h1-detection/protocol.json", "research/v018-h1-detection/h1lib.py",
              "research/v018-h1-detection/run_pilot.py", "src/rnglab14/models/pimc_harmonic.py",
              "src/rnglab14/sources/base.py", "src/rnglab14/sources/reference.py", "src/rnglab11/statistics.py",
              "src/rnglab12/uncertainty.py", "src/rnglab/seeds.py"]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_protocol() -> dict:
    if sha256(HERE / "protocol.json") != (HERE / "protocol.sha256").read_text().split()[0]:
        raise SystemExit("protocol.json does not match protocol.sha256; refusing to run")
    protocol = json.loads((HERE / "protocol.json").read_text())
    if protocol["status"] != "FROZEN_DESIGN":
        raise SystemExit("protocol is not frozen")
    return protocol


def pilot_cells(protocol: dict) -> list[dict]:
    pilot = protocol["data"]["pilot"]
    ranges = {"SHIFT": (1e-4, 0.1), "LAG1": (0.003, 0.5), "DUP": (0.001, 0.5)}  # protocol.json data.pilot.severities
    cells = [{"family": f, "index": k, "severity": float(s)}
             for f, (lo, hi) in ranges.items() for k, s in enumerate(np.geomspace(lo, hi, 10))]
    cells += [{"family": "GFSR", "index": k, "severity": list(pair)} for k, pair in enumerate(h1lib.GFSR_PAIRS)]
    for c in cells:
        c["replicas"] = pilot["replicas_per_severity"]
    return cells


def _write_jsonl(path: Path, records: list[dict]) -> None:
    tmp = path.with_suffix(".jsonl.tmp")
    tmp.write_text("".join(json.dumps(r) + "\n" for r in records))
    os.replace(tmp, path)


def run_cell(ns: str, cell: dict) -> str:
    path = RECORDS / f"{cell['family']}-{cell['index']:02d}.jsonl"
    if path.exists():
        return f"skip {path.name}"
    started = time.monotonic()
    severity = tuple(cell["severity"]) if cell["family"] == "GFSR" else cell["severity"]
    records = [h1lib.run_one(ns, "pilot", cell["family"], severity, i)[0] for i in range(cell["replicas"])]
    _write_jsonl(path, records)
    return f"{path.name}: {len(records)} runs in {time.monotonic() - started:.0f} s"


def run_calibration(ns: str, first: int, count: int) -> str:
    path = RECORDS / f"calibration-{first:03d}.jsonl"
    if path.exists():
        return f"skip {path.name}"
    started = time.monotonic()
    records = [h1lib.run_one(ns, "pilot-calibration", "healthy", 0, i)[0] for i in range(first, first + count)]
    _write_jsonl(path, records)
    return f"{path.name}: {count} runs in {time.monotonic() - started:.0f} s"


def run(workers: int) -> None:
    protocol = load_protocol()
    ns = protocol["data"]["pilot_namespace"]
    n_cal = protocol["data"]["pilot"]["healthy_replicas"]
    RECORDS.mkdir(parents=True, exist_ok=True)
    manifest = PILOT / "run_manifest.json"
    if not manifest.exists():
        manifest.write_text(json.dumps({"protocol_sha256": sha256(HERE / "protocol.json"),
                                        "code_sha256": {f: sha256(h1lib.REPO / f) for f in CODE_FILES},
                                        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                        "workers": workers, "python": platform.python_version(),
                                        "numpy": np.__version__}, indent=1) + "\n")
    context = multiprocessing.get_context("fork")
    with ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
        futures = [pool.submit(run_cell, ns, c) for c in pilot_cells(protocol)]
        futures += [pool.submit(run_calibration, ns, s, min(CAL_CHUNK, n_cal - s)) for s in range(0, n_cal, CAL_CHUNK)]
        for future in as_completed(futures):
            print(time.strftime("%H:%M:%S"), future.result(), flush=True)


def read(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines()]


def cp95(k: int, n: int) -> list[float]:
    return [0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1)), 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))]


def summarize() -> dict:
    protocol = load_protocol()
    order = protocol["data"]["pilot"]["threshold_order_statistic"]
    calibration = [r for p in sorted(RECORDS.glob("calibration-*.jsonl")) for r in read(p)]
    if len(calibration) != protocol["data"]["pilot"]["healthy_replicas"]:
        raise SystemExit("pilot calibration incomplete")
    constants = h1lib.application_constants(calibration)
    tau = {arm: h1lib.threshold([max(s for _, s in h1lib.arm_series(r, arm, constants)) for r in calibration], order)
           for arm in ("C", "A")}
    cal_wrong = sum(r["harm"]["wrong"] for r in calibration)
    out = {"status": "pilot; descriptive only; excluded from all endpoints", "thresholds": tau, "application_constants": constants,
           "calibration": {"runs": len(calibration), "wrong": cal_wrong, "wrong_cp95": cp95(cal_wrong, len(calibration))},
           "cells": []}
    for cell in pilot_cells(protocol):
        records = read(RECORDS / f"{cell['family']}-{cell['index']:02d}.jsonl")
        wrong = sum(r["harm"]["wrong"] for r in records)
        alarms = {arm: sum(h1lib.detection(h1lib.arm_series(r, arm, constants), tau[arm], r["onset_sweep"])["post_onset_alarm"]
                           for r in records) for arm in ("C", "A")}
        out["cells"].append({**cell, "runs": len(records), "wrong": wrong, "wrong_rate": wrong / len(records),
                             "wrong_cp95": cp95(wrong, len(records)),
                             "conventional_alarms": alarms["C"], "application_alarms": alarms["A"],
                             "onset_sweeps": [min(r["onset_sweep"] for r in records), max(r["onset_sweep"] for r in records)]})
    (PILOT / "summary.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("command", choices=["run", "summarize"])
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.command == "run":
        if os.environ.get("OPENBLAS_NUM_THREADS") != "1" or os.environ.get("OMP_NUM_THREADS") != "1":
            raise SystemExit("set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1")
        run(args.workers)
    else:
        s = summarize()
        for c in s["cells"]:
            print(c["family"], c["severity"], f"wrong {c['wrong']}/{c['runs']}", f"C {c['conventional_alarms']}", f"A {c['application_alarms']}")


if __name__ == "__main__":
    main()
