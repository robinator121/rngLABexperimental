#!/usr/bin/env python3
"""Generate the Stage C (H1) campaign (frozen protocol rnglab-v018-h1-detection-v1).

RNGLab v0.18 (written by Claude Code, 2026-11-06).

    PYTHONPATH=src python research/v018-h1-detection/run_campaign.py manifest
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python research/v018-h1-detection/run_campaign.py run --workers 8
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python research/v018-h1-detection/run_campaign.py replay

`run` refuses to start unless protocol.json matches protocol.sha256 and every
file hash in manifest_pregeneration.json (including grid.json) matches the
working tree. Each task is one (cell, block of 64 runs); it writes
traces/<cell>/block-XXXX.npz and .json (gitignored) through temporary files.
Completed blocks are skipped on a restart and never overwritten. When all
blocks exist, `run` merges them into results/cells/<cell>.npz and .jsonl and
results/run_manifest.json (committed). `replay` checks I-3 and I-4 on replica
0 of every cell and writes results/replay.json.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import multiprocessing
import os
import platform
import struct
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))  # appended, not prepended: other campaigns also have a module named run_campaign
import h1lib  # noqa: E402
from rnglab14.sources import PCG64Source  # noqa: E402

REPO = h1lib.REPO
BLOCK = 64  # protocol.json implementation_notes_added_at_freeze.storage
HASHED_FILES = [
    "research/v018-h1-detection/protocol.json", "research/v018-h1-detection/grid.json",
    "research/v018-h1-detection/h1lib.py", "research/v018-h1-detection/run_campaign.py",
    "research/v018-h1-detection/analyze.py", "research/v018-h1-detection/recount.py",
    "src/rnglab14/models/pimc_harmonic.py", "src/rnglab14/sources/__init__.py", "src/rnglab14/sources/base.py",
    "src/rnglab14/sources/reference.py", "src/rnglab11/statistics.py", "src/rnglab12/uncertainty.py",
    "src/rnglab/seeds.py", "research/phase2-20261030/environment/requirements-lock.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array_sha256(a) -> str:
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def load_protocol() -> dict:
    if sha256(HERE / "protocol.json") != (HERE / "protocol.sha256").read_text().split()[0]:
        raise SystemExit("protocol.json does not match protocol.sha256; refusing to run")
    protocol = json.loads((HERE / "protocol.json").read_text())
    if protocol["status"] != "FROZEN_DESIGN":
        raise SystemExit("protocol is not frozen")
    return protocol


def campaign_plan(protocol: dict, grid: dict) -> list[dict]:
    data = protocol["data"]
    cells = [{"cell": "calibration", "condition": "calibration", "family": "healthy", "severity": 0,
              "replicas": data["calibration"]["replicas"]},
             {"cell": "healthy", "condition": "test", "family": "healthy", "severity": 0,
              "replicas": data["test"]["healthy_replicas"]}]
    for family in ("SHIFT", "LAG1", "DUP", "GFSR"):
        for k, severity in enumerate(grid["families"][family]["severities"]):
            cells.append({"cell": f"{family}-{k}", "condition": "test", "family": family, "severity": severity,
                          "replicas": data["test"]["fault_cell_replicas"]})
    return cells


def severity_arg(cell: dict):
    return tuple(cell["severity"]) if cell["family"] == "GFSR" else cell["severity"]


def block_paths(out: Path, cell: str, block: int) -> tuple[Path, Path]:
    stem = out / "traces" / cell / f"block-{block:04d}"
    return stem.with_suffix(".npz"), stem.with_suffix(".json")


def pack_run(record: dict, arrays: dict) -> tuple[dict, dict]:
    """Split one run into scalar fields and fixed-shape arrays."""
    marks = arrays["marks"]
    conv, app = record.pop("conventional"), record.pop("application")
    packed = {
        "CV": arrays["CV"], "TD": arrays["TD"], "rg2": arrays["rg2"], "marks": marks,
        "conv_p": np.array([[np.nan if c[k] is None else c[k] for k in h1lib.CONV_P_KEYS] for c in conv]),
        "conv_end": np.array([c["end"] for c in conv]),
        "conv_bytes": np.array([marks[c["end"] - 1] for c in conv]),
        "app_stats": np.array([[c["mean_D"], c["uphill_fraction"], c["mean_rg2"]] for c in app]),
        "app_end": np.array([c["end"] for c in app]),
        "app_bytes": np.array([marks[c["end"] - 1] for c in app]),
    }
    record["cv_trace_sha256"] = array_sha256(arrays["CV"])
    return record, packed


def run_block(ns: str, cell: dict, block: int, out: Path, block_size: int = BLOCK, sweeps: int = h1lib.SWEEPS) -> str:
    npz, meta = block_paths(out, cell["cell"], block)
    if npz.exists() and meta.exists():
        return f"skip {cell['cell']} block {block} (complete)"
    npz.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    indices = list(range(block * block_size, min((block + 1) * block_size, cell["replicas"])))
    records, stacks = [], {}
    for i in indices:
        record, arrays = h1lib.run_one(ns, cell["condition"], cell["family"], severity_arg(cell), i, sweeps=sweeps)
        record, packed = pack_run(record, arrays)
        records.append(record)
        for key, value in packed.items():
            stacks.setdefault(key, []).append(value)
    tmp = npz.with_suffix(".npz.tmp")
    with tmp.open("wb") as handle:
        np.savez_compressed(handle, replica=np.asarray(indices), **{k: np.asarray(v) for k, v in stacks.items()})
    os.replace(tmp, npz)
    meta_tmp = meta.with_suffix(".json.tmp")
    meta_tmp.write_text(json.dumps({"cell": cell, "block": block, "npz_sha256": sha256(npz),
                                    "seconds": time.monotonic() - started, "records": records}) + "\n")
    os.replace(meta_tmp, meta)
    return f"{cell['cell']} block {block}: {len(indices)} runs in {time.monotonic() - started:.0f} s"


COMMITTED_ARRAYS = ("replica", "conv_p", "conv_end", "conv_bytes", "app_stats", "app_end", "app_bytes")


def merge(out: Path, plan: list[dict], block_size: int = BLOCK) -> dict:
    blocks = []
    (out / "results" / "cells").mkdir(parents=True, exist_ok=True)
    for cell in plan:
        records, parts = [], {k: [] for k in COMMITTED_ARRAYS}
        for block in range(-(-cell["replicas"] // block_size)):
            npz, meta = block_paths(out, cell["cell"], block)
            info = json.loads(meta.read_text())
            if info["npz_sha256"] != sha256(npz):
                raise SystemExit(f"{npz} does not match its block record")
            records += info["records"]
            with np.load(npz) as z:
                for k in COMMITTED_ARRAYS:
                    parts[k].append(z[k])
            blocks.append({"cell": cell["cell"], "block": block, "npz_sha256": info["npz_sha256"],
                           "json_sha256": sha256(meta), "seconds": info["seconds"]})
        if [r["replica"] for r in records] != list(range(cell["replicas"])):
            raise SystemExit(f"{cell['cell']}: runs incomplete or out of order")
        base = out / "results" / "cells" / cell["cell"]
        np.savez_compressed(base.with_suffix(".npz"), **{k: np.concatenate(v) for k, v in parts.items()})
        base.with_suffix(".jsonl").write_text("".join(json.dumps(r) + "\n" for r in records))
    return {"blocks": blocks}


def _task(ns, cell, block, out, block_size, sweeps):
    return run_block(ns, cell, block, Path(out), block_size, sweeps)


def run_all(workers: int, out: Path = HERE, ns: str | None = None, plan: list[dict] | None = None,
            block_size: int = BLOCK, sweeps: int = h1lib.SWEEPS, verify: bool = True) -> None:
    if verify:
        protocol = load_protocol()
        manifest = verify_manifest()
        ns = protocol["data"]["namespace"]
        plan = campaign_plan(protocol, json.loads((HERE / "grid.json").read_text()))
    log = out / "results" / "run.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    tasks = [(c, b) for c in sorted(plan, key=lambda c: c["cell"] == "calibration")
             for b in range(-(-c["replicas"] // block_size))]
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with log.open("a") as handle, ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context("fork")) as pool:
        handle.write(f"start {started} workers={workers}\n")
        futures = [pool.submit(_task, ns, c, b, str(out), block_size, sweeps) for c, b in tasks]
        for future in as_completed(futures):
            line = f"{time.strftime('%H:%M:%S')} {future.result()}"
            print(line, flush=True)
            handle.write(line + "\n")
            handle.flush()
        handle.write(f"end {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
    merged = merge(out, plan, block_size)
    if verify:
        import scipy
        record = {"protocol_sha256": sha256(HERE / "protocol.json"), "grid_sha256": sha256(HERE / "grid.json"),
                  "manifest_pregeneration_sha256": sha256(HERE / "manifest_pregeneration.json"),
                  "manifest_git_head": manifest["git_head_before_commit"],
                  "git_head_at_run": subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip(),
                  "first_start_utc": started, "workers": workers, "python": platform.python_version(),
                  "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform(), **merged,
                  "cell_sha256": {c["cell"]: {ext: sha256(out / "results" / "cells" / f"{c['cell']}.{ext}") for ext in ("npz", "jsonl")}
                                  for c in plan}}
        (out / "results" / "run_manifest.json").write_text(json.dumps(record, indent=1) + "\n")


# --------------------------------------------------------------------------
# I-3 and I-4 (PROTOCOL.md section 5.1)
# --------------------------------------------------------------------------
def expected_stream(family: str, severity, seeds: dict, n_bytes: int, onset: int = h1lib.ONSET_BYTES) -> bytes:
    """Independent re-application of the stated fault map to the healthy base words."""
    base = PCG64Source(seeds["local"])
    out = bytearray(base.read_bytes(onset))
    prev, z = out[-8:], None
    rng = np.random.default_rng(seeds["fault"])
    gfsr = h1lib.GFSR(*severity, seeds["gfsr"]) if family == "GFSR" else None
    while len(out) < n_bytes:
        if family == "GFSR":
            word = gfsr.word()
        else:
            raw = base.read_bytes(8)
            w = struct.unpack(">Q", raw)[0]
            u = (w >> 11) * 2.0**-53
            if family == "DUP":
                word = bytes(prev) if rng.random() < severity else raw
            elif family == "SHIFT":
                word = struct.pack(">Q", (int(min(u + severity, h1lib.TOP) * 2.0**53) << 11) | (w & 0x7FF))
            else:  # LAG1
                e = h1lib.STD.inv_cdf(u + 2.0**-54)
                z = e if z is None else severity * z + math.sqrt(1 - severity * severity) * e
                word = struct.pack(">Q", (int(min(h1lib.STD.cdf(z), h1lib.TOP) * 2.0**53) << 11) | (w & 0x7FF))
        out += word
        prev = word
    return bytes(out)


def replay(out: Path = HERE, ns: str | None = None, plan: list[dict] | None = None, sweeps: int = h1lib.SWEEPS) -> dict:
    if plan is None:
        protocol = load_protocol()
        ns = protocol["data"]["namespace"]
        plan = campaign_plan(protocol, json.loads((HERE / "grid.json").read_text()))
    result = {}
    for cell in plan:
        stored = json.loads((out / "results" / "cells" / f"{cell['cell']}.jsonl").read_text().splitlines()[0])
        record, arrays = h1lib.run_one(ns, cell["condition"], cell["family"], severity_arg(cell), 0, sweeps=sweeps, return_stream=True)
        entry = {"I-4": record["local_sha256"] == stored["local_sha256"] and array_sha256(arrays["CV"]) == stored["cv_trace_sha256"]}
        if cell["family"] != "healthy":
            entry["I-3"] = arrays["stream"] == expected_stream(cell["family"], severity_arg(cell), record["seeds"], len(arrays["stream"]))
        result[cell["cell"]] = entry
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "results" / "replay.json").write_text(json.dumps(result, indent=1) + "\n")
    return result


# --------------------------------------------------------------------------
# Pre-generation manifest (PROTOCOL.md section 0, step 4)
# --------------------------------------------------------------------------
def write_manifest() -> None:
    protocol = load_protocol()
    path = HERE / "manifest_pregeneration.json"
    if path.exists():
        raise SystemExit(f"{path} exists; a new manifest needs a new protocol version")
    if (HERE / "traces").exists():
        raise SystemExit("traces/ already exists; the manifest must precede every campaign draw")
    import scipy

    grid = json.loads((HERE / "grid.json").read_text())
    plan = campaign_plan(protocol, grid)
    ns = protocol["data"]["namespace"]
    table = [f"{ns}/{c['condition']}/{role}/{i:05d}\t{seed}"
             for c in plan if c["cell"] in ("calibration", "healthy")
             for i in range(c["replicas"]) for role, seed in h1lib.seeds_for(ns, c["condition"], i).items()]
    seeds = [line.split("\t")[1] for line in table]
    if len(set(seeds)) != len(seeds):
        raise SystemExit("seed collision inside the namespace")
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    manifest = {
        "protocol_id": protocol["id"], "protocol_sha256": sha256(HERE / "protocol.json"),
        "grid_sha256": sha256(HERE / "grid.json"), "grid": grid["families"],
        "cells": [{k: c[k] for k in ("cell", "condition", "family", "severity", "replicas")} for c in plan],
        "total_runs": sum(c["replicas"] for c in plan),
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "git_head_before_commit": head,
        "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(), "file_sha256": {f: sha256(REPO / f) for f in HASHED_FILES},
        "seed_units": len(table), "seed_collisions": 0,
        "seed_table_note": "fault cells reuse the 'test' seeds of the same replica index (pairing); calibration seeds are disjoint",
        "seed_table_sha256": hashlib.sha256("\n".join(table).encode()).hexdigest(),
        "seed_table_first_rows": table[:5],
        "statement": "Written before any campaign draw; run_campaign.py run verifies these hashes."}
    path.write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"wrote {path} ({manifest['total_runs']} runs, {len(table)} seed units)")


def verify_manifest() -> dict:
    path = HERE / "manifest_pregeneration.json"
    if not path.exists():
        raise SystemExit("no manifest_pregeneration.json; run the manifest command and commit it first")
    manifest = json.loads(path.read_text())
    for name, digest in manifest["file_sha256"].items():
        if sha256(REPO / name) != digest:
            raise SystemExit(f"{name} changed since the pre-generation manifest; refusing to run")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("command", choices=["manifest", "run", "replay"])
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    if args.command == "manifest":
        write_manifest()
        return
    if os.environ.get("OPENBLAS_NUM_THREADS") != "1" or os.environ.get("OMP_NUM_THREADS") != "1":
        raise SystemExit("set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1")
    if args.command == "run":
        run_all(args.workers)
    else:
        print(json.dumps(replay(), indent=1))


if __name__ == "__main__":
    main()
