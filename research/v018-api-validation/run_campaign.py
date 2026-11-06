#!/usr/bin/env python3
"""Generate the v0.18 Phase 2 API-validation campaign (frozen protocol rnglab-v018-api-validation-v1).

RNGLab v0.18 (written by Claude Code, 2026-11-05).

    PYTHONPATH=src python research/v018-api-validation/run_campaign.py manifest
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python research/v018-api-validation/run_campaign.py run --workers 4
    OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python research/v018-api-validation/run_campaign.py replay

`run` refuses to start unless protocol.json matches protocol.sha256 and every
file hash in manifest_pregeneration.json matches the working tree. Each task is
one (condition, block of 256 replicas); it writes traces/<condition>/block-XXXX.npz
and .json (gitignored) through temporary files. Completed blocks are skipped on a
restart and never overwritten. When every block exists, `run` merges them into
results/<condition>_replicas.csv and results/run_manifest.json (committed).
`replay` regenerates replica 0 of each condition (check I-4) into results/replay.json.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import multiprocessing
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))

from rnglab.seeds import namespace_seed  # noqa: E402
from rnglab12.uncertainty import summarize_series  # noqa: E402
from rnglab14.models.pimc_harmonic import (  # noqa: E402
    PathIntegralHarmonicOscillator, exact_finite_bead_energy, harmonic_energy_estimators,
    sample_exact_finite_bead_paths,
)
from rnglab14.sources import PCG64Source, StochasticChannels  # noqa: E402

M, BETA, DELTA, N, BATCH = 16, 2.0, 0.4, 8192, 512
CONDITIONS = ("local", "refresh", "iid")
BLOCK = 256  # protocol.json implementation_notes_added_at_freeze.storage
LEDGER_ROLES = ("proposal", "acceptance", "initialization", "resampling", "centroid")
START_SAMPLER = "sample_exact_finite_bead_paths(source, 16, 2.0, 2)[0]"
HASHED_FILES = [
    "research/v018-api-validation/protocol.json",
    "research/v018-api-validation/run_campaign.py",
    "research/v018-api-validation/analyze.py",
    "research/v018-api-validation/recount.py",
    "src/rnglab14/models/pimc_harmonic.py",
    "src/rnglab14/sources/__init__.py",
    "src/rnglab14/sources/base.py",
    "src/rnglab14/sources/reference.py",
    "src/rnglab12/uncertainty.py",
    "src/rnglab/longrun.py",
    "src/rnglab/seeds.py",
    "research/phase2-20261030/environment/requirements-lock.txt",
]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def array_sha256(a) -> str:
    """SHA-256 of little-endian float64 bytes in C order."""
    return hashlib.sha256(np.ascontiguousarray(a, dtype="<f8").tobytes()).hexdigest()


def load_protocol() -> dict:
    expected = (HERE / "protocol.sha256").read_text().split()[0]
    if sha256(HERE / "protocol.json") != expected:
        raise SystemExit("protocol.json does not match protocol.sha256; refusing to run")
    protocol = json.loads((HERE / "protocol.json").read_text())
    if protocol["status"] != "FROZEN_DESIGN":
        raise SystemExit("protocol is not frozen")
    model = protocol["model"]
    if (model["beads"], model["beta"], model["displacement_delta"], model["sweeps"], model["burn_in"]) != (M, BETA, DELTA, N, 0):
        raise SystemExit("script constants differ from the frozen protocol")
    if protocol["definitions"]["reference_value"] != exact_finite_bead_energy(M, BETA):
        raise SystemExit("reference value differs from exact_finite_bead_energy")
    return protocol


def seeds_for(ns: str, condition: str, i: int) -> dict:
    unit = f"{i:05d}"
    if condition == "iid":
        return {"paths": namespace_seed(ns, "iid", "paths", unit)}
    seeds = {role: namespace_seed(ns, "metropolis", role, unit) for role in ("start", "proposal", "acceptance")}
    if condition == "refresh":
        seeds["centroid"] = namespace_seed(ns, "refresh", "centroid", unit)
    return seeds


def _bm512_fields(prefix: str, summary: dict) -> dict:
    bm = summary["batch_means"][str(BATCH)]
    bounds = bm["interval_95"] or [None, None]
    return {f"{prefix}_mean": summary["mean"], f"{prefix}_mcse": bm["mcse"], f"{prefix}_lo": bounds[0],
            f"{prefix}_hi": bounds[1], f"{prefix}_covered": bm["covered"], f"{prefix}_status": bm["status"],
            f"{prefix}_geyer_tau": summary["geyer"]["tau"]}


def start_path(ns: str, i: int) -> tuple[np.ndarray, PCG64Source]:
    source = PCG64Source(seeds_for(ns, "local", i)["start"])
    return sample_exact_finite_bead_paths(source, M, BETA, 2)[0], source


def run_metropolis(ns: str, i: int, refresh: bool, n: int = N) -> tuple[dict, dict]:
    condition = "refresh" if refresh else "local"
    seeds = seeds_for(ns, condition, i)
    start, start_source = start_path(ns, i)
    channels = StochasticChannels(PCG64Source(seeds["proposal"]), PCG64Source(seeds["acceptance"]), start_source)
    options = {"centroid_source": PCG64Source(seeds["centroid"])} if refresh else {}
    provenance = {"namespace": ns, "condition": "metropolis", "role": "start", "replica": i,
                  "seed": seeds["start"], "sampler": START_SAMPLER}
    model = PathIntegralHarmonicOscillator(beads=M, temperature=1 / BETA, displacement_delta=DELTA)
    result = model.run_simulation(channels, sweeps=n, burn_in=0, initial_path=start, initial_path_provenance=provenance,
                                  uncertainty_batch_sizes=(BATCH,), record_trace=True, **options)
    trace = result.trace
    _, cv = harmonic_energy_estimators(trace["paths"], BETA)
    roles = result.source_provenance["roles"]
    row = {"replica": i, **{f"seed_{role}": seed for role, seed in seeds.items()},
           "start_sha256": result.initialization_protocol["initial_path_sha256"],
           "paths_sha256": array_sha256(trace["paths"]), "cv_trace_sha256": array_sha256(trace["CV"]),
           "acceptance_rate_4dp": result.acceptance_rate,
           "cv_rederive_max_abs": float(np.max(np.abs(cv - trace["CV"]))),
           **_bm512_fields("cv", result.virial_uncertainty), **_bm512_fields("td", result.thermodynamic_uncertainty),
           "start_bytes_before_run": roles["initialization"]["consumed_bytes_before"],
           **{f"bytes_{role}": (None if roles[role] is None else roles[role]["consumed_bytes_delta"]) for role in LEDGER_ROLES},
           **{f"digest_{role}": (None if roles[role] is None else roles[role]["stream_digest_after"]) for role in LEDGER_ROLES},
           "aliases": ";".join(sorted({a for entry in roles.values() if entry for a in entry["aliases"]}))}
    arrays = {"CV": trace["CV"], "TD": trace["TD"], "centroid": trace["centroid"], "rg2": trace["rg2"],
              "start": trace["initial_path"]}
    return row, arrays


def run_iid(ns: str, i: int, n: int = N) -> tuple[dict, dict]:
    seeds = seeds_for(ns, "iid", i)
    source = PCG64Source(seeds["paths"])
    paths = sample_exact_finite_bead_paths(source, M, BETA, n)
    td, cv = harmonic_energy_estimators(paths, BETA)
    centroid = paths.mean(axis=1)
    rg2 = ((paths - centroid[:, None]) ** 2).mean(axis=1)
    reference = exact_finite_bead_energy(M, BETA)
    row = {"replica": i, "seed_paths": seeds["paths"], "paths_sha256": array_sha256(paths),
           "cv_trace_sha256": array_sha256(cv),
           **_bm512_fields("cv", summarize_series(cv, reference=reference, batch_sizes=(BATCH,))),
           **_bm512_fields("td", summarize_series(td, reference=reference, batch_sizes=(BATCH,))),
           "bytes_paths": source.consumed_bytes, "digest_paths": source.stream_digest}
    return row, {"CV": cv, "TD": td, "centroid": centroid, "rg2": rg2}


def run_replica(ns: str, condition: str, i: int, n: int = N) -> tuple[dict, dict]:
    return run_iid(ns, i, n) if condition == "iid" else run_metropolis(ns, i, condition == "refresh", n)


def block_paths(out: Path, condition: str, block: int) -> tuple[Path, Path]:
    stem = out / "traces" / condition / f"block-{block:04d}"
    return stem.with_suffix(".npz"), stem.with_suffix(".json")


def run_block(ns: str, condition: str, block: int, replicas: int, out: Path, n: int = N, block_size: int = BLOCK) -> str:
    npz, meta = block_paths(out, condition, block)
    if npz.exists() and meta.exists():
        return f"skip {condition} block {block} (complete)"
    npz.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    indices = list(range(block * block_size, min((block + 1) * block_size, replicas)))
    rows, arrays = [], {}
    for i in indices:
        row, values = run_replica(ns, condition, i, n)
        rows.append(row)
        for key, value in values.items():
            arrays.setdefault(key, []).append(value)
    tmp = npz.with_suffix(".npz.tmp")
    with tmp.open("wb") as handle:
        np.savez_compressed(handle, replica=np.asarray(indices), **{k: np.asarray(v) for k, v in arrays.items()})
    os.replace(tmp, npz)
    record = {"namespace": ns, "condition": condition, "block": block, "n": n, "replicas": indices,
              "npz": npz.name, "npz_sha256": sha256(npz), "seconds": time.monotonic() - started, "rows": rows}
    tmp = meta.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(record) + "\n")
    os.replace(tmp, meta)
    return f"{condition} block {block}: {len(indices)} replicas in {record['seconds']:.1f} s"


def merge(out: Path, replicas: int, block_size: int = BLOCK) -> dict:
    """Write results/<condition>_replicas.csv; return per-block hashes."""
    blocks = []
    (out / "results").mkdir(parents=True, exist_ok=True)
    for condition in CONDITIONS:
        rows = []
        for block in range(-(-replicas // block_size)):
            npz, meta = block_paths(out, condition, block)
            record = json.loads(meta.read_text())
            if record["npz_sha256"] != sha256(npz):
                raise SystemExit(f"{npz} does not match its block record")
            rows += record["rows"]
            blocks.append({"condition": condition, "block": block, "npz_sha256": record["npz_sha256"],
                           "json_sha256": sha256(meta), "seconds": record["seconds"]})
        if [row["replica"] for row in rows] != list(range(replicas)):
            raise SystemExit(f"{condition}: replicas are incomplete or out of order")
        with (out / "results" / f"{condition}_replicas.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
    return {"blocks": blocks}


def _task(ns: str, condition: str, block: int, replicas: int, out: str, n: int, block_size: int) -> str:
    return run_block(ns, condition, block, replicas, Path(out), n, block_size)


def run_all(workers: int, out: Path = HERE, ns: str | None = None, replicas: int | None = None,
            n: int = N, block_size: int = BLOCK, verify: bool = True) -> None:
    if verify:
        protocol = load_protocol()
        manifest = verify_manifest()
        ns, replicas = protocol["data"]["namespace"], protocol["data"]["replicas_per_condition"]
    log = out / "results" / "run.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    order = ("refresh", "local", "iid")  # longest tasks first
    tasks = [(c, b) for c in order for b in range(-(-replicas // block_size))]
    context = multiprocessing.get_context("fork")
    started = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with log.open("a") as handle, ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
        handle.write(f"start {started} workers={workers}\n")
        futures = [pool.submit(_task, ns, c, b, replicas, str(out), n, block_size) for c, b in tasks]
        for future in as_completed(futures):
            line = f"{time.strftime('%H:%M:%S')} {future.result()}"
            print(line, flush=True)
            handle.write(line + "\n")
            handle.flush()
        handle.write(f"end {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")
    merged = merge(out, replicas, block_size)
    if verify:
        import scipy
        record = {"protocol_sha256": sha256(HERE / "protocol.json"),
                  "manifest_pregeneration_sha256": sha256(HERE / "manifest_pregeneration.json"),
                  "manifest_git_head": manifest["git_head_before_commit"],
                  "git_head_at_run": subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True,
                                                    text=True).stdout.strip(),
                  "first_start_utc": started, "workers": workers, "python": platform.python_version(),
                  "numpy": np.__version__, "scipy": scipy.__version__, "platform": platform.platform(), **merged,
                  "csv_sha256": {c: sha256(out / "results" / f"{c}_replicas.csv") for c in CONDITIONS}}
        (out / "results" / "run_manifest.json").write_text(json.dumps(record, indent=1) + "\n")


def replay(out: Path = HERE, ns: str | None = None, n: int = N) -> dict:
    """Check I-4: replica 0 of each condition regenerated gives the same CV trace."""
    if ns is None:
        ns = load_protocol()["data"]["namespace"]
    result = {}
    for condition in CONDITIONS:
        npz, _ = block_paths(out, condition, 0)
        with np.load(npz) as stored:
            if int(stored["replica"][0]) != 0:
                raise SystemExit(f"{npz} does not start at replica 0")
            stored_sha = array_sha256(stored["CV"][0])
        _, arrays = run_replica(ns, condition, 0, n)
        result[condition] = {"stored_cv_sha256": stored_sha, "replayed_cv_sha256": array_sha256(arrays["CV"]),
                             "identical": stored_sha == array_sha256(arrays["CV"])}
    (out / "results").mkdir(parents=True, exist_ok=True)
    (out / "results" / "replay.json").write_text(json.dumps(result, indent=1) + "\n")
    return result


def seed_table(ns: str, replicas: int) -> list[str]:
    lines = []
    for i in range(replicas):
        unit = f"{i:05d}"
        for condition, role in (("metropolis", "start"), ("metropolis", "proposal"), ("metropolis", "acceptance"),
                                ("refresh", "centroid"), ("iid", "paths")):
            lines.append(f"{ns}/{condition}/{role}/{unit}\t{namespace_seed(ns, condition, role, unit)}")
    return lines


def write_manifest() -> None:
    protocol = load_protocol()
    path = HERE / "manifest_pregeneration.json"
    if path.exists():
        raise SystemExit(f"{path} exists; a new manifest needs a new protocol version")
    if (HERE / "traces").exists():
        raise SystemExit("traces/ already exists; the manifest must precede every campaign draw")
    import scipy

    ns, replicas = protocol["data"]["namespace"], protocol["data"]["replicas_per_condition"]
    table = seed_table(ns, replicas)
    seeds = [line.split("\t")[1] for line in table]
    if len(set(seeds)) != len(seeds):
        raise SystemExit("seed collision inside the namespace")
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    manifest = {
        "protocol_id": protocol["id"], "protocol_sha256": sha256(HERE / "protocol.json"),
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "git_head_before_commit": head,
        "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(), "file_sha256": {f: sha256(REPO / f) for f in HASHED_FILES},
        "seed_units": len(table), "seed_collisions": 0,
        "seed_table_sha256": hashlib.sha256("\n".join(table).encode()).hexdigest(),
        "seed_table_format": "one line per unit, '<namespace>/<condition>/<role>/<replica:05d>\\t<seed>', replica-major, joined by newlines",
        "seed_table_first_rows": table[:5],
        "statement": "Written before any campaign draw; run_campaign.py run verifies these hashes."}
    path.write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"wrote {path} ({len(table)} seed units)")


def verify_manifest() -> dict:
    path = HERE / "manifest_pregeneration.json"
    if not path.exists():
        raise SystemExit("no manifest_pregeneration.json; run the manifest command and commit it first")
    manifest = json.loads(path.read_text())
    for name, digest in manifest["file_sha256"].items():
        if sha256(REPO / name) != digest:
            raise SystemExit(f"{name} changed since the pre-generation manifest; refusing to run")
    if manifest["protocol_sha256"] != sha256(HERE / "protocol.json"):
        raise SystemExit("protocol.json changed since the pre-generation manifest")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("command", choices=["manifest", "run", "replay"])
    parser.add_argument("--workers", type=int, default=4)
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
