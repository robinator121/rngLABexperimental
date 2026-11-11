"""Integrity of a stored campaign: atomic commits, partial-run detection, manifests.

Layout (research/v019-d-pilot/protocol.json, storage):
    <root>/traces/{condition}/{i:05d}.npz        arrays
    <root>/results/runs/{condition}/{i:05d}.json  run record, holding trace_sha256

The run record is the commit marker. A run counts as COMPLETE only if its
record parses, names the planned identity, and its trace file exists with the
recorded SHA-256. Everything else is classified, never repaired:

    MISSING    neither file exists
    PARTIAL    a trace without a record, or leftover temporary files
    CORRUPT    unreadable or empty record, identity mismatch, missing trace, hash mismatch
    COMPLETE   as above

`commit_run` writes both files through temporary names in the target directory
and renames the record last, so a process killed at any instant leaves either a
COMPLETE run or one that this module reports as MISSING or PARTIAL. Nothing in
this file generates data or decides whether an interrupted identity may be run
again; that is the campaign driver's rule.
"""
from __future__ import annotations

import hashlib
import io
import json
import os
import time
from pathlib import Path
from typing import Any, Iterable, Mapping

import numpy as np

STATES = ("COMPLETE", "MISSING", "PARTIAL", "CORRUPT")
TMP_SUFFIX = ".tmp"


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def run_paths(root: Path, condition: str, replica: int) -> tuple[Path, Path]:
    """(record path, trace path) of one planned run."""
    return (Path(root) / "results/runs" / condition / f"{replica:05d}.json",
            Path(root) / "traces" / condition / f"{replica:05d}.npz")


def canonical_json(value: Any) -> bytes:
    """protocol.json storage.json encoding."""
    return (json.dumps(value, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode("utf-8")


def _fsync_directory(directory: Path) -> None:
    descriptor = os.open(directory, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def atomic_write(path: Path, data: bytes) -> None:
    """Write through a temporary file in the same directory, flush to disk, then rename."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".{os.getpid()}{TMP_SUFFIX}")
    with open(temporary, "wb") as handle:
        handle.write(data)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)
    _fsync_directory(path.parent)


def commit_run(root: Path, record: Mapping[str, Any], arrays: Mapping[str, np.ndarray], *,
               overwrite: bool = False) -> dict[str, Any]:
    """Store one run atomically; the record is renamed into place last.

    Refuses to touch a run that already has a record unless overwrite=True: a
    completed draw is never silently replaced.
    """
    condition, replica = record["condition"], int(record["replica"])
    record_path, trace_path = run_paths(root, condition, replica)
    if record_path.exists() and not overwrite:
        raise FileExistsError(f"run {condition}/{replica:05d} already has a record")
    buffer = io.BytesIO()
    np.savez_compressed(buffer, **{key: np.ascontiguousarray(value) for key, value in arrays.items()})
    trace_bytes = buffer.getvalue()
    stored = {**record, "trace_path": trace_path.relative_to(Path(root)).as_posix(),
              "trace_sha256": sha256_bytes(trace_bytes)}
    record_bytes = canonical_json(stored)
    atomic_write(trace_path, trace_bytes)
    atomic_write(record_path, record_bytes)
    return {"condition": condition, "replica": replica,
            "files": {stored["trace_path"]: stored["trace_sha256"],
                      record_path.relative_to(Path(root)).as_posix(): sha256_bytes(record_bytes)}}


def classify_run(root: Path, condition: str, replica: int, *, verify_hash: bool = True) -> dict[str, Any]:
    """State of one planned run on disk, with the reason."""
    record_path, trace_path = run_paths(root, condition, replica)
    leftovers = sorted(p.name for directory in (record_path.parent, trace_path.parent) if directory.is_dir()
                       for p in directory.glob(f"{replica:05d}.*{TMP_SUFFIX}"))
    out: dict[str, Any] = {"condition": condition, "replica": replica, "leftover_temporary_files": leftovers}
    if not record_path.exists():
        if trace_path.exists() or leftovers:
            return {**out, "state": "PARTIAL", "reason": "trace or temporary file without a record"}
        return {**out, "state": "MISSING", "reason": None}
    try:
        raw = record_path.read_bytes()
        record = json.loads(raw.decode("utf-8"))
        if not isinstance(record, dict):
            raise ValueError("record is not a JSON object")
    except (ValueError, OSError) as error:
        return {**out, "state": "CORRUPT", "reason": f"unreadable record: {type(error).__name__}: {error}"}
    if record.get("condition") != condition or record.get("replica") != replica:
        return {**out, "state": "CORRUPT", "reason": "record identity differs from its path"}
    if not trace_path.exists():
        return {**out, "state": "CORRUPT", "reason": "record without a trace"}
    if trace_path.stat().st_size == 0:
        return {**out, "state": "CORRUPT", "reason": "empty trace"}
    out.update(record_sha256=sha256_bytes(raw), record_bytes=len(raw), trace_bytes=trace_path.stat().st_size,
               run_status=record.get("status"))
    if verify_hash:
        actual = sha256_file(trace_path)
        if record.get("trace_sha256") != actual:
            return {**out, "state": "CORRUPT", "reason": "trace hash differs from the record"}
        out["trace_sha256"] = actual
    return {**out, "state": "COMPLETE", "reason": None}


def scan_store(root: Path, plan: Mapping[str, int], *, verify_hash: bool = True) -> dict[str, Any]:
    """Classify every planned run and list files that belong to no planned run."""
    root = Path(root)
    runs = [classify_run(root, condition, replica, verify_hash=verify_hash)
            for condition, count in plan.items() for replica in range(count)]
    planned = set()
    for condition, count in plan.items():
        for replica in range(count):
            planned.update(run_paths(root, condition, replica))
    unplanned = sorted(p.relative_to(root).as_posix() for base in ("results/runs", "traces") if (root / base).is_dir()
                       for p in (root / base).rglob("*") if p.is_file() and p not in planned)
    counts = {state: sum(1 for r in runs if r["state"] == state) for state in STATES}
    return {"root": str(root), "planned": len(runs), "counts": counts, "unplanned_files": unplanned,
            "not_complete": [r for r in runs if r["state"] != "COMPLETE"], "runs": runs,
            "all_complete": counts["COMPLETE"] == len(runs) and not unplanned}


def build_manifest(root: Path, plan: Mapping[str, int], extra: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Every planned identity with its state and file hashes; missing files are explicit entries."""
    scan = scan_store(root, plan)
    entries = []
    for run in scan["runs"]:
        record_path, trace_path = run_paths(root, run["condition"], run["replica"])
        entries.append({"condition": run["condition"], "replica": run["replica"], "state": run["state"],
                        "reason": run["reason"], "run_status": run.get("run_status"),
                        "record": {"path": record_path.relative_to(Path(root)).as_posix(),
                                   "sha256": run.get("record_sha256"), "bytes": run.get("record_bytes")},
                        "trace": {"path": trace_path.relative_to(Path(root)).as_posix(),
                                  "sha256": run.get("trace_sha256"), "bytes": run.get("trace_bytes")}})
    body = {"schema": "rnglab-v019-audit-manifest-v1", "plan": dict(plan), "counts": scan["counts"],
            "unplanned_files": scan["unplanned_files"], "runs": entries, **(dict(extra) if extra else {})}
    body["manifest_sha256"] = sha256_bytes(canonical_json(body))
    return body


def verify_manifest(root: Path, manifest: Mapping[str, Any]) -> dict[str, Any]:
    """Compare the store with a manifest: changed, missing, added files, and the manifest's own hash."""
    body = {k: v for k, v in manifest.items() if k != "manifest_sha256"}
    problems: list[dict[str, Any]] = []
    if manifest.get("manifest_sha256") != sha256_bytes(canonical_json(body)):
        problems.append({"kind": "MANIFEST_ALTERED", "path": None})
    listed = set()
    for entry in manifest["runs"]:
        for part in ("record", "trace"):
            path = Path(root) / entry[part]["path"]
            listed.add(path)
            expected = entry[part]["sha256"]
            if expected is None:
                if path.exists():
                    problems.append({"kind": "APPEARED", "path": entry[part]["path"]})
            elif not path.exists():
                problems.append({"kind": "MISSING", "path": entry[part]["path"]})
            elif sha256_file(path) != expected:
                problems.append({"kind": "CHANGED", "path": entry[part]["path"]})
    for base in ("results/runs", "traces"):
        if (Path(root) / base).is_dir():
            for path in (Path(root) / base).rglob("*"):
                if path.is_file() and path not in listed:
                    problems.append({"kind": "UNLISTED", "path": path.relative_to(Path(root)).as_posix()})
    return {"ok": not problems, "problems": problems, "files_checked": 2 * len(manifest["runs"])}


def append_journal(root: Path, event: str, **fields: Any) -> None:
    """Append-only event log (start, commit, interruption found, resume); one JSON object per line."""
    path = Path(root) / "results" / "journal.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps({"time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "pid": os.getpid(), "event": event, **fields},
                      sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    with open(path, "a", encoding="utf-8") as handle:
        handle.write(line)
        handle.flush()
        os.fsync(handle.fileno())


def pending(root: Path, plan: Mapping[str, int]) -> list[tuple[str, int, str]]:
    """Planned runs that are not COMPLETE, with their state: the work list a restart must look at."""
    return [(r["condition"], r["replica"], r["state"]) for r in scan_store(root, plan, verify_hash=False)["runs"]
            if r["state"] != "COMPLETE"]


def namespace_collisions(seed_of, namespaces: Mapping[str, Mapping[str, int]], roles: Iterable[str]) -> list[Any]:
    """Seeds shared by two different (namespace, condition, role, replica) labels. Empty when disjoint."""
    seen: dict[int, tuple[str, str, str, int]] = {}
    clashes = []
    roles = tuple(roles)
    for namespace, plan in namespaces.items():
        for condition, count in plan.items():
            for replica in range(count):
                for role in roles:
                    seed = seed_of(namespace, condition, role, replica)
                    if seed in seen:
                        clashes.append((seen[seed], (namespace, condition, role, replica)))
                    seen[seed] = (namespace, condition, role, replica)
    return clashes
