#!/usr/bin/env python3
"""Pre-analysis record for Stage C Amendment 01.

Run once, after the amended scripts are installed and the test suite has passed, and before any analysis
of campaign data. Writes results/amendment01_preanalysis.json.

Records the approved source hashes, the superseded ones, checks that every other manifest-listed file and
every campaign data file is unchanged, checks that no frozen constant of the analysis changed, and checks
that no analysis output exists yet.
"""
import ast
import datetime
import hashlib
import json
import subprocess
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SUPERSEDED_AT = "e5c0140"  # commit holding the manifest and the superseded scripts (0aca45a before the 2026-11-08 rewrite)
AMENDED = ["research/v018-h1-detection/analyze.py", "research/v018-h1-detection/recount.py"]
ANALYSIS_OUTPUTS = ["replay.json", "gates.json", "secondary.json", "recount.json"]
TEST_EVIDENCE = ROOT / "docs/v018/evidence/pytest_stageC_amendment01.txt"


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def git(*args: str) -> bytes:
    return subprocess.run(["git", "-C", str(ROOT), *args], check=True, capture_output=True).stdout


def constants(source: str) -> dict:
    """Top-level NAME = literal assignments (single or tuple targets)."""
    out = {}
    for node in ast.parse(source).body:
        if not isinstance(node, ast.Assign) or len(node.targets) != 1:
            continue
        target = node.targets[0]
        try:
            value = ast.literal_eval(node.value)
        except ValueError:
            continue
        if isinstance(target, ast.Name):
            out[target.id] = value
        elif isinstance(target, ast.Tuple) and isinstance(value, tuple) and len(value) == len(target.elts):
            out.update({t.id: v for t, v in zip(target.elts, value) if isinstance(t, ast.Name)})
    return out


def main() -> None:
    manifest = json.loads((HERE / "manifest_pregeneration.json").read_text())
    run_manifest_bytes = (HERE / "results" / "run_manifest.json").read_bytes()
    run_manifest = json.loads(run_manifest_bytes)

    changed = sorted(p for p, h in manifest["file_sha256"].items() if sha((ROOT / p).read_bytes()) != h)
    cell_mismatches = sorted(f"{c}.{ext}" for c, h in run_manifest["cell_sha256"].items() for ext in ("npz", "jsonl")
                             if sha((HERE / "results" / "cells" / f"{c}.{ext}").read_bytes()) != h[ext])

    scripts = {}
    for p in AMENDED:
        old = git("show", f"{SUPERSEDED_AT}:{p}")
        new = (ROOT / p).read_bytes()
        c_old, c_new = constants(old.decode()), constants(new.decode())
        scripts[p] = {"superseded_sha256": sha(old), "manifest_sha256": manifest["file_sha256"].get(p),
                      "approved_sha256": sha(new),
                      "constants_added": {k: c_new[k] for k in sorted(set(c_new) - set(c_old))},
                      "constants_removed": sorted(set(c_old) - set(c_new)),
                      "constants_changed": sorted(k for k in set(c_old) & set(c_new) if c_old[k] != c_new[k])}

    tests = TEST_EVIDENCE.read_text().splitlines()
    record = {
        "written_utc": datetime.datetime.now(datetime.UTC).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "git_head": git("rev-parse", "HEAD").decode().strip(),
        "approval": "Author approval, 2026-11-08, text in AMENDMENT_01.md",
        "scripts": scripts,
        "manifest_files": len(manifest["file_sha256"]),
        "manifest_files_changed": changed,
        "only_amended_files_changed": changed == sorted(AMENDED),
        "protocol_sha256_matches": sha((HERE / "protocol.json").read_bytes()) == manifest["protocol_sha256"],
        "grid_sha256_matches": sha((HERE / "grid.json").read_bytes()) == manifest["grid_sha256"],
        "run_manifest_sha256": sha(run_manifest_bytes),
        "campaign_cells": len(run_manifest["cell_sha256"]),
        "campaign_runs": sum(len((HERE / "results" / "cells" / f"{c}.jsonl").read_text().splitlines())
                             for c in run_manifest["cell_sha256"]),
        "planned_runs": manifest["total_runs"],
        "campaign_file_mismatches": cell_mismatches,
        "analysis_outputs_present": [f for f in ANALYSIS_OUTPUTS if (HERE / "results" / f).exists()],
        "tests": {"evidence": str(TEST_EVIDENCE.relative_to(ROOT)), "command": tests[0].lstrip("# "), "result": tests[-1]},
    }
    (HERE / "results" / "amendment01_preanalysis.json").write_text(json.dumps(record, indent=1) + "\n")
    print(json.dumps(record, indent=1))


if __name__ == "__main__":
    main()
