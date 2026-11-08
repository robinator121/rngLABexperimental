#!/usr/bin/env python3
"""Stage C sensitivity analysis 01 (SENSITIVITY_01.md): the battery without the two targeted word tests.

    PYTHONPATH=src python research/v018-h1-detection/sensitivity01.py [--self-check]

Writes results/sensitivity01.json. With --self-check it only checks that the full battery reproduces
results/gates.json, and writes nothing.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from itertools import islice
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))
import analyze  # noqa: E402
import h1lib  # noqa: E402

SP_COLUMNS = tuple(range(8))
assert h1lib.CONV_P_KEYS[:8] == ("frequency", "block_frequency", "runs", "serial_1", "serial_2",
                                  "approximate_entropy", "cumulative_sums", "longest_run")
assert h1lib.CONV_P_KEYS[8:] == ("lag1_correlation", "adjacent_duplicates")
NAMESPACE = "rnglab-v018-h1-sensitivity01"


def conv_series(cell: dict, columns: tuple[int, ...]) -> np.ndarray:
    p = cell["conv_p"][..., list(columns)]
    p = np.where(np.isnan(p), np.inf, p)
    return -np.log10(np.maximum(p.min(axis=-1), h1lib.P_FLOOR))


def gates(columns: tuple[int, ...], ns: str | None = None) -> dict:
    original = analyze.conv_series
    analyze.conv_series = lambda cell: conv_series(cell, columns)
    try:
        return analyze.gates(HERE) if ns is None else analyze.gates(HERE, ns=ns)
    finally:
        analyze.conv_series = original


def differing_paths(a, b, path: str = "$"):
    if a == b:
        return
    if isinstance(a, dict) and isinstance(b, dict):
        for key in sorted(a.keys() | b.keys()):
            if key not in a or key not in b:
                yield f"{path}.{key}"
            else:
                yield from differing_paths(a[key], b[key], f"{path}.{key}")
    elif isinstance(a, list) and isinstance(b, list):
        for i in range(max(len(a), len(b))):
            if i >= len(a) or i >= len(b):
                yield f"{path}[{i}]"
            else:
                yield from differing_paths(a[i], b[i], f"{path}[{i}]")
    else:
        yield path


def self_check(stage_c: dict) -> dict:
    actual = json.loads(json.dumps(gates(tuple(range(len(h1lib.CONV_P_KEYS))))))
    return {"identical": actual == stage_c, "differing_paths": list(islice(differing_paths(actual, stage_c), 10))}


def freeze_record() -> dict:
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=HERE)

    root = Path(git("rev-parse", "--show-toplevel").decode().strip())
    paths = [HERE / "SENSITIVITY_01.md", Path(__file__).resolve()]
    names = [p.relative_to(root).as_posix() for p in paths]
    if git("status", "--porcelain", "--", *(str(p) for p in paths)).strip():
        raise RuntimeError("SENSITIVITY_01.md and sensitivity01.py must be tracked and unchanged against HEAD")
    head = git("rev-parse", "HEAD").decode().strip()
    files = {}
    for path, name in zip(paths, names):
        git("ls-files", "--error-unmatch", "--", str(path))
        commit = git("log", "-1", "--format=%H", "--", str(path)).decode().strip()
        sha = hashlib.sha256(path.read_bytes()).hexdigest()
        matches_head = sha == hashlib.sha256(git("show", f"{head}:{name}")).hexdigest()
        matches_commit = sha == hashlib.sha256(git("show", f"{commit}:{name}")).hexdigest()
        if not matches_head or not matches_commit:
            raise RuntimeError(f"{name} does not match its committed bytes")
        files[name] = {"last_changed_commit": commit, "sha256": sha,
                       "matches_head": matches_head, "matches_freezing_commit": matches_commit}
    return {"HEAD": head, "frozen": True, "files": files}


def first_firing(name: str, tau: float) -> dict:
    cell = analyze.load_cell(HERE, name)
    det = analyze.detect(conv_series(cell, SP_COLUMNS), cell["conv_end"][0], cell["conv_bytes"], cell["onset"], tau)
    firing = {}
    for i in np.flatnonzero(det["post"]):
        row = cell["conv_p"][i, det["first_index"][i], list(SP_COLUMNS)]
        key = h1lib.CONV_P_KEYS[SP_COLUMNS[int(np.nanargmin(row))]]
        firing[key] = firing.get(key, 0) + 1
    return firing


def sensitivity(stage_c: dict, freeze: dict, check: dict) -> dict:
    g = gates(SP_COLUMNS, ns=NAMESPACE)
    healthy = g["validity"]["V2"]["rate_C"]
    v2_pass = healthy["cp95"][1] <= 0.10
    implementation_pass = g["verdicts"]["implementation_checks_pass"]
    cells, changes = {}, []
    for name, c in g["cells"].items():
        old = stage_c["cells"][name]["class"]
        cells[name] = {"family": c["family"], "severity": c["severity"], "runs": c["runs"],
                       "onset_flagged": c["onset_flagged"], "C_prime_alarm": c["C_alarm"],
                       "wrong": c["wrong"], "A_alarm": {"k": c["A_alarm"]["k"]},
                       "A_minus_C_prime": c["A_minus_C"], "pass_but_broken": c["pass_but_broken"],
                       "class_c_prime": c["class"], "class_stage_c": old, "class_changed": old != c["class"],
                       "latency": c["latency"], "first_firing_sp80022": first_firing(name, g["thresholds"]["C"])}
        if old != c["class"]:
            changes.append({"cell": name, "stage_c": old, "c_prime": c["class"]})
    return {"freeze": freeze, "self_check": check, "c_prime_computed": True, "namespace": NAMESPACE,
            "thresholds": {"tau_C_prime": g["thresholds"]["C"], "tau_A": g["thresholds"]["A"],
                           "tau_A_matches_stage_c": g["thresholds"]["A"] == stage_c["thresholds"]["A"]},
            "healthy_C_prime": healthy, "v2_c_prime_pass": v2_pass, "V1": g["validity"]["V1"],
            "implementation_checks_pass": implementation_pass,
            "interpretable": g["validity"]["V1"]["result"] == "PASS" and v2_pass and implementation_pass,
            "cells": cells, "verdicts_c_prime": g["verdicts"]["H1-C"],
            "cells_h1_or_undetected_harm": [n for n, c in cells.items() if c["class_c_prime"] in ("H1", "undetected_harm")],
            "class_changes": changes}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-check", action="store_true")
    args = parser.parse_args()
    freeze = None
    if not args.self_check:
        try:
            freeze = freeze_record()
        except (RuntimeError, subprocess.CalledProcessError, OSError) as e:
            parser.exit(1, f"error: freeze check failed: {e}\n")
    stage_c = json.loads((HERE / "results" / "gates.json").read_text())
    check = self_check(stage_c)
    if args.self_check:
        print(json.dumps(check, indent=1))
        return 0 if check["identical"] else 1
    result = {"freeze": freeze, "self_check": check, "c_prime_computed": False}
    if check["identical"]:
        result = sensitivity(stage_c, freeze, check)
    (HERE / "results" / "sensitivity01.json").write_text(json.dumps(result, indent=1) + "\n")
    if not check["identical"]:
        print(json.dumps(result, indent=1))
        return 1
    print(json.dumps({"classes": {n: c["class_c_prime"] for n, c in result["cells"].items()},
                      "verdicts_c_prime": result["verdicts_c_prime"], "interpretable": result["interpretable"]}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
