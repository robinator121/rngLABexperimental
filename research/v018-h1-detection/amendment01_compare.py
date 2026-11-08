#!/usr/bin/env python3
"""Run the superseded analyze.py (commit 0aca45a) on the campaign data and compare it with the approved one.

Pre-declared in AMENDMENT_01.md. Run after `analyze.py` has written results/gates.json.

    PYTHONPATH=src python research/v018-h1-detection/amendment01_compare.py

Writes results/amendment01_superseded/gates.json (or the error, if the superseded version fails) and
results/amendment01_comparison.json: every path where the two gates.json differ, and the paths present in
only one of them.
"""
from __future__ import annotations

import json
import subprocess
import traceback
import types
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SUPERSEDED_AT = "0aca45a"


def superseded_module() -> types.ModuleType:
    source = subprocess.run(["git", "-C", str(ROOT), "show", f"{SUPERSEDED_AT}:research/v018-h1-detection/analyze.py"],
                            check=True, capture_output=True, text=True).stdout
    module = types.ModuleType("analyze_superseded")
    module.__file__ = str(HERE / "analyze.py")  # so its HERE is this folder and it reads the same data
    exec(compile(source, f"{SUPERSEDED_AT}:analyze.py", "exec"), module.__dict__)
    return module


def walk(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            p = f"{path}/{k}"
            if k not in b:
                out["only_superseded"].append(p)
            elif k not in a:
                out["only_approved"].append(p)
            else:
                walk(a[k], b[k], p, out)
    elif a != b:
        out["different"].append({"path": path, "superseded": a, "approved": b})


def main() -> None:
    target = HERE / "results" / "amendment01_superseded"
    target.mkdir(parents=True, exist_ok=True)
    approved = json.loads((HERE / "results" / "gates.json").read_text())
    try:
        old = json.loads(json.dumps(superseded_module().gates(HERE)))
    except Exception:  # noqa: BLE001 - a failure of the superseded version is itself a result to report
        error = traceback.format_exc()
        (target / "error.txt").write_text(error)
        result = {"superseded_ran": False, "error": error.strip().splitlines()[-1]}
    else:
        (target / "gates.json").write_text(json.dumps(old, indent=1) + "\n")
        diff = {"different": [], "only_approved": [], "only_superseded": []}
        walk(old, approved, "", diff)
        decision = {"class", "verdicts", "interpretable", "result", "k", "n", "rate", "runs", "bonferroni_level", "thresholds"}
        result = {"superseded_ran": True, **diff,
                  "decision_paths_different": [d["path"] for d in diff["different"] if decision & set(d["path"].split("/"))],
                  "identical_except_added_fields": not diff["different"] and not diff["only_superseded"]}
    (HERE / "results" / "amendment01_comparison.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({k: (v if not isinstance(v, list) else len(v)) for k, v in result.items()}, indent=1))


if __name__ == "__main__":
    main()
