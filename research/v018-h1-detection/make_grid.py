#!/usr/bin/env python3
"""Campaign severity grids from the pilot, by the frozen rule (PROTOCOL.md section 3.4).

RNGLab v0.18 (written by Claude Code, 2026-11-06).

    python research/v018-h1-detection/make_grid.py

Reads pilot/summary.json; writes grid.json. For SHIFT, LAG1 and DUP: s_lo is the
smallest pilot severity with wrong-rate >= 0.10, s_hi the smallest with >= 0.90
(else the largest); the grid is 6 log-spaced severities from s_lo/2 to s_hi.
If no pilot severity reaches 0.10, the grid is 6 log-spaced severities from
s_max/32 to s_max and the family is labelled "no harm observed in the pilot
range". GFSR uses all lag pairs that passed the irreducibility check.
"""
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import h1lib  # noqa: E402


def grid_for(cells: list[dict]) -> dict:
    cells = sorted(cells, key=lambda c: c["severity"])
    severities = [c["severity"] for c in cells]
    rates = [c["wrong_rate"] for c in cells]
    reached = [s for s, w in zip(severities, rates) if w >= 0.10]
    if not reached:
        top = severities[-1]
        return {"severities": np.geomspace(top / 32, top, 6).tolist(), "label": "no harm observed in the pilot range",
                "s_lo": None, "s_hi": None}
    s_lo = reached[0]
    high = [s for s, w in zip(severities, rates) if w >= 0.90]
    s_hi = high[0] if high else severities[-1]
    return {"severities": np.geomspace(s_lo / 2, s_hi, 6).tolist(), "label": "harm transition from the pilot",
            "s_lo": s_lo, "s_hi": s_hi, "s_hi_from": "first severity with wrong-rate >= 0.90" if high else "largest pilot severity"}


def main() -> None:
    summary_path = HERE / "pilot" / "summary.json"
    summary = json.loads(summary_path.read_text())
    grid = {"rule": "PROTOCOL.md section 3.4", "pilot_summary_sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
            "families": {}}
    for family in ("SHIFT", "LAG1", "DUP"):
        grid["families"][family] = grid_for([c for c in summary["cells"] if c["family"] == family])
    pairs = [list(p) for p in h1lib.GFSR_PAIRS if h1lib.trinomial_irreducible(*p)]
    grid["families"]["GFSR"] = {"severities": pairs, "label": "all irreducible lag pairs",
                                "dropped": [list(p) for p in h1lib.GFSR_PAIRS if list(p) not in pairs]}
    (HERE / "grid.json").write_text(json.dumps(grid, indent=1) + "\n")
    print(json.dumps(grid["families"], indent=1))


if __name__ == "__main__":
    main()
