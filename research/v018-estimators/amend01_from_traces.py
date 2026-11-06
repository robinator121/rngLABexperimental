#!/usr/bin/env python3
"""Amendment 01: recompute split-R-hat and add a burn-in diagnostic from saved Ising traces.

RNGLab v0.18 (written by Claude Code, 2026-11-05). See AMENDMENT_01.md. The
generation run stored split-R-hat values computed with the pre-amendment
normal-score denominator and NaN handling (results/ising/*.npz: rhat_m,
rhat_m2, rhat_e). Those files are left byte-identical. This script recomputes
R-hat with the corrected rnglab.convergence.split_rhat from the integer traces
and adds a split-half diagnostic: per-chain means of m^2 and e over the first
and the second half of the retained window. Output: results/amend01/<cell>.npz.

    PYTHONPATH=src python research/v018-estimators/amend01_from_traces.py
"""
import hashlib
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from rnglab.convergence import split_rhat  # noqa: E402

BLOCK = 4096
K = 16


def main() -> None:
    out_dir = HERE / "results" / "amend01"
    out_dir.mkdir(parents=True, exist_ok=True)
    for trace in sorted((HERE / "traces").glob("ising_*.npz")):
        out = out_dir / trace.name
        if out.exists():
            print(f"skip {out.name} (exists)")
            continue
        original = np.load(HERE / "results" / "ising" / trace.name)
        recorded = str(original["trace_sha256"])
        digest = hashlib.sha256(trace.read_bytes()).hexdigest()
        if digest != recorded:
            raise SystemExit(f"{trace.name}: trace SHA-256 differs from the one recorded at generation")
        z = np.load(trace)
        L, burn = int(z["L"]), int(z["burn_in"])
        N = L * L
        cols = {key: [] for key in ("rhat_m", "rhat_m2", "rhat_e", "m2_first", "m2_second", "e_first", "e_second")}
        for first in range(0, len(z["M"]), BLOCK):
            m = z["M"][first:first + BLOCK, burn:].astype(float) / N
            e = z["E"][first:first + BLOCK, burn:].astype(float) / N
            n = m.shape[1]
            for key, values in (("rhat_m", m), ("rhat_m2", m * m), ("rhat_e", e)):
                cols[key].append(split_rhat(values.reshape(-1, K, n)))
            half = n // 2
            cols["m2_first"].append((m[:, :half] ** 2).mean(axis=1))
            cols["m2_second"].append((m[:, half:] ** 2).mean(axis=1))
            cols["e_first"].append(e[:, :half].mean(axis=1))
            cols["e_second"].append(e[:, half:].mean(axis=1))
        np.savez_compressed(out, trace_sha256=digest, **{k: np.concatenate(v) for k, v in cols.items()})
        print(f"wrote {out.name}")


if __name__ == "__main__":
    main()
