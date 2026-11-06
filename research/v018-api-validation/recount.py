#!/usr/bin/env python3
"""Independent recount of the v0.18 Phase 2 primary endpoints (protocol section 5.4).

RNGLab v0.18 (written by Claude Code, 2026-11-05). Shares no code with
analyze.py, rnglab12.uncertainty or the package models: E_M comes from its closed
form, BM512 from a NumPy reshape of the saved CV traces, intervals from SciPy.
The bootstrap uses the protocol's seeds and chunked resampling scheme, so its
intervals must equal analyze.py's to rounding.

    PYTHONPATH=src python research/v018-api-validation/recount.py

Reads traces/, results/<condition>_replicas.csv and results/gates.json; writes
results/recount.json.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
M, BETA, BATCH = 16, 2.0, 512  # BM512: n // 512 batches (16 at n = 8192)
RESAMPLES, CHUNK = 10_000, 500


def reference_energy() -> float:
    z = BETA / (2 * M)
    return 0.5 / (math.tanh(M * math.asinh(z)) * math.sqrt(1 + z * z))


def seed(ns: str, *parts) -> int:
    return int.from_bytes(hashlib.sha256("/".join([ns, *parts]).encode()).digest()[:8], "big")


def load_cv(out: Path, condition: str) -> np.ndarray:
    blocks, replicas = [], []
    for path in sorted((out / "traces" / condition).glob("block-*.npz")):
        with np.load(path) as z:
            blocks.append(z["CV"])
            replicas.append(z["replica"])
    if not np.array_equal(np.concatenate(replicas), np.arange(sum(len(r) for r in replicas))):
        raise SystemExit(f"{condition}: trace replicas incomplete")
    return np.concatenate(blocks)


def bm512(cv: np.ndarray, reference: float) -> dict:
    r, n = cv.shape
    batches = n // BATCH
    batch_means = cv[:, :batches * BATCH].reshape(r, batches, BATCH).mean(axis=2)
    mean = batch_means.mean(axis=1)
    mcse = np.sqrt(batch_means.var(axis=1, ddof=1) / batches)
    q = stats.t.ppf(0.975, batches - 1)
    return {"mean": mean, "mcse": mcse, "covered": (mean - q * mcse <= reference) & (reference <= mean + q * mcse)}


def percentile_bootstrap(fn, arrays, s: int, resamples: int) -> list[float]:
    rng = np.random.default_rng(s)
    r = len(arrays[0])
    draws = []
    for start in range(0, resamples, CHUNK):
        idx = rng.integers(r, size=(min(CHUNK, resamples - start), r))
        draws.append(fn(*(a[idx] for a in arrays)))
    return np.quantile(np.concatenate(draws), [0.025, 0.975]).tolist()


def recount(out: Path = HERE, ns: str | None = None, resamples: int = RESAMPLES) -> dict:
    ns = ns or json.loads((HERE / "protocol.json").read_text())["data"]["namespace"]
    reference = reference_energy()
    gates = json.loads((out / "results" / "gates.json").read_text())
    comparisons, flags = [], {}

    def compare(name, mine, theirs, tolerance=0.0):
        ok = (mine == theirs) if tolerance == 0 else all(abs(a - b) <= tolerance * max(1.0, abs(b))
                                                          for a, b in zip(np.atleast_1d(mine), np.atleast_1d(theirs)))
        comparisons.append({"quantity": name, "recount": mine, "analysis": theirs, "agree": bool(ok)})

    for condition in ("local", "refresh", "iid"):
        cv = load_cv(out, condition)
        mine = bm512(cv, reference)
        with (out / "results" / f"{condition}_replicas.csv").open() as handle:
            api_flags = np.array([row["cv_covered"] == "True" for row in csv.DictReader(handle)])
        flags[condition] = mine["covered"]
        p = gates["primary"][condition]
        r, k = len(cv), int(mine["covered"].sum())
        cp = [float(stats.beta.ppf(0.025, k, r - k + 1)) if k else 0.0, float(stats.beta.ppf(0.975, k + 1, r - k)) if k < r else 1.0]
        dev = mine["mean"] - reference
        q = stats.t.ppf(0.975, r - 1)
        half = q * dev.std(ddof=1) / math.sqrt(r)
        ratio = lambda m2, x: m2.mean(axis=-1) / x.var(axis=-1, ddof=1)
        compare(f"{condition}: covered count", k, p["covered"])
        compare(f"{condition}: per-replica flags equal to API flags", int((mine["covered"] == api_flags).sum()), r)
        compare(f"{condition}: CP95", cp, p["cp95"], 1e-12)
        compare(f"{condition}: variance ratio", float(ratio(mine["mcse"] ** 2, mine["mean"])), p["variance_ratio"]["ratio"], 1e-9)
        compare(f"{condition}: variance-ratio bootstrap CI95",
                percentile_bootstrap(ratio, [mine["mcse"] ** 2, mine["mean"]], seed(ns, "analysis", "bootstrap", f"cv/{condition}"), resamples),
                p["variance_ratio"]["ci95"], 1e-9)
        compare(f"{condition}: bias CI95", [float(dev.mean() - half), float(dev.mean() + half)], p["bias"]["ci95"], 1e-9)
    d = flags["refresh"].astype(float) - flags["local"].astype(float)
    compare("paired: improved_only", int((d == 1).sum()), gates["P-RL"]["improved_only"])
    compare("paired: local_only", int((d == -1).sum()), gates["P-RL"]["local_only"])
    compare("paired: bootstrap CI95", percentile_bootstrap(lambda x: x.mean(axis=-1), [d], seed(ns, "analysis", "bootstrap", "paired"), resamples),
            gates["P-RL"]["bootstrap_ci95"], 1e-12)
    result = {"reference_E_M": reference, "comparisons": comparisons, "all_agree": all(c["agree"] for c in comparisons)}
    (out / "results" / "recount.json").write_text(json.dumps(result, indent=1) + "\n")
    return result


if __name__ == "__main__":
    r = recount()
    print(f"{sum(c['agree'] for c in r['comparisons'])}/{len(r['comparisons'])} agree; all_agree = {r['all_agree']}")
