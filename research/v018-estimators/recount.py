#!/usr/bin/env python3
"""Independent recount of the primary Phase 1 coverage counts.

RNGLab v0.18 (written by Claude Code, 2026-11-05). Recomputes, without using
rnglab.longrun, rnglab.susceptibility or rnglab16.core.enumeration:

- exact Ising references by its own brute-force enumeration;
- S coverage from the saved integer traces (own batch means, t15);
- G coverage from the saved integer traces (explicit leave-one-chain-out loop);
- AR(1) BM(32) and BM(512) coverage at n = 8192 for every phi, regenerating the
  series from the frozen seeds by an explicit time-step recursion (not lfilter).

Then compares with results/gates.json and results/ar1_table.csv. The two
implementations share NumPy, SciPy and the seed derivation.

    PYTHONPATH=src python research/v018-estimators/recount.py
"""
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

HERE = Path(__file__).resolve().parent
PROTOCOL = json.loads((HERE / "protocol.json").read_text())
NS = PROTOCOL["seed"]["namespace"]


def seed(*parts):
    return int.from_bytes(hashlib.sha256("/".join([NS, *map(str, parts)]).encode()).digest()[:8], "big")


def enumerate_ising(L, T):
    n = L * L
    states = np.array([[1 if (s >> i) & 1 else -1 for i in range(n)] for s in range(2 ** n)], dtype=np.int64)
    grid = states.reshape(-1, L, L)
    energy = np.zeros(len(grid))
    for i in range(L):
        for j in range(L):
            energy -= grid[:, i, j] * (grid[:, (i + 1) % L, j] + grid[:, i, (j + 1) % L])
    mag = grid.sum(axis=(1, 2)).astype(float)
    w = np.exp(-(energy - energy.min()) / T)
    w /= w.sum()
    return float((w @ mag ** 2 - (w @ mag) ** 2) / (n * T))


def recount_ising():
    out = {}
    q = student_t.ppf(0.975, 15)
    for path in sorted((HERE / "traces").glob("ising_*.npz")):
        z = np.load(path)
        L, T, burn = int(z["L"]), float(z["T"]), int(z["burn_in"])
        N = L * L
        chi = enumerate_ising(L, T)
        m = z["M"][:, burn:].astype(float) / N
        series = N / T * m * m
        means = series.reshape(len(series), 16, -1).mean(axis=2)
        centre = means.mean(axis=1)
        se = means.std(axis=1, ddof=1) / math.sqrt(16)
        s_cov = int(np.sum(np.abs(centre - chi) <= q * se))
        a2, a1 = (m * m).mean(axis=1), m.mean(axis=1)
        g_cov = 0
        for e in range(len(m) // 16):
            A, B = a2[16 * e:16 * e + 16], a1[16 * e:16 * e + 16]
            full = N / T * (A.mean() - B.mean() ** 2)
            loo = np.array([N / T * (np.delete(A, i).mean() - np.delete(B, i).mean() ** 2) for i in range(16)])
            est = 16 * full - 15 * loo.mean()
            half = q * math.sqrt(15 / 16 * np.sum((loo - loo.mean()) ** 2))
            g_cov += abs(est - chi) <= half
        out[path.stem] = {"chi_exact_recount": chi, "S_covered": s_cov, "S_units": len(m),
                          "G_covered": int(g_cov), "G_units": len(m) // 16}
    return out


def recount_ar1():
    out = {}
    spec = PROTOCOL["data"]["ar1"]
    n = 8192
    size = spec["chunk_size"][str(n)]
    for phi in spec["phi"]:
        cond = f"ar1_phi{float(phi)!r}_n{n}"
        xs = []
        for c in range(spec["replicas"] // size):
            children = np.random.SeedSequence(seed(cond, "innovations", f"chunk-{c:04d}")).spawn(size)
            eps = np.stack([np.random.default_rng(ch).standard_normal(n) for ch in children])
            x = np.empty_like(eps)
            x[:, 0] = eps[:, 0]
            scale = math.sqrt(1 - phi * phi)
            for t in range(1, n):
                x[:, t] = phi * x[:, t - 1] + scale * eps[:, t]
            xs.append(x)
        x = np.concatenate(xs)
        counts = {}
        for b in (32, 512):
            a = n // b
            means = x.reshape(len(x), a, b).mean(axis=2)
            se = means.std(axis=1, ddof=1) / math.sqrt(a)
            counts[f"BM{b}"] = int(np.sum(np.abs(means.mean(axis=1)) <= student_t.ppf(0.975, a - 1) * se))
        out[cond] = counts
    return out


def main():
    gates = json.loads((HERE / "results" / "gates.json").read_text())["gates"]
    with (HERE / "results" / "ar1_table.csv").open() as handle:
        table = {(r["setting"], r["method"], r["b"]): int(r["covered"]) for r in csv.DictReader(handle)}
    ising = recount_ising()
    ar1 = recount_ar1()
    comparisons = []
    for cell, r in ising.items():
        comparisons.append({"item": f"S-S[{cell}]", "recount": r["S_covered"], "analysis": gates[f"S-S[{cell}]"]["coverage"]["covered"]})
        comparisons.append({"item": f"S-G[{cell}]", "recount": r["G_covered"], "analysis": gates[f"S-G[{cell}]"]["coverage"]["covered"]})
    for cond, r in ar1.items():
        for key, value in r.items():
            comparisons.append({"item": f"{cond} {key}", "recount": value, "analysis": table[(cond, "BM", key[2:])]})
    for c in comparisons:
        c["agree"] = c["recount"] == c["analysis"]
    result = {"script": "research/v018-estimators/recount.py", "all_agree": all(c["agree"] for c in comparisons),
              "comparisons": comparisons, "ising_exact_recount": {k: v["chi_exact_recount"] for k, v in ising.items()}}
    (HERE / "results" / "recount.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({"all_agree": result["all_agree"], "items": len(comparisons),
                      "disagreements": [c for c in comparisons if not c["agree"]]}, indent=1))


if __name__ == "__main__":
    main()
