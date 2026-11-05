#!/usr/bin/env python3
"""Figures for the v0.18 Phase 1 estimator report.

RNGLab v0.18 (written by Claude Code, 2026-11-05). Reads results/gates.json,
results/ar1_table.csv and results/ising_table.csv (from analyze.py) and writes
figures/*.png. Every plotted number comes from those files.

    MPLBACKEND=Agg PYTHONPATH=src python research/v018-estimators/plot.py
"""
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
# Validated categorical order (dataviz reference palette, light mode), fixed per entity.
COLOR = {"BM": "#2a78d6", "OBM": "#eb6834", "SV-B": "#1baf7a", "SV-TH": "#eda100", "GEYER": "#e87ba4",
         "P": "#2a78d6", "S": "#eb6834", "G": "#1baf7a"}
MARKER = {"BM": "o", "OBM": "s", "SV-B": "^", "SV-TH": "D", "GEYER": "v", "P": "o", "S": "s", "G": "^"}
INK, MUTED, GRID, BAND = "#1f1f1e", "#5c5b55", "#e4e3de", "#ecebe6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "lines.linewidth": 2,
                     "legend.frameon": False, "savefig.dpi": 160, "figure.facecolor": "white"})


def rows(name):
    with (HERE / "results" / name).open() as handle:
        return list(csv.DictReader(handle))


def num(row, key):
    value = row.get(key, "")
    return float(value) if value not in ("", "None") else np.nan


def band(ax, lo, hi, label):
    ax.axhspan(lo, hi, color=BAND, zorder=0, label=label)


def ar1_coverage(table):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)
    for ax, b in zip(axes, ("32", "512")):
        band(ax, 0.92, 0.98, "gate band [0.92, 0.98]")
        for method in ("BM", "OBM", "SV-B", "SV-TH", "GEYER"):
            pts = sorted((num(r, "tau_n"), num(r, "coverage"), num(r, "cp95_lo"), num(r, "cp95_hi"))
                         for r in table if r["method"] == method and r["n"] == "8192" and (r["b"] == b or method == "GEYER"))
            if not pts:
                continue
            tau, cov, lo, hi = map(np.array, zip(*pts))
            ax.errorbar(tau, cov, yerr=[cov - lo, hi - cov], color=COLOR[method], marker=MARKER[method],
                        markersize=6, capsize=3, label=method if method != "GEYER" else "GEYER (no b)")
        ax.set_xscale("log")
        ax.set_title(f"n = 8192, b = {b}", color=INK)
        ax.set_xlabel("integrated autocorrelation time τ (finite-n)")
    axes[0].set_ylabel("coverage of nominal 95% interval")
    axes[0].set_ylim(0, 1.0)
    axes[1].legend(loc="lower left")
    fig.suptitle("AR(1) control: interval coverage against τ (Clopper–Pearson 95% bars)", color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "ar1_coverage.png")
    plt.close(fig)


def ar1_variance_ratio(table):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    ax.axhspan(0.75, 1.25, color=BAND, zorder=0, label="gate band [0.75, 1.25]")
    for method in ("BM", "OBM", "SV-B", "SV-TH"):
        sel = [r for r in table if r["method"] == method]
        x = np.array([num(r, "b") / num(r, "tau_n") for r in sel])
        y = np.array([num(r, "var_ratio") for r in sel])
        order = np.argsort(x)
        ax.plot(x[order], y[order], linestyle="none", color=COLOR[method], marker=MARKER[method],
                markersize=6, label=method)
    grid = np.logspace(-2, 3.5, 200)
    ax.plot(grid, np.clip(1 - 1 / (2 * grid), 0, None), color=MUTED, linewidth=1.2, linestyle="--",
            label="BM first-order AR(1): 1 − τ/(2b)")
    ax.set_xscale("log")
    ax.set_ylim(-0.05, 1.4)
    ax.set_xlabel("batch or truncation length / τ")
    ax.set_ylabel("mean estimated / exact n·Var(x̄)")
    ax.set_title("AR(1): variance-estimator bias against b/τ (all φ, n, b)", color=INK)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "ar1_variance_ratio.png")
    plt.close(fig)


def diagnostic(gates):
    entries = [(k, g) for k, g in gates.items() if k.startswith("L-D[") and k.endswith(",BM]")]
    labels, abstain, cov, lo, hi = [], [], [], [], []
    for key, g in sorted(entries, key=lambda kv: (int(kv[0].split(",")[1]), float(kv[0][4:].split(",")[0]))):
        phi, n = key[4:].split(",")[:2]
        counts = g["status_counts"]
        total = sum(counts.values())
        labels.append(f"φ={phi}\nn={n}")
        abstain.append(1 - g["non_abstaining"] / total)
        c = g["coverage"]
        cov.append(c["coverage"] if c["coverage"] is not None else np.nan)
        lo.append(c["cp95"][0] if c["cp95"] else np.nan)
        hi.append(c["cp95"][1] if c["cp95"] else np.nan)
    x = np.arange(len(labels))
    fig, axes = plt.subplots(2, 1, figsize=(9, 5.5), sharex=True)
    axes[0].bar(x, abstain, color=MUTED, width=0.6)
    axes[0].set_ylabel("fraction abstaining\n(RUN_TOO_SHORT)")
    axes[0].set_ylim(0, 1.05)
    band(axes[1], 0.92, 0.98, "gate band")
    cov, lo, hi = map(np.array, (cov, lo, hi))
    axes[1].errorbar(x, cov, yerr=[cov - lo, hi - cov], color=COLOR["BM"], marker="o", linestyle="none", capsize=3,
                     label="adaptive BM, non-abstaining replicas")
    axes[1].set_ylabel("coverage")
    axes[1].set_ylim(0.85, 1.0)
    axes[1].set_xticks(x, labels)
    axes[1].legend(loc="lower left")
    fig.suptitle("Batch-length diagnostic (b ≥ 10·τ̂, ≥ 16 batches): abstention and coverage", color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "diagnostic.png")
    plt.close(fig)


def susceptibility(table, gates):
    cells = sorted(table, key=lambda r: (int(r["L"]), float(r["T"])))
    labels = [f"L{r['L']} T{float(r['T']):.2f}" for r in cells]
    x = np.arange(len(cells))
    fig, axes = plt.subplots(2, 1, figsize=(9.5, 6), sharex=True)
    for offset, (name, key) in zip((-0.2, 0.0, 0.2), (("P", "P_bias"), ("S", "S_bias"), ("G", "G_bias"))):
        axes[0].plot(x + offset, [num(r, key) for r in cells], linestyle="none", marker=MARKER[name],
                     color=COLOR[name], markersize=7, label={"P": "P plug-in", "S": "S symmetry-aware", "G": "G 16-chain jackknife"}[name])
    axes[0].axhline(0, color=MUTED, linewidth=1)
    axes[0].set_yscale("symlog", linthresh=0.01)
    axes[0].set_ylabel("mean relative bias")
    axes[0].legend(loc="lower left")
    band(axes[1], 0.92, 0.98, "gate band")
    for offset, name in ((-0.1, "S"), (0.1, "G")):
        cov = []
        for r in cells:
            c = gates[f"S-{name}[{r['cell']}]"]["coverage"]
            cov.append((c["coverage"], c["cp95"][0], c["cp95"][1]))
        c, lo, hi = map(np.array, zip(*cov))
        axes[1].errorbar(x + offset, c, yerr=[c - lo, hi - c], linestyle="none", marker=MARKER[name],
                         color=COLOR[name], capsize=3, markersize=6, label=name)
    axes[1].set_ylabel("coverage of nominal 95%")
    axes[1].set_xticks(x, labels, rotation=30)
    axes[1].legend(loc="lower left")
    fig.suptitle("Ising susceptibility: bias (top, symlog) and interval coverage (bottom)", color=INK)
    fig.tight_layout()
    fig.savefig(FIG / "susceptibility.png")
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    gates = json.loads((HERE / "results" / "gates.json").read_text())["gates"]
    ar1 = rows("ar1_table.csv")
    ar1_coverage(ar1)
    ar1_variance_ratio(ar1)
    diagnostic(gates)
    susceptibility(rows("ising_table.csv"), gates)
    print("wrote", sorted(p.name for p in FIG.glob("*.png")))


if __name__ == "__main__":
    main()
