#!/usr/bin/env python3
"""Figures for the v0.18 Phase 2 API-validation report.

RNGLab v0.18 (written by Claude Code, 2026-11-06). Reads results/gates.json and
results/secondary.json (from analyze.py) and the Phase II confirmation decision
file; writes figures/*.png. Every plotted number comes from those files.

    MPLBACKEND=Agg python research/v018-api-validation/plot.py
"""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
PHASE2 = HERE.parents[1] / "research" / "phase2-20261030" / "experiment" / "confirmation-v1" / "primary_decision.json"
CONDITIONS = ("local", "refresh", "iid")
LABEL = {"local": "local Metropolis", "refresh": "local + centroid refresh", "iid": "IID Gaussian paths"}
PHASE2_KEY = {"local": "local_stationary", "refresh": "centroid_refresh", "iid": "iid_gaussian"}
# Validated categorical order (dataviz reference palette, light mode), fixed per entity, as in Phase 1.
COLOR = {"local": "#2a78d6", "refresh": "#eb6834", "iid": "#1baf7a", "api": "#2a78d6", "phase2": "#eb6834",
         "centroid": "#2a78d6", "shape": "#eb6834"}
MARKER = {"local": "o", "refresh": "s", "iid": "^", "api": "o", "phase2": "s"}
INK, MUTED, GRID, BAND = "#1f1f1e", "#5c5b55", "#e4e3de", "#ecebe6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "lines.linewidth": 2,
                     "legend.frameon": False, "savefig.dpi": 160, "figure.facecolor": "white"})


def coverage(g, phase2):
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.axhspan(0.92, 0.98, color=BAND, zorder=0, label="gate band [0.92, 0.98]")
    ax.axhline(0.95, color=MUTED, linewidth=1, linestyle=":", zorder=1)
    x = np.arange(len(CONDITIONS))
    for offset, key in ((-0.12, "api"), (0.12, "phase2")):
        if key == "api":
            pts = [(g["primary"][c]["coverage"], *g["primary"][c]["cp95"]) for c in CONDITIONS]
            label = "package API, R = 4,096 (this campaign)"
        else:
            pts = [(phase2[PHASE2_KEY[c]]["row"]["coverage"], *phase2[PHASE2_KEY[c]]["row"]["coverage_ci95"]) for c in CONDITIONS]
            label = "Phase II standalone, R = 1,024"
        cov, lo, hi = map(np.array, zip(*pts))
        ax.errorbar(x + offset, cov, yerr=[cov - lo, hi - cov], color=COLOR[key], marker=MARKER[key], markersize=7,
                    linestyle="none", capsize=4, label=label, zorder=3)
    ax.set_xticks(x, [LABEL[c] for c in CONDITIONS])
    ax.set_ylim(0.74, 1.0)
    ax.set_ylabel("coverage of CV BM512 95% interval\n(Clopper–Pearson 95% CI)")
    ax.set_title("Coverage through the package API matches Phase II", color=INK)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "fig1_coverage_vs_phase2.png")
    plt.close(fig)


def batch_ratio(s):
    fig, ax = plt.subplots(figsize=(7, 4.6))
    ax.axhspan(0.75, 1.25, color=BAND, zorder=0, label="variance-ratio gate band")
    ax.axhline(1.0, color=MUTED, linewidth=1, linestyle=":", zorder=1)
    batches = (32, 128, 512)
    for c in CONDITIONS:
        ratios = [s[c]["cv_methods"][f"BM{b}"]["variance_ratio"] for b in batches]
        ax.plot(batches, ratios, color=COLOR[c], marker=MARKER[c], markersize=7, label=LABEL[c], zorder=3)
        lift = {"local": -3, "refresh": -11, "iid": 5}[c]  # keep the two end labels near 1.0 apart
        ax.annotate(f"{ratios[-1]:.2f}", (batches[-1], ratios[-1]), textcoords="offset points", xytext=(8, lift),
                    color=INK, fontsize=9)
    ax.set_xscale("log", base=2)
    ax.set_xticks(batches, [str(b) for b in batches])
    ax.set_xlim(24, 800)
    ax.set_ylim(0, 1.3)
    ax.set_xlabel("batch length b (n = 8192)")
    ax.set_ylabel("mean BM variance / empirical variance\nof replica means")
    ax.set_title("After the refresh, short batches still miss a slow component", color=INK)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.17), ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "fig2_variance_ratio_vs_batch.png")
    plt.close(fig)


def decomposition(s):
    fig, ax = plt.subplots(figsize=(7, 4))
    x = np.arange(len(CONDITIONS))
    width = 0.34
    for offset, part, label in ((-width / 2 - 0.01, "centroid", "centroid part  Var(mean c̄²/2)"),
                                (width / 2 + 0.01, "shape", "shape part  Var(mean r_g²)")):
        values = [s[c]["variance_decomposition"][f"{part}_part"] for c in CONDITIONS]
        bars = ax.bar(x + offset, values, width, color=COLOR[part], label=label, zorder=3)
        for bar, value in zip(bars, values):
            ax.annotate(f"{value:.1e}", (bar.get_x() + bar.get_width() / 2, value), textcoords="offset points",
                        xytext=(0, 3), ha="center", color=INK, fontsize=8)
    ax.set_yscale("log")
    ax.set_ylim(5e-7, 3e-2)
    ax.set_xticks(x, [LABEL[c] for c in CONDITIONS])
    ax.set_ylabel("between-replica variance of the time average")
    ax.set_title("The refresh removes the centroid variance; the shape part remains", color=INK)
    ax.legend(loc="upper right")
    fig.tight_layout()
    fig.savefig(FIG / "fig3_variance_decomposition.png")
    plt.close(fig)


def main():
    FIG.mkdir(exist_ok=True)
    g = json.loads((HERE / "results" / "gates.json").read_text())
    s = json.loads((HERE / "results" / "secondary.json").read_text())
    phase2 = json.loads(PHASE2.read_text())["gates"]
    coverage(g, phase2)
    batch_ratio(s)
    decomposition(s)
    print("wrote", ", ".join(sorted(p.name for p in FIG.glob("*.png"))))


if __name__ == "__main__":
    main()
