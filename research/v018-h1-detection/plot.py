"""Figures for the Stage C (H1) report, from the analysis outputs.

RNGLab v0.18 (written by OpenAI Codex, gpt-6-astra, 2026-11-07; reviewed by Claude Code)
Reads results/gates.json, results/secondary.json, results/recount.json and,
if present, results/run_manifest.json. Writes figures/fig1_rates_by_severity.png
and figures/fig2_latency.png. Every plotted interval comes from gates.json.

    MPLBACKEND=Agg python research/v018-h1-detection/plot.py
"""
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.transforms import ScaledTranslation
import numpy as np

HERE = Path(__file__).resolve().parent
FIG = HERE / "figures"
FAMILIES = ("SHIFT", "LAG1", "DUP", "GFSR")
# Validated categorical order, as in the API-validation campaign.
COLOR = {"wrong": "#2a78d6", "C_alarm": "#eb6834", "A_alarm": "#1baf7a"}
MARKER = {"wrong": "o", "C_alarm": "s", "A_alarm": "^"}
LABEL = {"wrong": "Wrong rate", "C_alarm": "C-alarm rate", "A_alarm": "A-alarm rate"}
INK, MUTED, GRID, BAND = "#1f1f1e", "#5c5b55", "#e4e3de", "#ecebe6"
plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "lines.linewidth": 2,
                     "legend.frameon": False, "savefig.dpi": 160, "figure.facecolor": "white"})


def cell_key(item):
    name, cell = item
    value = cell["severity"]
    return FAMILIES.index(cell["family"]), tuple(value) if isinstance(value, list) else (value,), name


def rates(g):
    fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharey=True)
    cells = sorted(g["cells"].items(), key=cell_key)
    for ax, family in zip(axes.flat, FAMILIES):
        group = [c for _, c in cells if c["family"] == family]
        ax.set_title(family, color=INK)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Lag pair (p,q)" if family == "GFSR" else "Severity")
        ax.set_ylabel("Fraction of runs (95% interval)")
        if family != "GFSR":
            ax.set_xscale("log")
        if not group:
            ax.text(0.5, 0.5, "No cells in input", transform=ax.transAxes, ha="center", color=MUTED)
            continue
        x = np.arange(len(group)) if family == "GFSR" else np.array([c["severity"] for c in group])
        for offset, key in zip((-5, 0, 5), COLOR):
            y = np.array([c[key]["k"] / (c["wrong"]["n"] if key == "wrong" else c["runs"]) for c in group])
            lo, hi = np.array([c[key]["cp95"] for c in group]).T
            # Small display offsets keep coincident intervals visible at the same severity.
            transform = ax.transData + ScaledTranslation(offset / 72, 0, fig.dpi_scale_trans)
            ax.errorbar(x, y, yerr=[y - lo, hi - y], color=COLOR[key], marker=MARKER[key],
                        markersize=6, linestyle="none", capsize=3, label=LABEL[key], zorder=3,
                        clip_on=False, transform=transform)
        labels = [",".join(map(str, c["severity"])) if family == "GFSR" else f"{c['severity']:.3g}" for c in group]
        ax.set_xticks(x, labels, rotation=25 if len(group) > 4 else 0)
        ax.minorticks_off()
        if family == "GFSR":
            ax.set_xlim(-0.5, len(group) - 0.5)
        else:
            pad = (x[-1] / x[0]) ** 0.08 if x[-1] > x[0] else 1.5
            ax.set_xlim(x[0] / pad, x[-1] * pad)
    handles = [plt.Line2D([], [], color=COLOR[key], marker=MARKER[key], linestyle="none", label=LABEL[key]) for key in COLOR]
    fig.legend(handles=handles, loc="lower center", ncol=3)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(FIG / "fig1_rates_by_severity.png")
    plt.close(fig)


def latency(g):
    cells = [(name, c["latency"]["sweeps"]) for name, c in sorted(g["cells"].items(), key=cell_key)
             if c["latency"]["sweeps"]["both"] >= 2 and c["latency"]["sweeps"]["ci95"] is not None
             and c["latency"]["sweeps"]["median_A_minus_C"] is not None]
    fig, ax = plt.subplots(figsize=(max(7, len(cells) * 0.45), 4.8))
    ax.axhline(0, color=MUTED, linewidth=1, linestyle="--", zorder=2)
    if cells:
        x = np.arange(len(cells))
        y = np.array([c["median_A_minus_C"] for _, c in cells])
        lo, hi = np.array([c["ci95"] for _, c in cells]).T
        ax.errorbar(x, y, yerr=[y - lo, hi - y], color=COLOR["wrong"], marker="o", linestyle="none",
                    capsize=4, zorder=3)
        ax.set_xticks(x, [name for name, _ in cells], rotation=45, ha="right")
        ax.set_xlim(-0.5, len(cells) - 0.5)
    else:
        ax.set_xticks([])
        ax.text(0.5, 0.6, "No cells with at least two runs where both arms detect", transform=ax.transAxes,
                ha="center", color=MUTED)
    ax.set_xlabel("Fault cell")
    ax.set_ylabel("Median A-minus-C detection delay\n(sweeps, 95% interval)")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_latency.png")
    plt.close(fig)


def main():
    inputs = {name: json.loads((HERE / "results" / f"{name}.json").read_text())
              for name in ("gates", "secondary", "recount")}
    manifest = HERE / "results" / "run_manifest.json"
    if manifest.exists():
        inputs["run_manifest"] = json.loads(manifest.read_text())
    FIG.mkdir(exist_ok=True)
    rates(inputs["gates"])
    latency(inputs["gates"])
    print("wrote figures/fig1_rates_by_severity.png, figures/fig2_latency.png")


if __name__ == "__main__":
    main()
