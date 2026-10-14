"""Publication-readable PNG summaries; raw trial data remain available."""
from collections import defaultdict
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import t


def plot_all(config, output, summary, diagnostics, windows, histograms):
    output = Path(output) / "plots"
    output.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False,
                         "axes.grid": True, "grid.alpha": 0.2, "savefig.dpi": 160})
    names = [c["name"] for c in config["scenarios"]]
    preferred = ["synthetic_baseline", "prng", "csprng_control", "bias_severe", "autocorrelation_severe",
                 "periodic_severe", "hybrid_frozen", "csprng_initial_frozen"]
    chosen = [s for s in preferred if s in names] or names[:8]
    workloads = config["workloads"]
    colors = {s: plt.get_cmap("tab10")(i % 10) for i, s in enumerate(chosen)}
    fixed = {r["scenario"] for r in summary if r.get("distinct_initial_keys") == 1}
    label = lambda name: name + (" †" if name in fixed else "")
    number = lambda value: np.nan if value is None else value

    def finish(fig, name):
        fig.savefig(output / name, bbox_inches="tight")
        plt.close(fig)

    fig, axes = plt.subplots(2, len(workloads), figsize=(7*len(workloads), 8), squeeze=False, layout="constrained")
    for col, workload in enumerate(workloads):
        for name in chosen:
            data = sorted([r for r in summary if r["scenario"] == name and r["workload"] == workload], key=lambda r:r["n"])
            n = [r["n"] for r in data]
            axes[0, col].plot(n, [r["mean_estimate"] for r in data], label=label(name), color=colors[name])
            axes[0, col].fill_between(n, [number(r["mean_ci95_low"]) for r in data], [number(r["mean_ci95_high"]) for r in data], alpha=.12, color=colors[name])
            axes[1, col].plot(n, [max(r["rmse"], 1e-12) for r in data], marker=".", color=colors[name])
        target = next(r["target"] for r in summary if r["workload"] == workload)
        axes[0, col].axhline(target, color="black", ls="--", lw=1, label="analytic target")
        base = next((r for r in summary if r["scenario"] == "synthetic_baseline" and r["workload"] == workload), None)
        if base:
            n = np.array(config["checkpoints"])
            axes[1, col].plot(n, base["rmse"]*np.sqrt(base["n"]/n), "k:", label="N⁻¹ᐟ² reference")
        axes[0, col].set(title=f"{workload}: mean estimate · 95% trial-mean interval", xscale="log", ylabel="Estimate")
        axes[1, col].set(title="Root mean squared error across trials", xscale="log", yscale="log", xlabel="Monte Carlo samples", ylabel="RMSE")
        if base:
            axes[1, col].legend(fontsize=8)
    axes[0, 0].legend(fontsize=8, loc="best")
    fig.suptitle("Source defects and Monte Carlo convergence" + ("\n† One initial key across trials; no confidence interval" if fixed else ""), fontsize=14)
    finish(fig, "convergence.png")

    fig, axes = plt.subplots(1, len(workloads), figsize=(7*len(workloads), max(6, len(names)*.29)), squeeze=False, layout="constrained")
    for col, workload in enumerate(workloads):
        data = [next(r for r in summary if r["scenario"] == name and r["workload"] == workload and r["n"] == config["samples"]) for name in names]
        ax = axes[0, col]
        means = np.array([r["mean_abs_error"] for r in data])
        half = np.array([number(r["mae_ci95_half"]) for r in data])
        ax.errorbar(means, np.arange(len(names)), xerr=np.vstack((np.minimum(half, means*.999), half)), fmt="o", capsize=2)
        ax.set(yticks=np.arange(len(names)), yticklabels=[label(s) for s in names], xscale="log", xlabel="Mean absolute error · approximate 95% interval", title=f"{workload} · N={config['samples']:,}")
        ax.invert_yaxis()
    fig.suptitle(f"All providers · {config['trials']} trials per configuration" + ("\n† One initial key; no confidence interval" if fixed else ""), fontsize=14)
    finish(fig, "provider_comparison.png")

    faults = ["bias", "freeze", "periodic", "autocorrelation"]
    if any(c.get("mechanism") for c in config["scenarios"]):
        fig, axes = plt.subplots(1, 4, figsize=(16, 4.5), layout="constrained")
        for ax, fault in zip(axes, faults):
            for workload in workloads:
                configs = [c for c in config["scenarios"] if c.get("mechanism") == fault]
                configs = sorted(configs, key=lambda c:c["strength"])
                baseline = next((r for r in summary if r["scenario"] == "synthetic_baseline" and r["n"] == config["samples"] and r["workload"] == workload), None)
                data = [next(r for r in summary if r["scenario"] == c["name"] and r["n"] == config["samples"] and r["workload"] == workload) for c in configs]
                x = [c["strength"] for c in configs]
                if baseline:
                    data, x = [baseline]+data, [0]+x
                ax.errorbar(x, [r["mean_abs_error"] for r in data], yerr=[r["mae_ci95_half"] for r in data], marker="o", capsize=3, label=workload)
            ax.set(title=fault, xlabel="Mechanism-specific strength", ylabel="Mean absolute error")
            ax.legend()
        fig.suptitle("Severity parameters differ by mechanism; bars are approximate 95% trial intervals")
        finish(fig, "severity.png")

    # Single final window and one workload prevent pseudo-replication across windows.
    final_diag = [r for r in diagnostics if r["stream"] == "output" and r["workload"] == workloads[0] and r["end_sample"] == config["samples"]]
    fig, axes = plt.subplots(1, 3, figsize=(15, max(6, .29*len(names))), sharey=True, layout="constrained")
    for col, (field, title, reject) in enumerate([("bit_bias", "Mean output bit bias", False), ("lag1_rho", "Mean lag-1 correlation", False), ("monobit_p_iid", "Monobit rejection fraction (α=.05)", True)]):
        means, errs = [], []
        for name in names:
            vals = [r[field] for r in final_diag if r["scenario"] == name and r[field] is not None]
            vals = np.array([v < .05 for v in vals] if reject else vals, dtype=float)
            means.append(float(vals.mean()) if len(vals) else np.nan)
            errs.append(float(t.ppf(.975, len(vals)-1)*vals.std(ddof=1)/np.sqrt(len(vals))) if len(vals)>1 and not reject else 0)
        axes[col].barh(np.arange(len(names)), means, xerr=errs, height=.6)
        axes[col].set(title=title, yticks=np.arange(len(names)), yticklabels=[label(s) for s in names])
        if reject:
            axes[col].axvline(.05, color="black", ls="--", lw=1)
            axes[col].set_xlim(0, 1.05)
    axes[0].invert_yaxis()
    fig.suptitle("Final output window · IID-null diagnostics are not entropy/security tests" + ("\n† Repeated stream; rejection fraction has no independent-trial interpretation" if fixed else ""), fontsize=13)
    finish(fig, "diagnostics.png")

    fig, axes = plt.subplots(1, len(workloads), figsize=(7*len(workloads), 4.5), squeeze=False, layout="constrained")
    hist_group = defaultdict(lambda:np.zeros(16))
    for row in histograms:
        hist_group[(row["scenario"], row["workload"])][row["bin"]] += row["count"]
    for col, workload in enumerate(workloads):
        for name in chosen:
            hist = hist_group[(name, workload)]
            axes[0, col].stairs(hist/hist.sum()*16, np.linspace(0, 1, 17), color=colors[name], label=name)
        axes[0, col].axhline(1, color="black", ls="--", lw=1)
        axes[0, col].set(title=workload, xlabel="Uniform variate", ylabel="Density (all trials and windows)")
    axes[0, 0].legend(fontsize=8)
    fig.suptitle("Distribution of values supplied to the simulations")
    finish(fig, "histograms.png")

    timed = [c for c in config["scenarios"] if c.get("timeline")]
    if timed:
        fig, axes = plt.subplots(2, len(workloads), figsize=(7*len(workloads), 8), squeeze=False, layout="constrained")
        for col, workload in enumerate(workloads):
            for row_index, provider_kind in enumerate(["synthetic", "hybrid"]):
                ax = axes[row_index, col]
                for c in timed:
                    if c["provider"] != provider_kind:
                        continue
                    points = defaultdict(list)
                    for row in windows:
                        if row["scenario"] == c["name"] and row["workload"] == workload:
                            points[(row["start_sample"]+row["end_sample"])/2].append(row["window_estimate"])
                    x = sorted(points)
                    means = np.array([np.mean(points[n]) for n in x])
                    half = np.array([t.ppf(.975, len(points[n])-1)*np.std(points[n], ddof=1)/np.sqrt(len(points[n])) for n in x])
                    line, = ax.plot(x, means, marker=".", label=c["name"])
                    ax.fill_between(x, means-half, means+half, color=line.get_color(), alpha=.15)
                    f = c["faults"][0]
                    dim = 2 if workload == "pi" else 1
                    onset = f["onset_bit"]/(32*dim) if provider_kind == "synthetic" else (f["onset_bit"]//256+1)*c["reseed_interval"]
                    ax.axvline(onset, color=line.get_color(), ls="--", alpha=.6)
                if provider_kind == "hybrid":
                    interval = next(c["reseed_interval"] for c in timed if c["provider"] == "hybrid")
                    for n in range(interval, config["samples"], interval):
                        ax.axvline(n, color="gray", lw=.5, alpha=.25)
                target = next(r["target"] for r in summary if r["workload"] == workload)
                ax.axhline(target, color="black", ls=":")
                ax.set(title=f"{workload} · {provider_kind}", xlabel="Sample index (window midpoint)", ylabel="Window estimate · 95% trial-mean interval")
                ax.legend(fontsize=8)
        fig.suptitle("Fault onset (dashed) and reseeds (thin gray) · window estimates reveal time variation")
        finish(fig, "onset_reseeding.png")
