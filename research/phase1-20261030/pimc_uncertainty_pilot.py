#!/usr/bin/env python3
"""Reproducible, exploratory scalar replay of RNGLab's harmonic PIMC control.

Uses only the Python standard library. It independently reimplements the
repository's M=16, beta=2, delta=0.4 local Metropolis kernel and energy
estimators so the pilot can run where NumPy/SciPy are unavailable. It is not
the RNGLab NumPy implementation and is not confirmatory evidence.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import platform
import random
import statistics
import sys
import time
from pathlib import Path


OUT = Path(__file__).resolve().parent
REPLICAS = 64
SWEEPS = 12_288
BEADS = 16
BETA = 2.0
DELTA = 0.4
BATCHES = (32, 128, 512)
SEED_NAMESPACE = "rnglab-pimc-pilot-20261030-scalar-v1"
NOMINAL = 0.95


def _beta_fraction(a: float, b: float, x: float) -> float:
    """Continued fraction for the incomplete beta (Lentz/Numerical Recipes)."""
    max_iter, eps, floor = 300, 3.0e-14, 1.0e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < floor:
        d = floor
    d = 1.0 / d
    h = d
    for m in range(1, max_iter + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if abs(d) < floor:
            d = floor
        c = 1.0 + aa / c
        if abs(c) < floor:
            c = floor
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if abs(d) < floor:
            d = floor
        c = 1.0 + aa / c
        if abs(c) < floor:
            c = floor
        d = 1.0 / d
        change = d * c
        h *= change
        if abs(change - 1.0) < eps:
            return h
    raise ArithmeticError("incomplete-beta continued fraction did not converge")


def regularized_beta(x: float, a: float, b: float) -> float:
    if x <= 0.0:
        return 0.0
    if x >= 1.0:
        return 1.0
    front = math.exp(
        math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
        + a * math.log(x) + b * math.log1p(-x)
    )
    if x < (a + 1.0) / (a + b + 2.0):
        return front * _beta_fraction(a, b, x) / a
    return 1.0 - front * _beta_fraction(b, a, 1.0 - x) / b


def student_t_cdf(value: float, df: int) -> float:
    x = df / (df + value * value)
    tail = 0.5 * regularized_beta(x, df / 2.0, 0.5)
    return 1.0 - tail if value >= 0.0 else tail


def student_t_quantile(probability: float, df: int) -> float:
    if not 0.0 < probability < 1.0 or df < 1:
        raise ValueError("probability must be inside (0,1) and df positive")
    if probability == 0.5:
        return 0.0
    if probability < 0.5:
        return -student_t_quantile(1.0 - probability, df)
    lo, hi = 0.0, 64.0
    for _ in range(100):
        mid = (lo + hi) / 2.0
        if student_t_cdf(mid, df) < probability:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def finite_bead_energy(beads: int, beta: float) -> float:
    z = beta / (2.0 * beads)
    return 0.5 / (math.tanh(beads * math.asinh(z)) * math.sqrt(1.0 + z * z))


def log_partition(beads: int, beta: float) -> float:
    return -math.log(2.0 * math.sinh(beads * math.asinh(beta / (2.0 * beads))))


def wilson_interval(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = k / n
    den = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / den
    radius = z / den * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n))
    return max(0.0, center - radius), min(1.0, center + radius)


def run_chain(replica: int) -> tuple[list[float], list[float], dict]:
    seed_text = f"{SEED_NAMESPACE}/{replica:03d}"
    seed_digest = hashlib.sha256(seed_text.encode("ascii")).hexdigest()
    rng = random.Random(int(seed_digest[:32], 16))
    x = [rng.random() - 0.5 for _ in range(BEADS)]
    dt = BETA / BEADS
    td_trace: list[float] = []
    cv_trace: list[float] = []
    accepted = 0

    for _sweep in range(SWEEPS):
        for j in range(BEADS):
            up, ua = rng.random(), rng.random()
            old = x[j]
            new = old + DELTA * (2.0 * up - 1.0)
            left, right = x[(j - 1) % BEADS], x[(j + 1) % BEADS]
            ds = (
                (new - left) ** 2 + (right - new) ** 2
                - (old - left) ** 2 - (right - old) ** 2
            ) / (2.0 * dt) + dt * (new * new - old * old) / 2.0
            if ds <= 0.0 or ua < math.exp(-ds):
                x[j] = new
                accepted += 1

        spring_sum = sum((x[(j + 1) % BEADS] - x[j]) ** 2 for j in range(BEADS))
        x2 = sum(value * value for value in x) / BEADS
        centroid = sum(x) / BEADS
        td_trace.append(BEADS / (2.0 * BETA) - BEADS * spring_sum / (2.0 * BETA**2) + 0.5 * x2)
        cv_trace.append(1.0 / (2.0 * BETA) + x2 - 0.5 * centroid * centroid)

    return td_trace, cv_trace, {
        "replica": replica,
        "seed_namespace": SEED_NAMESPACE,
        "seed_sha256": seed_digest,
        "acceptance_rate": accepted / (SWEEPS * BEADS),
    }


def batch_interval(values: list[float], reference: float, batch: int) -> dict:
    count = len(values) // batch
    if count < 4:
        return {"batches": count, "covered": None, "status": "too_few_batches"}
    used = values[:count * batch]
    means = [statistics.fmean(used[i * batch:(i + 1) * batch]) for i in range(count)]
    center = statistics.fmean(means)
    se = statistics.stdev(means) / math.sqrt(count)
    critical = student_t_quantile(0.975, count - 1)
    lo, hi = center - critical * se, center + critical * se
    return {
        "batches": count,
        "degrees_of_freedom": count - 1,
        "t_critical_975": critical,
        "estimate": center,
        "mcse": se,
        "lower": lo,
        "upper": hi,
        "covered": lo <= reference <= hi,
        "status": "estimated",
    }


def make_svg(rows: list[dict], path: Path) -> None:
    width, height = 1100, 760
    panel_w, panel_h = 470, 285
    starts = [(85, 90), (600, 90), (85, 415), (600, 415)]
    colors = {32: "#275DAD", 128: "#D97706", 512: "#188977"}
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#18212b}.small{font-size:13px}.axis{stroke:#6b7280;stroke-width:1}.grid{stroke:#d8dee8;stroke-width:1}.nom{stroke:#a33;stroke-width:2;stroke-dasharray:7 5}</style>',
        '<text x="550" y="34" text-anchor="middle" font-size="22" font-weight="bold">PIMC healthy-chain interval coverage — exploratory pilot</text>',
        '<text x="550" y="58" text-anchor="middle" font-size="14">64 seeded replicas per panel; points are observed coverage, bars are Wilson 95% intervals</text>',
    ]
    panels = [("short", "TD"), ("short", "CV"), ("long", "TD"), ("long", "CV")]
    for (window, observable), (x0, y0) in zip(panels, starts):
        left, top = x0 + 54, y0 + 30
        plot_w, plot_h = panel_w - 90, panel_h - 78
        parts.append(f'<text x="{x0 + panel_w/2}" y="{y0 + 4}" text-anchor="middle" font-size="17" font-weight="bold">{window.title()} window — {observable}</text>')
        for tick in (0.0, 0.25, 0.5, 0.75, 1.0):
            yy = top + (1.0 - tick) * plot_h
            parts.append(f'<line class="grid" x1="{left}" y1="{yy:.1f}" x2="{left+plot_w}" y2="{yy:.1f}"/>')
            parts.append(f'<text class="small" x="{left-9}" y="{yy+4:.1f}" text-anchor="end">{tick:.2f}</text>')
        nom_y = top + (1.0 - NOMINAL) * plot_h
        parts.append(f'<line class="nom" x1="{left}" y1="{nom_y:.1f}" x2="{left+plot_w}" y2="{nom_y:.1f}"/>')
        parts.append(f'<line class="axis" x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}"/>')
        parts.append(f'<line class="axis" x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}"/>')
        for idx, batch in enumerate(BATCHES):
            row = next(r for r in rows if r["window"] == window and r["observable"] == observable and r["batch_size"] == batch)
            xx = left + (idx + 1) * plot_w / 4
            p = row["coverage"]
            low, high = row["wilson_95"]
            y = top + (1.0 - p) * plot_h
            ylow, yhigh = top + (1.0 - low) * plot_h, top + (1.0 - high) * plot_h
            color = colors[batch]
            parts.extend([
                f'<line x1="{xx:.1f}" y1="{yhigh:.1f}" x2="{xx:.1f}" y2="{ylow:.1f}" stroke="{color}" stroke-width="3"/>',
                f'<line x1="{xx-6:.1f}" y1="{yhigh:.1f}" x2="{xx+6:.1f}" y2="{yhigh:.1f}" stroke="{color}" stroke-width="2"/>',
                f'<line x1="{xx-6:.1f}" y1="{ylow:.1f}" x2="{xx+6:.1f}" y2="{ylow:.1f}" stroke="{color}" stroke-width="2"/>',
                f'<circle cx="{xx:.1f}" cy="{y:.1f}" r="6" fill="{color}"/>',
                f'<text class="small" x="{xx:.1f}" y="{top+plot_h+22}" text-anchor="middle">{batch}</text>',
                f'<text class="small" x="{xx:.1f}" y="{y-10:.1f}" text-anchor="middle">{row["covered"]}/64</text>',
            ])
        parts.append(f'<text class="small" x="{left+plot_w/2}" y="{top+plot_h+43}" text-anchor="middle">Batch size (sweeps)</text>')
    parts.append('<line class="nom" x1="110" y1="735" x2="150" y2="735"/><text class="small" x="158" y="740">Nominal 95%</text>')
    for idx, batch in enumerate(BATCHES):
        x = 360 + idx * 150
        parts.append(f'<circle cx="{x}" cy="735" r="6" fill="{colors[batch]}"/><text class="small" x="{x+13}" y="740">batch {batch}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    started = time.monotonic()
    reference = finite_bead_energy(BEADS, BETA)
    h = 1.0e-5
    derivative_reference = -(log_partition(BEADS, BETA + h) - log_partition(BEADS, BETA - h)) / (2.0 * h)
    if abs(reference - derivative_reference) > 1.0e-9:
        raise AssertionError("finite-bead reference failed its partition-derivative check")

    # Check the numerical t quantile solver at the nominal probability before
    # using it in any interval calculations.
    quantile_check = {}
    for df in (5, 15, 23, 63, 95, 255):
        q = student_t_quantile(0.975, df)
        error = abs(student_t_cdf(q, df) - 0.975)
        quantile_check[str(df)] = {"quantile": q, "cdf_error": error}
        if error > 1.0e-12:
            raise AssertionError(f"Student-t quantile check failed for df={df}")

    chains: dict[int, dict[str, list[float]]] = {}
    chain_meta = []
    with (OUT / "pimc-uncertainty-pilot-traces.csv.gz").open("wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", filename="", mtime=0) as compressed:
            with io.TextIOWrapper(compressed, newline="", encoding="utf-8") as stream:
                writer = csv.writer(stream)
                writer.writerow(("replica", "sweep_1based", "thermodynamic_energy", "centroid_virial_energy"))
                for replica in range(REPLICAS):
                    td, cv, metadata = run_chain(replica)
                    chains[replica] = {"TD": td, "CV": cv}
                    chain_meta.append(metadata)
                    writer.writerows((replica, i + 1, repr(td[i]), repr(cv[i])) for i in range(SWEEPS))

    windows = {
        "short": (1024, 4096),  # 3072 retained sweeps
        "long": (4096, 12288),  # 8192 retained sweeps
    }
    rows = []
    for window, (start, stop) in windows.items():
        n_retained = stop - start
        for observable in ("TD", "CV"):
            for batch in BATCHES:
                covers, estimates, ses = [], [], []
                per_chain = []
                interval_result = None
                for replica in range(REPLICAS):
                    values = chains[replica][observable][start:stop]
                    result = batch_interval(values, reference, batch)
                    if result["status"] != "estimated":
                        raise AssertionError("unexpected insufficient batch count")
                    covers.append(bool(result["covered"]))
                    estimates.append(result["estimate"])
                    ses.append(result["mcse"])
                    per_chain.append(result)
                    interval_result = result
                covered = sum(covers)
                rows.append({
                    "window": window,
                    "sweeps_1based": [start + 1, stop],
                    "retained_sweeps_per_chain": n_retained,
                    "observable": observable,
                    "batch_size": batch,
                    "batch_count": interval_result["batches"],
                    "degrees_of_freedom": interval_result["degrees_of_freedom"],
                    "covered": covered,
                    "replicas": REPLICAS,
                    "coverage": covered / REPLICAS,
                    "wilson_95": list(wilson_interval(covered, REPLICAS)),
                    "mean_chain_estimate": statistics.fmean(estimates),
                    "median_chain_mcse": statistics.median(ses),
                    "reference": reference,
                    "scope": "same 64 trajectories across all windows/methods; Wilson interval describes replica coverage variability only",
                })

    with (OUT / "pimc-uncertainty-pilot-chain-summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=("replica", "seed_namespace", "seed_sha256", "acceptance_rate"))
        writer.writeheader()
        writer.writerows(chain_meta)

    make_svg(rows, OUT / "pimc-uncertainty-pilot-coverage.svg")
    report = {
        "study": "Exploratory independent scalar replay of RNGLab harmonic PIMC uncertainty control",
        "generated_utc_date": "2026-10-30",
        "classification": "exploratory pilot; not confirmatory; not an exact replay of NumPy implementation",
        "source_revision_reference": "rnglabexperimental main dd96adba23b8100aba9147d71fffc63932fc6b0d",
        "environment": {"python": sys.version, "platform": platform.platform(), "dependencies": "Python standard library only"},
        "configuration": {
            "replicas": REPLICAS,
            "sweeps_per_replica": SWEEPS,
            "beads": BEADS,
            "beta": BETA,
            "proposal_delta": DELTA,
            "initializer": "independent uniform bead coordinates on [-0.5, 0.5)",
            "proposal": "one bead per sweep in fixed order; delta*(2u-1), u~Python random.Random",
            "acceptance": "Metropolis with an independent sequential uniform draw from the same replica stream",
            "seed_rule": "first 128 bits of SHA256(f'{namespace}/{replica:03d}') fed to random.Random",
            "seed_namespace": SEED_NAMESPACE,
            "windows": {name: {"sweeps_1based": [a + 1, b], "retained": b - a} for name, (a, b) in windows.items()},
            "estimators": {
                "TD": "M/(2 beta) - M*sum((x[j+1]-x[j])^2)/(2 beta^2) + 0.5*mean(x^2)",
                "CV": "1/(2 beta) + mean(x^2) - 0.5*mean(x)^2",
            },
            "interval": "non-overlapping batch means; Student-t critical value; batches consume the window prefix",
            "batches": list(BATCHES),
            "minimum_batches": 4,
            "replica_coverage_interval": "Wilson 95% interval; descriptive variability, not a proof of coverage",
        },
        "reference": {
            "finite_bead_energy": reference,
            "continuum_energy": 0.5 / math.tanh(BETA / 2.0),
            "partition_function_derivative_check": derivative_reference,
            "absolute_difference": abs(reference - derivative_reference),
        },
        "student_t_quantile_checks": quantile_check,
        "coverage": rows,
        "trace_file": "pimc-uncertainty-pilot-traces.csv.gz",
        "elapsed_seconds": time.monotonic() - started,
        "limitations": [
            "Scalar Python implementation is an independent reimplementation, not an execution of RNGLab's NumPy/SciPy code.",
            "Python random.Random streams are deterministic PRNG streams, not physical entropy and not certified independent in a mathematical sense.",
            "Only 64 replicas; coverage estimates have wide intervals.",
            "All methods and windows reuse the same 64 trajectories; comparisons are correlated.",
            "Burn-in and retained length differ between the short and long windows, so their effects are confounded.",
            "No stationarity or mixing proof; pilot cannot establish which interval procedure is valid asymptotically.",
            "No held-out confirmatory cohort or preregistered acceptance threshold was used.",
        ],
    }
    (OUT / "pimc-uncertainty-pilot-results.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"reference": reference, "coverage": rows, "elapsed_seconds": report["elapsed_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
