#!/usr/bin/env python3
"""Analyse the v0.18 Phase 2 campaign (frozen protocol rnglab-v018-api-validation-v1).

RNGLab v0.18 (written by Claude Code, 2026-11-05).

    PYTHONPATH=src python research/v018-api-validation/analyze.py

Reads results/<condition>_replicas.csv, results/replay.json, results/run_manifest.json
and the traces. Writes results/gates.json (implementation checks I-1..I-4, the
required gates G-R and G-I, the replication gates F-L and P-RL, and the
verdicts) and results/secondary.json (descriptive endpoints, no gates).
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_dist, t as student_t

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
for _path in (REPO / "src", HERE):
    if str(_path) not in sys.path:
        sys.path.insert(0, str(_path))

from rnglab import longrun  # noqa: E402
from rnglab.seeds import namespace_seed  # noqa: E402
from rnglab14.models.pimc_harmonic import exact_finite_bead_energy  # noqa: E402

M, BETA, N, BATCH = 16, 2.0, 8192, 512
CONDITIONS = ("local", "refresh", "iid")
RESAMPLES, CHUNK = 10_000, 500
COVERAGE_GATE, RATIO_GATE, BIAS_GATE, MARGIN, FL_UPPER = (0.92, 0.98), (0.75, 1.25), (-0.01, 0.01), 0.10, 0.92
I3_TOLERANCE = 1e-12
SECONDARY_METHODS = (("BM", 32), ("BM", 128), ("BM", 512), ("OBM", 512), ("SV-B", 512), ("SV-TH", 512), ("GEYER", None))
PHASE2 = REPO / "research" / "phase2-20261030" / "experiment" / "confirmation-v1" / "primary_decision.json"
PHASE2_KEYS = {"local": "local_stationary", "refresh": "centroid_refresh", "iid": "iid_gaussian"}


def load_rows(results: Path, condition: str) -> list[dict]:
    with (results / f"{condition}_replicas.csv").open() as handle:
        return list(csv.DictReader(handle))


def column(rows, name, kind=float):
    return np.array([kind(r[name]) if r[name] != "" else np.nan for r in rows])


def flag(rows, name):
    return np.array([r[name] == "True" for r in rows])


def cp95(k: int, n: int) -> list[float]:
    return [0.0 if k == 0 else float(beta_dist.ppf(0.025, k, n - k + 1)),
            1.0 if k == n else float(beta_dist.ppf(0.975, k + 1, n - k))]


def t_interval(values) -> dict:
    v = np.asarray(values, dtype=float)
    mean, se = float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v)))
    q = float(student_t.ppf(0.975, len(v) - 1))
    return {"mean": mean, "ci95": [mean - q * se, mean + q * se]}


def bootstrap(statistic, arrays, seed: int, resamples: int = RESAMPLES) -> list[float]:
    """Percentile 95% interval over replicas; draws in chunks of CHUNK rows."""
    rng = np.random.default_rng(seed)
    r = len(arrays[0])
    values = []
    for start in range(0, resamples, CHUNK):
        idx = rng.integers(r, size=(min(CHUNK, resamples - start), r))
        values.append(statistic(*(a[idx] for a in arrays)))
    return np.quantile(np.concatenate(values), [0.025, 0.975]).tolist()


def ratio_statistic(mcse2, means):
    return np.nanmean(mcse2, axis=-1) / means.var(axis=-1, ddof=1)


def inside(interval, bounds) -> bool:
    return bool(bounds[0] <= interval[0] and interval[1] <= bounds[1])


def primary(rows, ns: str, reference: float, prefix: str = "cv", resamples: int = RESAMPLES) -> dict:
    estimated = np.array([r[f"{prefix}_status"] == "ESTIMATED" for r in rows])
    covered = estimated & flag(rows, f"{prefix}_covered")
    means = column(rows, f"{prefix}_mean")
    mcse2 = np.where(estimated, column(rows, f"{prefix}_mcse") ** 2, np.nan)
    k, r = int(covered.sum()), len(rows)
    condition = rows[0].get("condition", "")
    return {"covered": k, "units": r, "coverage": k / r, "cp95": cp95(k, r), "not_estimated": int((~estimated).sum()),
            "variance_ratio": {"ratio": float(ratio_statistic(mcse2, means)),
                               "ci95": bootstrap(ratio_statistic, [mcse2, means],
                                                 namespace_seed(ns, "analysis", "bootstrap", f"{prefix}/{condition}"), resamples)},
            "bias": t_interval(means - reference), "covered_flags": covered}


def required_gate(p: dict) -> dict:
    parts = {"coverage": inside(p["cp95"], COVERAGE_GATE), "variance_ratio": inside(p["variance_ratio"]["ci95"], RATIO_GATE),
             "bias": inside(p["bias"]["ci95"], BIAS_GATE)}
    return {**parts, "result": "PASS" if all(parts.values()) else "FAIL"}


def implementation_checks(rows: dict, n: int, ns: str, out: Path) -> dict:
    from run_campaign import start_path, array_sha256

    expect = {"proposal": 8 * M * n, "initialization": 0, "centroid": 16 * n, "iid_paths": 16 * ((M * n + 1) // 2), "start": 256}
    bad = []
    for condition in ("local", "refresh"):
        for r in rows[condition]:
            want = {"bytes_proposal": str(expect["proposal"]), "bytes_initialization": "0", "bytes_resampling": "",
                    "bytes_centroid": str(expect["centroid"]) if condition == "refresh" else "", "aliases": "",
                    "start_bytes_before_run": str(expect["start"])}
            bad += [(condition, r["replica"], k) for k, v in want.items() if r[k] != v]
    bad += [("iid", r["replica"], "bytes_paths") for r in rows["iid"] if r["bytes_paths"] != str(expect["iid_paths"])]
    i2 = [r["replica"] for r, s in zip(rows["local"], rows["refresh"])
          if not (r["start_sha256"] == s["start_sha256"] == array_sha256(start_path(ns, int(r["replica"]))[0]))]
    i3 = max(float(r["cv_rederive_max_abs"]) for c in ("local", "refresh") for r in rows[c])
    replay = json.loads((out / "results" / "replay.json").read_text())
    checks = {"I-1": {"mismatches": len(bad), "first": bad[:5], "expected_bytes": expect},
              "I-2": {"mismatched_replicas": len(i2), "first": i2[:5]},
              "I-3": {"max_abs_difference": i3, "tolerance": I3_TOLERANCE},
              "I-4": replay}
    checks["I-1"]["result"] = "PASS" if not bad else "FAIL"
    checks["I-2"]["result"] = "PASS" if not i2 else "FAIL"
    checks["I-3"]["result"] = "PASS" if i3 <= I3_TOLERANCE else "FAIL"
    checks["I-4"]["result"] = "PASS" if all(v["identical"] for v in replay.values()) else "FAIL"
    return checks


def gates(out: Path = HERE, ns: str | None = None, n: int = N, resamples: int = RESAMPLES) -> dict:
    protocol = json.loads((HERE / "protocol.json").read_text())
    ns = ns or protocol["data"]["namespace"]
    reference = exact_finite_bead_energy(M, BETA)
    rows = {c: load_rows(out / "results", c) for c in CONDITIONS}
    for c in CONDITIONS:
        for r in rows[c]:
            r["condition"] = c
    p = {c: primary(rows[c], ns, reference, resamples=resamples) for c in CONDITIONS}
    if [r["replica"] for r in rows["local"]] != [r["replica"] for r in rows["refresh"]]:
        raise SystemExit("local and refresh replicas are not aligned")
    d = p["refresh"]["covered_flags"].astype(float) - p["local"]["covered_flags"].astype(float)
    paired_ci = bootstrap(lambda x: x.mean(axis=-1), [d], namespace_seed(ns, "analysis", "bootstrap", "paired"), resamples)
    result = {
        "implementation": implementation_checks(rows, n, ns, out),
        "G-R": {**required_gate(p["refresh"]), "condition": "refresh"},
        "G-I": {**required_gate(p["iid"]), "condition": "iid"},
        "F-L": {"cp95_upper": p["local"]["cp95"][1], "criterion": f"< {FL_UPPER}",
                "result": "PASS" if p["local"]["cp95"][1] < FL_UPPER else "FAIL"},
        "P-RL": {"mean_difference": float(d.mean()), "bootstrap_ci95": paired_ci, "improved_only": int((d == 1).sum()),
                 "local_only": int((d == -1).sum()), "criterion": f"lower > {MARGIN}",
                 "result": "PASS" if paired_ci[0] > MARGIN else "FAIL"},
        "primary": {c: {k: v for k, v in p[c].items() if k != "covered_flags"} for c in CONDITIONS},
        "reference_E_M": reference, "resamples": resamples,
    }
    implementation_ok = all(v["result"] == "PASS" for v in result["implementation"].values())
    h_api = result["G-R"]["result"] == result["G-I"]["result"] == "PASS"
    h_rep = result["F-L"]["result"] == result["P-RL"]["result"] == "PASS"
    result["verdicts"] = {"implementation_checks_pass": implementation_ok,
                          "H-API": "SUPPORTED" if h_api else "FAILED", "H-REP": "SUPPORTED" if h_rep else "FAILED",
                          "core_result_reproduced": h_api and h_rep,
                          "interpretable": implementation_ok}
    return result


def trace_blocks(out: Path, condition: str, run_manifest: dict | None):
    hashes = {(b["condition"], b["block"]): b["npz_sha256"] for b in run_manifest["blocks"]} if run_manifest else {}
    block = 0
    while (path := out / "traces" / condition / f"block-{block:04d}.npz").exists():
        if hashes and hashlib.sha256(path.read_bytes()).hexdigest() != hashes[(condition, block)]:
            raise SystemExit(f"{path} does not match run_manifest.json")
        with np.load(path) as z:
            yield {k: z[k] for k in ("replica", "CV", "centroid", "rg2")}
        block += 1


def secondary(out: Path = HERE, primary_result: dict | None = None, n: int = N) -> dict:
    reference = exact_finite_bead_energy(M, BETA)
    manifest_path = out / "results" / "run_manifest.json"
    run_manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else None
    phase2 = json.loads(PHASE2.read_text())["gates"]
    result = {}
    for condition in CONDITIONS:
        rows = load_rows(out / "results", condition)
        stats = {f"{m}{b or ''}": {"covered": [], "sigma2_over_n": [], "mean": []} for m, b in SECONDARY_METHODS}
        tau, adaptive_cov, statuses, cmean, smean, means = [], [], [], [], [], []
        for z in trace_blocks(out, condition, run_manifest):
            x = z["CV"]
            for method, b in SECONDARY_METHODS:
                e = longrun.estimate(x, method, b)
                lo, hi = longrun.interval(e["mean"], e["sigma2"], e["n_used"], e["df"])
                s = stats[f"{method}{b or ''}"]
                s["covered"].append(longrun.covers(lo, hi, reference))
                s["sigma2_over_n"].append(e["sigma2"] / e["n_used"])
                s["mean"].append(e["mean"])
                if method == "GEYER":
                    tau.append(e["tau"])
            t_hat = tau[-1]
            statuses.append(longrun.adaptive_status(n, t_hat))
            chosen = longrun.adaptive_batch(n, t_hat)
            for b in sorted(set(chosen[chosen > 0].tolist())):
                rows_b = chosen == b
                e = longrun.estimate(x[rows_b], "BM", int(b))
                lo, hi = longrun.interval(e["mean"], e["sigma2"], e["n_used"], e["df"])
                adaptive_cov.append(longrun.covers(lo, hi, reference))
            cmean.append((0.5 * z["centroid"] ** 2).mean(axis=1))
            smean.append(z["rg2"].mean(axis=1))
            means.append(x.mean(axis=1))
        all_means = np.concatenate(means)
        empirical = float(all_means.var(ddof=1))
        methods = {}
        for name, s in stats.items():
            cov = np.concatenate(s["covered"])
            v = np.concatenate(s["sigma2_over_n"])
            methods[name] = {"covered": int(cov.sum()), "units": len(cov), "cp95": cp95(int(cov.sum()), len(cov)),
                             "variance_ratio": float(np.nanmean(v) / empirical), "undefined": int(np.isnan(v).sum())}
        status = np.concatenate(statuses)
        adaptive = np.concatenate(adaptive_cov) if adaptive_cov else np.zeros(0, dtype=bool)
        t_all = np.concatenate(tau)
        c, s = np.concatenate(cmean), np.concatenate(smean)
        cov_cs = float(np.cov(c, s, ddof=1)[0, 1])
        td_cov = np.array([r["td_status"] == "ESTIMATED" and r["td_covered"] == "True" for r in rows])
        td_means, td_mcse2 = column(rows, "td_mean"), column(rows, "td_mcse") ** 2
        k2 = phase2[PHASE2_KEYS[condition]]["row"]
        p_api, p_ph2 = methods["BM512"]["covered"] / len(rows), k2["covered"] / k2["replicas"]
        se = math.sqrt(p_api * (1 - p_api) / len(rows) + p_ph2 * (1 - p_ph2) / k2["replicas"])
        result[condition] = {
            "cv_methods": methods,
            "cv_BM512_longrun_matches_api_flags": int((np.concatenate(stats["BM512"]["covered"])
                                                       == flag(rows, "cv_covered")).sum()),
            "td_BM512": {"covered": int(td_cov.sum()), "units": len(rows), "cp95": cp95(int(td_cov.sum()), len(rows)),
                         "variance_ratio": float(np.nanmean(td_mcse2) / td_means.var(ddof=1)),
                         "bias": t_interval(td_means - reference)},
            "diagnostic": {"status_counts": {k: int((status == k).sum()) for k in ("ADEQUATE", "RUN_TOO_SHORT", "TAU_UNDEFINED")},
                           "adaptive_covered": int(adaptive.sum()), "adaptive_units": len(adaptive),
                           "adaptive_cp95": cp95(int(adaptive.sum()), len(adaptive)) if len(adaptive) else None},
            "geyer_tau_quantiles": dict(zip(("q05", "q50", "q95"), np.nanquantile(t_all, [0.05, 0.5, 0.95]).tolist())),
            "ess_quantiles": dict(zip(("q05", "q50", "q95"), np.nanquantile(n / t_all, [0.05, 0.5, 0.95]).tolist())),
            "variance_decomposition": {"centroid_part": float(c.var(ddof=1)), "shape_part": float(s.var(ddof=1)),
                                       "twice_covariance": 2 * cov_cs, "sum": float(c.var(ddof=1) + s.var(ddof=1) + 2 * cov_cs),
                                       "cv_means_variance": empirical},
            "phase2_comparison": {"api_coverage": p_api, "phase2_coverage": p_ph2, "difference": p_api - p_ph2,
                                  "wald_ci95": [p_api - p_ph2 - 1.959963984540054 * se, p_api - p_ph2 + 1.959963984540054 * se]},
        }
        if condition != "iid":
            acc = column(rows, "acceptance_rate_4dp")
            result[condition]["acceptance_rate_4dp_quantiles"] = dict(zip(("q05", "q50", "q95"), np.quantile(acc, [0.05, 0.5, 0.95]).tolist()))
    return result


def main() -> None:
    g = gates()
    (HERE / "results" / "gates.json").write_text(json.dumps(g, indent=1) + "\n")
    s = secondary()
    (HERE / "results" / "secondary.json").write_text(json.dumps(s, indent=1) + "\n")
    print(json.dumps({k: g[k]["result"] for k in ("G-R", "G-I", "F-L", "P-RL")} | {"verdicts": g["verdicts"]}, indent=1))


if __name__ == "__main__":
    main()
