#!/usr/bin/env python3
"""Analyse the Stage C (H1) campaign (frozen protocol rnglab-v018-h1-detection-v1).

RNGLab v0.18 (written by Claude Code, 2026-11-06).

    PYTHONPATH=src python research/v018-h1-detection/analyze.py

Reads results/cells/*.npz and *.jsonl, results/replay.json and, for the
trace-based secondaries, traces/. Writes results/gates.json (implementation
checks, validity gates, thresholds, cell classes, verdicts, latencies) and
results/secondary.json (descriptive, no gates).
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import beta as beta_dist

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))  # appended, not prepended: other campaigns also have a module named run_campaign
import h1lib  # noqa: E402
from rnglab.seeds import namespace_seed  # noqa: E402
from rnglab12.uncertainty import summarize_series  # noqa: E402

RESAMPLES, CHUNK = 10_000, 500
BONFERRONI_CELLS = 24  # PROTOCOL.md sections 4 and 5.3: frozen multiplicity.
V1_BAND, V2_MAX = (0.92, 0.98), 0.10
HIGH, LOW = 0.5, 0.2
CUSUM_K = 0.5


def cp(k: int, n: int, level: float) -> list[float] | None:
    if n == 0:
        return None
    a = 1 - level
    return [0.0 if k == 0 else float(beta_dist.ppf(a / 2, k, n - k + 1)),
            1.0 if k == n else float(beta_dist.ppf(1 - a / 2, k + 1, n - k))]


def bootstrap_mean(d: np.ndarray, seed: int, level: float, resamples: int = RESAMPLES) -> list[float]:
    rng = np.random.default_rng(seed)
    means = []
    for start in range(0, resamples, CHUNK):
        means.append(d[rng.integers(len(d), size=(min(CHUNK, resamples - start), len(d)))].mean(axis=1))
    a = 1 - level
    return np.quantile(np.concatenate(means), [a / 2, 1 - a / 2]).tolist()


def load_cell(out: Path, name: str) -> dict:
    with np.load(out / "results" / "cells" / f"{name}.npz") as z:
        cell = {k: z[k] for k in z.files}
    cell["records"] = [json.loads(line) for line in (out / "results" / "cells" / f"{name}.jsonl").read_text().splitlines()]
    cell["onset"] = np.array([r["onset_sweep"] for r in cell["records"]])
    return cell


def conv_series(cell: dict) -> np.ndarray:
    p = np.where(np.isnan(cell["conv_p"]), np.inf, cell["conv_p"])
    return -np.log10(np.maximum(p.min(axis=2), h1lib.P_FLOOR))


def app_series(cell: dict, c: dict, parts: str = "A1A2") -> np.ndarray:
    z1 = np.abs(cell["app_stats"][:, :, 0]) / c["sigma_D"]
    z2 = np.abs(cell["app_stats"][:, :, 1] - c["mu_f"]) / c["sigma_f"]
    return {"A1A2": np.maximum(z1, z2), "A1": z1, "A2": z2}[parts]


def detect(series: np.ndarray, ends: np.ndarray, cp_bytes: np.ndarray, onset: np.ndarray, tau: float) -> dict:
    over = series > tau
    post = over & (ends[None, :] > onset[:, None])
    first = np.where(post.any(axis=1), post.argmax(axis=1), -1)
    rows = np.arange(len(series))
    return {"whole_run": over.any(axis=1), "post": post.any(axis=1),
            "latency_sweeps": np.where(first >= 0, ends[np.maximum(first, 0)] - onset, -1),
            "latency_bytes": np.where(first >= 0, cp_bytes[rows, np.maximum(first, 0)] - h1lib.ONSET_BYTES, -1),
            "first_index": first}


def classify(wrong_cp, c_cp, a_cp, gap_ci) -> str:
    if wrong_cp is None:
        return "indeterminate"
    if wrong_cp[0] > HIGH and c_cp[1] < LOW and gap_ci[0] > 0:
        return "H1"
    if wrong_cp[0] > HIGH and c_cp[1] < LOW and a_cp[1] < LOW:
        return "undetected_harm"
    if c_cp[0] > HIGH and wrong_cp[1] < LOW:
        return "conventional_first"
    if wrong_cp[0] > HIGH and c_cp[0] > HIGH and a_cp[0] > HIGH:
        return "both_detect_harmed"
    if wrong_cp[1] < LOW and c_cp[1] < LOW and a_cp[1] < LOW:
        return "harmless_quiet"
    return "indeterminate"


def cell_names(out: Path) -> list[str]:
    return sorted(p.stem for p in (out / "results" / "cells").glob("*.npz"))


def gates(out: Path = HERE, ns: str | None = None, order_statistic: int | None = None, resamples: int = RESAMPLES) -> dict:
    protocol = json.loads((HERE / "protocol.json").read_text())
    ns = ns or protocol["data"]["namespace"]
    names = cell_names(out)
    fault_names = [n for n in names if n not in ("calibration", "healthy")]
    level = 1 - 0.05 / BONFERRONI_CELLS
    cal, healthy = load_cell(out, "calibration"), load_cell(out, "healthy")
    order_statistic = order_statistic or math.ceil((len(cal["records"]) + 1) * 0.95)
    constants = {"sigma_D": float(cal["app_stats"][:, :, 0].std(ddof=1)), "mu_f": float(cal["app_stats"][:, :, 1].mean()),
                 "sigma_f": float(cal["app_stats"][:, :, 1].std(ddof=1))}
    tau = {"C": float(np.sort(conv_series(cal).max(axis=1))[order_statistic - 1]),
           "A": float(np.sort(app_series(cal, constants).max(axis=1))[order_statistic - 1])}

    def arms(cell):
        return {"C": detect(conv_series(cell), cell["conv_end"][0], cell["conv_bytes"], cell["onset"], tau["C"]),
                "A": detect(app_series(cell, constants), cell["app_end"][0], cell["app_bytes"], cell["onset"], tau["A"])}

    # Implementation checks
    every = [load_cell(out, n) for n in names]
    i1 = sum(not r["per_sweep_bytes_ok"] for c in every for r in c["records"])
    hp = [r["prefix_sha256"] for r in healthy["records"]]
    i2 = sum(r["prefix_sha256"] != hp[r["replica"]] for n in fault_names for r in load_cell(out, n)["records"])
    replay = json.loads((out / "results" / "replay.json").read_text())
    implementation = {
        "I-1": {"failing_runs": i1, "result": "PASS" if i1 == 0 else "FAIL"},
        "I-2": {"mismatched_prefixes": i2, "result": "PASS" if i2 == 0 else "FAIL"},
        "I-3": {"cells": {k: v["I-3"] for k, v in replay.items() if "I-3" in v},
                "result": "PASS" if all(replay.get(n, {}).get("I-3", False) for n in fault_names)
                and all(v["I-3"] for v in replay.values() if "I-3" in v) else "FAIL"},
        "I-4": {"cells": {k: v.get("I-4", False) for k, v in replay.items()},
                "result": "PASS" if all(replay.get(n, {}).get("I-4", False) for n in names)
                and all(v.get("I-4", False) for v in replay.values()) else "FAIL"}}

    # Validity gates
    def harm(cell):
        ok = cell["onset"] <= h1lib.HARM_START
        wrong = np.array([r["harm"]["wrong"] for r in cell["records"]])
        return wrong[ok], int((~ok).sum())

    hw, h_flagged = harm(healthy)
    covered = int((~hw).sum())
    v1_cp = cp(covered, len(hw), 0.95)
    ha = arms(healthy)
    fpr = {a: {"alarms": int(ha[a]["whole_run"].sum()), "runs": len(hw) + h_flagged,
               "cp95": cp(int(ha[a]["whole_run"].sum()), len(healthy["records"]), 0.95)} for a in ("C", "A")}
    d = ha["A"]["whole_run"].astype(float) - ha["C"]["whole_run"].astype(float)
    validity = {
        "V1": {"covered": covered, "runs": len(hw), "onset_flagged": h_flagged, "cp95": v1_cp,
               "result": "PASS" if v1_cp is not None and V1_BAND[0] <= v1_cp[0] and v1_cp[1] <= V1_BAND[1] else "FAIL"},
        "V2": {**{f"rate_{a}": fpr[a] for a in fpr}, "paired_A_minus_C": float(d.mean()),
               "paired_ci95": bootstrap_mean(d, namespace_seed(ns, "analysis", "bootstrap", "V2"), 0.95, resamples),
               "result": "PASS" if all(fpr[a]["cp95"][1] <= V2_MAX for a in fpr) else "FAIL"}}

    # Cells
    cells = {}
    for name in fault_names:
        cell = load_cell(out, name)
        w, flagged = harm(cell)
        ok = cell["onset"] <= h1lib.HARM_START
        det = arms(cell)
        n = len(cell["records"])
        kc, ka, kw = int(det["C"]["post"].sum()), int(det["A"]["post"].sum()), int(w.sum())
        gap = det["A"]["post"].astype(float) - det["C"]["post"].astype(float)
        entry = {"family": cell["records"][0]["family"], "severity": cell["records"][0]["severity"], "runs": n,
                 "onset_flagged": flagged,
                 "wrong": {"k": kw, "n": len(w), "rate": kw / len(w) if len(w) else None,
                           "cp_bonferroni": cp(kw, len(w), level), "cp95": cp(kw, len(w), 0.95)},
                 "C_alarm": {"k": kc, "rate": kc / n, "cp_bonferroni": cp(kc, n, level), "cp95": cp(kc, n, 0.95)},
                 "A_alarm": {"k": ka, "rate": ka / n, "cp_bonferroni": cp(ka, n, level), "cp95": cp(ka, n, 0.95)},
                 "A_minus_C": {"mean": float(gap.mean()),
                               "ci_bonferroni": bootstrap_mean(gap, namespace_seed(ns, "analysis", "bootstrap", f"gap/{name}"), level, resamples),
                               "ci95": bootstrap_mean(gap, namespace_seed(ns, "analysis", "bootstrap", f"gap95/{name}"), 0.95, resamples)},
                 "pass_but_broken": int((~det["C"]["post"][ok] & w).sum())}
        entry["class"] = classify(entry["wrong"]["cp_bonferroni"], entry["C_alarm"]["cp_bonferroni"],
                                  entry["A_alarm"]["cp_bonferroni"], entry["A_minus_C"]["ci_bonferroni"])
        both = det["C"]["post"] & det["A"]["post"]
        lat = {}
        for unit in ("sweeps", "bytes"):
            diff = (det["A"][f"latency_{unit}"] - det["C"][f"latency_{unit}"])[both].astype(float)
            lat[unit] = {"both": int(both.sum()), "median_A_minus_C": float(np.median(diff)) if len(diff) else None,
                         "ci95": (np.quantile([np.median(diff[idx]) for idx in np.random.default_rng(
                             namespace_seed(ns, "analysis", "bootstrap", f"latency/{unit}/{name}")).integers(len(diff), size=(2000, len(diff)))],
                             [0.025, 0.975]).tolist() if len(diff) > 1 else None)}
        entry["latency"] = {**lat, "only_C": int((det["C"]["post"] & ~det["A"]["post"]).sum()),
                            "only_A": int((det["A"]["post"] & ~det["C"]["post"]).sum()),
                            "neither": int((~det["C"]["post"] & ~det["A"]["post"]).sum())}
        cells[name] = entry
    families = sorted({c["family"] for c in cells.values()})
    verdicts = {f: "SUPPORTED" if any(c["class"] == "H1" for c in cells.values() if c["family"] == f) else "NOT SUPPORTED"
                for f in families}
    interpretable = (validity["V1"]["result"] == validity["V2"]["result"] == "PASS"
                     and all(v["result"] == "PASS" for v in implementation.values()))
    return {"thresholds": tau, "order_statistic": order_statistic, "application_constants": constants,
            "bonferroni_level": level, "resamples": resamples, "implementation": implementation, "validity": validity,
            "cells": cells, "verdicts": {"H1-C": verdicts, "interpretable": interpretable,
                                         "implementation_checks_pass": all(v["result"] == "PASS" for v in implementation.values())}}


def secondary(out: Path = HERE, g: dict | None = None) -> dict:
    g = g or json.loads((out / "results" / "gates.json").read_text())
    names = cell_names(out)
    cal = load_cell(out, "calibration")
    c = g["application_constants"]
    k = g["order_statistic"]
    e_m = h1lib.exact_finite_bead_energy(h1lib.M, h1lib.BETA)
    eye = np.eye(h1lib.M)
    lap = (np.roll(eye, -1, 0) - eye).T @ (np.roll(eye, -1, 0) - eye)
    sigma = np.linalg.inv(h1lib.BETA / h1lib.M * eye + h1lib.M / h1lib.BETA * lap)
    e_rg2 = float(np.trace((eye - 1 / h1lib.M) @ sigma) / h1lib.M)
    sd_rg2 = float(cal["app_stats"][:, :, 2].std(ddof=1))

    def shape_series(cell):
        return np.abs(cell["app_stats"][:, :, 2] - e_rg2) / sd_rg2

    def trace_blocks(name):
        """(replica, CV, TD) per 64-run block, to bound memory."""
        for path in sorted((out / "traces" / name).glob("block-*.npz")):
            with np.load(path) as z:
                yield z["replica"], z["CV"], z["TD"]

    acc = [0.0, 0.0, 0]
    for _, cv, td in trace_blocks("calibration"):
        d = td - cv
        acc[0] += float(d.sum())
        acc[1] += float((d * d).sum())
        acc[2] += d.size
    sd_d = math.sqrt((acc[1] - acc[0]**2 / acc[2]) / (acc[2] - 1))

    def cusum(cv, td, onset=None, h=None):
        """Two-sided Page CUSUM of standardized per-sweep D: whole-run maxima, or post-onset alarms at h."""
        z = (td - cv) / sd_d
        sp, sm = np.zeros(len(z)), np.zeros(len(z))
        best, alarm = np.zeros(len(z)), np.zeros(len(z), dtype=bool)
        for t in range(z.shape[1]):
            sp = np.maximum(0, sp + z[:, t] - CUSUM_K)
            sm = np.maximum(0, sm - z[:, t] - CUSUM_K)
            stat = np.maximum(sp, sm)
            best = np.maximum(best, stat)
            if h is not None:
                alarm |= (stat > h) & (t >= onset)
        return alarm if h is not None else best

    cal_cusum = np.concatenate([cusum(cv, td) for _, cv, td in trace_blocks("calibration")])
    taus = {"A1": float(np.sort(app_series(cal, c, "A1").max(axis=1))[k - 1]),
            "A2": float(np.sort(app_series(cal, c, "A2").max(axis=1))[k - 1]),
            "shape": float(np.sort(shape_series(cal).max(axis=1))[k - 1]),
            "cusum_D": float(np.sort(cal_cusum)[k - 1])}
    out_cells = {}
    for name in names:
        if name == "calibration":
            continue
        cell = load_cell(out, name)
        ends = cell["app_end"][0]
        res = {key: int(detect(series, ends, cell["app_bytes"], cell["onset"], taus[key])["post"].sum())
               for key, series in (("A1", app_series(cell, c, "A1")), ("A2", app_series(cell, c, "A2")), ("shape", shape_series(cell)))}
        res["cusum_D"] = 0
        td_wrong = 0
        for replicas, cv, td in trace_blocks(name):
            res["cusum_D"] += int(cusum(cv, td, cell["onset"][replicas], taus["cusum_D"]).sum())
            for series in td:
                bm = summarize_series(series[h1lib.HARM_START:h1lib.HARM_END], reference=e_m,
                                      batch_sizes=(h1lib.HARM_BATCH,))["batch_means"][str(h1lib.HARM_BATCH)]
                td_wrong += not (bm["status"] == "ESTIMATED" and bm["covered"])
        conv = detect(conv_series(cell), cell["conv_end"][0], cell["conv_bytes"], cell["onset"], g["thresholds"]["C"])
        firing = {}
        for i in np.flatnonzero(conv["post"]):
            row = cell["conv_p"][i, conv["first_index"][i]]
            name_k = h1lib.CONV_P_KEYS[int(np.nanargmin(row))]
            firing[name_k] = firing.get(name_k, 0) + 1
        ht =[r["health_tests"] for r in cell["records"]]
        out_cells[name] = {"runs": len(cell["records"]), "post_onset_alarms": res, "first_firing_conventional_test": firing,
                           "health_test_post_onset_alarm_runs": sum(h["rct_post"] + h["apt_post"] > 0 for h in ht),
                           "mean_cv_bias": float(np.mean([r["harm"]["mean"] for r in cell["records"]]) - e_m),
                           "td_wrong": td_wrong}
    return {"thresholds": taus, "E_rg2": e_rg2, "cusum_k": CUSUM_K, "cells": out_cells,
            "status": "descriptive, no gates"}


def main() -> None:
    g = gates()
    (HERE / "results" / "gates.json").write_text(json.dumps(g, indent=1) + "\n")
    s = secondary(HERE, g)
    (HERE / "results" / "secondary.json").write_text(json.dumps(s, indent=1) + "\n")
    print(json.dumps({"validity": {k: v["result"] for k, v in g["validity"].items()},
                      "implementation": {k: v["result"] for k, v in g["implementation"].items()},
                      "classes": {k: v["class"] for k, v in g["cells"].items()}, "verdicts": g["verdicts"]}, indent=1))


if __name__ == "__main__":
    main()
