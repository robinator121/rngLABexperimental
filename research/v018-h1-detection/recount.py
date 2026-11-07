#!/usr/bin/env python3
"""Independent recount of the Stage C (H1) endpoints (PROTOCOL.md section 5.4).

RNGLab v0.18 (written by Claude Code, 2026-11-06). Shares no code with
analyze.py, h1lib.py or rnglab12: E_M comes from its closed form; the harm
interval is a NumPy reshape of the CV traces (16 batches of 480, t15); the arm
statistics, thresholds, detections, rates, Clopper-Pearson bounds and cell
classes are recomputed here. The paired bootstrap uses the protocol's seeds and
chunked resampling, so its interval must equal analyze.py's.

    PYTHONPATH=src python research/v018-h1-detection/recount.py

Reads results/cells/, traces/ and results/gates.json; writes results/recount.json.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
M, BETA, ONSET_BYTES = 16, 2.0, 2**19
HARM = (2560, 10240, 480)
P_FLOOR = 1e-300


def reference_energy() -> float:
    z = BETA / (2 * M)
    return 0.5 / (math.tanh(M * math.asinh(z)) * math.sqrt(1 + z * z))


def seed(ns: str, *parts) -> int:
    return int.from_bytes(hashlib.sha256("/".join([ns, *parts]).encode()).digest()[:8], "big")


def clopper_pearson(k: int, n: int, level: float) -> tuple[float, float]:
    a = 1 - level
    lo = stats.beta.ppf(a / 2, k, n - k + 1) if k > 0 else 0.0
    hi = stats.beta.ppf(1 - a / 2, k + 1, n - k) if k < n else 1.0
    return float(lo), float(hi)


def paired_interval(d: np.ndarray, s: int, level: float, resamples: int, chunk: int = 500) -> list[float]:
    rng = np.random.default_rng(s)
    draws = [d[rng.integers(len(d), size=(min(chunk, resamples - i), len(d)))].mean(axis=1) for i in range(0, resamples, chunk)]
    return np.quantile(np.concatenate(draws), [(1 - level) / 2, 1 - (1 - level) / 2]).tolist()


def load(out: Path, name: str) -> dict:
    with np.load(out / "results" / "cells" / f"{name}.npz") as z:
        cell = {k: z[k] for k in z.files}
    records = [json.loads(line) for line in (out / "results" / "cells" / f"{name}.jsonl").read_text().splitlines()]
    cell["onset"] = np.array([r["onset_sweep"] for r in records])
    cell["records"] = records
    return cell


def wrong_from_traces(out: Path, name: str, reference: float) -> np.ndarray:
    start, end, b = HARM
    flags = []
    for path in sorted((out / "traces" / name).glob("block-*.npz")):
        with np.load(path) as z:
            x = z["CV"][:, start:end]
        means = x.reshape(len(x), (end - start) // b, b).mean(axis=2)
        k = means.shape[1]
        centre = means.mean(axis=1)
        half = stats.t.ppf(0.975, k - 1) * np.sqrt(means.var(axis=1, ddof=1) / k)
        flags.append(~((centre - half <= reference) & (reference <= centre + half)))
    return np.concatenate(flags)


def recount(out: Path = HERE, ns: str | None = None, resamples: int = 10_000) -> dict:
    ns = ns or json.loads((HERE / "protocol.json").read_text())["data"]["namespace"]
    g = json.loads((out / "results" / "gates.json").read_text())
    reference = reference_energy()
    comparisons = []

    def compare(name, mine, theirs, tol=0.0):
        if tol:
            ok = all(abs(a - b) <= tol * max(1.0, abs(b)) for a, b in zip(np.atleast_1d(mine), np.atleast_1d(theirs)))
        else:
            ok = mine == theirs
        comparisons.append({"quantity": name, "recount": mine, "analysis": theirs, "agree": bool(ok)})

    names = sorted(p.stem for p in (out / "results" / "cells").glob("*.npz"))
    faults = [n for n in names if n not in ("calibration", "healthy")]
    level = 1 - 0.05 / len(faults)
    cal = load(out, "calibration")
    order = g["order_statistic"]
    stats_d, stats_f = cal["app_stats"][:, :, 0], cal["app_stats"][:, :, 1]
    sd, mu, sf = stats_d.std(ddof=1), stats_f.mean(), stats_f.std(ddof=1)

    def s_conv(cell):
        p = cell["conv_p"].copy()
        p[np.isnan(p)] = 1.0
        return -np.log10(np.clip(p.min(axis=2), P_FLOOR, None))

    def s_app(cell):
        return np.maximum(np.abs(cell["app_stats"][:, :, 0]) / sd, np.abs(cell["app_stats"][:, :, 1] - mu) / sf)

    tau = {"C": float(np.sort(s_conv(cal).max(axis=1))[order - 1]), "A": float(np.sort(s_app(cal).max(axis=1))[order - 1])}
    compare("threshold C", tau["C"], g["thresholds"]["C"], 1e-12)
    compare("threshold A", tau["A"], g["thresholds"]["A"], 1e-12)

    def post_alarm(series, ends, onset, t):
        return ((series > t) & (ends[None, :] > onset[:, None])).any(axis=1)

    for name in ["healthy", *faults]:
        cell = load(out, name)
        ok = cell["onset"] <= HARM[0]
        wrong = wrong_from_traces(out, name, reference)
        stored = np.array([r["harm"]["wrong"] for r in cell["records"]])
        compare(f"{name}: wrong flags equal stored", int((wrong == stored).sum()), len(stored))
        c_post = post_alarm(s_conv(cell), cell["conv_end"][0], cell["onset"], tau["C"])
        a_post = post_alarm(s_app(cell), cell["app_end"][0], cell["onset"], tau["A"])
        if name == "healthy":
            v = g["validity"]
            compare("V1 covered", int((~wrong[ok]).sum()), v["V1"]["covered"])
            compare("V1 CP95", list(clopper_pearson(int((~wrong[ok]).sum()), int(ok.sum()), 0.95)), v["V1"]["cp95"], 1e-12)
            for arm, s in (("C", s_conv(cell)), ("A", s_app(cell))):
                compare(f"V2 whole-run alarms {arm}", int((s > tau[arm]).any(axis=1).sum()), v["V2"][f"rate_{arm}"]["alarms"])
            continue
        e = g["cells"][name]
        kw, kc, ka = int(wrong[ok].sum()), int(c_post.sum()), int(a_post.sum())
        compare(f"{name}: wrong", kw, e["wrong"]["k"])
        compare(f"{name}: C alarms", kc, e["C_alarm"]["k"])
        compare(f"{name}: A alarms", ka, e["A_alarm"]["k"])
        compare(f"{name}: pass but broken", int((~c_post[ok] & wrong[ok]).sum()), e["pass_but_broken"])
        n, nw = len(cell["records"]), int(ok.sum())
        w_cp, c_cp, a_cp = clopper_pearson(kw, nw, level), clopper_pearson(kc, n, level), clopper_pearson(ka, n, level)
        gap = paired_interval(a_post.astype(float) - c_post.astype(float), seed(ns, "analysis", "bootstrap", f"gap/{name}"), level, resamples)
        compare(f"{name}: A-C interval", gap, e["A_minus_C"]["ci_bonferroni"], 1e-12)
        high, low = 0.5, 0.2
        if w_cp[0] > high and c_cp[1] < low and gap[0] > 0:
            cls = "H1"
        elif w_cp[0] > high and c_cp[1] < low and a_cp[1] < low:
            cls = "undetected_harm"
        elif c_cp[0] > high and w_cp[1] < low:
            cls = "conventional_first"
        elif w_cp[0] > high and c_cp[0] > high and a_cp[0] > high:
            cls = "both_detect_harmed"
        elif w_cp[1] < low and c_cp[1] < low and a_cp[1] < low:
            cls = "harmless_quiet"
        else:
            cls = "indeterminate"
        compare(f"{name}: class", cls, e["class"])
    result = {"reference_E_M": reference, "comparisons": comparisons, "all_agree": all(c["agree"] for c in comparisons)}
    (out / "results" / "recount.json").write_text(json.dumps(result, indent=1) + "\n")
    return result


if __name__ == "__main__":
    r = recount()
    print(f"{sum(c['agree'] for c in r['comparisons'])}/{len(r['comparisons'])} agree; all_agree = {r['all_agree']}")
