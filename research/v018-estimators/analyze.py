#!/usr/bin/env python3
"""Evaluate the frozen gates of protocol rnglab-v018-estimators-v1.

RNGLab v0.18 (written by Claude Code, 2026-11-05; hashed in
manifest_pregeneration.json before any campaign draw).

    PYTHONPATH=src python research/v018-estimators/analyze.py

Reads results/ar1/*.npz and results/ising/*.npz; writes results/gates.json,
results/ar1_table.csv and results/ising_table.csv. Deterministic: bootstrap
seeds come from the protocol namespace.

Amendment 01 (AMENDMENT_01.md, before any campaign output was analyzed):
verdicts require every registered cell, setting, unit count and gate ID;
R-hat comes from results/amend01 (recomputed from traces), with the original
values reported beside it; the control gates get a secondary
multiplicity-adjusted view and the Ising cells a split-half burn-in
diagnostic. No frozen gate definition or criterion changes.
"""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import beta, chi2, t as student_t

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

from rnglab import longrun  # noqa: E402
from rnglab.seeds import namespace_seed  # noqa: E402

PROTOCOL = json.loads((HERE / "protocol.json").read_text())
NS = PROTOCOL["seed"]["namespace"]
GATES = PROTOCOL["gates"]["validation"]
COVERAGE_BOUNDS = tuple(GATES["S-S"]["coverage_cp95_within"])
BIAS_BOUNDS = tuple(GATES["S-S"]["relative_bias_ci_within"])
RATIO_BOUNDS = tuple(GATES["S-S"]["variance_ratio_ci_within"])
BOOTSTRAP = 2000
RHAT_THRESHOLD = 1.01
ISING = PROTOCOL["data"]["ising"]
AR1 = PROTOCOL["data"]["ar1"]
ISING_CELLS = [f"ising_L{L}_{label}" for L in ISING["L"] for label in ISING["T"]]
AR1_SETTINGS = [f"ar1_phi{float(phi)!r}_n{n}" for phi in AR1["phi"] for n in AR1["n"]]


AMENDMENT: dict = {"L-C_adjusted": {}, "rhat_original_vs_amended": {}, "split_half_burn_in": {}}


def expected_gate_ids() -> dict:
    """Every gate ID the frozen protocol registers, grouped by verdict prefix (Amendment 01)."""
    return {
        "L-C[": [f"L-C[{s}]" for s in AR1_SETTINGS],
        "S-C[": [f"S-C[{c}]" for c in ISING_CELLS],
        "S-S[": [f"S-S[{c}]" for c in ISING_CELLS],
        "S-G[": [f"S-G[{c}]" for c in ISING_CELLS],
        "L-A[": [f"L-A[{float(p)},{m}]" for p in GATES["L-A"]["phi"] for m in GATES["L-A"]["methods"]],
        "L-G[": [f"L-G[{float(p)}]" for p in GATES["L-G"]["phi"]],
        "L-D[": [f"L-D[{float(p)},{n},{f}]" for p in AR1["phi"] for n in AR1["n"] for f in GATES["L-D"]["families"]],
        "single": ["S-F1", "D-F1", "L-F1", "L-F2"],
    }


def completeness() -> dict:
    """Check that every registered input exists with its registered unit count (Amendment 01)."""
    problems = []
    chains = ISING["chains_per_cell"]
    for cell in ISING_CELLS:
        for folder, checks in (("ising", {"S": chains, "G": chains // 16}), ("amend01", {"rhat_m": chains // 16})):
            path = HERE / "results" / folder / f"{cell}.npz"
            if not path.exists():
                problems.append(f"missing results/{folder}/{cell}.npz")
                continue
            z = np.load(path)
            problems += [f"{folder}/{cell}: {key} has {len(z[key])} units, expected {count}"
                         for key, count in checks.items() if len(z[key]) != count]
    for setting in AR1_SETTINGS:
        path = HERE / "results" / "ar1" / f"{setting}.npz"
        if not path.exists():
            problems.append(f"missing results/ar1/{setting}.npz")
        elif len(np.load(path)["xbar"]) != AR1["replicas"]:
            problems.append(f"ar1/{setting}: wrong replica count")
    return {"complete": not problems, "problems": problems}


def cp(k: int, n: int) -> list[float]:
    lo = 0.0 if k == 0 else float(beta.ppf(0.025, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(0.975, k + 1, n - k))
    return [lo, hi]


def t_interval(values) -> dict:
    v = np.asarray(values, dtype=float)
    mean = float(v.mean())
    half = float(student_t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v)))
    return {"mean": mean, "ci95": [mean - half, mean + half], "n": int(len(v))}


def bootstrap(statistic, arrays: tuple, endpoint: str) -> list[float]:
    rng = np.random.default_rng(namespace_seed(NS, "analysis", "bootstrap", endpoint))
    n = len(arrays[0])
    draws = np.empty(BOOTSTRAP)
    for i in range(BOOTSTRAP):
        idx = rng.integers(0, n, n)
        draws[i] = statistic(*(a[idx] for a in arrays))
    return [float(q) for q in np.quantile(draws, [0.025, 0.975])]


def inside(interval, bounds) -> bool:
    return bounds[0] <= interval[0] and interval[1] <= bounds[1]


def coverage(covered) -> dict:
    covered = np.asarray(covered, dtype=bool)
    k, n = int(covered.sum()), int(len(covered))
    return {"covered": k, "units": n, "coverage": k / n if n else None, "cp95": cp(k, n) if n else None}


# ---------------------------------------------------------------- AR(1)

def ar1_rows():
    for path in sorted((HERE / "results" / "ar1").glob("*.npz")):
        z = np.load(path)
        yield path.stem, {key: z[key] for key in z.files}


def ar1_interval(z, sigma2, df, rows=None):
    xbar = z["xbar"] if rows is None else z["xbar"][rows]
    lo, hi = longrun.interval(xbar, sigma2, int(z["n"]), df)
    return longrun.covers(lo, hi, 0.0)


def ar1_ratio_vs_exact(z, sigma2, endpoint):
    exact = float(z["n"] * z["exact_var_mean"])
    finite = np.where(np.isfinite(sigma2), sigma2, np.nan)
    ratio = float(np.nanmean(finite) / exact)
    ci = bootstrap(lambda s: np.nanmean(s) / exact, (finite,), endpoint)
    return {"ratio": ratio, "ci95": ci}


def analyze_ar1(gates: dict, table: list) -> None:
    for name, z in ar1_rows():
        n, phi = int(z["n"]), float(z["phi"])
        exact = float(z["exact_var_mean"])
        tau_n = n * exact
        if name in ("ar1_phi0.99_n8192", "ar1_phi0.997_n8192"):
            gid = "L-F1" if "0.99_" in name else "L-F2"
            cov = coverage(ar1_interval(z, z["sigma2_BM_b32"], z["df_BM_b32"]))
            gates[gid] = {"data": name, "BM(32)": cov, "criterion": "CP95 upper < 0.92",
                          "result": "REPRODUCED" if cov["cp95"][1] < 0.92 else "NOT_REPRODUCED"}
        var_ratio = float(np.var(z["xbar"], ddof=1) / exact)
        var_ci = bootstrap(lambda x: np.var(x, ddof=1) / exact, (z["xbar"],), f"{name}/L-C/var")
        mean = t_interval(z["xbar"])
        gates[f"L-C[{name}]"] = {"empirical_over_exact_var": var_ratio, "ci95": var_ci, "mean_xbar": mean,
                                 "result": "PASS" if var_ci[0] <= 1 <= var_ci[1] and mean["ci95"][0] <= 0 <= mean["ci95"][1] else "FAIL"}
        # Amendment 01 secondary view: the frozen L-C gate is 20 separate 95% intervals, so a perfect
        # generator misses at least one with probability about 1 - 0.95**20. Bonferroni over those 20,
        # using exact chi-square and t intervals (the AR(1) replica mean is exactly Gaussian).
        level = 0.05 / (2 * len(AR1_SETTINGS))
        r = len(z["xbar"])
        s2 = float(np.var(z["xbar"], ddof=1))
        adj_var = [(r - 1) * s2 / chi2.ppf(1 - level / 2, r - 1) / exact, (r - 1) * s2 / chi2.ppf(level / 2, r - 1) / exact]
        half = float(student_t.ppf(1 - level / 2, r - 1) * np.sqrt(s2 / r))
        adj_mean = [float(z["xbar"].mean()) - half, float(z["xbar"].mean()) + half]
        AMENDMENT["L-C_adjusted"][name] = {"familywise_level": 0.05, "per_interval_level": level,
                                           "var_ratio_interval": adj_var, "mean_interval": adj_mean,
                                           "flagged": not (adj_var[0] <= 1 <= adj_var[1] and adj_mean[0] <= 0 <= adj_mean[1])}
        for method in ("BM", "OBM", "SV-B", "SV-TH"):
            for b in PROTOCOL["definitions"]["batch_sizes"]:
                sigma2, df = z[f"sigma2_{method}_b{b}"], z[f"df_{method}_b{b}"]
                cov = coverage(ar1_interval(z, sigma2, df))
                ratio = ar1_ratio_vs_exact(z, sigma2, f"{name}/{method}/b{b}")
                status = longrun.batch_length_status(b, n, z["tau_GEYER"])
                ess_ratio = float(np.nanmedian(n * z["gamma0"] / sigma2) * exact)
                table.append({"setting": name, "phi": phi, "n": n, "tau_n": tau_n, "method": method, "b": b,
                              **_flat(cov), "var_ratio": ratio["ratio"], "var_ratio_lo": ratio["ci95"][0],
                              "var_ratio_hi": ratio["ci95"][1], "median_ess_over_true": ess_ratio,
                              "adequate_fraction": float(np.mean(status == "ADEQUATE"))})
                if n == 8192 and b == 512 and phi in GATES["L-A"]["phi"]:
                    ok = inside(cov["cp95"], COVERAGE_BOUNDS) and inside(ratio["ci95"], RATIO_BOUNDS)
                    gates[f"L-A[{phi},{method}]"] = {"coverage": cov, "variance_ratio_vs_exact": ratio,
                                                     "result": "PASS" if ok else "FAIL"}
        sigma2 = np.where(z["valid_GEYER"], z["sigma2_GEYER"], np.nan)
        cov = coverage(ar1_interval(z, sigma2, np.full(len(sigma2), np.inf)))
        ratio = ar1_ratio_vs_exact(z, sigma2, f"{name}/GEYER")
        tau_q = np.nanquantile(z["tau_GEYER"] / tau_n, [0.1, 0.5, 0.9])
        table.append({"setting": name, "phi": phi, "n": n, "tau_n": tau_n, "method": "GEYER", "b": "",
                      **_flat(cov), "var_ratio": ratio["ratio"], "var_ratio_lo": ratio["ci95"][0],
                      "var_ratio_hi": ratio["ci95"][1],
                      "median_ess_over_true": float(np.nanmedian(n * z["gamma0"] / sigma2) * exact),
                      "adequate_fraction": "", "tau_hat_over_tau_q10": tau_q[0],
                      "tau_hat_over_tau_q50": tau_q[1], "tau_hat_over_tau_q90": tau_q[2],
                      "geyer_invalid": int((~z["valid_GEYER"]).sum())})
        if n == 8192 and phi in GATES["L-G"]["phi"]:
            ok = inside(cov["cp95"], COVERAGE_BOUNDS) and inside(ratio["ci95"], RATIO_BOUNDS)
            gates[f"L-G[{phi}]"] = {"coverage": cov, "variance_ratio_vs_exact": ratio, "result": "PASS" if ok else "FAIL"}
        status = longrun.adaptive_status(n, z["tau_GEYER"])
        for family in GATES["L-D"]["families"]:
            rows = z["bstar"] > 0
            k = int(rows.sum())
            cov = coverage(ar1_interval(z, z[f"sigma2_adaptive_{family}"][rows], z[f"df_adaptive_{family}"][rows], rows)) if k else coverage([])
            if k >= 2048:
                result = "PASS" if inside(cov["cp95"], COVERAGE_BOUNDS) else "FAIL"
            elif k >= 200:
                result = "FAIL" if cov["cp95"][1] < COVERAGE_BOUNDS[0] else "INCONCLUSIVE"
            else:
                result = "NOT_EVALUABLE"
            gates[f"L-D[{phi},{n},{family}]"] = {
                "non_abstaining": k, "coverage": cov, "result": result,
                "status_counts": {s: int(np.sum(status == s)) for s in ("ADEQUATE", "RUN_TOO_SHORT", "TAU_UNDEFINED")},
                "bstar_counts": {str(int(b)): int(np.sum(z["bstar"] == b)) for b in np.unique(z["bstar"])}}
            table.append({"setting": name, "phi": phi, "n": n, "tau_n": tau_n, "method": f"adaptive-{family}",
                          "b": "b*", **_flat(cov), "non_abstaining": k})


def _flat(cov: dict) -> dict:
    return {"covered": cov["covered"], "units": cov["units"], "coverage": cov["coverage"],
            "cp95_lo": cov["cp95"][0] if cov["cp95"] else None, "cp95_hi": cov["cp95"][1] if cov["cp95"] else None}


# ---------------------------------------------------------------- Ising

def analyze_ising(gates: dict, table: list) -> None:
    for path in sorted((HERE / "results" / "ising").glob("*.npz")):
        cell = path.stem
        z = np.load(path)
        chi, n_ret = float(z["chi_exact"]), int(z["n_retained"])

        def s_cov(sigma2, df):
            lo, hi = longrun.interval(z["S"], sigma2, n_ret, df)
            return coverage(longrun.covers(lo, hi, chi))

        cov_s = s_cov(z["S_sigma2_BM16"], z["S_df_BM16"])
        bias_s = t_interval(z["S"] / chi - 1)
        var_s = z["S_sigma2_BM16"] / n_ret
        ratio_s = float(var_s.mean() / z["S"].var(ddof=1))
        ratio_s_ci = bootstrap(lambda v, e: v.mean() / e.var(ddof=1), (var_s, z["S"]), f"{cell}/S-S/variance_ratio")
        ok = inside(cov_s["cp95"], COVERAGE_BOUNDS) and inside(bias_s["ci95"], BIAS_BOUNDS) and inside(ratio_s_ci, RATIO_BOUNDS)
        gates[f"S-S[{cell}]"] = {"coverage": cov_s, "relative_bias": bias_s,
                                 "variance_ratio": {"ratio": ratio_s, "ci95": ratio_s_ci}, "result": "PASS" if ok else "FAIL"}

        half = student_t.ppf(0.975, 15) * np.sqrt(z["G_variance"])
        cov_g = coverage((z["G"] - half <= chi) & (chi <= z["G"] + half))
        bias_g = t_interval(z["G"] / chi - 1)
        ratio_g = float(z["G_variance"].mean() / z["G"].var(ddof=1))
        ratio_g_ci = bootstrap(lambda v, e: v.mean() / e.var(ddof=1), (z["G_variance"], z["G"]), f"{cell}/S-G/variance_ratio")
        ok = inside(cov_g["cp95"], COVERAGE_BOUNDS) and inside(bias_g["ci95"], BIAS_BOUNDS) and inside(ratio_g_ci, RATIO_BOUNDS)
        gates[f"S-G[{cell}]"] = {"coverage": cov_g, "relative_bias": bias_g,
                                 "variance_ratio": {"ratio": ratio_g, "ci95": ratio_g_ci}, "result": "PASS" if ok else "FAIL"}

        energy = t_interval(z["energy_mean"] / float(z["energy_exact"]) - 1)
        gates[f"S-C[{cell}]"] = {"relative_bias_energy": energy,
                                 "result": "PASS" if inside(energy["ci95"], BIAS_BOUNDS) else "FAIL"}
        bias_p = t_interval(z["P"] / chi - 1)
        # Amendment 01: R-hat from results/amend01 (corrected Blom scores; +inf for separated constant
        # chains, counted as flagged; NaN only when all draws are tied, counted as not flagged).
        a = np.load(HERE / "results" / "amend01" / f"{cell}.npz")
        flagged = {k: coverage(np.nan_to_num(a[f"rhat_{k}"], nan=0.0, posinf=np.inf) > RHAT_THRESHOLD)
                   for k in ("m", "m2", "e")}
        comparison = {}
        for k in ("m", "m2", "e"):
            old, new = z[f"rhat_{k}"], a[f"rhat_{k}"]
            both = np.isfinite(old) & np.isfinite(new)
            comparison[k] = {
                "max_abs_difference_finite": float(np.max(np.abs(old[both] - new[both]))) if both.any() else None,
                "original_flagged": int(np.sum(np.nan_to_num(old, nan=0.0) > RHAT_THRESHOLD)),
                "amended_flagged": flagged[k]["covered"],
                "original_nan": int(np.isnan(old).sum()), "amended_nan": int(np.isnan(new).sum()),
                "amended_inf": int(np.isposinf(new).sum())}
        AMENDMENT["rhat_original_vs_amended"][cell] = comparison
        N, T = int(z["N"]), float(z["T"])
        AMENDMENT["split_half_burn_in"][cell] = {
            "second_minus_first_m2_relative_to_chi": t_interval(N / T * (a["m2_second"] - a["m2_first"]) / chi),
            "second_minus_first_energy_relative": t_interval((a["e_second"] - a["e_first"]) / abs(float(z["energy_exact"]))),
            "reading": "descriptive; an interval excluding 0 indicates drift within the retained window (burn-in or slow mixing)"}
        if cell == "ising_L4_T1.5":
            gates["S-F1"] = {"relative_bias_P": bias_p, "criterion": "CI upper < -0.5",
                             "result": "REPRODUCED" if bias_p["ci95"][1] < -0.5 else "NOT_REPRODUCED"}
            gates["D-F1"] = {"rhat_m_not_converged": flagged["m"], "criterion": "CP95 lower >= 0.90",
                             "result": "HOLDS" if flagged["m"]["cp95"][0] >= 0.90 else "DOES_NOT_HOLD"}
        row = {"cell": cell, "L": int(z["L"]), "T": float(z["T"]), "chi_exact": chi,
               "S_coverage": cov_s["coverage"], "S_bias": bias_s["mean"], "S_var_ratio": ratio_s,
               "G_coverage": cov_g["coverage"], "G_bias": bias_g["mean"], "G_var_ratio": ratio_g,
               "P_bias": bias_p["mean"], "chi_abs_exact": float(z["chi_abs_exact"]),
               "chi_abs_plugin_bias": float(np.mean(z["chi_abs_plugin"] / float(z["chi_abs_exact"]) - 1)),
               "energy_bias": energy["mean"],
               "S_OBM_coverage": s_cov(z["S_sigma2_OBM"], z["S_df_OBM"])["coverage"],
               "S_SVTH_coverage": s_cov(z["S_sigma2_SVTH"], z["S_df_SVTH"])["coverage"],
               "S_GEYER_coverage": s_cov(z["S_sigma2_GEYER"], np.full(len(z["S"]), np.inf))["coverage"],
               "S_status_b96_adequate": float(np.mean(z["S_status_b"] == "ADEQUATE")),
               "S_tau_median": float(np.nanmedian(z["S_tau_GEYER"])),
               "S_adaptive_abstain": float(np.mean(z["S_bstar"] == 0)),
               "rhat_m_flagged": flagged["m"]["coverage"], "rhat_m2_flagged": flagged["m2"]["coverage"],
               "rhat_e_flagged": flagged["e"]["coverage"]}
        rows = z["S_bstar"] > 0
        if rows.any():
            lo, hi = longrun.interval(z["S"][rows], z["S_sigma2_adaptive_BM"][rows], n_ret, z["S_df_adaptive_BM"][rows])
            row["S_adaptive_coverage"] = float(np.mean(longrun.covers(lo, hi, chi)))
        table.append(row)


def verdicts(gates: dict, complete: dict) -> dict:
    expected = expected_gate_ids()
    missing = [g for ids in expected.values() for g in ids if g not in gates]
    if not complete["complete"] or missing:
        reason = {"problems": complete["problems"], "missing_gates": missing}
        return {h: "NOT_EVALUABLE_INCOMPLETE" for h in ("controls_valid", "H-S", "H-L", "H-D")} | {"incomplete": reason}

    def all_pass(prefix):
        return all(gates[g]["result"] == "PASS" for g in expected[prefix])

    l_d = {g: gates[g] for g in expected["L-D["]}
    families_ok = all(any(g["result"] == "PASS" for k, g in l_d.items() if k.endswith(f",{f}]"))
                      for f in GATES["L-D"]["families"])
    return {
        "controls_valid": all_pass("L-C[") and all_pass("S-C["),
        "H-S": gates.get("S-F1", {}).get("result") == "REPRODUCED" and all_pass("S-S[") and all_pass("S-G["),
        "H-L": all(gates.get(g, {}).get("result") == "REPRODUCED" for g in ("L-F1", "L-F2"))
               and all_pass("L-A[") and all_pass("L-G["),
        "H-D": gates.get("D-F1", {}).get("result") == "HOLDS"
               and not any(g["result"] == "FAIL" for g in l_d.values()) and families_ok,
    }


def write_csv(path: Path, rows: list) -> None:
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    complete = completeness()
    if not complete["complete"]:
        print(json.dumps({"complete": False, "problems": complete["problems"]}, indent=1))
        raise SystemExit("campaign outputs incomplete; no gates evaluated (Amendment 01)")
    gates: dict = {}
    ar1_table: list = []
    ising_table: list = []
    analyze_ar1(gates, ar1_table)
    analyze_ising(gates, ising_table)
    out = {"protocol_id": PROTOCOL["id"], "script": "research/v018-estimators/analyze.py",
           "amendment": "AMENDMENT_01.md", "completeness": complete,
           "verdicts": verdicts(gates, complete), "gates": gates, "amendment_01_secondary": AMENDMENT}
    (HERE / "results" / "gates.json").write_text(json.dumps(out, indent=1) + "\n")
    write_csv(HERE / "results" / "ar1_table.csv", ar1_table)
    write_csv(HERE / "results" / "ising_table.csv", ising_table)
    counts: dict = {}
    for gate in gates.values():
        counts[gate["result"]] = counts.get(gate["result"], 0) + 1
    print(json.dumps({"verdicts": out["verdicts"], "gate_result_counts": counts}, indent=1))


if __name__ == "__main__":
    main()
