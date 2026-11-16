"""Design numbers of the B2 verification, from analytic references only (no random numbers).

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/design.py

Writes design.json: scheme digests (S1-S3 must equal the frozen Layer A
schemes), the eligible states and pairs of the ensemble tests, the planned
name of every test in each family, the expected applicability of the path
tests, and the approximate power of the positive controls.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

from rnglab.nuclear import DecayScheme, decay_probabilities, state_probabilities
from rnglab.nuclear.aggregate import AggregatePlan
from rnglab.nuclear.checks import multinomial_check
from rnglab.nuclear.comparison import planned_path_cells
from rnglab.nuclear.decay import Mode

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ENS_TIMES = (2.0, 8.0)
PATH_MIN_EXPECTED = 10.0


def naive_probabilities(scheme: DecayScheme, start: str, times) -> np.ndarray:
    """State law of the deliberately wrong one-step aggregation (at most one decay per interval)."""
    size = len(scheme)
    p = np.zeros(size)
    p[scheme.index[start]] = 1.0
    out = [p.copy()]
    for delta in np.diff(times):
        step = np.zeros((size, size))
        for i, nuclide in enumerate(scheme.nuclides):
            if nuclide.stable:
                step[i, i] = 1.0
                continue
            stay = math.exp(-math.log(2.0) / nuclide.half_life * delta)
            step[i, i] = stay
            for channel, probability in zip(nuclide.channels, nuclide.channel_probabilities()):
                step[scheme.index[channel.daughter], i] += (1.0 - stay) * probability
        p = step @ p
        out.append(p.copy())
    return np.array(out)


def chi2_applicable(total: int, probabilities) -> bool:
    """Whether the pooling rule leaves two or more cells: it depends only on the total and the probabilities."""
    probabilities = np.asarray(probabilities, dtype=np.float64)
    counts = np.zeros(probabilities.size, dtype=np.int64)
    counts[int(np.argmax(probabilities))] = total
    try:
        multinomial_check(counts, probabilities)
    except ValueError as error:
        if "fewer than two cells" in str(error):
            return False
        raise
    return True


def chi2_power(noncentrality: float, dof: int, level: float) -> float:
    return float(stats.ncx2.sf(stats.chi2.isf(level, dof), dof, noncentrality))


def single_start(case: dict) -> str | None:
    return next(iter(case["initial"])) if len(case["initial"]) == 1 else None


def planned_battery(key: str, scheme: DecayScheme, case: dict, prefix: str) -> list[str]:
    names = []
    times = case["times"]
    groups = groups_of(scheme, case)
    for start, members in groups:
        n = case["initial"][start]
        p = state_probabilities(scheme, start, times)
        for it, t in enumerate(times[1:], start=1):
            if chi2_applicable(n, p[it, members]):
                names.append(f"{prefix} {key} BT-STATE {start} t={t!r}")
    for start, members in groups:
        for j in members:
            if not scheme.nuclides[j].stable:
                names.append(f"{prefix} {key} BT-DECAY {scheme.nuclides[j].name}")
    for start, _ in groups:
        n = case["initial"][start]
        q = decay_probabilities(scheme, start, times)[:, scheme.index[start]]
        if chi2_applicable(n, np.append(np.diff(q), 1.0 - q[-1])):
            names.append(f"{prefix} {key} BT-BINS {start}")
    for nuclide in scheme.nuclides:
        if len(nuclide.channels) >= 2:
            names.append(f"{prefix} {key} BT-CHANNEL {nuclide.name}")
        for channel in nuclide.channels:
            if 0.0 < channel.photon_yield < 1.0:
                names.append(f"{prefix} {key} BT-PHOTON {nuclide.name}/{channel.mode.value}->{channel.daughter}")
    if key == "S4":
        for start, _ in groups:
            i = scheme.index[start]
            n = case["initial"][start]
            k = next(c for c, ch in enumerate(scheme.nuclides[i].channels) if ch.mode is Mode.BETA_DELAYED_NEUTRON)
            pn = scheme.nuclides[i].channel_probabilities()[k]
            q = decay_probabilities(scheme, start, times)[:, i]
            if chi2_applicable(n, np.append(pn * np.diff(q), 1.0 - pn * q[-1])):
                names.append(f"{prefix} {key} BT-NEUTRON {start}")
    return names


def groups_of(scheme: DecayScheme, case: dict) -> list[tuple[str, list[int]]]:
    """(start nuclide, the states it can reach) for every initial nuclide."""
    out = []
    for start in case["initial"]:
        reach, frontier = {scheme.index[start]}, [scheme.index[start]]
        while frontier:
            i = frontier.pop()
            for channel in scheme.nuclides[i].channels:
                j = scheme.index[channel.daughter]
                if j not in reach:
                    reach.add(j)
                    frontier.append(j)
        out.append((start, sorted(reach)))
    return out


def path_plan(scheme: DecayScheme, case: dict) -> dict:
    """Cells of the BT-PATH tests, fixed before generation: {source name: {interval index: cells}}.

    Expected start counts come from the exact law; cells follow comparison.pooled_cells with a
    minimum expected count of 10 (stricter than the 5 used elsewhere, because the realised start
    counts fluctuate around the expected ones).
    """
    start = single_start(case)
    n0 = case["initial"][start]
    times = case["times"]
    plan = AggregatePlan(scheme, times, n0)
    expected = [state_probabilities(scheme, start, t)[0] * n0 for t in times[:-1]]
    cells = planned_path_cells(plan, expected, min_expected=PATH_MIN_EXPECTED)
    return {scheme.nuclides[i].name: {str(k): c for k, c in per_interval.items()} for i, per_interval in cells.items()}


def main() -> None:
    protocol = json.loads((HERE / "protocol.json").read_text())
    layer_a = json.loads((REPO / "research/v020-nuclear-a-verify/design.json").read_text())
    schemes = {k: DecayScheme.from_dict(v) for k, v in protocol["schemes"].items()}
    digests = {k: s.sha256() for k, s in schemes.items()}
    for key in ("S1", "S2", "S3"):
        if digests[key] != layer_a["scheme_sha256"][key]:
            raise SystemExit(f"{key} differs from the frozen Layer A scheme")
    cases = protocol["cases"]
    s2 = schemes["S2"]
    n_ens = protocol["runs"]["ENS-B2"]["initial"]["Xa-400"]
    names = [n.name for n in s2.nuclides]
    eligible = {}
    for t in ENS_TIMES:
        p = state_probabilities(s2, "Xa-400", t)[0]
        eligible[repr(t)] = [names[j] for j in range(len(s2)) if n_ens * p[j] * (1 - p[j]) >= 20.0]
    at2, at8 = eligible["2.0"], eligible["8.0"]
    pairs_t2 = [(a, b) for i, a in enumerate(at2) for b in at2[i + 1:]]
    pairs_two_time = [(j, j) for j in at2 if j in at8] + [("Xa-400", k) for k in at8 if k != "Xa-400"]
    unstable = [n.name for n in s2.nuclides if not n.stable]
    yield_channels = [f"{n.name}/{c.mode.value}->{c.daughter}" for n in s2.nuclides for c in n.channels
                      if 0.0 < c.photon_yield < 1.0]

    path_cells = path_plan(s2, cases["S2"])
    path_detail = {"min_expected": PATH_MIN_EXPECTED, "cells": path_cells,
                   "dof": {name: sum(len(c) - 1 for c in per.values()) for name, per in path_cells.items()}}
    families = {"NC1": [], "NC2": []}
    for family in families:
        for key in ("S1", "S2", "S3", "S4", "S5"):
            families[family] += planned_battery(key, schemes[key], cases[key], family)
        families[family] += [f"{family} S2 BT-PATH {name}" for name in path_cells]
    ens = ([f"ENS ET-MEAN {j} t={t}" for t in ("2.0", "8.0") for j in eligible[t]]
           + [f"ENS ET-VAR {j} t={t}" for t in ("2.0", "8.0") for j in eligible[t]]
           + [f"ENS ET-VARSUM {a}+{b} t=2.0" for a, b in pairs_t2]
           + [f"ENS ET-VARSUM2 {a}@2+{b}@8" for a, b in pairs_two_time])
    comparison = ([f"COMP CT-HOMO t={t!r}" for t in cases["S2"]["times"][1:]]
                  + [f"COMP CT-VAR {j} t=2.0" for j in at2]
                  + [f"COMP CT-VARSUM {a}+{b} t=2.0" for a, b in pairs_t2]
                  + [f"COMP CT-DECAY {j}" for j in unstable]
                  + [f"COMP CT-PHOTON {c}" for c in yield_channels])
    big = []
    for size in protocol["runs"]["BIG"]["initial_sizes"]:
        p_big = state_probabilities(s2, "Xa-400", cases["S2"]["times"])
        big += [f"BIG {size} BIG-STATE t={t!r}" for it, t in enumerate(cases["S2"]["times"][1:], start=1)
                if chi2_applicable(size, p_big[it])]
        big += [f"BIG {size} BIG-DECAY {j}" for j in unstable]
    families["NC1"] += ens + comparison + big
    alpha = protocol["families"]["alpha"]
    level = alpha / len(families["NC1"])

    # Positive controls.
    s1 = schemes["S1"]
    n_pc1 = protocol["runs"]["PC-B1"]["initial"]["Sa-400"]
    p_ref = 2.0 ** -2.0
    p_true = 2.0 ** (-2.0 / protocol["runs"]["PC-B1"]["perturbation"]["Sa-400"]["half_life"])
    z_pc1 = n_pc1 * (p_true - p_ref) / math.sqrt(n_pc1 * p_ref * (1 - p_ref))
    critical = stats.norm.isf(level / 2)
    power_pc1 = float(stats.norm.sf(critical - z_pc1) + stats.norm.cdf(-critical - z_pc1))
    times2 = cases["S2"]["times"]
    exact12 = state_probabilities(s2, "Xa-400", 12.0)[0]
    naive12 = naive_probabilities(s2, "Xa-400", times2)[-1]
    keep = exact12 > 0
    n_pc2 = protocol["runs"]["PC-B2"]["initial"]["Xa-400"]
    lam_single = float(n_pc2 * np.sum((naive12[keep] - exact12[keep]) ** 2 / exact12[keep]))
    pooled = protocol["runs"]["PC-B2"]["replicates"] * protocol["runs"]["PC-B2"]["replicate_initial"]["Xa-400"]
    # Homogeneity of two samples of M nuclei each with laws a and b, m = (a + b) / 2:
    # each row contributes M sum (a - m)^2 / m = (M / 4) sum (a - b)^2 / m, so the total is (M / 2) sum (a - b)^2 / m.
    mid = (naive12 + exact12) / 2.0
    keep2 = mid > 0
    lam_comp = float(pooled / 2.0 * np.sum((naive12[keep2] - exact12[keep2]) ** 2 / mid[keep2]))
    dof = int(keep.sum()) - 1

    plans = {}
    for key, scheme in schemes.items():
        plan = AggregatePlan(scheme, cases[key]["times"], sum(cases[key]["initial"].values()))
        plans[key] = {"layout": plan.layout(), "diagnostics": plan.diagnostics}

    design = {
        "generated_by": "research/v020-nuclear-b2-verify/design.py (analytic references only, no random numbers)",
        "scheme_sha256": digests,
        "layer_a_schemes_identical": True,
        "ens": {"times": list(ENS_TIMES), "per_replicate": n_ens, "eligible": eligible, "pairs_t2": pairs_t2,
                "pairs_two_time": pairs_two_time},
        "comparison": {"unstable": unstable, "yield_channels": yield_channels},
        "path_tests": path_detail,
        "planned_tests": families,
        "family_sizes": {k: len(v) for k, v in families.items()},
        "positive_control_level": level,
        "positive_controls": {
            "PC-B1": {"p_reference": p_ref, "p_true": p_true, "expected_z": z_pc1, "approximate_power": power_pc1},
            "PC-B2": {"exact_p12": exact12.tolist(), "naive_p12": naive12.tolist(), "dof": dof,
                      "noncentrality_single_run": lam_single,
                      "approximate_power_single_run": chi2_power(lam_single, dof, level),
                      "noncentrality_comparison": lam_comp,
                      "approximate_power_comparison": chi2_power(lam_comp, dof, level)},
            "method": "normal (PC-B1) and noncentral chi-square (PC-B2) approximations at the level alpha / NC1 family size",
        },
        "plans": plans,
    }
    (HERE / "design.json").write_text(json.dumps(design, indent=2) + "\n")
    print(json.dumps({"family_sizes": design["family_sizes"], "level": level, "eligible": eligible,
                      "pairs_t2": len(pairs_t2), "pairs_two_time": len(pairs_two_time),
                      "path_tests": {k: (len(v), path_detail["dof"][k]) for k, v in path_cells.items()},
                      "PC-B1 power": power_pc1, "PC-B2 power single": chi2_power(lam_single, dof, level),
                      "PC-B2 power comparison": chi2_power(lam_comp, dof, level),
                      "plans": {k: v["layout"]["words_per_interval"] for k, v in plans.items()}}, indent=1))


if __name__ == "__main__":
    main()
