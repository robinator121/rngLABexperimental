"""Design numbers of the B2 verification, from analytic references only (no random numbers).

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/design.py

Writes design.json: scheme digests (S1-S3 must equal the frozen Layer A
schemes), the eligible states and pairs of the ensemble tests, the planned
name of every test in each family, the cells fixed before generation (path,
per-interval photon and Layer A / B2 homogeneity tests, taken from the same
functions run_verification.py uses), the approximate power of the positive
controls, and the sensitivity of the variance tests (the variance ratios they
detect), so that a non-rejection is read for what it is.
"""
from __future__ import annotations

import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

from rnglab.nuclear import DecayScheme, decay_probabilities, state_probabilities
from rnglab.nuclear.aggregate import AggregatePlan
from rnglab.nuclear.binomial import PART_MAX
from rnglab.nuclear.checks import multinomial_check
from rnglab.nuclear.comparison import planned_path_cells
from rnglab.nuclear.decay import Mode
from rnglab.nuclear.moments import state_count_moments, state_sum_moments, two_time_sum_moments

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
ENS_TIMES = (2.0, 8.0)
PATH_MIN_EXPECTED = 10.0
SENSITIVITY_RATIOS = (1.1, 1.25, 1.5, 2.0)


def load_verification():
    """run_verification.py as a module (its main() is not run): the planned cells come from its functions."""
    spec = importlib.util.spec_from_file_location("b2_verify_run_verification", HERE / "run_verification.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def variance_test_power(nu: float, level: float, ratio: float, kind: str) -> float:
    """Approximate power of a two-sided variance test at ``level`` when the true variance is ``ratio`` times
    the null one (scaled chi-square: ET-*) or B2's variance is ``ratio`` times Layer A's (F ratio: CT-*)."""
    if kind == "scaled_chi2":
        low, high = stats.chi2.ppf(level / 2, nu), stats.chi2.isf(level / 2, nu)
        return float(stats.chi2.cdf(low / ratio, nu) + stats.chi2.sf(high / ratio, nu))
    low, high = stats.f.ppf(level / 2, nu, nu), stats.f.isf(level / 2, nu, nu)
    return float(stats.f.cdf(low * ratio, nu, nu) + stats.f.sf(high * ratio, nu, nu))


def detectable_ratio(nu: float, level: float, power: float, kind: str) -> float:
    """Smallest variance ratio above 1 detected with the given power."""
    from scipy.optimize import brentq
    return float(brentq(lambda r: variance_test_power(nu, level, r, kind) - power, 1.0 + 1e-9, 1e3))


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


def planned_battery(key: str, scheme: DecayScheme, case: dict, prefix: str, verify) -> list[str]:
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
    for j, k, _ in verify.photon_interval_cells(scheme, case):
        names.append(f"{prefix} {key} BT-PHOTON-T {verify.channel_label(scheme.nuclides[j], scheme.nuclides[j].channels[k])}")
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
    verify = load_verification()
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
            families[family] += planned_battery(key, schemes[key], cases[key], family, verify)
        families[family] += [f"{family} S2 BT-PATH {name}" for name in path_cells]
    photon_cells = {key: {verify.channel_label(schemes[key].nuclides[j], schemes[key].nuclides[j].channels[k]): cells
                          for j, k, cells in verify.photon_interval_cells(schemes[key], cases[key])}
                    for key in ("S1", "S2", "S3", "S4", "S5")}
    reps = protocol["runs"]["ENS-B2"]["replicates"]
    times2 = cases["S2"]["times"]
    comp_cells = verify.comparison_cells(s2, times2, reps, n_ens)
    ens = ([f"ENS ET-MEAN {j} t={t}" for t in ("2.0", "8.0") for j in eligible[t]]
           + [f"ENS ET-VAR {j} t={t}" for t in ("2.0", "8.0") for j in eligible[t]]
           + [f"ENS ET-VARSUM {a}+{b} t=2.0" for a, b in pairs_t2]
           + [f"ENS ET-VARSUM2 {a}@2+{b}@8" for a, b in pairs_two_time])
    comparison = ([f"COMP CT-HOMO t={t!r}" for t in times2[1:] if f"COMP CT-HOMO t={t!r}" in comp_cells]
                  + [f"COMP CT-VAR {j} t=2.0" for j in at2]
                  + [f"COMP CT-VARSUM {a}+{b} t=2.0" for a, b in pairs_t2]
                  + [f"COMP CT-DECAY {j}" for j in unstable if f"COMP CT-DECAY {j}" in comp_cells]
                  + [f"COMP CT-PHOTON {c}" for c in yield_channels if f"COMP CT-PHOTON {c}" in comp_cells])
    big = []
    for size in protocol["runs"]["BIG"]["initial_sizes"]:
        p_big = state_probabilities(s2, "Xa-400", cases["S2"]["times"])
        big += [f"BIG {size} BIG-STATE t={t!r}" for it, t in enumerate(cases["S2"]["times"][1:], start=1)
                if chi2_applicable(size, p_big[it])]
        big += [f"BIG {size} BIG-DECAY {j}" for j in unstable]
    families["NC1"] += ens + comparison + big
    alpha = protocol["families"]["alpha"]
    level = alpha / len(families["NC1"])
    big_parts = [-(-size // PART_MAX) for size in protocol["runs"]["BIG"]["initial_sizes"]]
    if big_parts != protocol["runs"]["BIG"]["parts"]:
        raise SystemExit(f"BIG parts {big_parts} differ from protocol.json")
    # calibration.py draws the BIG null with numpy's binomial and multinomial, which return doubles cast to
    # integers: at counts above 2**53 the null lies on a lattice whose spacing is at most the float spacing at
    # twice the count. Relative to the standard deviation of the state counts the BIG tests use:
    null_lattice = []
    for size in protocol["runs"]["BIG"]["initial_sizes"]:
        p_big = state_probabilities(s2, "Xa-400", cases["S2"]["times"])
        for it, t in enumerate(cases["S2"]["times"][1:], start=1):
            for j, q in enumerate(np.asarray(p_big[it], dtype=np.float64).tolist()):
                mean = size * q
                if 0.0 < q < 1.0 and mean >= 2.0 ** 53:
                    spacing = float(np.spacing(2.0 * mean))
                    null_lattice.append({"size": size, "t": t, "nuclide": names[j], "mean": mean, "spacing": spacing,
                                         "spacing_over_sd": spacing / math.sqrt(size * q * (1.0 - q))})
    null_lattice_worst = max(null_lattice, key=lambda row: row["spacing_over_sd"]) if null_lattice else None
    # numpy's binomial is also wrong at n above 2**53 when the expected count is small (numerics.py
    # numpy_binomial: at n = 2**60 + 12345 the variance is 1.2 times too large at mean 115, the mean is biased at
    # 1,000 and 10,000, and the variance is still slightly inflated at 1e6). The BIG null calls it inside multinomial (categories in order, each a binomial in the count
    # still unassigned with the conditional probability) and for photons (a binomial in the decay count). With
    # expected counts in place of draws, the smallest expected count min(n p, n (1 - p)) of a call with n > 2**53:
    numpy_calls = []
    for size in protocol["runs"]["BIG"]["initial_sizes"]:
        plan = AggregatePlan(s2, cases["S2"]["times"], size)
        p_big = np.asarray(state_probabilities(s2, "Xa-400", cases["S2"]["times"]), dtype=np.float64)
        steps = plan.category_steps.astype(np.float64)
        for k, delta in enumerate(plan.deltas.tolist()):
            expected = np.zeros(plan.n_categories)
            for i, categories in enumerate(plan.categories):
                start, first = size * p_big[k, i], int(plan.first_category[i])
                if len(categories) == 1:
                    expected[first] = start
                    continue
                q = plan.interval_plans[delta][i].law.as_float()
                law = q / q.sum()
                expected[first:first + len(law)] = start * law
                left, remaining_p = start, 1.0
                for c in range(len(law) - 1):
                    p = min(float(law[c] / remaining_p), 1.0)
                    if left > 2.0 ** 53:
                        numpy_calls.append({"size": size, "interval": k + 1, "nuclide": names[i], "call": "multinomial",
                                            "category": c, "n": left, "p": p, "expected": left * min(p, 1.0 - p)})
                    left -= start * law[c]
                    remaining_p -= law[c]
            decays = (expected @ steps).reshape(plan.n_nuclides, plan.n_channels)
            for j, channel, y in plan.yield_channels:
                if decays[j, channel] > 2.0 ** 53:
                    numpy_calls.append({"size": size, "interval": k + 1, "nuclide": names[j], "call": "photons",
                                        "channel": channel, "n": float(decays[j, channel]), "p": y,
                                        "expected": float(decays[j, channel]) * min(y, 1.0 - y)})
    numpy_smallest = min(numpy_calls, key=lambda row: row["expected"]) if numpy_calls else None
    if numpy_smallest is not None and numpy_smallest["expected"] < 1e12:
        raise SystemExit(f"the BIG null calls numpy's binomial above 2**53 with an expected count of "
                         f"{numpy_smallest['expected']:.3g}, where numpy is not reliable; see numerics.py")

    # Sensitivity of the variance tests at the Holm level (approximate: nu taken from the null moments).
    index = {name: j for j, name in enumerate(names)}
    nus = {"ET-VAR": [], "ET-VARSUM": [], "ET-VARSUM2": [], "CT-VAR": [], "CT-VARSUM": []}
    for t in ENS_TIMES:
        for j in eligible[repr(t)]:
            m = state_count_moments(s2, "Xa-400", n_ens, t, index[j], index[j])
            nu = 2.0 * m.sigma_ab ** 2 / m.sample_covariance_variance(reps)
            nus["ET-VAR"].append(nu)
            if t == 2.0:
                nus["CT-VAR"].append(nu)
    for a, b in pairs_t2:
        m = state_sum_moments(s2, "Xa-400", n_ens, 2.0, index[a], index[b])
        nu = 2.0 * m.sigma_ab ** 2 / m.sample_covariance_variance(reps)
        nus["ET-VARSUM"].append(nu)
        nus["CT-VARSUM"].append(nu)
    for a, b in pairs_two_time:
        m = two_time_sum_moments(s2, "Xa-400", n_ens, 2.0, 8.0, index[a], index[b])
        nus["ET-VARSUM2"].append(2.0 * m.sigma_ab ** 2 / m.sample_covariance_variance(reps))
    sensitivity = {"level": level, "replicates": reps, "method": (
        "two-sided tests at the Holm level alpha / NC1 size; nu = 2 sigma^4 / Var(s^2) from the exact null moments, "
        "kept under the alternative. ET-*: true variance = ratio x null variance (scaled chi-square). CT-*: B2's "
        "variance = ratio x Layer A's (F(nu, nu)). Non-rejection does not establish equivalence: ratios below "
        "'ratio_detected_with_power_0.5' are more likely missed than found."), "tests": {}}
    for kind_name, values in nus.items():
        kind = "f_ratio" if kind_name.startswith("CT") else "scaled_chi2"
        sensitivity["tests"][kind_name] = {
            "instances": len(values), "nu_min": min(values), "nu_max": max(values),
            "power_at_ratio": {repr(r): [min(variance_test_power(nu, level, r, kind) for nu in values),
                                         max(variance_test_power(nu, level, r, kind) for nu in values)]
                               for r in SENSITIVITY_RATIOS},
            "ratio_detected_with_power_0.5": [min(detectable_ratio(nu, level, 0.5, kind) for nu in values),
                                              max(detectable_ratio(nu, level, 0.5, kind) for nu in values)],
            "ratio_detected_with_power_0.9": [min(detectable_ratio(nu, level, 0.9, kind) for nu in values),
                                              max(detectable_ratio(nu, level, 0.9, kind) for nu in values)]}

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
    # Homogeneity of two samples of M nuclei each with laws a and b, m = (a + b) / 2, on the cells fixed for
    # CT-HOMO at t = 12: each row contributes M sum (a - m)^2 / m = (M / 4) sum (a - b)^2 / m over cells,
    # so the total is (M / 2) sum (a - b)^2 / m.
    cells12 = comp_cells[f"COMP CT-HOMO t={times2[-1]!r}"]
    naive_cells = np.array([naive12[c].sum() for c in cells12])
    exact_cells = np.array([exact12[c].sum() for c in cells12])
    mid = (naive_cells + exact_cells) / 2.0
    lam_comp = float(pooled / 2.0 * np.sum((naive_cells - exact_cells) ** 2 / mid))
    dof_comp = len(cells12) - 1
    dof = int(keep.sum()) - 1

    plans = {}
    for key, scheme in schemes.items():
        plan = AggregatePlan(scheme, cases[key]["times"], sum(cases[key]["initial"].values()))
        plans[key] = {"layout": plan.layout(), "diagnostics": dict(plan.diagnostics)}

    design = {
        "generated_by": "research/v020-nuclear-b2-verify/design.py (analytic references only, no random numbers)",
        "scheme_sha256": digests,
        "layer_a_schemes_identical": True,
        "ens": {"times": list(ENS_TIMES), "per_replicate": n_ens, "eligible": eligible, "pairs_t2": pairs_t2,
                "pairs_two_time": pairs_two_time},
        "comparison": {"unstable": unstable, "yield_channels": yield_channels,
                       "min_expected": verify.COMPARISON_MIN_EXPECTED, "cells": comp_cells},
        "path_tests": path_detail,
        "photon_interval_tests": {"min_expected": verify.PHOTON_MIN_EXPECTED, "cells": photon_cells},
        "big_parts": big_parts,
        "calibration_null_lattice": {
            "note": "numpy's binomial is not exact above 2**53: its results lie on a lattice of spacing at most "
                    "np.spacing(2 * count), and for small expected counts its mean and variance are wrong "
                    "(numerics.py numpy_binomial). The BIG null of calibration.py stays where neither matters: the "
                    "lattice spacing is at most 'worst' standard deviations of a BIG state count, and every numpy "
                    "binomial call above 2**53 has an expected count of at least numpy_calls_above_2_53 "
                    "(asserted here to be at least 1e12); numerics.py measures numpy at that count",
            "cells_above_2_53": len(null_lattice), "worst": null_lattice_worst,
            "numpy_calls_above_2_53": {"calls": len(numpy_calls),
                                       "smallest_expected_count": numpy_smallest["expected"] if numpy_smallest else None,
                                       "smallest": numpy_smallest}},
        "sensitivity": sensitivity,
        "planned_tests": families,
        "family_sizes": {k: len(v) for k, v in families.items()},
        "positive_control_level": level,
        "positive_controls": {
            "PC-B1": {"p_reference": p_ref, "p_true": p_true, "expected_z": z_pc1, "approximate_power": power_pc1},
            "PC-B2": {"exact_p12": exact12.tolist(), "naive_p12": naive12.tolist(), "dof": dof,
                      "noncentrality_single_run": lam_single,
                      "approximate_power_single_run": chi2_power(lam_single, dof, level),
                      "noncentrality_comparison": lam_comp, "dof_comparison": dof_comp,
                      "approximate_power_comparison": chi2_power(lam_comp, dof_comp, level)},
            "method": "normal (PC-B1) and noncentral chi-square (PC-B2) approximations at the level alpha / NC1 family size",
        },
        "plans": plans,
    }
    (HERE / "design.json").write_text(json.dumps(design, indent=2) + "\n")
    print(json.dumps({"family_sizes": design["family_sizes"], "level": level, "eligible": eligible,
                      "pairs_t2": len(pairs_t2), "pairs_two_time": len(pairs_two_time),
                      "path_tests": {k: (len(v), path_detail["dof"][k]) for k, v in path_cells.items()},
                      "PC-B1 power": power_pc1, "PC-B2 power single": chi2_power(lam_single, dof, level),
                      "PC-B2 power comparison": chi2_power(lam_comp, dof_comp, level),
                      "photon_interval_tests": sum(len(v) for v in photon_cells.values()),
                      "big_parts": big_parts,
                      "sensitivity_ratio_power_0.5": {k: v["ratio_detected_with_power_0.5"]
                                                      for k, v in sensitivity["tests"].items()},
                      "plans": {k: v["layout"]["words_per_interval"] for k, v in plans.items()}}, indent=1))


if __name__ == "__main__":
    main()
