"""Design numbers of the Layer A verification, from the analytic references only.

Draws no random numbers. Writes design.json next to this file: scheme digests,
the ENS eligibility list, the size of each test family, the smallest expected
cell of every chi-square test before pooling, and the approximate power of the
two positive controls at the Bonferroni level of the NC1 family.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-a-verify/design.py
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy import stats

from rnglab.nuclear import DecayScheme, decay_probabilities, state_probabilities

HERE = Path(__file__).resolve().parent


def unit_bin_probabilities(scheme: DecayScheme, start: str, nuclide: str, horizon: float) -> np.ndarray:
    edges = np.arange(0.0, horizon + 0.5, 1.0)
    q = decay_probabilities(scheme, start, edges)[:, scheme.index[nuclide]]
    return np.append(np.diff(q), 1.0 - q[-1])


def two_sided_power(effect_z: float, level: float) -> float:
    critical = stats.norm.isf(level / 2.0)
    return float(stats.norm.sf(critical - effect_z) + stats.norm.cdf(-critical - effect_z))


def main() -> None:
    protocol = json.loads((HERE / "protocol.json").read_text())
    schemes = {key: DecayScheme.from_dict(value) for key, value in protocol["schemes"].items()}
    cases = protocol["cases"]
    tests = protocol["tests"]
    s2 = schemes["S2"]
    n2 = cases["S2"]["initial"]["Xa-400"]

    ens = protocol["runs"]["ENS"]
    per_replicate = ens["initial"]["Xa-400"]
    p_ens = state_probabilities(s2, "Xa-400", ens["time"])[0]
    binomial_variance = per_replicate * p_ens * (1.0 - p_ens)
    eligible = [n.name for n, v in zip(s2.nuclides, binomial_variance) if v >= 20.0]

    replicates = ens["replicates"]
    eligible_mask = binomial_variance >= 20.0
    pq = p_ens * (1.0 - p_ens)
    excess_kurtosis = (1.0 - 6.0 * pq[eligible_mask]) / (per_replicate * pq[eligible_mask])
    # Var(s^2) = sigma^4 (2/(R-1) + kappa/R); normal theory assumes kappa = 0.
    variance_inflation = 1.0 + excess_kurtosis * (replicates - 1) / (2.0 * replicates)
    t10b = {"max_excess_kurtosis": float(excess_kurtosis.max()),
            "max_relative_variance_change": float(variance_inflation.max() - 1.0),
            "max_relative_sd_change": float(np.sqrt(variance_inflation.max()) - 1.0)}

    unstable_s2 = [n.name for n in s2.nuclides if not n.stable]
    yield_channels = [f"{n.name}/{c.mode.value}" for n in s2.nuclides for c in n.channels if 0.0 < c.photon_yield < 1.0]
    sizes = {"T1": 1, "T2": len(tests["T2"]["times"]), "T3": 1, "T4": len(tests["T4"]["times"]),
             "T5": len(unstable_s2), "T6": len(tests["T6"]["nuclides"]), "T7": len(tests["T7"]["nuclides"]),
             "T8": len(yield_channels), "T9": len(tests["T9"]["times"])}
    nc2 = sum(sizes.values())
    nc1 = nc2 + 2 * len(eligible)
    level = protocol["families"]["alpha"] / nc1

    smallest_expected = {}
    for key, scheme_key, start in (("T4", "S2", "Xa-400"), ("T9", "S3", "Ya-410")):
        count = cases[scheme_key]["initial"][start]
        probabilities = state_probabilities(schemes[scheme_key], start, tests[key]["times"])
        smallest_expected[key] = {repr(t): float(count * row[row > 0].min())
                                  for t, row in zip(tests[key]["times"], probabilities)}
    for nuclide in tests["T7"]["nuclides"]:
        probabilities = unit_bin_probabilities(s2, "Xa-400", nuclide, cases["S2"]["horizon"])
        smallest_expected[f"T7 {nuclide}"] = float(n2 * probabilities.min())

    n1 = cases["S1"]["initial"]["Sa-400"]
    p_ref = 2.0 ** -2.0
    p_true = 2.0 ** (-2.0 / protocol["runs"]["PC1"]["perturbation"]["Sa-400"]["half_life"])
    z_pc1 = n1 * (p_true - p_ref) / math.sqrt(n1 * p_ref * (1.0 - p_ref))
    q_xa = float(decay_probabilities(s2, "Xa-400", cases["S2"]["horizon"])[0][s2.index["Xa-400"]])
    expected_decays = n2 * q_xa
    b_ref, b_true = 0.64, protocol["runs"]["PC2"]["perturbation"]["Xa-400"]["branching"][0]
    z_pc2 = expected_decays * (b_true - b_ref) / math.sqrt(expected_decays * b_ref * (1.0 - b_ref))

    # Rounding-error bound of the closed-form Bateman route (the quantity bateman_chain refuses above
    # ROUNDING_LIMIT) over every path of every scheme, at the A1 times of that scheme and at t = 0.
    worst_bound = 0.0
    for key, scheme in schemes.items():
        times = np.array(sorted({0.0, *tests[{"S1": "T2", "S2": "T4", "S3": "T9"}[key]]["times"],
                                 cases[key]["horizon"]}))
        rates = [0.0 if n.stable else math.log(2.0) / n.half_life for n in scheme.nuclides]
        for nuclide in scheme.nuclides:
            for path, _ in scheme.paths(nuclide.name):
                r = [rates[i] for i in path]
                terms = np.array([np.exp(-x * times) / math.prod(r[k] - x for k in range(len(r)) if k != q)
                                  for q, x in enumerate(r)])
                bound = 4.0 * len(r) * np.finfo(float).eps * math.prod(r[:-1]) * np.abs(terms).sum(axis=0)
                worst_bound = max(worst_bound, float(np.max(bound)))

    design = {
        "generated_by": "research/v020-nuclear-a-verify/design.py (analytic references only, no random numbers)",
        "scheme_sha256": {key: scheme.sha256() for key, scheme in schemes.items()},
        "bateman_rounding_bound_max": {"value": worst_bound, "limit": 1e-10,
                                       "scope": "every decay path of S1, S2, S3 from every start, at t = 0, "
                                                "the test times and the horizon"},
        "ens": {"time": ens["time"], "per_replicate": per_replicate,
                "p": dict(zip([n.name for n in s2.nuclides], p_ens.tolist())),
                "binomial_variance": dict(zip([n.name for n in s2.nuclides], binomial_variance.tolist())),
                "eligible": eligible, "T10b_normal_theory_error": t10b},
        "family_sizes": {"NC1": nc1, "NC2": nc2, "per_test": sizes, "T8_channels": yield_channels,
                         "T5_nuclides": unstable_s2},
        "positive_control_level": level,
        "smallest_expected_cell_before_pooling": smallest_expected,
        "positive_controls": {
            "PC1": {"target": "T2 at t = 2", "p_reference": p_ref, "p_true": p_true, "expected_z": z_pc1,
                    "approximate_power": two_sided_power(z_pc1, level)},
            "PC2": {"target": "T6 Xa-400", "expected_decays": expected_decays, "b_reference": b_ref,
                    "b_true": b_true, "expected_z": z_pc2, "approximate_power": two_sided_power(z_pc2, level)},
            "method": "normal approximation to the exact binomial test at the Bonferroni level alpha / NC1 family size",
        },
    }
    (HERE / "design.json").write_text(json.dumps(design, indent=2) + "\n")
    print(json.dumps({"eligible": eligible, "NC1": nc1, "NC2": nc2, "level": level,
                      "PC1_power": design["positive_controls"]["PC1"]["approximate_power"],
                      "PC2_power": design["positive_controls"]["PC2"]["approximate_power"],
                      "smallest_expected": smallest_expected}, indent=1))


if __name__ == "__main__":
    main()
