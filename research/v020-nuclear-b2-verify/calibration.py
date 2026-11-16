"""Null calibration of the approximate tests of the B2 verification, by direct sampling of the exact laws.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/calibration.py [--scale 1.0]

No RNGLab engine is involved. Counts are drawn from the exact laws with numpy's own Generator and
multinomial sampler (namespace rnglab-v020-nuclear-b2-calibration), so the calibration is independent of
the samplers under test; the exact laws themselves (state, two-time and interval laws) come from the
module's analytic references. For each test type the null rejection rate is estimated at the nominal
levels 1e-2, 1e-3, 1e-4 and at the Holm level alpha / (NC1 family size), with exact (Clopper-Pearson)
95% intervals. Writes calibration.json; it is design evidence produced before the freeze, and the frozen
run does not read it.

Tests calibrated (names as in protocol.json):
  ET-VAR, ET-VARSUM, ET-VARSUM2   scaled chi-square test of a sample variance (R = 400 replicates of 2,000 nuclei)
  CT-VAR, CT-VARSUM               F-test of the ratio of two sample variances, degrees of freedom from exact moments
  CT-VAR z-difference (not used)  the alternative that was tried first and rejected: z-test of s_a^2 - s_b^2
  BT-PATH                         summed fixed-cell Pearson statistics along the population chain of S2 (1e5 nuclei)
  BT-STATE                        chi-square with the pooling rule, S2 at t = 0.5 and t = 12 (1e5 nuclei)
  CT-HOMO                         2 x K homogeneity chi-square of two pooled ensembles, t = 0.5 and t = 12
"""
from __future__ import annotations

import argparse
import datetime
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

from rnglab.nuclear import DecayScheme, state_probabilities
from rnglab.nuclear.aggregate import AggregatePlan
from rnglab.nuclear.checks import multinomial_check
from rnglab.nuclear.comparison import (
    fixed_cells_chisquare, homogeneity_check, null_difference_check, planned_path_cells,
    scaled_chisquare_variance_check, variance_ratio_f_check,
)
from rnglab.nuclear.manifest import code_identity, environment
from rnglab.nuclear.moments import state_count_moments, state_sum_moments, two_time_law, two_time_sum_moments
from rnglab.seeds import namespace_seed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
NAMESPACE = "rnglab-v020-nuclear-b2-calibration"


def scaled_chi2_p(s2: np.ndarray, sigma2: float, var_s2: float) -> np.ndarray:
    nu = 2.0 * sigma2 ** 2 / var_s2
    statistic = nu * s2 / sigma2
    return np.minimum(1.0, 2.0 * np.minimum(stats.chi2.cdf(statistic, nu), stats.chi2.sf(statistic, nu)))


def difference_p(a: np.ndarray, b: np.ndarray, var_each: float) -> np.ndarray:
    return 2.0 * stats.norm.sf(np.abs(a - b) / np.sqrt(2.0 * var_each))


def f_ratio_p(s2_a: np.ndarray, s2_b: np.ndarray, sigma2: float, var_s2: float) -> np.ndarray:
    nu = 2.0 * sigma2 ** 2 / var_s2
    ratio = s2_a / s2_b
    return np.minimum(1.0, 2.0 * np.minimum(stats.f.cdf(ratio, nu, nu), stats.f.sf(ratio, nu, nu)))


def variance(x: np.ndarray) -> np.ndarray:
    return x.var(axis=1, ddof=1)


class Tally:
    def __init__(self, levels):
        self.levels = levels
        self.rejections: dict[str, np.ndarray] = {}
        self.trials: dict[str, int] = {}

    def add(self, name: str, p_values: np.ndarray) -> None:
        counts = np.array([int(np.sum(p_values < level)) for level in self.levels])
        self.rejections[name] = self.rejections.get(name, 0) + counts
        self.trials[name] = self.trials.get(name, 0) + int(p_values.size)

    def summary(self, prefix: str) -> dict:
        names = [n for n in self.rejections if n.startswith(prefix)]
        total = sum(self.rejections[n] for n in names)
        trials = sum(self.trials[n] for n in names)
        rows = []
        for level, k in zip(self.levels, total):
            interval = stats.binomtest(int(k), trials).proportion_ci(confidence_level=0.95, method="exact")
            rows.append({"nominal_level": level, "rejections": int(k), "trials": trials, "rate": int(k) / trials,
                         "rate_over_nominal": int(k) / trials / level, "ci95": [interval.low, interval.high],
                         "ci95_over_nominal": [interval.low / level, interval.high / level]})
        worst = max(names, key=lambda n: self.rejections[n][-1] / self.trials[n])
        return {"instances": len(names), "pooled": rows,
                "worst_instance_at_last_level": {"name": worst, "rejections": int(self.rejections[worst][-1]),
                                                 "trials": self.trials[worst]}}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scale", type=float, default=1.0, help="fraction of the null replications (development)")
    parser.add_argument("--out", type=Path, default=HERE / "calibration.json")
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"{args.out} exists; results are never overwritten")
    protocol = json.loads((HERE / "protocol.json").read_text())
    design = json.loads((HERE / "design.json").read_text())
    s2 = DecayScheme.from_dict(protocol["schemes"]["S2"])
    names = [n.name for n in s2.nuclides]
    index = {n: i for i, n in enumerate(names)}
    n_rep, reps = design["ens"]["per_replicate"], protocol["runs"]["ENS-B2"]["replicates"]
    holm_level = design["positive_control_level"]
    levels = [1e-2, 1e-3, 1e-4, holm_level]
    tally = Tally(levels)
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    clock = time.perf_counter()
    timings = {}

    # 1. Single-time ensembles at t = 2 (two independent ensembles) and t = 8.
    t0 = time.perf_counter()
    p2 = state_probabilities(s2, "Xa-400", 2.0)[0]
    p8 = state_probabilities(s2, "Xa-400", 8.0)[0]
    eligible2, eligible8 = design["ens"]["eligible"]["2.0"], design["ens"]["eligible"]["8.0"]
    pairs2 = [tuple(pair) for pair in design["ens"]["pairs_t2"]]
    var2 = {j: state_count_moments(s2, "Xa-400", n_rep, 2.0, index[j], index[j]) for j in eligible2}
    var8 = {j: state_count_moments(s2, "Xa-400", n_rep, 8.0, index[j], index[j]) for j in eligible8}
    sum2 = {pair: state_sum_moments(s2, "Xa-400", n_rep, 2.0, index[pair[0]], index[pair[1]]) for pair in pairs2}
    rng = np.random.default_rng(namespace_seed(NAMESPACE, "single-time"))
    total, batch, checked = int(100000 * args.scale), 1000, False
    for start in range(0, total, batch):
        size = min(batch, total - start)
        a = rng.multinomial(n_rep, p2, size=(size, reps)).astype(np.float64)
        b = rng.multinomial(n_rep, p2, size=(size, reps)).astype(np.float64)
        c = rng.multinomial(n_rep, p8, size=(size, reps)).astype(np.float64)
        for j in eligible2:
            m = var2[j]
            s_a, s_b = variance(a[:, :, index[j]]), variance(b[:, :, index[j]])
            p = scaled_chi2_p(s_a, m.sigma_ab, m.sample_covariance_variance(reps))
            tally.add(f"ET-VAR {j} t=2.0", p)
            p_ratio = f_ratio_p(s_a, s_b, m.sigma_ab, m.sample_covariance_variance(reps))
            tally.add(f"CT-VAR {j} t=2.0", p_ratio)
            tally.add(f"ZDIFF-VAR {j} t=2.0", difference_p(s_a, s_b, m.sample_covariance_variance(reps)))
            if not checked:  # the vectorised formulas must equal the functions the study calls
                want = scaled_chisquare_variance_check(s_a[0], m.sigma_ab, m.sample_covariance_variance(reps))["p_value"]
                ratio = variance_ratio_f_check(s_a[0], s_b[0], m.sigma_ab, m.sample_covariance_variance(reps))["p_value"]
                diff = null_difference_check(s_a[0], s_b[0], m.sample_covariance_variance(reps))["p_value"]
                assert abs(p[0] - want) <= 1e-12 * max(want, 1e-300)
                assert abs(p_ratio[0] - ratio) <= 1e-12 * max(ratio, 1e-300)
                assert abs(difference_p(s_a[:1], s_b[:1], m.sample_covariance_variance(reps))[0] - diff) <= 1e-12 * diff
        checked = True
        for j in eligible8:
            m = var8[j]
            tally.add(f"ET-VAR {j} t=8.0", scaled_chi2_p(variance(c[:, :, index[j]]), m.sigma_ab,
                                                         m.sample_covariance_variance(reps)))
        for pair in pairs2:
            ja, jb = index[pair[0]], index[pair[1]]
            m = sum2[pair]
            sum_a, sum_b = variance(a[:, :, ja] + a[:, :, jb]), variance(b[:, :, ja] + b[:, :, jb])
            tally.add(f"ET-VARSUM {pair[0]}+{pair[1]}", scaled_chi2_p(sum_a, m.sigma_ab,
                                                                     m.sample_covariance_variance(reps)))
            tally.add(f"CT-VARSUM {pair[0]}+{pair[1]}", f_ratio_p(sum_a, sum_b, m.sigma_ab,
                                                                 m.sample_covariance_variance(reps)))
    timings["single_time"] = time.perf_counter() - t0

    # 2. Two-time sums between t = 2 and t = 8.
    t0 = time.perf_counter()
    joint = two_time_law(s2, "Xa-400", 2.0, 8.0)
    pairs28 = [tuple(pair) for pair in design["ens"]["pairs_two_time"]]
    sum28 = {pair: two_time_sum_moments(s2, "Xa-400", n_rep, 2.0, 8.0, index[pair[0]], index[pair[1]]) for pair in pairs28}
    rng = np.random.default_rng(namespace_seed(NAMESPACE, "two-time"))
    total, batch = int(60000 * args.scale), 250
    size9 = len(names)
    for start in range(0, total, batch):
        size = min(batch, total - start)
        draws = rng.multinomial(n_rep, joint.ravel() / joint.sum(), size=(size, reps)).reshape(size, reps, size9, size9)
        at2 = draws.sum(axis=3).astype(np.float64)
        at8 = draws.sum(axis=2).astype(np.float64)
        for pair in pairs28:
            m = sum28[pair]
            tally.add(f"ET-VARSUM2 {pair[0]}@2+{pair[1]}@8",
                      scaled_chi2_p(variance(at2[:, :, index[pair[0]]] + at8[:, :, index[pair[1]]]), m.sigma_ab,
                                    m.sample_covariance_variance(reps)))
    timings["two_time"] = time.perf_counter() - t0

    # 3. BT-PATH: the population chain of S2 with numpy multinomials, summed fixed-cell statistics.
    t0 = time.perf_counter()
    case = protocol["cases"]["S2"]
    n0, times = case["initial"]["Xa-400"], case["times"]
    plan = AggregatePlan(s2, times, n0)
    expected = [state_probabilities(s2, "Xa-400", t)[0] * n0 for t in times[:-1]]
    cells = planned_path_cells(plan, expected, min_expected=design["path_tests"]["min_expected"])
    rng = np.random.default_rng(namespace_seed(NAMESPACE, "path-chain"))
    total, batch, zero_start_counts, checked = int(1000000 * args.scale), 10000, 0, False
    for start in range(0, total, batch):
        size = min(batch, total - start)
        populations = np.zeros((size, len(names)), dtype=np.int64)
        populations[:, index["Xa-400"]] = n0
        statistic = {i: np.zeros(size) for i in cells}
        for k, delta in enumerate(plan.deltas.tolist()):
            new = np.zeros_like(populations)
            for i, categories in enumerate(plan.categories):
                if len(categories) == 1:
                    new[:, categories[0].end] += populations[:, i]
                    continue
                q = plan.interval_plans[delta][i].law.as_float()
                counts = rng.multinomial(populations[:, i], q / q.sum())
                for c, category in enumerate(categories):
                    new[:, category.end] += counts[:, c]
                if i in cells and k in cells[i]:
                    n_start = populations[:, i].astype(np.float64)
                    zero_start_counts += int(np.sum(n_start == 0))
                    part = np.zeros(size)
                    for cell in cells[i][k]:
                        expected_cell = n_start * q[cell].sum()
                        with np.errstate(divide="ignore", invalid="ignore"):
                            part += (counts[:, cell].sum(axis=1) - expected_cell) ** 2 / expected_cell
                    statistic[i] += part
                    if not checked and i == index["Xa-400"]:
                        want = fixed_cells_chisquare(counts[0], q, cells[i][k])["statistic"]
                        assert abs(part[0] - want) <= 1e-9 * max(want, 1.0)
                        checked = True
            populations = new
        for i, value in statistic.items():
            dof = sum(len(c) - 1 for c in cells[i].values())
            tally.add(f"BT-PATH {names[i]}", stats.chi2.sf(value, dof))
    timings["path_chain"] = time.perf_counter() - t0

    # 4. BT-STATE and CT-HOMO through the functions the study calls.
    t0 = time.perf_counter()
    rng = np.random.default_rng(namespace_seed(NAMESPACE, "state-and-homogeneity"))
    pooled = reps * n_rep
    for t in (0.5, 12.0):
        p = state_probabilities(s2, "Xa-400", t)[0]
        draws = rng.multinomial(n0, p, size=int(400000 * args.scale))
        tally.add(f"BT-STATE t={t}", np.array([multinomial_check(row, p)["p_value"] for row in draws]))
        a = rng.multinomial(pooled, p, size=int(200000 * args.scale))
        b = rng.multinomial(pooled, p, size=int(200000 * args.scale))
        tally.add(f"CT-HOMO t={t}", np.array([homogeneity_check(x, y)["p_value"] for x, y in zip(a, b)]))
    timings["state_and_homogeneity"] = time.perf_counter() - t0

    result = {"generated_by": "research/v020-nuclear-b2-verify/calibration.py", "started_utc": started,
              "code": code_identity(REPO), "environment": environment(), "namespace": NAMESPACE,
              "scale": args.scale, "levels": levels, "holm_level": holm_level,
              "sampler": "numpy Generator.multinomial on the exact laws; no RNGLab engine",
              "zero_start_counts_in_planned_path_intervals": zero_start_counts,
              "tests": {label: tally.summary(prefix) for prefix, label in
                        (("ET-VAR ", "ET-VAR"), ("ET-VARSUM ", "ET-VARSUM"), ("ET-VARSUM2 ", "ET-VARSUM2"),
                         ("CT-VAR ", "CT-VAR"), ("CT-VARSUM ", "CT-VARSUM"),
                         ("ZDIFF-VAR ", "CT-VAR z-difference (not used)"), ("BT-PATH ", "BT-PATH"),
                         ("BT-STATE ", "BT-STATE"), ("CT-HOMO ", "CT-HOMO"))},
              "timings_seconds": timings, "elapsed_seconds": time.perf_counter() - clock}
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    for prefix, summary in result["tests"].items():
        rows = summary["pooled"]
        print(f"{prefix:<31} trials {rows[0]['trials']:>9}: rate/nominal "
              + "  ".join(f"{r['nominal_level']:.1e}: {r['rate_over_nominal']:.2f} "
                          f"[{r['ci95_over_nominal'][0]:.2f}, {r['ci95_over_nominal'][1]:.2f}]" for r in rows))
    print(f"wrote {args.out} in {result['elapsed_seconds']:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
