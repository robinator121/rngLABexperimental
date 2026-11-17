"""Null calibration of the complete B2 verification batteries (protocol.json: calibration, criterion C1).

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/calibration.py --development [--batteries 2000]
    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/calibration.py   # the record: once, after the freeze

A null battery draws every run of NC1 and NC2 from the exact laws with numpy's own Generator: the S1-S5
single runs, the ENS-B2 and ENS-A ensembles (400 replicates of 2,000 nuclei each) and both BIG runs.
Each interval's path-category counts are Multinomial(start count, decimal interval law) for every
nuclide, decays follow from the paths, photons are Binomial(decays, yield). No RNGLab sampler is
involved, so the calibration does not depend on the samplers under test. Every planned test is then
evaluated through the functions run_verification.py uses (battery, path_tests, ensemble_tests,
comparison_tests, big_tests), so the calibration measures the study's own tests.

Batteries are independent (battery b uses default_rng(namespace_seed(namespace, 'battery', b))), so:
- per family, the number of batteries in which Holm rejects at least one test (or a planned test cannot
  be computed: R1 would FAIL) is binomial, and its Clopper-Pearson interval is valid. Criterion C1: the
  95% upper bound is at most 0.02 for NC1 and for NC2;
- per test, rejections at 1e-2, 1e-3 and the Holm level have valid per-test Clopper-Pearson intervals;
- per test type, the tests of one battery are dependent, so no pooled binomial interval is used. The
  number of batteries with at least one rejection of the type is binomial, with a valid Clopper-Pearson
  interval; divided by instances x level it is the type's share of the union bound on the family-wise
  error (1 at nominal calibration if the type's tests rarely reject together). The mean number of
  rejections per battery is reported beside it, with only an approximate (normal) interval.

The record run refuses unless protocol.sha256 matches the frozen files, the Git tree is clean and
results/calibration.json does not exist. --development uses the development namespace and writes
calibration_development.json next to this script (design evidence before the freeze, never results).
"""
from __future__ import annotations

import argparse
import datetime
import importlib.util
import json
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy import stats

from rnglab.nuclear import DecayScheme
from rnglab.nuclear.aggregate import AggregatePlan
from rnglab.nuclear.checks import holm
from rnglab.nuclear.decay import Mode
from rnglab.nuclear.manifest import code_identity, environment
from rnglab.seeds import namespace_seed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
LEVELS = (1e-2, 1e-3)  # the family's Holm level alpha / m is added per family
TYPE_PREFIXES = ("BT-", "ET-", "CT-", "BIG-")


def load_verification():
    spec = importlib.util.spec_from_file_location("b2_verify_run_verification", HERE / "run_verification.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_type(name: str) -> str:
    return next(token for token in name.split() if token.startswith(TYPE_PREFIXES))


class NullModel:
    """Exact-law sampler of one scheme on one grid with numpy (no RNGLab sampler)."""

    def __init__(self, scheme: DecayScheme, times, total: int):
        self.plan = AggregatePlan(scheme, times, total)
        self.scheme = scheme
        width = max(scheme.max_channels, 1)
        self.yields = np.zeros((len(scheme), width))
        self.neutron = np.zeros((len(scheme), width), dtype=bool)
        for j, nuclide in enumerate(scheme.nuclides):
            for k, channel in enumerate(nuclide.channels):
                self.yields[j, k] = channel.photon_yield
                self.neutron[j, k] = channel.mode is Mode.BETA_DELAYED_NEUTRON
        self.laws = {}
        for delta, plans in self.plan.interval_plans.items():
            for i, source in plans.items():
                q = source.law.as_float()
                self.laws[(delta, i)] = q / q.sum()

    def draw(self, initial: dict, rng: np.random.Generator, replicates: int) -> SimpleNamespace:
        plan, scheme = self.plan, self.scheme
        size, intervals = len(scheme), plan.deltas.size
        pops = np.zeros((replicates, intervals + 1, size), dtype=np.int64)
        for name, count in initial.items():
            pops[:, 0, scheme.index[name]] = count
        decays = np.zeros((replicates, intervals, size, self.yields.shape[1]), dtype=np.int64)
        paths = np.zeros((replicates, intervals, plan.n_categories), dtype=np.int64)
        for k, delta in enumerate(plan.deltas.tolist()):
            for i, categories in enumerate(plan.categories):
                start = pops[:, k, i]
                first = int(plan.first_category[i])
                if len(categories) == 1:
                    pops[:, k + 1, categories[0].end] += start
                    paths[:, k, first] = start
                    continue
                counts = rng.multinomial(start, self.laws[(delta, i)])
                paths[:, k, first:first + len(categories)] = counts
                for c, category in enumerate(categories):
                    pops[:, k + 1, category.end] += counts[:, c]
                    for j, channel in category.steps:
                        decays[:, k, j, channel] += counts[:, c]
        photons = rng.binomial(decays, np.broadcast_to(self.yields, decays.shape))
        neutrons = (decays * self.neutron).sum(axis=(2, 3))
        return SimpleNamespace(populations=pops, decays=decays, photons=photons, neutrons=neutrons, paths=paths)

    def single(self, initial: dict, rng: np.random.Generator) -> SimpleNamespace:
        d = self.draw(initial, rng, 1)
        return SimpleNamespace(plan=self.plan, initial=dict(initial), populations=d.populations[0],
                               decays=d.decays[0], photons=d.photons[0], neutrons=d.neutrons[0],
                               path_counts=d.paths[0])


class Context:
    def __init__(self):
        self.verify = load_verification()
        self.protocol = json.loads((HERE / "protocol.json").read_text())
        self.design = json.loads((HERE / "design.json").read_text())
        p = self.protocol
        self.schemes = {k: DecayScheme.from_dict(v) for k, v in p["schemes"].items()}
        self.cases = p["cases"]
        self.models = {k: NullModel(self.schemes[k], self.cases[k]["times"], sum(self.cases[k]["initial"].values()))
                       for k in self.schemes}
        self.s2 = self.schemes["S2"]
        self.times2 = self.cases["S2"]["times"]
        self.reps = p["runs"]["ENS-B2"]["replicates"]
        self.n_ens = p["runs"]["ENS-B2"]["initial"]["Xa-400"]
        self.ensemble = NullModel(self.s2, self.times2, self.n_ens)
        self.big = {size: NullModel(self.s2, self.times2, size) for size in p["runs"]["BIG"]["initial_sizes"]}
        self.planned = {f: set(self.design["planned_tests"][f]) for f in ("NC1", "NC2")}
        self.alpha = p["families"]["alpha"]

    def battery(self, rng: np.random.Generator) -> dict[str, dict]:
        v, found = self.verify, {"NC1": {}, "NC2": {}}
        for family in ("NC1", "NC2"):
            for key in self.protocol["runs"][family]["schemes"]:
                run = self.models[key].single(self.cases[key]["initial"], rng)
                found[family].update(v.battery(key, run, self.schemes[key], self.cases[key], family))
                if key == "S2":
                    found[family].update(v.path_tests(run, family, self.cases[key],
                                                      self.design["path_tests"]["min_expected"])[0])
        start = {"Xa-400": self.n_ens}
        b = self.ensemble.draw(start, rng, self.reps)
        a = self.ensemble.draw(start, rng, self.reps)
        found["NC1"].update(v.ensemble_tests(self.s2, self.design, self.reps, self.n_ens, self.times2, b.populations))
        found["NC1"].update(v.comparison_tests(self.s2, self.design, self.reps, self.n_ens, self.times2,
                                               a.populations, b.populations, a.decays, b.decays, a.photons, b.photons))
        for size, model in self.big.items():
            run = model.single({"Xa-400": size}, rng)
            found["NC1"].update(v.big_tests(self.s2, self.times2, size, run.populations, run.decays))
        return found


_CONTEXT: Context | None = None


def _init_worker() -> None:
    global _CONTEXT
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    _CONTEXT = Context()


def run_chunk(namespace: str, batteries: list[int]) -> dict:
    ctx = _CONTEXT
    out = {"batteries": 0, "families": {}}
    for family in ("NC1", "NC2"):
        holm_level = ctx.alpha / len(ctx.planned[family])
        levels = (*LEVELS, holm_level)
        types = sorted({test_type(n) for n in ctx.planned[family]})
        out["families"][family] = {
            "levels": levels, "failed_batteries": 0, "batteries_with_errors": 0, "errors": {},
            "per_test": {n: [0] * len(levels) for n in sorted(ctx.planned[family])},
            "type_sum": {t: [0] * len(levels) for t in types}, "type_sum_squares": {t: [0] * len(levels) for t in types},
            "type_batteries": {t: [0] * len(levels) for t in types}}
    for b in batteries:
        rng = np.random.default_rng(namespace_seed(namespace, "battery", b))
        found = ctx.battery(rng)
        out["batteries"] += 1
        for family, results in found.items():
            agg = out["families"][family]
            members = {n: r["p_value"] for n, r in results.items() if "p_value" in r}
            errors = sorted(n for n, r in results.items() if "error" in r)
            if set(members) | set(errors) != ctx.planned[family]:
                raise RuntimeError(f"battery {b}: the tests of {family} differ from design.json "
                                   f"(missing {sorted(ctx.planned[family] - set(members) - set(errors))[:5]}, "
                                   f"extra {sorted(set(members) - ctx.planned[family])[:5]})")
            for n in errors:
                agg["errors"][n] = agg["errors"].get(n, 0) + 1
            names = sorted(members)
            rejected_any = bool(holm([members[n] for n in names], ctx.alpha)["rejected"]) if names else False
            agg["failed_batteries"] += int(rejected_any or bool(errors))
            agg["batteries_with_errors"] += int(bool(errors))
            per_type = {t: [0] * len(agg["levels"]) for t in agg["type_sum"]}
            for n in names:
                for li, level in enumerate(agg["levels"]):
                    if members[n] < level:
                        agg["per_test"][n][li] += 1
                        per_type[test_type(n)][li] += 1
            for t, counts in per_type.items():
                for li, c in enumerate(counts):
                    agg["type_sum"][t][li] += c
                    agg["type_sum_squares"][t][li] += c * c
                    agg["type_batteries"][t][li] += int(c > 0)
    return out


def merge(total: dict | None, part: dict) -> dict:
    if total is None:
        return part
    total["batteries"] += part["batteries"]
    for family, agg in part["families"].items():
        t = total["families"][family]
        t["failed_batteries"] += agg["failed_batteries"]
        t["batteries_with_errors"] += agg["batteries_with_errors"]
        for n, c in agg["errors"].items():
            t["errors"][n] = t["errors"].get(n, 0) + c
        for key in ("per_test", "type_sum", "type_sum_squares", "type_batteries"):
            for n, counts in agg[key].items():
                t[key][n] = [x + y for x, y in zip(t[key][n], counts)]
    return total


def clopper_pearson(k: int, n: int) -> list[float]:
    interval = stats.binomtest(int(k), int(n)).proportion_ci(confidence_level=0.95, method="exact")
    return [float(interval.low), float(interval.high)]


def summarise(raw: dict, criterion_limit: float) -> dict:
    m = raw["batteries"]
    out = {}
    for family, agg in raw["families"].items():
        levels = agg["levels"]
        ci = clopper_pearson(agg["failed_batteries"], m)
        instances = {t: 0 for t in agg["type_sum"]}
        for n in agg["per_test"]:
            instances[test_type(n)] += 1
        types = {}
        for t in agg["type_sum"]:
            rows = []
            for li, level in enumerate(levels):
                mean = agg["type_sum"][t][li] / m
                var = max(agg["type_sum_squares"][t][li] / m - mean ** 2, 0.0)
                se = (var / m) ** 0.5
                scale = instances[t] * level
                hit = agg["type_batteries"][t][li]
                low, high = clopper_pearson(hit, m)
                rows.append({"level": level, "batteries_with_rejection": hit, "batteries": m,
                             "battery_rate_over_union_bound": hit / m / scale,
                             "battery_rate_ci95_over_union_bound": [low / scale, high / scale],
                             "rejections": agg["type_sum"][t][li], "trials": instances[t] * m,
                             "mean_rejections_over_nominal": mean / scale,
                             "mean_rejections_approx_ci95_over_nominal": [(mean - 1.96 * se) / scale,
                                                                          (mean + 1.96 * se) / scale]})
            types[t] = {"instances": instances[t], "levels": rows}
        per_test = {n: {"rejections": counts,
                        "ci95": [clopper_pearson(c, m) for c in counts]} for n, counts in agg["per_test"].items()}
        worst = max(agg["per_test"], key=lambda n: agg["per_test"][n][-1])
        out[family] = {"tests": len(agg["per_test"]), "levels": list(levels), "batteries": m,
                       "failed_batteries": agg["failed_batteries"], "false_fail_rate": agg["failed_batteries"] / m,
                       "false_fail_ci95": ci, "C1_holds": ci[1] <= criterion_limit,
                       "batteries_with_errors": agg["batteries_with_errors"], "errors": agg["errors"],
                       "types": types, "per_test": per_test,
                       "most_rejections_at_holm_level": {"test": worst, "rejections": agg["per_test"][worst][-1]}}
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--development", action="store_true", help="development namespace; design evidence only")
    parser.add_argument("--batteries", type=int, default=None, help="development only (the record uses protocol.json)")
    parser.add_argument("--workers", type=int, default=min(16, os.cpu_count() or 1))
    parser.add_argument("--chunk", type=int, default=10)
    parser.add_argument("--out", type=Path, default=None, help="development only: another output file")
    args = parser.parse_args()
    if args.out is not None and not args.development:
        raise SystemExit("--out is for development runs only")
    protocol = json.loads((HERE / "protocol.json").read_text())
    spec = protocol["calibration"]
    verify = load_verification()
    if args.development:
        namespace, batteries = spec["development_namespace"], args.batteries or 2000
        out = args.out or HERE / "calibration_development.json"
        frozen = verify.frozen_hashes(smoke=True)
    else:
        if args.batteries is not None:
            raise SystemExit("the record run uses the number of batteries in protocol.json")
        frozen = verify.frozen_hashes()
        namespace, batteries = spec["namespace"], spec["batteries"]
        if code_identity(REPO)["git_dirty"] is not False:
            raise SystemExit("the record run needs a clean, committed Git tree")
        out = HERE / "results" / "calibration.json"
    if out.exists():
        raise SystemExit(f"{out} exists; results are never overwritten")
    started = datetime.datetime.now(datetime.timezone.utc).isoformat()
    clock = time.perf_counter()
    chunks = [list(range(s, min(s + args.chunk, batteries))) for s in range(0, batteries, args.chunk)]
    raw = None
    with ProcessPoolExecutor(max_workers=args.workers, initializer=_init_worker) as pool:
        for done, part in enumerate(pool.map(run_chunk, [namespace] * len(chunks), chunks), start=1):
            raw = merge(raw, part)
            if done % max(1, len(chunks) // 20) == 0:
                print(f"{raw['batteries']} / {batteries} batteries, {time.perf_counter() - clock:.0f} s", flush=True)
    summary = summarise(raw, 0.02)
    result = {"generated_by": "research/v020-nuclear-b2-verify/calibration.py", "development": args.development,
              "started_utc": started, "namespace": namespace, "batteries": batteries, "protocol_files": frozen,
              "code": code_identity(REPO), "environment": environment(), "criterion": "C1: false_fail_ci95[1] <= 0.02",
              "sampler": "numpy Generator.multinomial / binomial on the exact laws; no RNGLab sampler",
              "families": summary, "C1": all(f["C1_holds"] for f in summary.values()),
              "elapsed_seconds": time.perf_counter() - clock}
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=1) + "\n")
    for family, f in summary.items():
        print(f"{family}: {f['failed_batteries']} / {batteries} batteries fail under the null "
              f"(rate {f['false_fail_rate']:.4f}, 95% CI [{f['false_fail_ci95'][0]:.4f}, {f['false_fail_ci95'][1]:.4f}]); "
              f"C1 {'holds' if f['C1_holds'] else 'FAILS'}; batteries with errors {f['batteries_with_errors']}")
        for t, row in sorted(f["types"].items()):
            r = row["levels"][-1]
            print(f"  {t:<12} x{row['instances']:<3} Holm level: {r['batteries_with_rejection']:>4} batteries, "
                  f"rate / union bound {r['battery_rate_over_union_bound']:.2f} "
                  f"[{r['battery_rate_ci95_over_union_bound'][0]:.2f}, {r['battery_rate_ci95_over_union_bound'][1]:.2f}]"
                  f"   1e-2 mean/nominal {row['levels'][0]['mean_rejections_over_nominal']:.2f}")
    print(f"wrote {out} in {result['elapsed_seconds']:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
