"""B2 verification against analytic laws and the frozen Layer A engine: the frozen study in protocol.json.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/run_verification.py [--smoke]

The script refuses to draw if a frozen file differs from protocol.sha256. The
main run also refuses on a dirty Git tree or when results/verification.json
exists. --smoke uses the smoke namespace and reduced sizes, writes to smoke/
(not in Git) and is excluded from every result.
"""
from __future__ import annotations

import argparse
import copy
import datetime
import json
import math
import sys
import time
from decimal import Context, Decimal, localcontext
from pathlib import Path

import numpy as np

from rnglab.nuclear import (
    DecayScheme, NumpyPCG64Source, WordStream, decay_probabilities, simulate_decays, state_probabilities,
)
from rnglab.nuclear.aggregate import aggregate_manifest, simulate_aggregate, simulate_aggregate_reference
from rnglab.nuclear.binomial import binomial_inverse
from rnglab.nuclear.checks import binomial_check, decay_invariants, holm, multinomial_check
from rnglab.nuclear.comparison import (
    fixed_cells_chisquare, homogeneity_check, planned_path_cells, sample_covariance,
    scaled_chisquare_variance_check, summed_chisquare, variance_ratio_f_check, z_check,
)
from rnglab.nuclear.decay import Mode
from rnglab.nuclear.manifest import code_identity, environment, file_sha256
from rnglab.nuclear.moments import state_count_moments, state_sum_moments, two_time_sum_moments
from rnglab.nuclear.paths import bounds_hold, interval_law, tree_probabilities
from rnglab.seeds import namespace_seed
from rnglab14.sources import Xoshiro256PlusPlusSource

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FROZEN = ("protocol.json", "design.json", "design.py", "run_verification.py")


def frozen_hashes() -> dict[str, str]:
    if not (HERE / "protocol.sha256").exists():
        raise SystemExit("protocol.sha256 is missing: the protocol is not frozen; refusing to draw")
    listed = {}
    for line in (HERE / "protocol.sha256").read_text().splitlines():
        if line.strip():
            digest, name = line.split()
            listed[name] = digest
    if set(listed) != set(FROZEN):
        raise SystemExit("protocol.sha256 must list exactly the frozen files")
    for name in FROZEN:
        if file_sha256(HERE / name) != listed[name]:
            raise SystemExit(f"{name} differs from protocol.sha256; refusing to draw")
    return listed


def make_source(generator: str, seed: int, label: str):
    if generator == "numpy-pcg64":
        return NumpyPCG64Source(seed, label=label)
    if generator == "xoshiro256pp":
        return Xoshiro256PlusPlusSource(seed=seed)
    raise ValueError(generator)


def groups_of(scheme: DecayScheme, initial: dict) -> list[tuple[str, list[int]]]:
    out = []
    for start in initial:
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


def bookkeeping(run) -> dict[str, bool]:
    scheme = run.plan.scheme
    pops, decays, photons = run.populations, run.decays, run.photons
    total = sum(run.initial.values())
    inflow = np.zeros_like(pops[1:])
    yields = np.zeros(decays.shape[1:])
    neutron = np.zeros(decays.shape[1:], dtype=bool)
    for j, nuclide in enumerate(scheme.nuclides):
        for k, channel in enumerate(nuclide.channels):
            inflow[:, scheme.index[channel.daughter]] += decays[:, j, k]
            yields[j, k] = channel.photon_yield
            neutron[j, k] = channel.mode is Mode.BETA_DELAYED_NEUTRON
    outflow = decays.sum(axis=2)
    return {"conserved_and_nonnegative": bool(np.all(pops >= 0) and set(pops.sum(axis=1).tolist()) == {total}),
            "flow_balance": bool(np.array_equal(np.diff(pops, axis=0), inflow - outflow)),
            "outflow_bounded": bool(np.all(outflow <= pops[:-1] + inflow)),
            "photons_bounded": bool(np.all((photons >= 0) & (photons <= decays))),
            "photon_yield_extremes": bool(np.all(photons[:, yields == 0.0] == 0)
                                          and np.array_equal(photons[:, yields == 1.0], decays[:, yields == 1.0])),
            "neutrons_equal_beta_n_decays": bool(np.array_equal(run.neutrons, decays[:, neutron].sum(axis=1)))}


def b2_consumption(run, stream) -> dict:
    words = run.plan.deltas.size * run.plan.words_per_interval
    return {"words_expected": words, "words_read": run.words_read,
            "ok": run.words_read == words == stream.words_read and stream.source.consumed_bytes == 8 * words}


def chi2(counts, probabilities) -> dict:
    """multinomial_check, recording the pre-declared 'not applicable' case (fewer than two cells after pooling)."""
    try:
        return multinomial_check(counts, probabilities)
    except ValueError as error:
        if "fewer than two cells" in str(error):
            return {"not_applicable": str(error)}
        raise


def battery(key: str, run, ref: DecayScheme, case: dict, prefix: str) -> dict:
    out = {}
    times = case["times"]
    pops, decays, photons = run.populations, run.decays, run.photons
    for start, members in groups_of(ref, case["initial"]):
        n = case["initial"][start]
        p = state_probabilities(ref, start, times)
        for it, t in enumerate(times[1:], start=1):
            out[f"{prefix} {key} BT-STATE {start} t={t!r}"] = chi2(pops[it, members], p[it, members])
        q_end = decay_probabilities(ref, start, times[-1])[0]
        for j in members:
            if not ref.nuclides[j].stable:
                out[f"{prefix} {key} BT-DECAY {ref.nuclides[j].name}"] = binomial_check(
                    decays[:, j, :].sum(), n, q_end[j])
        i = ref.index[start]
        q = decay_probabilities(ref, start, times)[:, i]
        per_interval = decays[:, i, :].sum(axis=1)
        out[f"{prefix} {key} BT-BINS {start}"] = chi2(
            np.append(per_interval, n - per_interval.sum()), np.append(np.diff(q), 1.0 - q[-1]))
        if key == "S4":
            channels = ref.nuclides[i].channels
            k = next(c for c, ch in enumerate(channels) if ch.mode is Mode.BETA_DELAYED_NEUTRON)
            pn = ref.nuclides[i].channel_probabilities()[k]
            neutrons = decays[:, i, k]
            probabilities = np.append(pn * np.diff(q), 1.0 - pn * q[-1])
            out[f"{prefix} {key} BT-NEUTRON {start}"] = chi2(np.append(neutrons, n - neutrons.sum()), probabilities)
    for j, nuclide in enumerate(ref.nuclides):
        total_j = int(decays[:, j, :].sum())
        if len(nuclide.channels) >= 2:
            out[f"{prefix} {key} BT-CHANNEL {nuclide.name}"] = binomial_check(
                decays[:, j, 0].sum(), total_j, nuclide.channel_probabilities()[0])
        for k, channel in enumerate(nuclide.channels):
            if 0.0 < channel.photon_yield < 1.0:
                out[f"{prefix} {key} BT-PHOTON {nuclide.name}/{channel.mode.value}->{channel.daughter}"] = \
                    binomial_check(photons[:, j, k].sum(), decays[:, j, k].sum(), channel.photon_yield)
    return out


def path_tests(run, prefix: str, case: dict, min_expected: float) -> tuple[dict, dict]:
    """One test per source: Pearson statistics of the planned intervals on pre-declared cells, summed.

    Given the start count, the category counts of an interval are Multinomial(start count, q(delta))
    whatever happened before; with cells fixed in advance each interval's statistic is asymptotically
    chi-square with constant degrees of freedom given the past, so the sum is asymptotically
    chi-square(summed degrees of freedom). A statistic that cannot be computed is recorded as an
    error, which makes check D10 fail; it never reads as a passing test.
    """
    plan = run.plan
    scheme = plan.scheme
    (start, n0), = case["initial"].items()
    expected = [state_probabilities(scheme, start, t)[0] * n0 for t in case["times"][:-1]]
    cells = planned_path_cells(plan, expected, min_expected=min_expected)
    out = {}
    for i, per_interval in cells.items():
        first, width = int(plan.first_category[i]), len(plan.categories[i])
        name = f"{prefix} S2 BT-PATH {scheme.nuclides[i].name}"
        try:
            parts = [fixed_cells_chisquare(run.path_counts[k, first:first + width],
                                           plan.interval_plans[float(plan.deltas[k])][i].law.as_float(),
                                           interval_cells)
                     for k, interval_cells in per_interval.items()]
            out[name] = summed_chisquare(parts)
            out[name]["intervals"] = sorted(per_interval)
            out[name]["interval_statistics"] = [part["statistic"] for part in parts]
        except ValueError as error:
            out[name] = {"error": str(error)}
    names = {scheme.nuclides[i].name: {str(k): c for k, c in per.items()} for i, per in cells.items()}
    return out, names


def naive_aggregate(scheme: DecayScheme, initial: dict, times, stream: WordStream) -> np.ndarray:
    """The deliberately wrong one-step aggregation (positive control): at most one decay per interval."""
    pops = np.array([initial.get(n.name, 0) for n in scheme.nuclides], dtype=np.int64)
    out = [pops.copy()]
    for delta in np.diff(times):
        new = np.zeros_like(pops)
        for i, nuclide in enumerate(scheme.nuclides):
            if nuclide.stable:
                new[i] += pops[i]
                continue
            stay = math.exp(-math.log(2.0) / nuclide.half_life * delta)
            probabilities = [stay] + [(1.0 - stay) * b for b in nuclide.channel_probabilities()]
            targets = [i] + [scheme.index[c.daughter] for c in nuclide.channels]
            remaining, tail = int(pops[i]), 1.0
            for k in range(len(probabilities) - 1):
                u = stream.uniforms(1)
                p = probabilities[k] / tail if tail > 0 else 0.0
                p = min(max(p, 0.0), 1.0)
                if p > 0.5:
                    x = remaining - int(binomial_inverse([remaining], [1.0 - p], u)[0])
                else:
                    x = int(binomial_inverse([remaining], [p], u)[0])
                new[targets[k]] += x
                remaining -= x
                tail -= probabilities[k]
            new[targets[-1]] += remaining
        pops = new
        out.append(pops.copy())
    return np.array(out)


def layer_a_observations(run, times):
    times = np.asarray(times)
    scheme = run.scheme
    pops = run.populations(times)
    ev = run.events
    interval = np.searchsorted(times, ev["time"], side="left") - 1
    at_zero = int(np.sum(ev["time"] <= times[0]))
    shape = (times.size - 1, len(scheme), max(scheme.max_channels, 1))
    decays = np.zeros(shape, dtype=np.int64)
    photons = np.zeros(shape, dtype=np.int64)
    keep = interval >= 0
    np.add.at(decays, (interval[keep], ev["parent"][keep], ev["channel"][keep]), 1)
    flag = keep & ev["photon"]
    np.add.at(photons, (interval[flag], ev["parent"][flag], ev["channel"][flag]), 1)
    return pops, decays, photons, at_zero


def law_checks(schemes: dict, cases: dict) -> dict:
    worst_tree = worst_layer_a = 0.0
    worst_sum = Decimal(0)
    brackets, digits, laws, error = True, 0, 0, None
    try:
        for key, scheme in schemes.items():
            for delta in sorted(set(np.diff(cases[key]["times"]).tolist())):
                for i, nuclide in enumerate(scheme.nuclides):
                    law = interval_law(scheme, i, delta)
                    laws += 1
                    digits = max(digits, law.digits_used)
                    q = law.as_float()
                    worst_tree = max(worst_tree, float(np.max(np.abs(q - tree_probabilities(scheme, i, delta)))))
                    ends = np.zeros(len(scheme))
                    for category, value in zip(law.categories, q):
                        ends[category.end] += value
                    worst_layer_a = max(worst_layer_a, float(np.max(np.abs(
                        ends - state_probabilities(scheme, nuclide.name, delta)[0]))))
                    with localcontext(Context(prec=300)):
                        worst_sum = max(worst_sum, abs(sum(law.probabilities, Decimal(0)) - 1))
                    brackets &= bounds_hold(scheme, law)
    except ValueError as caught:
        error = str(caught)
    ok = (error is None and worst_tree <= 1e-12 and worst_layer_a <= 1e-12 and worst_sum <= Decimal("1e-30")
          and brackets)
    return {"laws_checked": laws, "max_abs_decimal_minus_tree": worst_tree,
            "max_abs_end_state_minus_layer_a": worst_layer_a, "max_abs_sum_minus_one": float(worst_sum),
            "brackets_hold": brackets, "max_digits": digits, "error": error, "ok": bool(ok)}


def scaled(initial: dict, scale: float) -> dict:
    return {k: max(1, round(v * scale)) for k, v in initial.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--smoke", action="store_true")
    args = parser.parse_args()
    frozen = frozen_hashes()
    protocol = json.loads((HERE / "protocol.json").read_text())
    design = json.loads((HERE / "design.json").read_text())
    code = code_identity(REPO)
    out_dir = HERE / ("smoke" if args.smoke else "results")
    if not args.smoke:
        if code["git_dirty"] is not False:
            raise SystemExit("the main run needs a clean, committed Git tree")
        if (out_dir / "verification.json").exists():
            raise SystemExit("results/verification.json exists; results are never overwritten")
    (out_dir / "manifests").mkdir(parents=True, exist_ok=True)
    ns = protocol["smoke_namespace"] if args.smoke else protocol["namespace"]
    smoke = protocol["smoke"]
    scale = smoke["scale"] if args.smoke else 1.0
    runs = protocol["runs"]
    experiment = protocol["id"] + ("-smoke" if args.smoke else "")
    schemes = {k: DecayScheme.from_dict(v) for k, v in protocol["schemes"].items()}
    if {k: s.sha256() for k, s in schemes.items()} != design["scheme_sha256"]:
        raise SystemExit("schemes differ from design.json")
    cases = {k: {"initial": scaled(c["initial"], scale), "times": c["times"]} for k, c in protocol["cases"].items()}
    clock = time.perf_counter()
    result = {"protocol": {"id": protocol["id"], "frozen_sha256": frozen}, "namespace": ns, "smoke": args.smoke,
              "code": code, "environment": environment(),
              "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "runs": {},
              "timings_seconds": {}}

    def seed(*parts):
        return namespace_seed(ns, *parts), "/".join([ns, *map(str, parts)])

    def b2(scheme, initial, times, generator, label_parts, **kw):
        s, label = seed(*label_parts)
        stream = WordStream(make_source(generator, s, label))
        t0 = time.perf_counter()
        run = simulate_aggregate(scheme, initial, times, stream, **kw)
        return run, stream, time.perf_counter() - t0, s, label

    def manifest(run_id, run, stream, label):
        (out_dir / "manifests" / f"{run_id}.json").write_text(json.dumps(
            aggregate_manifest(experiment=experiment, run_id=run_id, run=run, stream=stream, code=code,
                               extra={"seed_label": label}), indent=1) + "\n")

    deterministic: dict[str, bool] = {}
    result["laws"] = law_checks(schemes, protocol["cases"])
    deterministic["D8 interval laws"] = result["laws"]["ok"]

    tests = {"NC1": {}, "NC2": {}}
    d1, d2 = True, True
    for family in ("NC1", "NC2"):
        generator = runs[family]["generator"]
        result["runs"][family] = {}
        for key in runs[family]["schemes"]:
            run, stream, seconds, s, label = b2(schemes[key], cases[key]["initial"], cases[key]["times"], generator,
                                                (family, key, "physics"))
            entry = {"seed": s, "seed_label": label, "seconds": seconds, "summary": run.summary(),
                     "stream_sha256": stream.source.stream_digest, "consumption": b2_consumption(run, stream),
                     "bookkeeping": bookkeeping(run),
                     "populations": run.populations.tolist(), "decays": run.decays.tolist(),
                     "photons": run.photons.tolist(), "path_counts": run.path_counts.tolist()}
            found = battery(key, run, schemes[key], cases[key], family)
            if key == "S2":
                path_found, path_cells = path_tests(run, family, cases[key], design["path_tests"]["min_expected"])
                found.update(path_found)
                if not args.smoke and path_cells != design["path_tests"]["cells"]:
                    raise SystemExit("BT-PATH cells differ from design.json")
            entry["tests"] = found
            tests[family].update({k: v for k, v in found.items() if "p_value" in v})
            d1 &= entry["consumption"]["ok"]
            d2 &= all(entry["bookkeeping"].values())
            result["runs"][family][key] = entry
            result["timings_seconds"][f"{family} {key}"] = seconds
            manifest(f"{family}_{key}", run, stream, label)
            if family == "NC1" and key == "S2":
                nc1_s2 = (run, stream, s, label)

    # Ensembles of S2: B2, Layer A, and the naive aggregation (positive control).
    s2, times2 = schemes["S2"], protocol["cases"]["S2"]["times"]
    reps = smoke["ens_replicates"] if args.smoke else runs["ENS-B2"]["replicates"]
    n_ens = runs["ENS-B2"]["initial"]["Xa-400"]
    pop_b, dec_b, pho_b, pop_a, dec_a, pho_a, pop_n = [], [], [], [], [], [], []
    seeds, digests = {"ENS-B2": [], "ENS-A": [], "ENS-NAIVE": []}, {"ENS-B2": [], "ENS-A": [], "ENS-NAIVE": []}
    a_invariants, a_consumption, a_zero = True, True, 0
    t0 = time.perf_counter()
    for r in range(reps):
        run, stream, _, s, _ = b2(s2, {"Xa-400": n_ens}, times2, "numpy-pcg64", ("ENS-B2", "S2", "physics", f"{r:04d}"))
        pop_b.append(run.populations)
        dec_b.append(run.decays)
        pho_b.append(run.photons)
        seeds["ENS-B2"].append(s)
        digests["ENS-B2"].append(stream.source.stream_digest)
        d1 &= b2_consumption(run, stream)["ok"]
        d2 &= all(bookkeeping(run).values())
        s, label = seed("ENS-A", "S2", "physics", f"{r:04d}")
        stream = WordStream(make_source("numpy-pcg64", s, label))
        layer_a = simulate_decays(s2, {"Xa-400": n_ens}, times2[-1], stream)
        pops, decays, photons, at_zero = layer_a_observations(layer_a, times2)
        pop_a.append(pops)
        dec_a.append(decays)
        pho_a.append(photons)
        a_zero += at_zero
        a_invariants &= all(decay_invariants(layer_a).values())
        a_consumption &= stream.words_read == n_ens * s2.words_per_nucleus
        seeds["ENS-A"].append(s)
        digests["ENS-A"].append(stream.source.stream_digest)
        s, label = seed("ENS-NAIVE", "S2", "physics", f"{r:04d}")
        stream = WordStream(make_source("numpy-pcg64", s, label))
        pop_n.append(naive_aggregate(s2, {"Xa-400": n_ens}, times2, stream))
        seeds["ENS-NAIVE"].append(s)
        digests["ENS-NAIVE"].append(stream.source.stream_digest)
    result["timings_seconds"]["ensembles"] = time.perf_counter() - t0
    pop_b, dec_b, pho_b = np.array(pop_b), np.array(dec_b), np.array(pho_b)
    pop_a, dec_a, pho_a, pop_n = np.array(pop_a), np.array(dec_a), np.array(pho_a), np.array(pop_n)
    d1 &= a_consumption
    d2 &= a_invariants
    deterministic["D7 ensemble streams"] = (all(len(set(v)) == reps for v in seeds.values())
                                            and all(len(set(v)) == reps for v in digests.values()) and a_zero == 0)
    names = [n.name for n in s2.nuclides]
    index_of_time = {repr(t): k for k, t in enumerate(times2)}
    ens_tests, comp_tests = {}, {}
    for t in ("2.0", "8.0"):
        it = index_of_time[t]
        p = state_probabilities(s2, "Xa-400", float(t))[0]
        for name in design["ens"]["eligible"][t]:
            j = s2.index[name]
            ens_tests[f"ENS ET-MEAN {name} t={t}"] = binomial_check(pop_b[:, it, j].sum(), reps * n_ens, p[j])
            moments = state_count_moments(s2, "Xa-400", n_ens, float(t), j, j)
            ens_tests[f"ENS ET-VAR {name} t={t}"] = scaled_chisquare_variance_check(
                sample_covariance(pop_b[:, it, j], pop_b[:, it, j]), moments.sigma_ab,
                moments.sample_covariance_variance(reps))
    i2, i8 = index_of_time["2.0"], index_of_time["8.0"]
    for a, b in design["ens"]["pairs_t2"]:
        moments = state_sum_moments(s2, "Xa-400", n_ens, 2.0, s2.index[a], s2.index[b])
        total = pop_b[:, i2, s2.index[a]] + pop_b[:, i2, s2.index[b]]
        ens_tests[f"ENS ET-VARSUM {a}+{b} t=2.0"] = scaled_chisquare_variance_check(
            sample_covariance(total, total), moments.sigma_ab, moments.sample_covariance_variance(reps))
    for a, b in design["ens"]["pairs_two_time"]:
        moments = two_time_sum_moments(s2, "Xa-400", n_ens, 2.0, 8.0, s2.index[a], s2.index[b])
        total = pop_b[:, i2, s2.index[a]] + pop_b[:, i8, s2.index[b]]
        ens_tests[f"ENS ET-VARSUM2 {a}@2+{b}@8"] = scaled_chisquare_variance_check(
            sample_covariance(total, total), moments.sigma_ab, moments.sample_covariance_variance(reps))
    for t in times2[1:]:
        it = index_of_time[repr(t)]
        comp_tests[f"COMP CT-HOMO t={t!r}"] = homogeneity_check(pop_a[:, it, :].sum(axis=0), pop_b[:, it, :].sum(axis=0))
    for name in design["ens"]["eligible"]["2.0"]:
        j = s2.index[name]
        moments = state_count_moments(s2, "Xa-400", n_ens, 2.0, j, j)
        comp_tests[f"COMP CT-VAR {name} t=2.0"] = variance_ratio_f_check(
            sample_covariance(pop_a[:, i2, j], pop_a[:, i2, j]), sample_covariance(pop_b[:, i2, j], pop_b[:, i2, j]),
            moments.sigma_ab, moments.sample_covariance_variance(reps))
    for a, b in design["ens"]["pairs_t2"]:
        ja, jb = s2.index[a], s2.index[b]
        moments = state_sum_moments(s2, "Xa-400", n_ens, 2.0, ja, jb)
        total_a, total_b = pop_a[:, i2, ja] + pop_a[:, i2, jb], pop_b[:, i2, ja] + pop_b[:, i2, jb]
        comp_tests[f"COMP CT-VARSUM {a}+{b} t=2.0"] = variance_ratio_f_check(
            sample_covariance(total_a, total_a), sample_covariance(total_b, total_b), moments.sigma_ab,
            moments.sample_covariance_variance(reps))
    for name in design["comparison"]["unstable"]:
        j = s2.index[name]
        va = dec_a[:, :, j, :].sum(axis=(0, 2))
        vb = dec_b[:, :, j, :].sum(axis=(0, 2))
        comp_tests[f"COMP CT-DECAY {name}"] = homogeneity_check(np.append(va, reps * n_ens - va.sum()),
                                                                np.append(vb, reps * n_ens - vb.sum()))
    for j, nuclide in enumerate(s2.nuclides):
        for k, channel in enumerate(nuclide.channels):
            if 0.0 < channel.photon_yield < 1.0:
                pa, da = pho_a[:, :, j, k].sum(), dec_a[:, :, j, k].sum()
                pb, db = pho_b[:, :, j, k].sum(), dec_b[:, :, j, k].sum()
                comp_tests[f"COMP CT-PHOTON {nuclide.name}/{channel.mode.value}->{channel.daughter}"] = \
                    homogeneity_check([pa, da - pa], [pb, db - pb])
    result["ensemble"] = {"replicates": reps, "per_replicate": n_ens, "nuclides": names, "times": times2,
                          "populations_b2": pop_b.tolist(), "populations_layer_a": pop_a.tolist(),
                          "populations_naive": pop_n.tolist(), "layer_a_decays_at_time_zero": a_zero,
                          "tests": ens_tests, "comparison_tests": comp_tests}
    (out_dir / "manifests" / "ENSEMBLES_S2.json").write_text(json.dumps(
        {"schema": "rnglab-nuclear-ensemble-manifest-v1", "experiment": experiment, "code": code,
         "seed_recipe": protocol["seed_recipe"], "seeds": seeds, "stream_sha256": digests}, indent=1) + "\n")
    tests["NC1"].update(ens_tests)
    tests["NC1"].update(comp_tests)

    # Huge populations.
    big_tests, result["big"] = {}, {}
    q_end = decay_probabilities(s2, "Xa-400", times2[-1])[0]
    p_grid = state_probabilities(s2, "Xa-400", times2)
    d9 = True
    for size in (smoke["big_sizes"] if args.smoke else runs["BIG"]["initial_sizes"]):
        run, stream, seconds, s, label = b2(s2, {"Xa-400": size}, times2, "numpy-pcg64", ("BIG", "S2", "physics", size))
        d1 &= b2_consumption(run, stream)["ok"]
        d2 &= all(bookkeeping(run).values())
        d9 &= run.plan.parts == -(-size // 2 ** 53)
        for it, t in enumerate(times2[1:], start=1):
            big_tests[f"BIG {size} BIG-STATE t={t!r}"] = chi2(run.populations[it], p_grid[it])
        for j, nuclide in enumerate(s2.nuclides):
            if not nuclide.stable:
                total_j = int(run.decays[:, j, :].sum())
                big_tests[f"BIG {size} BIG-DECAY {nuclide.name}"] = z_check(
                    total_j, size * q_end[j], size * q_end[j] * (1 - q_end[j]), "binomial_normal_approximation")
        result["big"][str(size)] = {"parts": run.plan.parts, "seconds": seconds, "summary": run.summary(),
                                    "populations": [[int(x) for x in row] for row in run.populations]}
        result["timings_seconds"][f"BIG {size}"] = seconds
        manifest(f"BIG_{size}", run, stream, label)
    if not args.smoke:
        d9 &= [result["big"][str(s)]["parts"] for s in runs["BIG"]["initial_sizes"]] == [1, 129]
    deterministic["D9 big populations"] = d9
    tests["NC1"].update({k: v for k, v in big_tests.items() if "p_value" in v})
    result["big_tests"] = big_tests

    # Positive controls.
    pc = runs["PC-B1"]
    s1_perturbed = copy.deepcopy(protocol["schemes"]["S1"])
    s1_perturbed["nuclides"][0]["half_life"] = pc["perturbation"]["Sa-400"]["half_life"]
    s1_perturbed["name"] += "-perturbed"
    pc1_case = {"initial": scaled(pc["initial"], scale), "times": protocol["cases"]["S1"]["times"]}
    run, stream, _, _, label = b2(DecayScheme.from_dict(s1_perturbed), pc1_case["initial"], pc1_case["times"],
                                  pc["generator"], ("PC-B1", "S1", "physics"))
    pc1_tests = battery("S1", run, schemes["S1"], pc1_case, "PC-B1")
    manifest("PC-B1_S1", run, stream, label)
    s_pc2, label_pc2 = seed("PC-B2", "S2", "physics")
    naive = naive_aggregate(s2, scaled(runs["PC-B2"]["initial"], scale), times2,
                            WordStream(make_source("numpy-pcg64", s_pc2, label_pc2)))
    p12 = state_probabilities(s2, "Xa-400", times2)
    pc2_tests = {f"PC-B2 S2 BT-STATE Xa-400 t={t!r}": chi2(naive[it], p12[it])
                 for it, t in enumerate(times2[1:], start=1)}
    pc2_comp = homogeneity_check(pop_a[:, -1, :].sum(axis=0), pop_n[:, -1, :].sum(axis=0))

    # Software checks.
    xref = []
    for case in runs["XREF"]["cases"]:
        key = case["scheme"]
        initial = scaled(case["initial"], smoke["xref_scale"] if args.smoke else 1.0)
        s, label = seed("XREF", key, "physics", case["part_max"])
        fast_stream = WordStream(make_source("numpy-pcg64", s, label))
        slow_stream = WordStream(make_source("numpy-pcg64", s, label))
        fast = simulate_aggregate(schemes[key], initial, protocol["cases"][key]["times"], fast_stream,
                                  part_max=case["part_max"])
        t0 = time.perf_counter()
        slow = simulate_aggregate_reference(schemes[key], initial, protocol["cases"][key]["times"], slow_stream,
                                            part_max=case["part_max"])
        same = all(np.array_equal(getattr(fast, f), getattr(slow, f))
                   for f in ("populations", "path_counts", "decays", "photons", "neutrons"))
        xref.append({"scheme": key, "initial": initial, "part_max": case["part_max"], "parts": fast.plan.parts,
                     "digests_equal": fast.digest == slow.digest, "arrays_equal": bool(same),
                     "stream_digests_equal": fast_stream.source.stream_digest == slow_stream.source.stream_digest,
                     "reference_seconds": time.perf_counter() - t0})
        d1 &= b2_consumption(fast, fast_stream)["ok"] and b2_consumption(slow, slow_stream)["ok"]
    result["xref"] = xref
    deterministic["D3 xref"] = all(x["digests_equal"] and x["arrays_equal"] and x["stream_digests_equal"] for x in xref)
    run0, stream0, s0, label0 = nc1_s2
    replay, replay_stream, _, _, _ = b2(s2, cases["S2"]["initial"], times2, "numpy-pcg64", ("NC1", "S2", "physics"))
    deterministic["D4 replay"] = (replay.digest == run0.digest
                                  and replay_stream.source.stream_digest == stream0.source.stream_digest)
    checkpoints = []
    full, _, _, _, _ = b2(s2, cases["S2"]["initial"], times2, "numpy-pcg64", ("NC1", "S2", "physics"),
                          checkpoint_every=1, on_checkpoint=checkpoints.append)
    cut = runs["CHECKPOINT"]["resume_after_interval"]
    checkpoint = json.loads(json.dumps(checkpoints[cut - 1]))
    resumed = simulate_aggregate(s2, cases["S2"]["initial"], times2,
                                 WordStream(NumpyPCG64Source(s0, label=label0, start_word=checkpoint["words_consumed"])),
                                 resume_from=checkpoint)
    deterministic["D5 checkpoint"] = (resumed.digest == full.digest == run0.digest
                                      and np.array_equal(resumed.populations, run0.populations[cut:]))
    records = []
    lean, _, _, _, _ = b2(s2, cases["S2"]["initial"], times2, "numpy-pcg64", ("NC1", "S2", "physics"),
                          keep_outputs=False, sink=lambda k, record: records.append(record["populations"]))
    deterministic["D6 sink"] = lean.digest == run0.digest and np.array_equal(np.array(records), run0.populations[1:])
    deterministic["D1 consumption"] = bool(d1)
    deterministic["D2 bookkeeping and Layer A invariants"] = bool(d2)

    # Families, Holm, verdicts.
    alpha = protocol["families"]["alpha"]
    result["families"] = {}
    for family, members in tests.items():
        names_f = list(members)
        correction = holm([members[n]["p_value"] for n in names_f], alpha)
        result["families"][family] = {"tests": names_f, "p_values": [members[n]["p_value"] for n in names_f],
                                      "holm": correction,
                                      "rejected": [n for n, r in zip(names_f, correction["reject"]) if r]}
    if not args.smoke:
        deterministic["D10 planned tests ran"] = all(
            sorted(tests[f]) == sorted(design["planned_tests"][f]) for f in ("NC1", "NC2"))
    level = alpha / len(tests["NC1"])
    t1 = pc1_tests["PC-B1 S1 BT-STATE Sa-400 t=2.0"]["p_value"]
    t2 = pc2_tests["PC-B2 S2 BT-STATE Xa-400 t=12.0"]["p_value"]
    t3 = pc2_comp["p_value"]
    result["positive_controls"] = {
        "level": level,
        "PC-B1": {"target": "BT-STATE Sa-400 t=2.0", "p_value": t1, "rejected": t1 < level, "all_tests": pc1_tests},
        "PC-B2 single run": {"target": "BT-STATE Xa-400 t=12.0", "p_value": t2, "rejected": t2 < level,
                             "all_tests": pc2_tests},
        "PC-B2 comparison": {"target": "CT-HOMO t=12.0, ENS-A against ENS-NAIVE", "p_value": t3,
                             "rejected": t3 < level, "test": pc2_comp}}
    failed = [k for k, ok in deterministic.items() if not ok]
    failed += [f"Holm rejection in {f}: {n}" for f, v in result["families"].items() for n in v["rejected"]]
    result["verdicts"] = {"deterministic": deterministic, "R1": "PASS" if not failed else "FAIL",
                          "R1_failed_items": failed,
                          "R2": "MET" if all(result["positive_controls"][k]["rejected"]
                                             for k in ("PC-B1", "PC-B2 single run", "PC-B2 comparison")) else "NOT MET"}
    result["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    result["elapsed_seconds"] = time.perf_counter() - clock
    (out_dir / "verification.json").write_text(json.dumps(result, indent=1, default=str) + "\n")
    print(json.dumps({"namespace": ns, "R1": result["verdicts"]["R1"], "failed": failed,
                      "R2": result["verdicts"]["R2"], "PC p-values": [t1, t2, t3], "level": level,
                      "family sizes": {f: len(v) for f, v in tests.items()},
                      "smallest p": {f: min(v["p_values"]) for f, v in result["families"].items()},
                      "elapsed_seconds": round(result["elapsed_seconds"], 1)}, indent=1))
    return 0 if result["verdicts"]["R1"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
