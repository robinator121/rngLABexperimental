"""Layer A verification against analytic decay laws: the frozen study in protocol.json.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-a-verify/run_verification.py [--smoke]

The script refuses to draw if a frozen file differs from protocol.sha256. The
main run also refuses on a dirty Git tree or when results/verification.json
exists (results are never overwritten). --smoke uses the smoke namespace and
reduced sizes, writes to smoke/ (not in Git) and is excluded from every result.
"""
from __future__ import annotations

import argparse
import copy
import datetime
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from rnglab.nuclear import (
    DecayScheme, NumpyPCG64Source, WordStream, decay_probabilities, generator_matrix, path_state_probabilities,
    simulate_decays, simulate_decays_reference, state_probabilities, survival_probability,
    truncated_exponential_cdf,
)
from rnglab.nuclear.checks import (
    binomial_check, decay_invariants, holm, ks_check, multinomial_check, variance_ratio_check,
)
from rnglab.nuclear.decay import EVENT_COLUMNS
from rnglab.nuclear.manifest import code_identity, environment, file_sha256, run_manifest
from rnglab.seeds import namespace_seed
from rnglab14.sources import Xoshiro256PlusPlusSource

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FROZEN = ("protocol.json", "design.json", "design.py", "run_verification.py")


def frozen_hashes() -> dict[str, str]:
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
    raise ValueError(f"unknown generator {generator!r}")


def perturbed(scheme_dict: dict, perturbation: dict) -> DecayScheme:
    data = copy.deepcopy(scheme_dict)
    for nuclide in data["nuclides"]:
        change = perturbation.get(nuclide["name"], {})
        if "half_life" in change:
            nuclide["half_life"] = change["half_life"]
        if "branching" in change:
            for channel, branching in zip(nuclide["channels"], change["branching"], strict=True):
                channel["branching"] = branching
    data["name"] += "-perturbed"
    return DecayScheme.from_dict(data)


def unit_bins(scheme: DecayScheme, start: str, nuclide: str, horizon: float) -> tuple[np.ndarray, np.ndarray]:
    edges = np.arange(0.0, horizon + 0.5, 1.0)
    q = decay_probabilities(scheme, start, edges)[:, scheme.index[nuclide]]
    return edges, np.append(np.diff(q), 1.0 - q[-1])


def battery(key: str, run, ref: DecayScheme, case: dict, tests: dict) -> dict:
    """The tests of protocol.json for one scheme; the run is tested against the reference scheme."""
    out = {}
    horizon = case["horizon"]
    (start, count), = case["initial"].items()
    if key == "S1":
        j = ref.index["Sa-400"]
        half_life = ref.nuclides[j].half_life
        out["T1"] = ks_check(run.decay_times("Sa-400"), truncated_exponential_cdf(math.log(2.0) / half_life, horizon))
        times = tests["T2"]["times"]
        for t, observed, p in zip(times, run.populations(times)[:, j], survival_probability(half_life, times)):
            out[f"T2 t={t!r}"] = binomial_check(observed, count, p)
        out["T3"] = binomial_check(run.photon_counts()[j, 0], run.channel_counts()[j, 0],
                                   ref.nuclides[j].channels[0].photon_yield)
    elif key == "S2":
        times = tests["T4"]["times"]
        for t, observed, p in zip(times, run.populations(times), state_probabilities(ref, start, times)):
            out[f"T4 t={t!r}"] = multinomial_check(observed, p)
        q = decay_probabilities(ref, start, horizon)[0]
        decays, channels, photons = run.decay_counts(), run.channel_counts(), run.photon_counts()
        for j, nuclide in enumerate(ref.nuclides):
            if not nuclide.stable:
                out[f"T5 {nuclide.name}"] = binomial_check(decays[j], count, q[j])
        for name in tests["T6"]["nuclides"]:
            j = ref.index[name]
            out[f"T6 {name}"] = binomial_check(channels[j, 0], decays[j], ref.nuclides[j].channel_probabilities()[0])
        for name in tests["T7"]["nuclides"]:
            edges, probabilities = unit_bins(ref, start, name, horizon)
            when = run.decay_times(name)
            binned = np.bincount(np.searchsorted(edges[1:-1], when, side="left"), minlength=edges.size - 1)
            out[f"T7 {name}"] = multinomial_check(np.append(binned, count - when.size), probabilities)
        for j, nuclide in enumerate(ref.nuclides):
            for k, channel in enumerate(nuclide.channels):
                if 0.0 < channel.photon_yield < 1.0:
                    out[f"T8 {nuclide.name}/{channel.mode.value}"] = binomial_check(
                        photons[j, k], channels[j, k], channel.photon_yield)
    elif key == "S3":
        times = tests["T9"]["times"]
        for t, observed, p in zip(times, run.populations(times), state_probabilities(ref, start, times)):
            out[f"T9 t={t!r}"] = multinomial_check(observed, p)
    return out


def consumption(run, stream) -> dict:
    words = run.n_nuclei * run.words_per_nucleus
    return {"words_expected": words, "words_read": stream.words_read,
            "source_bytes": int(stream.source.consumed_bytes),
            "ok": stream.words_read == words and stream.source.consumed_bytes == 8 * words}


def simulate(scheme, case, generator, seed, label, *, chunk=65536, engine=simulate_decays):
    stream = WordStream(make_source(generator, seed, label))
    started = time.perf_counter()
    if engine is simulate_decays:
        run = engine(scheme, case["initial"], case["horizon"], stream, chunk=chunk)
    else:
        run = engine(scheme, case["initial"], case["horizon"], stream)
    return run, stream, time.perf_counter() - started


def scaled(case: dict, scale: float) -> dict:
    return {"initial": {k: max(1, round(v * scale)) for k, v in case["initial"].items()}, "horizon": case["horizon"]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--smoke", action="store_true", help="smoke namespace, reduced sizes, output in smoke/")
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
    namespace = protocol["smoke_namespace"] if args.smoke else protocol["namespace"]
    scale = protocol["smoke"]["scale"] if args.smoke else 1.0
    started_utc = datetime.datetime.now(datetime.timezone.utc).isoformat()
    clock = time.perf_counter()

    schemes = {key: DecayScheme.from_dict(value) for key, value in protocol["schemes"].items()}
    if {k: s.sha256() for k, s in schemes.items()} != design["scheme_sha256"]:
        raise SystemExit("schemes differ from design.json")
    cases = {key: scaled(case, scale) for key, case in protocol["cases"].items()}
    tests, runs_spec = protocol["tests"], protocol["runs"]
    experiment = protocol["id"] + ("-smoke" if args.smoke else "")
    result: dict = {"protocol": {"id": protocol["id"], "frozen_sha256": frozen}, "namespace": namespace,
                    "smoke": args.smoke, "code": code, "environment": environment(), "started_utc": started_utc,
                    "runs": {}, "timings_seconds": {}}

    def save_manifest(run_id, run, stream, initial, extra):
        manifest = run_manifest(experiment=experiment, run_id=run_id, run=run, stream=stream, initial=initial,
                                code=code, extra=extra)
        (out_dir / "manifests" / f"{run_id}.json").write_text(json.dumps(manifest, indent=1) + "\n")

    # A1: analytic cross-check, no random numbers. NaN propagates; a refused evaluation is a failure.
    differences, column_sums, a1_error = [], [], None
    try:
        for key, scheme in schemes.items():
            times = sorted({*tests[{"S1": "T2", "S2": "T4", "S3": "T9"}[key]]["times"],
                            protocol["cases"][key]["horizon"]})
            column_sums.append(np.max(np.abs(generator_matrix(scheme).sum(axis=0))))
            for nuclide in scheme.nuclides:
                differences.append(np.max(np.abs(state_probabilities(scheme, nuclide.name, times)
                                                 - path_state_probabilities(scheme, nuclide.name, times))))
    except ValueError as error:
        a1_error = str(error)
    worst = float(np.max(differences)) if differences else float("nan")
    columns = float(np.max(column_sums)) if column_sums else float("nan")
    result["analytic_crosscheck"] = {"max_abs_expm_minus_paths": worst, "max_abs_generator_column_sum": columns,
                                     "error": a1_error,
                                     "ok": a1_error is None and bool(worst <= 1e-11) and bool(columns <= 1e-15)}

    # NC1, NC2, PC1, PC2: simulate, check invariants and consumption, run the battery.
    stored = {}
    for run_name in ("NC1", "NC2", "PC1", "PC2"):
        spec = runs_spec[run_name]
        keys = spec.get("schemes") or [spec["scheme"]]
        result["runs"][run_name] = {}
        for key in keys:
            scheme = perturbed(protocol["schemes"][key], spec["perturbation"]) if "perturbation" in spec else schemes[key]
            seed = namespace_seed(namespace, run_name, key, "physics")
            label = f"{namespace}/{run_name}/{key}/physics"
            run, stream, seconds = simulate(scheme, cases[key], spec["generator"], seed, label)
            entry = {"generator": spec["generator"], "seed": seed, "seed_label": label, "seconds": seconds,
                     "summary": run.summary(), "stream_sha256": stream.source.stream_digest,
                     "consumption": consumption(run, stream), "invariants": decay_invariants(run),
                     "tests": battery(key, run, schemes[key], cases[key], tests)}
            if run_name == "NC1" and key == "S3":
                ya, yb = schemes["S3"].index["Ya-410"], schemes["S3"].index["Yb-410"]
                ev = run.events
                in_window = (ev["time"] > 20.0) & (ev["time"] <= 40.0)
                # Per-nucleus indicators: each nucleus decays out of each nuclide at most once.
                a = np.zeros(run.n_nuclei)
                b = np.zeros(run.n_nuclei)
                a[ev["nucleus"][in_window & (ev["parent"] == ya)]] = 1.0
                b[ev["nucleus"][in_window & (ev["parent"] == yb)]] = 1.0
                d_a, d_b = int(a.sum()), int(b.sum())
                mean_a, mean_b = a.mean(), b.mean()
                covariance = np.cov(a, b, ddof=1)
                log_variance = (covariance[0, 0] / mean_a**2 + covariance[1, 1] / mean_b**2
                                - 2.0 * covariance[0, 1] / (mean_a * mean_b)) / run.n_nuclei
                half_width = 1.96 * math.sqrt(log_variance)
                q = decay_probabilities(schemes["S3"], "Ya-410", [20.0, 40.0])
                exact = (q[1, yb] - q[0, yb]) / (q[1, ya] - q[0, ya])
                rates = schemes["S3"].decay_constants
                entry["R3"] = {"decays_Ya_20_40": d_a, "decays_Yb_20_40": d_b, "ratio": d_b / d_a,
                               "log_ratio_se": math.sqrt(log_variance),
                               "indicator_correlation": float(covariance[0, 1] / math.sqrt(covariance[0, 0] * covariance[1, 1])),
                               "interval95": [d_b / d_a * math.exp(-half_width), d_b / d_a * math.exp(half_width)],
                               "exact_expected_ratio": exact,
                               "transient_equilibrium_limit": rates[yb] / (rates[yb] - rates[ya])}
            result["runs"][run_name][key] = entry
            result["timings_seconds"][f"{run_name} {key}"] = seconds
            save_manifest(f"{run_name}_{key}", run, stream, cases[key]["initial"],
                          {"seed_label": label, "perturbation": spec.get("perturbation")})
            stored[(run_name, key)] = (run, stream)

    # ENS: replicate ensemble at t = 2.
    ens = runs_spec["ENS"]
    replicates = protocol["smoke"]["ens_replicates"] if args.smoke else ens["replicates"]
    ens_case = {"initial": ens["initial"], "horizon": ens["horizon"]}
    counts, seeds, digests, invariant_ok, consumption_ok = [], [], [], True, True
    ens_clock = time.perf_counter()
    for r in range(replicates):
        seed = namespace_seed(namespace, "ENS", "S2", "physics", f"{r:04d}")
        run, stream, _ = simulate(schemes["S2"], ens_case, ens["generator"], seed,
                                  f"{namespace}/ENS/S2/physics/{r:04d}")
        counts.append(run.populations([ens["time"]])[0])
        seeds.append(seed)
        digests.append(stream.source.stream_digest)
        invariant_ok &= all(decay_invariants(run).values())
        consumption_ok &= consumption(run, stream)["ok"]
    result["timings_seconds"]["ENS"] = time.perf_counter() - ens_clock
    counts = np.array(counts)
    p_ens = state_probabilities(schemes["S2"], "Xa-400", ens["time"])[0]
    per_replicate = ens["initial"]["Xa-400"]
    eligible = [n.name for n, p in zip(schemes["S2"].nuclides, p_ens) if per_replicate * p * (1 - p) >= 20.0]
    if eligible != design["ens"]["eligible"]:
        raise SystemExit("ENS eligibility differs from design.json")
    ens_tests = {}
    for name in eligible:
        j = schemes["S2"].index[name]
        ens_tests[f"T10a {name}"] = binomial_check(counts[:, j].sum(), replicates * per_replicate, p_ens[j])
        ens_tests[f"T10b {name}"] = variance_ratio_check(counts[:, j], per_replicate * p_ens[j] * (1 - p_ens[j]))
    result["ensemble"] = {"replicates": replicates, "counts_at_t": counts.tolist(),
                          "nuclides": [n.name for n in schemes["S2"].nuclides], "seeds": seeds,
                          "stream_sha256": digests, "invariants_ok": bool(invariant_ok),
                          "consumption_ok": bool(consumption_ok),
                          "distinct_seeds": len(set(seeds)) == replicates,
                          "distinct_streams": len(set(digests)) == replicates, "tests": ens_tests}
    (out_dir / "manifests" / "ENS_S2.json").write_text(json.dumps(
        {"schema": "rnglab-nuclear-ensemble-manifest-v1", "experiment": experiment, "code": code,
         "seed_recipe": protocol["seed_recipe"], "generator": ens["generator"], "case": ens_case,
         "seeds": seeds, "stream_sha256": digests}, indent=1) + "\n")

    # XREF and REPLAY: software checks.
    xref = runs_spec["XREF"]
    xref_case = {"initial": {"Xa-400": protocol["smoke"]["xref_nuclei"] if args.smoke else xref["initial"]["Xa-400"]},
                 "horizon": xref["horizon"]}
    seed = namespace_seed(namespace, "XREF", "S2", "physics")
    label = f"{namespace}/XREF/S2/physics"
    fast, fast_stream, _ = simulate(schemes["S2"], xref_case, xref["generator"], seed, label, chunk=997)
    slow, slow_stream, slow_seconds = simulate(schemes["S2"], xref_case, xref["generator"], seed, label,
                                               engine=simulate_decays_reference)
    discrete = [c for c in EVENT_COLUMNS if c != "time"]
    same_shape = fast.n_events == slow.n_events
    relative = (float(np.max(np.abs(fast.events["time"] - slow.events["time"])
                             / np.maximum(np.abs(slow.events["time"]), np.finfo(float).tiny)))
                if same_shape and fast.n_events else 0.0)
    result["xref"] = {"nuclei": xref_case["initial"]["Xa-400"], "events": [fast.n_events, slow.n_events],
                      "discrete_columns_equal": same_shape and all(np.array_equal(fast.events[c], slow.events[c])
                                                                   for c in discrete),
                      "final_states_equal": bool(np.array_equal(fast.final_state, slow.final_state)),
                      "stream_digests_equal": fast_stream.source.stream_digest == slow_stream.source.stream_digest,
                      "max_relative_time_difference": relative, "digests_equal": fast.digest() == slow.digest(),
                      "reference_engine_seconds": slow_seconds}
    result["xref"]["ok"] = bool(result["xref"]["discrete_columns_equal"] and result["xref"]["final_states_equal"]
                                and result["xref"]["stream_digests_equal"] and relative <= 1e-12)
    result["xref"]["consumption"] = {"vectorised": consumption(fast, fast_stream),
                                     "reference": consumption(slow, slow_stream)}
    save_manifest("XREF_S2_vectorised", fast, fast_stream, xref_case["initial"], {"seed_label": label, "chunk": 997})
    save_manifest("XREF_S2_reference", slow, slow_stream, xref_case["initial"], {"seed_label": label})

    original, original_stream = stored[("NC1", "S2")]
    replay, replay_stream, _ = simulate(schemes["S2"], cases["S2"], runs_spec["NC1"]["generator"],
                                        namespace_seed(namespace, "NC1", "S2", "physics"),
                                        f"{namespace}/NC1/S2/physics", chunk=runs_spec["REPLAY"]["chunk"])
    result["replay"] = {"chunk": runs_spec["REPLAY"]["chunk"], "events_sha256": replay.digest(),
                        "original_events_sha256": original.digest(),
                        "stream_sha256": replay_stream.source.stream_digest,
                        "original_stream_sha256": original_stream.source.stream_digest}
    result["replay"]["ok"] = (replay.digest() == original.digest()
                              and replay_stream.source.stream_digest == original_stream.source.stream_digest)
    result["replay"]["consumption"] = consumption(replay, replay_stream)

    # Families, Holm, positive controls, verdicts.
    alpha = protocol["families"]["alpha"]
    families = {"NC1": {**{f"NC1 {k} {name}": t for k in ("S1", "S2", "S3")
                           for name, t in result["runs"]["NC1"][k]["tests"].items()},
                        **{f"ENS {name}": t for name, t in ens_tests.items()}},
                "NC2": {f"NC2 {k} {name}": t for k in ("S1", "S2", "S3")
                        for name, t in result["runs"]["NC2"][k]["tests"].items()}}
    result["families"] = {}
    for family, members in families.items():
        names = list(members)
        correction = holm([members[n]["p_value"] for n in names], alpha)
        result["families"][family] = {"tests": names, "p_values": [members[n]["p_value"] for n in names],
                                      "holm": correction,
                                      "rejected": [n for n, r in zip(names, correction["reject"]) if r]}
    sizes = {k: len(v) for k, v in families.items()}
    if not args.smoke and sizes != {"NC1": design["family_sizes"]["NC1"], "NC2": design["family_sizes"]["NC2"]}:
        raise SystemExit(f"family sizes {sizes} differ from design.json")
    level = alpha / sizes["NC1"]
    pc1 = result["runs"]["PC1"]["S1"]["tests"]["T2 t=2.0"]["p_value"]
    pc2 = result["runs"]["PC2"]["S2"]["tests"]["T6 Xa-400"]["p_value"]
    result["positive_controls"] = {"level": level, "PC1": {"target": "T2 t=2.0", "p_value": pc1, "rejected": pc1 < level},
                                   "PC2": {"target": "T6 Xa-400", "p_value": pc2, "rejected": pc2 < level}}

    deterministic = {
        "D1 consumption": all(e["consumption"]["ok"] for r in ("NC1", "NC2", "PC1", "PC2")
                              for e in result["runs"][r].values()) and result["ensemble"]["consumption_ok"]
        and all(c["ok"] for c in result["xref"]["consumption"].values()) and result["replay"]["consumption"]["ok"],
        "D2 invariants": all(all(e["invariants"].values()) for r in ("NC1", "NC2", "PC1", "PC2")
                             for e in result["runs"][r].values()) and result["ensemble"]["invariants_ok"]
        and all(decay_invariants(fast).values()) and all(decay_invariants(slow).values()),
        "D3 xref": result["xref"]["ok"],
        "D4 replay": result["replay"]["ok"],
        "D5 ensemble streams": result["ensemble"]["distinct_seeds"] and result["ensemble"]["distinct_streams"],
        "A1 analytic cross-check": result["analytic_crosscheck"]["ok"],
    }
    failed = [k for k, ok in deterministic.items() if not ok]
    failed += [f"Holm rejection in {f}: {n}" for f, v in result["families"].items() for n in v["rejected"]]
    result["verdicts"] = {
        "deterministic": deterministic,
        "R1": "PASS" if not failed else "FAIL",
        "R1_failed_items": failed,
        "R2": "MET" if result["positive_controls"]["PC1"]["rejected"] and result["positive_controls"]["PC2"]["rejected"]
        else "NOT MET",
    }
    result["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    result["elapsed_seconds"] = time.perf_counter() - clock
    (out_dir / "verification.json").write_text(json.dumps(result, indent=1) + "\n")
    print(json.dumps({"namespace": namespace, "R1": result["verdicts"]["R1"], "failed": failed,
                      "R2": result["verdicts"]["R2"], "PC p-values": [pc1, pc2], "level": level,
                      "smallest NC1 p": min(result["families"]["NC1"]["p_values"]),
                      "smallest NC2 p": min(result["families"]["NC2"]["p_values"]),
                      "elapsed_seconds": round(result["elapsed_seconds"], 1)}, indent=1))
    return 0 if result["verdicts"]["R1"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
