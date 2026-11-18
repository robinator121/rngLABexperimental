"""Null tail rates of single tests: BT-NEUTRON and BT-PHOTON-T at 1e-2, 1e-3 and the families' Holm levels.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/tail_probe.py [--workers N] [--scale S]

Design evidence before the freeze: tail_probe_namespace of protocol.json, output tail_probe_development.json
next to this script (never results). The development calibration counts, per test type, the batteries in
which a type rejects at the Holm level; at 4,000 batteries those are a few batteries per type, and the
largest of 27 type ratios is a selected one. This probe instead draws each test's own statistic under its
exact null law many times:

- BT-NEUTRON (S4, one test per start nuclide): the start nuclide's neutrons per interval, and the nuclei
  without one, are Multinomial(n, cells) under the null (run_verification.neutron_law); drawn directly.
- BT-PHOTON-T (one test per channel with pre-declared interval cells): decays drawn by calibration.NullModel
  (numpy on the exact laws) and photons binomial in them, as the calibration draws them.

The statistics are run_verification's own. BT-PHOTON-T calls photon_timing_test for every replicate.
BT-NEUTRON computes the Pearson statistic of a chunk at once from the cells and expected counts that
evaluate(chi2, ...) returns for the chunk's first draw (they depend only on n and the cells), and checks the
first CHECKED draws of every chunk against evaluate(chi2, ...) itself; a difference stops the probe. Each
rejection rate comes with its Clopper-Pearson interval, both divided by the level (1 at nominal
calibration). Draws that give no p-value (not applicable, or an error) are counted apart.
"""
from __future__ import annotations

import argparse
import datetime
import importlib.util
import json
import sys
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy import stats

from rnglab.nuclear import DecayScheme
from rnglab.nuclear.manifest import environment
from rnglab.seeds import namespace_seed

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
NEUTRON_DRAWS = 20_000_000   # per BT-NEUTRON test
PHOTON_DRAWS = 1_000_000     # NullModel replicates per scheme with BT-PHOTON-T tests
NEUTRON_CHUNK, PHOTON_CHUNK, CHECKED = 200_000, 10_000, 100

_STATE: dict = {}


def _module(name: str, file: str):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, HERE / file)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def _init_worker() -> None:
    verify = _module("b2_verify_run_verification", "run_verification.py")
    calibration = _module("b2_verify_calibration", "calibration.py")
    protocol = json.loads((HERE / "protocol.json").read_text())
    design = json.loads((HERE / "design.json").read_text())
    schemes = {k: DecayScheme.from_dict(v) for k, v in protocol["schemes"].items()}
    alpha = protocol["families"]["alpha"]
    holm = {f: alpha / len(design["planned_tests"][f]) for f in ("NC1", "NC2")}
    _STATE.update(verify=verify, calibration=calibration, protocol=protocol, schemes=schemes,
                  cases=protocol["cases"], levels=(1e-2, 1e-3, holm["NC2"], holm["NC1"]), holm=holm, models={})


def _count(results: list[dict], levels) -> dict:
    hits, missing = [0] * len(levels), 0
    for result in results:
        if "p_value" not in result:
            missing += 1
            continue
        for li, level in enumerate(levels):
            hits[li] += int(result["p_value"] < level)
    return {"hits": hits, "missing": missing, "draws": len(results), "checked": len(results)}


def neutron_chunk(namespace: str, start: str, index: int, size: int) -> dict:
    verify, levels = _STATE["verify"], _STATE["levels"]
    i, k, n, cells = verify.neutron_law(_STATE["schemes"]["S4"], _STATE["cases"]["S4"], start)
    rng = np.random.default_rng(namespace_seed(namespace, "neutron", start, index))
    counts = rng.multinomial(n, cells, size=size)
    first = verify.evaluate(verify.chi2, counts[0], cells)
    if "p_value" not in first:          # n and the cells fix the pooling: no draw of this test has a p-value
        return {"hits": [0] * len(levels), "missing": size, "draws": size, "checked": 1}
    groups, expected = first["cells"], np.asarray(first["cell_expected"], dtype=np.float64)
    observed = np.stack([counts[:, group].sum(axis=1) for group in groups], axis=1).astype(np.float64)
    p = stats.chi2.sf(((observed - expected) ** 2 / expected).sum(axis=1), len(groups) - 1)
    for r in range(min(size, CHECKED)):
        study = verify.evaluate(verify.chi2, counts[r], cells)["p_value"]
        if abs(p[r] - study) > 1e-12 * study:
            raise RuntimeError(f"BT-NEUTRON {start}: the chunk statistic differs from the study's at draw {r} "
                               f"({p[r]!r} against {study!r})")
    return {"hits": [int(np.sum(p < level)) for level in levels], "missing": 0, "draws": size,
            "checked": min(size, CHECKED)}


def photon_chunk(namespace: str, key: str, index: int, size: int) -> dict:
    verify, levels = _STATE["verify"], _STATE["levels"]
    ref, case = _STATE["schemes"][key], _STATE["cases"][key]
    if key not in _STATE["models"]:
        _STATE["models"][key] = _STATE["calibration"].NullModel(ref, case["times"], sum(case["initial"].values()))
    rng = np.random.default_rng(namespace_seed(namespace, "photons", key, index))
    drawn = _STATE["models"][key].draw(case["initial"], rng, size)
    out = {}
    for j, k, cells in verify.photon_interval_cells(ref, case):
        name = verify.channel_label(ref.nuclides[j], ref.nuclides[j].channels[k])
        out[name] = _count([verify.photon_timing_test(ref, j, k, cells, drawn.decays[r], drawn.photons[r])
                            for r in range(size)], levels)
    return out


def _merge(total: dict | None, part: dict) -> dict:
    if total is None:
        return {**part, "hits": list(part["hits"])}
    total["hits"] = [a + b for a, b in zip(total["hits"], part["hits"])]
    total["missing"] += part["missing"]
    total["draws"] += part["draws"]
    total["checked"] += part["checked"]
    return total


def summarise(merged: dict, levels, holm: dict) -> dict:
    calibration = _STATE["calibration"]
    m = merged["draws"] - merged["missing"]
    rows = []
    for level, hits in zip(levels, merged["hits"]):
        low, high = calibration.clopper_pearson(hits, m)
        name = next((f"{f} Holm level" for f, h in holm.items() if h == level), None)
        rows.append({"level": level, "level_name": name, "rejections": hits, "draws_with_p_value": m,
                     "rate_over_level": hits / m / level, "ci95_over_level": [low / level, high / level]})
    return {"draws": merged["draws"], "without_p_value": merged["missing"],
            "draws_evaluated_by_the_study_function": merged["checked"], "levels": rows}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workers", type=int, default=None)
    parser.add_argument("--scale", type=float, default=1.0, help="fraction of the draws (development checks)")
    parser.add_argument("--out", type=Path, default=HERE / "tail_probe_development.json")
    args = parser.parse_args()
    if args.out.exists():
        raise SystemExit(f"{args.out} exists; outputs are never overwritten")
    _init_worker()
    calibration, verify = _STATE["calibration"], _STATE["verify"]
    namespace = _STATE["protocol"]["calibration"]["tail_probe_namespace"]
    inherited, threads = calibration.single_thread_environment()
    workers = args.workers or calibration.allocated_cpus()
    code = verify.record_code(REPO)
    neutron_draws, photon_draws = int(NEUTRON_DRAWS * args.scale), int(PHOTON_DRAWS * args.scale)
    starts = [start for start, _ in verify.groups_of(_STATE["schemes"]["S4"], _STATE["cases"]["S4"]["initial"])]
    keys = [key for key in _STATE["protocol"]["runs"]["NC1"]["schemes"]
            if verify.photon_interval_cells(_STATE["schemes"][key], _STATE["cases"][key])]
    jobs = [("neutron", start, i, min(NEUTRON_CHUNK, neutron_draws - s)) for start in starts
            for i, s in enumerate(range(0, neutron_draws, NEUTRON_CHUNK))]
    jobs += [("photons", key, i, min(PHOTON_CHUNK, photon_draws - s)) for key in keys
             for i, s in enumerate(range(0, photon_draws, PHOTON_CHUNK))]
    clock = time.perf_counter()
    merged: dict[str, dict] = {}
    with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker) as pool:
        futures = [(kind, label, pool.submit(neutron_chunk if kind == "neutron" else photon_chunk, namespace, label,
                                             index, size)) for kind, label, index, size in jobs]
        for kind, label, future in futures:
            part = future.result()
            if kind == "neutron":
                merged[f"S4 BT-NEUTRON {label}"] = _merge(merged.get(f"S4 BT-NEUTRON {label}"), part)
            else:
                for channel, counts in part.items():
                    name = f"{label} BT-PHOTON-T {channel}"
                    merged[name] = _merge(merged.get(name), counts)
    levels, holm = _STATE["levels"], _STATE["holm"]
    result = {"generated_by": "research/v020-nuclear-b2-verify/tail_probe.py", "development": True,
              "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "namespace": namespace,
              "code": code, "environment": environment(), "workers": workers, "thread_environment": threads,
              "inherited_thread_environment": inherited, "scale": args.scale,
              "draws": {"BT-NEUTRON per test": neutron_draws, "BT-PHOTON-T NullModel replicates per scheme":
                        photon_draws}, "levels": list(levels), "holm_levels": holm,
              "tests": {name: summarise(counts, levels, holm) for name, counts in sorted(merged.items())},
              "elapsed_seconds": time.perf_counter() - clock}
    args.out.write_text(json.dumps(result, indent=1) + "\n")
    for name, test in result["tests"].items():
        cells = "  ".join(f"{row['rate_over_level']:.2f} [{row['ci95_over_level'][0]:.2f}, "
                          f"{row['ci95_over_level'][1]:.2f}]" for row in test["levels"])
        print(f"{name:<44} {cells}  (no p-value: {test['without_p_value']})")
    print(f"levels {[f'{level:.3g}' for level in levels]}; wrote {args.out} in {result['elapsed_seconds']:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
