#!/usr/bin/env python3
"""Generate the v0.18 Phase 1 estimator campaign (frozen protocol rnglab-v018-estimators-v1).

RNGLab v0.18 (written by Claude Code, 2026-11-05).

    PYTHONPATH=src python research/v018-estimators/run_campaign.py manifest   # pre-generation manifest
    PYTHONPATH=src python research/v018-estimators/run_campaign.py run        # campaign

`run` refuses to start unless protocol.json matches protocol.sha256 and every
file hash in manifest_pregeneration.json matches the working tree. Outputs:
results/ar1/*.npz and results/ising/*.npz (per-unit endpoints, committed) and
traces/*.npz (Ising integer M and E traces, gitignored). Existing outputs are
never overwritten.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import os
import platform
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
RESULTS = HERE / "results"
TRACES = HERE / "traces"
HASHED_FILES = [
    "research/v018-estimators/protocol.json",
    "research/v018-estimators/run_campaign.py",
    "research/v018-estimators/analyze.py",
    "src/rnglab/longrun.py",
    "src/rnglab/convergence.py",
    "src/rnglab/susceptibility.py",
    "src/rnglab/seeds.py",
    "src/rnglab12/uncertainty.py",
    "src/rnglab16/core.py",
    "research/phase2-20261030/environment/requirements-lock.txt",
]
METHODS = ("BM", "OBM", "SV-B", "SV-TH")
ADAPTIVE_FAMILIES = ("BM", "OBM", "SV-TH")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_protocol() -> dict:
    expected = (HERE / "protocol.sha256").read_text().split()[0]
    if sha256(HERE / "protocol.json") != expected:
        raise SystemExit("protocol.json does not match protocol.sha256; refusing to run")
    protocol = json.loads((HERE / "protocol.json").read_text())
    if protocol["status"] != "FROZEN_DESIGN":
        raise SystemExit("protocol is not frozen")
    return protocol


def phi_label(phi: float) -> str:
    return repr(float(phi))


def ising_conditions(protocol: dict) -> list[dict]:
    from rnglab16.core import TC

    temps = {"T1.5": 1.5, "Tc": float(TC), "T3.5": 3.5}
    return [{"condition": f"ising_L{L}_{label}", "L": L, "T": T}
            for L in protocol["data"]["ising"]["L"] for label, T in temps.items()]


def ar1_conditions(protocol: dict) -> list[dict]:
    return [{"condition": f"ar1_phi{phi_label(phi)}_n{n}", "phi": float(phi), "n": int(n)}
            for phi in protocol["data"]["ar1"]["phi"] for n in protocol["data"]["ar1"]["n"]]


def seed_table(protocol: dict) -> list[dict]:
    from rnglab.seeds import namespace_seed

    ns = protocol["seed"]["namespace"]
    rows = []
    ising = protocol["data"]["ising"]
    for cell in ising_conditions(protocol):
        for c in range(ising["chunks"]):
            unit = f"chunk-{c:04d}"
            rows.append({"label": f"{ns}/{cell['condition']}/chains/{unit}", "size": ising["chunk_size"],
                         "seed": namespace_seed(ns, cell["condition"], "chains", unit)})
    ar1 = protocol["data"]["ar1"]
    for cell in ar1_conditions(protocol):
        size = ar1["chunk_size"][str(cell["n"])]
        for c in range(ar1["replicas"] // size):
            unit = f"chunk-{c:04d}"
            rows.append({"label": f"{ns}/{cell['condition']}/innovations/{unit}", "size": size,
                         "seed": namespace_seed(ns, cell["condition"], "innovations", unit)})
    return rows


def write_manifest() -> None:
    protocol = load_protocol()
    path = HERE / "manifest_pregeneration.json"
    if path.exists():
        raise SystemExit(f"{path} exists; a new manifest needs a new protocol version")
    import numpy, scipy

    seeds = seed_table(protocol)
    if len({row["seed"] for row in seeds}) != len(seeds):
        raise SystemExit("seed collision inside the namespace")
    head = subprocess.run(["git", "rev-parse", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    manifest = {
        "protocol_id": protocol["id"], "protocol_sha256": sha256(HERE / "protocol.json"),
        "written_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "git_head_before_commit": head,
        "python": platform.python_version(), "numpy": numpy.__version__, "scipy": scipy.__version__,
        "platform": platform.platform(), "file_sha256": {f: sha256(REPO / f) for f in HASHED_FILES},
        "seed_units": len(seeds), "seeds": seeds,
        "statement": "Written before any campaign draw; run_campaign.py run verifies these hashes."}
    path.write_text(json.dumps(manifest, indent=1) + "\n")
    print(f"wrote {path} ({len(seeds)} seed units)")


def verify_manifest() -> dict:
    manifest = json.loads((HERE / "manifest_pregeneration.json").read_text())
    for name, digest in manifest["file_sha256"].items():
        if sha256(REPO / name) != digest:
            raise SystemExit(f"{name} changed since the pre-generation manifest; refusing to run")
    return manifest


def run_ar1(protocol: dict, cell: dict) -> str:
    import numpy as np

    from rnglab import longrun
    from rnglab.seeds import namespace_seed, replica_generators

    out = RESULTS / "ar1" / f"{cell['condition']}.npz"
    if out.exists():
        return f"skip {out.name} (exists)"
    start = time.monotonic()
    ns = protocol["seed"]["namespace"]
    spec = protocol["data"]["ar1"]
    n, phi = cell["n"], cell["phi"]
    size = spec["chunk_size"][str(n)]
    batches = protocol["definitions"]["batch_sizes"]
    rows: dict[str, list] = {}

    def add(key, value):
        rows.setdefault(key, []).append(np.asarray(value))

    for c in range(spec["replicas"] // size):
        seed = namespace_seed(ns, cell["condition"], "innovations", f"chunk-{c:04d}")
        x = longrun.ar1_series(replica_generators(seed, size), n, phi)
        add("chunk", np.full(size, c))
        add("xbar", x.mean(axis=1))
        add("gamma0", x.var(axis=1))
        for method in METHODS:
            for b in batches:
                est = longrun.estimate(x, method, b)
                add(f"sigma2_{method}_b{b}", est["sigma2"])
                add(f"df_{method}_b{b}", est["df"])
        g = longrun.geyer(x)
        add("sigma2_GEYER", g["sigma2"])
        add("tau_GEYER", g["tau"])
        add("valid_GEYER", g["valid"])
        bstar = longrun.adaptive_batch(n, g["tau"])
        add("bstar", bstar)
        for family in ADAPTIVE_FAMILIES:
            sigma2 = np.full(size, np.nan)
            df = np.full(size, np.nan)
            for b in np.unique(bstar[bstar > 0]):
                rows_b = bstar == b
                est = longrun.estimate(x[rows_b], family, int(b))
                sigma2[rows_b] = est["sigma2"]
                df[rows_b] = est["df"]
            add(f"sigma2_adaptive_{family}", sigma2)
            add(f"df_adaptive_{family}", df)
    arrays = {key: np.concatenate(value) for key, value in rows.items()}
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, n=n, phi=phi, exact_var_mean=longrun.ar1_exact_var_mean(phi, n),
                        wall_seconds=time.monotonic() - start, **arrays)
    return f"{out.name} {time.monotonic() - start:.1f}s"


def run_ising(protocol: dict, cell: dict) -> str:
    import numpy as np

    from rnglab import convergence, longrun, susceptibility
    from rnglab.seeds import namespace_seed
    from rnglab16.core import ising

    out = RESULTS / "ising" / f"{cell['condition']}.npz"
    trace_out = TRACES / f"{cell['condition']}.npz"
    if out.exists() or trace_out.exists():
        return f"skip {out.name} (exists)"
    start = time.monotonic()
    ns = protocol["seed"]["namespace"]
    spec = protocol["data"]["ising"]
    L, T = cell["L"], cell["T"]
    N = L * L
    size, chunks, sweeps, burn = spec["chunk_size"], spec["chunks"], spec["sweeps"], spec["burn_in"]
    M = np.empty((size * chunks, sweeps), dtype=np.int8)
    E = np.empty((size * chunks, sweeps), dtype=np.int16)
    for c in range(chunks):
        seed = namespace_seed(ns, cell["condition"], "chains", f"chunk-{c:04d}")
        trace, _, _ = ising(L, T, size, sweeps, seed)
        m_total = trace[:, :, 1] * N
        e_total = trace[:, :, 0] * N
        if np.abs(m_total - np.rint(m_total)).max() > 1e-9 or np.abs(e_total - np.rint(e_total)).max() > 1e-9:
            raise RuntimeError("non-integer total magnetization or energy")
        M[c * size:(c + 1) * size] = np.rint(m_total).astype(np.int8)
        E[c * size:(c + 1) * size] = np.rint(e_total).astype(np.int16)
    TRACES.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(trace_out, M=M, E=E, L=L, T=T, burn_in=burn)

    # Post-process in blocks of whole 16-chain ensembles to bound memory.
    cols: dict[str, list] = {}

    def add(key, value):
        cols.setdefault(key, []).append(np.asarray(value))

    block = 4096
    n_ret = sweeps - burn
    for first in range(0, len(M), block):
        m = M[first:first + block, burn:].astype(float) / N
        e = E[first:first + block, burn:].astype(float) / N
        sym = susceptibility.symmetric(m, N, T)
        series = sym.pop("series")
        g = longrun.geyer(series)
        bstar = longrun.adaptive_batch(n_ret, g["tau"])
        adaptive = np.full(len(m), np.nan)
        adaptive_df = np.full(len(m), np.nan)
        for b in np.unique(bstar[bstar > 0]):
            rows_b = bstar == b
            est = longrun.batch_means(series[rows_b], int(b))
            adaptive[rows_b] = est["sigma2"]
            adaptive_df[rows_b] = est["df"]
        obm = longrun.overlapping_batch_means(series, sym["batch"])
        svth = longrun.spectral_variance(series, sym["batch"], "tukey-hanning")
        jack = susceptibility.ensemble_jackknife(m, N, T, 16)
        add("P", susceptibility.plugin(m, N, T))
        add("chi_abs_plugin", susceptibility.abs_plugin(m, N, T))
        add("S", sym["estimate"])
        add("S_sigma2_BM16", sym["variance"] * n_ret)
        add("S_df_BM16", np.full(len(m), sym["batches"] - 1.0))
        add("S_sigma2_OBM", obm["sigma2"])
        add("S_df_OBM", obm["df"])
        add("S_sigma2_SVTH", svth["sigma2"])
        add("S_df_SVTH", svth["df"])
        add("S_sigma2_GEYER", g["sigma2"])
        add("S_tau_GEYER", g["tau"])
        add("S_bstar", bstar)
        add("S_sigma2_adaptive_BM", adaptive)
        add("S_df_adaptive_BM", adaptive_df)
        add("S_status_b", longrun.batch_length_status(sym["batch"], n_ret, g["tau"]).astype("U16"))
        add("energy_mean", e.mean(axis=1))
        add("G", jack["estimate"])
        add("G_uncorrected", jack["uncorrected"])
        add("G_variance", jack["variance"])
        for key, values in (("rhat_m", m), ("rhat_m2", m ** 2), ("rhat_e", e)):
            add(key, convergence.split_rhat(values.reshape(-1, 16, n_ret)))
    exact = susceptibility.exact_references(L, T)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        out, L=L, T=T, N=N, n_retained=n_ret, chi_exact=exact["chi"], chi_abs_exact=exact["chi_abs"],
        energy_exact=exact["energy_per_spin"], S_batch=n_ret // susceptibility.S_BATCHES,
        trace_sha256=sha256(trace_out), wall_seconds=time.monotonic() - start,
        **{key: np.concatenate(value) for key, value in cols.items()})
    return f"{out.name} {time.monotonic() - start:.1f}s"


def _task(kind: str, cell: dict) -> str:
    protocol = load_protocol()
    return run_ising(protocol, cell) if kind == "ising" else run_ar1(protocol, cell)


def run_all(workers: int) -> None:
    protocol = load_protocol()
    verify_manifest()
    tasks = [("ising", c) for c in sorted(ising_conditions(protocol), key=lambda c: -c["L"])]
    tasks += [("ar1", c) for c in ar1_conditions(protocol)]
    log = HERE / "results" / "run.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    context = multiprocessing.get_context("fork")
    with log.open("a") as handle, ProcessPoolExecutor(max_workers=workers, mp_context=context) as pool:
        handle.write(f"start {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} workers={workers}\n")
        futures = {pool.submit(_task, kind, cell): cell["condition"] for kind, cell in tasks}
        for future in as_completed(futures):
            line = f"{time.strftime('%H:%M:%S')} {future.result()}"
            print(line, flush=True)
            handle.write(line + "\n")
            handle.flush()
        handle.write(f"end {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("command", choices=["manifest", "run"])
    parser.add_argument("--workers", type=int, default=3)
    args = parser.parse_args()
    if os.environ.get("OPENBLAS_NUM_THREADS") != "1":
        raise SystemExit("set OPENBLAS_NUM_THREADS=1 and OMP_NUM_THREADS=1")
    write_manifest() if args.command == "manifest" else run_all(args.workers)


if __name__ == "__main__":
    sys.path.insert(0, str(REPO / "src"))
    main()
