#!/usr/bin/env python3
"""Runtime, byte-accounting and determinism smoke for the Phase 2 design.

RNGLab v0.18 (written by Claude Code, 2026-11-05). Uses the separate smoke
namespace rnglab-v018-api-validation-smoke; nothing produced here counts
toward the campaign. Runs a few replicas of each planned condition through the
package API (PathIntegralHarmonicOscillator.run_simulation, the package
exact-path sampler and rnglab12 summarize_series) at the planned size:
M=16, beta=2, delta=0.4, n=8192 retained sweeps, zero burn-in.

Usage: OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python research/v018-api-validation/smoke/smoke_timing.py > out.json
"""
import hashlib
import json
import platform
import resource
import sys
import time

import numpy as np

from rnglab.seeds import namespace_seed
from rnglab12.uncertainty import summarize_series
from rnglab14.models.pimc_harmonic import (
    PathIntegralHarmonicOscillator, exact_finite_bead_energy, harmonic_energy_estimators,
    sample_exact_finite_bead_paths,
)
from rnglab14.sources import PCG64Source, StochasticChannels

NAMESPACE = "rnglab-v018-api-validation-smoke"
M, BETA, DELTA, N = 16, 2.0, 0.4, 8192
REPLICAS = 4


def sha(array):
    return hashlib.sha256(np.ascontiguousarray(array, dtype="<f8").tobytes()).hexdigest()


def metropolis(replica, refresh):
    start_source = PCG64Source(namespace_seed(NAMESPACE, "metropolis", "start", replica))
    start = sample_exact_finite_bead_paths(start_source, M, BETA, 2)[0]
    channels = StochasticChannels(PCG64Source(namespace_seed(NAMESPACE, "metropolis", "proposal", replica)),
                                  PCG64Source(namespace_seed(NAMESPACE, "metropolis", "acceptance", replica)),
                                  start_source)
    options = {"centroid_source": PCG64Source(namespace_seed(NAMESPACE, "refresh", "centroid", replica))} if refresh else {}
    model = PathIntegralHarmonicOscillator(beads=M, temperature=1/BETA, displacement_delta=DELTA)
    result = model.run_simulation(channels, sweeps=N, burn_in=0, initial_path=start,
                                  initial_path_provenance={"namespace": NAMESPACE, "replica": replica},
                                  uncertainty_batch_sizes=(512,), record_trace=True, **options)
    bm = result.virial_uncertainty["batch_means"]["512"]
    td, cv = harmonic_energy_estimators(result.trace["paths"], BETA)
    roles = result.source_provenance["roles"]
    return {"start_sha256": result.initialization_protocol["initial_path_sha256"],
            "start_matches_stored": result.initialization_protocol["initial_path_sha256"] == sha(start),
            "acceptance_rate": result.acceptance_rate, "cv_mean": bm["prefix_mean"], "bm512_mcse": bm["mcse"],
            "bm512_covered": bm["covered"],
            "cv_rederived_max_abs_diff": float(np.max(np.abs(cv - result.trace["CV"]))),
            "cv_trace_sha256": sha(result.trace["CV"]), "paths_sha256": sha(result.trace["paths"]),
            "bytes": {role: (None if entry is None else entry["consumed_bytes_delta"]) for role, entry in roles.items()},
            "aliases": sorted({a for entry in roles.values() if entry for a in entry["aliases"]})}


def iid(replica):
    source = PCG64Source(namespace_seed(NAMESPACE, "iid", "paths", replica))
    paths = sample_exact_finite_bead_paths(source, M, BETA, N)
    td, cv = harmonic_energy_estimators(paths, BETA)
    bm = summarize_series(cv, reference=exact_finite_bead_energy(M, BETA), batch_sizes=(512,))["batch_means"]["512"]
    return {"cv_mean": bm["prefix_mean"], "bm512_mcse": bm["mcse"], "bm512_covered": bm["covered"],
            "cv_trace_sha256": sha(cv), "bytes": {"paths": source.consumed_bytes}}


def timed(function, *args):
    started = time.perf_counter()
    out = function(*args)
    return {**out, "seconds": time.perf_counter() - started}


def main():
    runs = {"local": [timed(metropolis, i, False) for i in range(REPLICAS)],
            "refresh": [timed(metropolis, i, True) for i in range(REPLICAS)],
            "iid": [timed(iid, i) for i in range(REPLICAS)]}
    replay = {"local_0": timed(metropolis, 0, False)["cv_trace_sha256"] == runs["local"][0]["cv_trace_sha256"],
              "refresh_0": timed(metropolis, 0, True)["cv_trace_sha256"] == runs["refresh"][0]["cv_trace_sha256"],
              "iid_0": timed(iid, 0)["cv_trace_sha256"] == runs["iid"][0]["cv_trace_sha256"]}
    print(json.dumps({
        "script": "research/v018-api-validation/smoke/smoke_timing.py", "namespace": NAMESPACE,
        "design": {"M": M, "beta": BETA, "delta": DELTA, "n": N, "burn_in": 0, "replicas_per_condition": REPLICAS},
        "python": platform.python_version(), "numpy": np.__version__, "machine": platform.machine(),
        "mean_seconds_per_replica": {k: float(np.mean([r["seconds"] for r in v])) for k, v in runs.items()},
        "deterministic_replay_identical": replay,
        "peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
        "runs": runs,
    }, indent=1))


if __name__ == "__main__":
    sys.exit(main())
