#!/usr/bin/env python3
"""Runtime and memory smoke for the Phase 1 estimator campaign design.

RNGLab v0.18 (written by Claude Code, 2026-11-05). Uses the separate smoke
namespace rnglab-v018-estimators-smoke; nothing produced here counts toward
the campaign. Measures (a) the v0.16 vectorized Ising kernel at L=4, T=1.5,
2048 sweeps, and (b) AR(1) generation plus FFT/cumsum work comparable to the
planned long-run variance estimators at n=8192 and n=65536.

Usage: OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python research/v018-estimators/smoke/smoke_timing.py > out.json
"""
import hashlib
import json
import platform
import resource
import time

import numpy as np
from scipy.signal import lfilter

from rnglab16.core import TC, ising

NAMESPACE = "rnglab-v018-estimators-smoke"


def seed(*parts) -> int:
    return int.from_bytes(hashlib.sha256("/".join([NAMESPACE, *map(str, parts)]).encode()).digest()[:8], "big")


def peak_mib() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def time_ising(L: int, T: float, reps: int, steps: int) -> dict:
    start = time.perf_counter()
    trace, _, _ = ising(L, T, reps, steps, seed("ising", L, T, reps))
    elapsed = time.perf_counter() - start
    return {"L": L, "T": T, "reps": reps, "steps": steps, "seconds": elapsed,
            "seconds_per_chain": elapsed / reps, "peak_rss_MiB_after": peak_mib(),
            "trace_shape": list(trace.shape)}


def time_ar1(reps: int, n: int, phi: float) -> dict:
    rng = np.random.default_rng(seed("ar1", reps, n, phi))
    start = time.perf_counter()
    eps = rng.standard_normal((reps, n))
    eps[:, 0] /= np.sqrt(1 - phi * phi) if phi < 1 else 1.0
    x = lfilter([np.sqrt(1 - phi * phi)], [1, -phi], eps, axis=1)
    generated = time.perf_counter()
    centered = x - x.mean(axis=1, keepdims=True)
    size = 1 << (2 * n - 1).bit_length()
    spectrum = np.fft.rfft(centered, n=size, axis=1)
    acov = np.fft.irfft(spectrum * spectrum.conj(), n=size, axis=1)[:, :n] / n
    csum = np.concatenate([np.zeros((reps, 1)), np.cumsum(x, axis=1)], axis=1)
    for b in (32, 128, 512, 1024):
        window_means = (csum[:, b:] - csum[:, :-b]) / b
        _ = ((window_means - x.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    end = time.perf_counter()
    return {"reps": reps, "n": n, "phi": phi, "generate_seconds": generated - start,
            "estimator_work_seconds": end - generated, "acov_lag1_mean": float(acov[:, 1].mean() / acov[:, 0].mean()),
            "peak_rss_MiB_after": peak_mib()}


if __name__ == "__main__":
    out = {"script": "research/v018-estimators/smoke/smoke_timing.py", "namespace": NAMESPACE,
           "python": platform.python_version(), "numpy": np.__version__,
           "classification": "smoke only; excluded from every campaign endpoint"}
    out["ising"] = [time_ising(4, 1.5, 256, 2048), time_ising(4, float(TC), 256, 2048),
                    time_ising(2, 1.5, 256, 2048)]
    out["ar1"] = [time_ar1(512, 8192, 0.99), time_ar1(64, 65536, 0.997)]
    print(json.dumps(out, indent=2))
