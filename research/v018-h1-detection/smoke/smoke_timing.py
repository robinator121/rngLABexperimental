#!/usr/bin/env python3
"""Runtime smoke for the draft Stage C (H1) design. No campaign data.

RNGLab v0.18 (written by Claude Code, 2026-11-06). Separate namespace
rnglab-v018-h1-detection-smoke; nothing here counts toward any campaign.
Prototypes, for timing only:
- one local source (proposal and acceptance share it) wrapped so that, after
  an onset byte offset, every 8-byte word is transformed (shift, lag-1 Gaussian
  copula or word duplication), with all emitted bytes recorded;
- a "sweep clock": the healthy centroid source is read exactly once per sweep,
  after the local updates, so wrapping it records the local byte offset of
  every sweep boundary;
- the in-repo SP 800-22 three-test subset on rolling byte windows aligned to
  512-sweep windows (about 0.87 Mbit; the in-repo implementation caps one
  call at 1,000,000 bits), stride 256 sweeps.

Usage: OPENBLAS_NUM_THREADS=1 PYTHONPATH=src python research/v018-h1-detection/smoke/smoke_timing.py > out.json
"""
import hashlib
import json
import math
import platform
import resource
import struct
import time
from statistics import NormalDist

import numpy as np

from rnglab11.statistics import evaluate_bits
from rnglab14.models.pimc_harmonic import PathIntegralHarmonicOscillator, sample_exact_finite_bead_paths
from rnglab14.sources import PCG64Source
from rnglab14.sources.base import RandomSource, SourceManifest

NS = "rnglab-v018-h1-detection-smoke"
M, BETA, DELTA, SWEEPS, ONSET_BYTES = 16, 2.0, 0.4, 10240, 2**19
WINDOW, STRIDE = 512, 256  # battery window: rnglab11.evaluate_bits accepts at most 1,000,000 bits per call
STD = NormalDist()


def seed(*parts):
    return int.from_bytes(hashlib.sha256("/".join([NS, *map(str, parts)]).encode()).digest()[:8], "big")


class FaultyWords(RandomSource):
    def __init__(self, base, family, severity, fault_seed):
        super().__init__(SourceManifest(source_id=f"smoke_{family}", family="synthetic_fault",
                                        is_deterministic=True, is_physical=False))
        self.base, self.family, self.s = base, family, severity
        self.emitted = bytearray()
        self.rng = np.random.default_rng(fault_seed)
        self.prev, self.z = None, 0.0

    def _read_bytes(self, n):
        raw = self.base.read_bytes(n)
        if n == 8 and self.family != "healthy" and len(self.emitted) >= ONSET_BYTES:
            word = struct.unpack(">Q", raw)[0]
            u = (word >> 11) * 2.0**-53
            if self.family == "shift":
                u = min(u + self.s, 1 - 2.0**-53)
            elif self.family == "lag1":
                self.z = self.s * self.z + math.sqrt(1 - self.s * self.s) * STD.inv_cdf(u + 2.0**-54)
                u = min(STD.cdf(self.z), 1 - 2.0**-53)
            if self.family == "dup" and self.prev is not None and self.rng.random() < self.s:
                raw = self.prev
            else:
                raw = struct.pack(">Q", (int(u * 2.0**53) << 11) | (word & 0x7FF))
        self.prev = raw
        self.emitted += raw
        return raw


class SweepClock(RandomSource):
    def __init__(self, base, local):
        super().__init__(SourceManifest(source_id="clock", family="reference_prng", is_deterministic=True,
                                        is_physical=False))
        self.base, self.local, self.marks = base, local, []

    def _read_bytes(self, n):
        self.marks.append(self.local.consumed_bytes)
        return self.base.read_bytes(n)


def one_run(family, severity, i):
    local = FaultyWords(PCG64Source(seed("local", i)), family, severity, seed("fault", i))
    clock = SweepClock(PCG64Source(seed("centroid", i)), local)
    start = sample_exact_finite_bead_paths(PCG64Source(seed("start", i)), M, BETA, 2)[0]
    model = PathIntegralHarmonicOscillator(beads=M, temperature=1 / BETA, displacement_delta=DELTA)
    t0 = time.perf_counter()
    result = model.run_simulation(local, sweeps=SWEEPS, burn_in=0, initial_path=start, centroid_source=clock,
                                  uncertainty_batch_sizes=(512,), record_trace=True)
    t_run = time.perf_counter() - t0
    marks = np.array(clock.marks + [local.consumed_bytes])  # bytes consumed before each sweep's refresh
    bits = np.unpackbits(np.frombuffer(bytes(local.emitted), dtype=np.uint8))
    t0 = time.perf_counter()
    minp = []
    for end in range(WINDOW, SWEEPS + 1, STRIDE):
        lo, hi = marks[end - WINDOW], marks[end - 1]
        tests = evaluate_bits(bits[8 * lo:8 * hi], block_size=max(128, (8 * (hi - lo)) // 99 + 1))["tests"]
        minp.append(min((t["p_value"] for t in tests.values() if t.get("p_value") is not None), default=1.0))
    t_battery = time.perf_counter() - t0
    d = result.trace["TD"] - result.trace["CV"]
    onset_sweep = int(np.searchsorted(marks, ONSET_BYTES))
    return {"family": family, "severity": severity, "replica": i, "run_seconds": t_run, "battery_seconds": t_battery,
            "local_bytes": local.consumed_bytes, "bytes_per_sweep": local.consumed_bytes / SWEEPS,
            "onset_sweep": onset_sweep, "windows": len(minp), "min_window_p": float(min(minp)),
            "td_minus_cv_mean_post_onset": float(d[onset_sweep:].mean()), "acceptance_rate_4dp": result.acceptance_rate}


def main():
    runs = [one_run(f, s, i) for f, s in (("healthy", 0.0), ("shift", 0.01), ("lag1", 0.1), ("dup", 0.05)) for i in range(2)]
    print(json.dumps({"script": "research/v018-h1-detection/smoke/smoke_timing.py", "namespace": NS,
                      "design": {"M": M, "beta": BETA, "delta": DELTA, "sweeps": SWEEPS, "onset_bytes": ONSET_BYTES,
                                 "window_sweeps": WINDOW, "stride_sweeps": STRIDE},
                      "python": platform.python_version(), "numpy": np.__version__,
                      "peak_rss_mib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, "runs": runs}, indent=1))


if __name__ == "__main__":
    main()
