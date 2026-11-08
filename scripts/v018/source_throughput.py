#!/usr/bin/env python3
"""Measure uint64 throughput of the rnglab14 reference sources on this host.

RNGLab v0.18 audit helper (written by Claude Code, 2026-11-05). Checks the
throughput column of docs/validation/SOURCE_COMPARISON.md against the
repository implementations. Each source reads 8-byte words through the audited
RandomSource.read_bytes boundary, the path every model uses. Bulk throughput
(one large read) is reported separately because it bypasses per-call overhead.

Usage: PYTHONPATH=src python scripts/v018/source_throughput.py > out.json
"""
import json
import platform
import statistics
import time

from rnglab14.sources import (ChaCha20Source, MersenneTwisterSource, PCG64Source,
                              Xoshiro256PlusPlusSource)

WORDS = 20_000
BULK_WORDS = 200_000
REPEATS = 3


def rate(factory, per_call: bool) -> float:
    times = []
    for _ in range(REPEATS):
        source = factory()
        start = time.perf_counter()
        if per_call:
            for _ in range(WORDS):
                source.read_bytes(8)
            count = WORDS
        else:
            source.read_bytes(8 * BULK_WORDS)
            count = BULK_WORDS
        times.append(count / (time.perf_counter() - start))
    return statistics.median(times)


def main() -> int:
    sources = {
        "PCG64Source": lambda: PCG64Source(seed=20261105),
        "Xoshiro256PlusPlusSource": lambda: Xoshiro256PlusPlusSource(seed=20261105),
        "MersenneTwisterSource": lambda: MersenneTwisterSource(seed=20261105),
        "ChaCha20Source(rounds=20)": lambda: ChaCha20Source(key=bytes(32), rounds=20),
        "ChaCha20Source(rounds=8)": lambda: ChaCha20Source(key=bytes(32), rounds=8),
    }
    result = {"script": "scripts/v018/source_throughput.py", "python": platform.python_version(),
              "machine": platform.machine(), "processor": platform.processor(),
              "repeats": REPEATS, "per_call_words": WORDS, "bulk_words": BULK_WORDS,
              "unit": "million uint64 words per second (median of repeats)", "sources": {}}
    for name, factory in sources.items():
        result["sources"][name] = {"per_call_8_byte_reads": rate(factory, True) / 1e6,
                                   "single_bulk_read": rate(factory, False) / 1e6}
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
