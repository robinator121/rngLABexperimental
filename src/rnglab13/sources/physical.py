"""Synthetic TRNG-like sources and recorded physical binary dataset streams."""
from __future__ import annotations

import hashlib
import math
import os
from pathlib import Path
import struct
import time
from typing import Any

from .base import RandomSource, SourceExhausted, SourceHealth, SourceKind


class SyntheticTRNGSource(RandomSource):
    """Simulates physical entropy processes (e.g. thermal noise, phase jitter).

    IMPORTANT: This source is strictly synthetic and MUST NEVER be labeled
    as genuine physical randomness in scientific reports.
    """

    def __init__(
        self,
        seed: int = 42,
        jitter_std: float = 1.0,
        thermal_drift_rate: float = 0.0,
        source_id: str = "synthetic-physical-trng",
    ):
        super().__init__(
            source_id=source_id,
            kind=SourceKind.SYNTHETIC_TRNG,
            is_physical=False,
            is_synthetic=True,
        )
        self.seed = seed
        self.jitter_std = jitter_std
        self.thermal_drift_rate = thermal_drift_rate
        self._state = seed
        self._step = 0
        self.metadata.extra.update({
            "seed": seed,
            "jitter_std": jitter_std,
            "thermal_drift_rate": thermal_drift_rate,
            "scientific_notice": "SYNTHETIC SIMULATION ONLY; NOT A PHYSICAL SOURCE",
        })

    def _next_sample(self) -> int:
        """Simulate single noise voltage sample quantized to an 8-bit byte."""
        self._step += 1
        # LCG core for underlying noise
        self._state = (self._state * 6364136223846793005 + 1442695040888963407) & 0xFFFFFFFFFFFFFFFF
        # Gaussian-like jitter via central limit sum of uniforms
        norm = sum(((self._state >> (i * 8)) & 0xFF) for i in range(4)) / 1020.0 - 0.5
        # Thermal drift
        drift = self.thermal_drift_rate * math.sin(self._step * 0.001)
        raw_val = int(128 + (norm * self.jitter_std * 64) + drift)
        return max(0, min(255, raw_val))

    def _read_bytes(self, n: int) -> bytes:
        return bytes(self._next_sample() for _ in range(n))


class BinaryFileSource(RandomSource):
    """Streams random bytes from an offline recorded binary file (e.g. lab dataset)."""

    def __init__(self, file_path: str | Path, source_id: str = "binary-file-source", chunk_size: int = 65536):
        super().__init__(
            source_id=source_id,
            kind=SourceKind.RECORDED_PHYSICAL,
            is_physical=True,
            is_synthetic=False,
        )
        self.file_path = Path(file_path).resolve(strict=True)
        self.chunk_size = chunk_size
        self._file = open(self.file_path, "rb")
        self.file_size = self.file_path.stat().st_size
        self.metadata.extra["file_path"] = str(self.file_path)
        self.metadata.extra["file_size_bytes"] = self.file_size

    def __del__(self):
        if hasattr(self, "_file") and not self._file.closed:
            self._file.close()

    def _read_bytes(self, n: int) -> bytes:
        data = self._file.read(n)
        if len(data) < n:
            self.metadata.health = SourceHealth.EXHAUSTED
            raise SourceExhausted(f"End of binary file reached: needed {n} bytes, got {len(data)}")
        return data
