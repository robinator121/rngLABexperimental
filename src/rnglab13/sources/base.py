"""Common Randomness Source Abstraction for RNGLab v0.13+.

Provides a unified, streaming and offline interface for:
- Deterministic PRNGs (PCG64, MT19937, Xoshiro)
- Cryptographic OS entropy
- Synthetic TRNG-like streams (explicitly marked as synthetic)
- Recorded physical binary datasets
- Future hardware adapters (Serial framing, Pulse events)
- Conditioned and reseeded streams

Includes byte counters, health state, provenance metadata, and reproducibility manifests.
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
import hashlib
import json
import struct
import time
from typing import Any, Mapping


class SourceKind(str, Enum):
    DETERMINISTIC_PRNG = "deterministic_prng"
    OS_CRYPTOGRAPHIC = "os_cryptographic"
    SYNTHETIC_TRNG = "synthetic_trng"
    RECORDED_PHYSICAL = "recorded_physical"
    HARDWARE_SERIAL = "hardware_serial"
    HARDWARE_PULSE = "hardware_pulse"
    CONDITIONED = "conditioned"
    RESEEDING_CSPRNG = "reseeding_csprng"
    MIXED = "mixed"


class SourceHealth(str, Enum):
    HEALTHY = "healthy"
    DEGRADED = "degraded"
    FAILED = "failed"
    EXHAUSTED = "exhausted"


class SourceExhausted(Exception):
    """Raised when an offline or bounded source runs out of bytes."""
    pass


@dataclass
class SourceMetadata:
    source_id: str
    kind: SourceKind
    created_at: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    total_bytes_read: int = 0
    read_calls: int = 0
    reseeding_events: int = 0
    health: SourceHealth = SourceHealth.HEALTHY
    is_physical: bool = False
    is_synthetic: bool = False
    extra: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d["kind"] = self.kind.value
        d["health"] = self.health.value
        return d


class RandomSource(ABC):
    """Abstract base class for all RNGLab randomness sources."""

    def __init__(self, source_id: str, kind: SourceKind, is_physical: bool = False, is_synthetic: bool = False):
        self.metadata = SourceMetadata(
            source_id=source_id,
            kind=kind,
            is_physical=is_physical,
            is_synthetic=is_synthetic,
        )
        self._hasher = hashlib.sha256()

    @property
    def source_id(self) -> str:
        return self.metadata.source_id

    @property
    def kind(self) -> SourceKind:
        return self.metadata.kind

    @property
    def total_bytes_read(self) -> int:
        return self.metadata.total_bytes_read

    @property
    def health(self) -> SourceHealth:
        return self.metadata.health

    @abstractmethod
    def _read_bytes(self, n: int) -> bytes:
        """Internal implementation to produce exactly n bytes."""
        raise NotImplementedError

    def read_bytes(self, n: int) -> bytes:
        """Read exactly n random bytes, tracking accounting and cumulative hash."""
        if n < 0:
            raise ValueError(f"Byte count must be non-negative: {n}")
        if n == 0:
            return b""
        if self.metadata.health == SourceHealth.EXHAUSTED:
            raise SourceExhausted(f"Source {self.source_id!r} is exhausted")

        data = self._read_bytes(n)
        if len(data) != n:
            self.metadata.health = SourceHealth.EXHAUSTED
            raise SourceExhausted(f"Requested {n} bytes but got {len(data)} from {self.source_id!r}")

        self.metadata.total_bytes_read += n
        self.metadata.read_calls += 1
        self._hasher.update(data)
        return data

    def read_u32(self) -> int:
        """Read a 32-bit unsigned integer (big-endian)."""
        data = self.read_bytes(4)
        return struct.unpack(">I", data)[0]

    def read_u64(self) -> int:
        """Read a 64-bit unsigned integer (big-endian)."""
        data = self.read_bytes(8)
        return struct.unpack(">Q", data)[0]

    def random_float(self) -> float:
        """Read a standard uniform float in [0, 1) using 53 bits of precision."""
        val = self.read_u64() >> 11  # 53 bits
        return val * (1.0 / (1 << 53))

    def reseed(self, entropy: bytes) -> None:
        """Reseed the source if supported."""
        raise NotImplementedError(f"Reseeding is not supported by {self.__class__.__name__}")

    def reproducibility_manifest(self) -> dict[str, Any]:
        """Return a reproducibility manifest with exact byte count and SHA-256 digest."""
        return {
            "source_id": self.source_id,
            "kind": self.kind.value,
            "total_bytes_read": self.total_bytes_read,
            "read_calls": self.metadata.read_calls,
            "stream_sha256": self._hasher.hexdigest(),
            "health": self.health.value,
            "metadata": self.metadata.to_dict(),
        }
