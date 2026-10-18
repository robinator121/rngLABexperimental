"""RNGLab v0.14: Core Randomness Source Abstractions & Stochastic Channels.

Provides the interchangeable base interface for all PRNGs, CSPRNGs, physical adapters,
and conditioners, along with multi-channel stochastic decoupling (proposal vs acceptance).
"""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass, field
import datetime
import hashlib
import struct
from typing import Any, Dict, List, Optional, Sequence


@dataclass(frozen=True)
class SourceManifest:
    source_id: str
    family: str  # e.g., 'reference_prng', 'csprng', 'physical_gm', 'physical_optical', 'conditioned'
    is_deterministic: bool
    is_physical: bool
    seed: Optional[int] = None
    created_at: str = field(default_factory=lambda: datetime.datetime.now(datetime.timezone.utc).isoformat())
    metadata: Dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class SourceBudgetExceeded(Exception):
    """Raised when a randomness consumer attempts to read beyond allocated byte budget."""
    pass


class RandomSource(ABC):
    """Abstract base class for all interchangeable randomness sources in RNGLab v0.14."""

    def __init__(self, manifest: SourceManifest, max_bytes: Optional[int] = None):
        if max_bytes is not None and (isinstance(max_bytes,bool) or not isinstance(max_bytes,int) or max_bytes<0):
            raise ValueError("max_bytes must be None or a nonnegative integer")
        self.manifest = manifest
        self.max_bytes = max_bytes
        self.consumed_bytes = 0
        self._sha256 = hashlib.sha256()

    def read_bytes(self, n: int) -> bytes:
        """Single audited consumption boundary (provenance format v2).

        Adapters implement _read_bytes. Every successful direct/float/bit/int read
        hashes the exact emitted bytes once at this layer. This is a digest of
        consumed output, not a claim about entropy or independent source origin.
        """
        if isinstance(n, bool) or not isinstance(n, int) or n < 0:
            raise ValueError("byte count must be a nonnegative integer")
        if self.max_bytes is not None and self.consumed_bytes + n > self.max_bytes:
            raise SourceBudgetExceeded("requested bytes exceed source budget")
        if n == 0:
            return b""
        raw = self._read_bytes(n)
        if not isinstance(raw, bytes) or len(raw) != n:
            raise ValueError("source adapter must return exactly n bytes")
        self.consumed_bytes += n
        self._sha256.update(raw)
        return raw

    @abstractmethod
    def _read_bytes(self, n: int) -> bytes:
        """Adapter implementation; account/hash only in read_bytes."""
        raise NotImplementedError

    def _check_and_account_budget(self, n: int):
        if self.max_bytes is not None and (self.consumed_bytes + n) > self.max_bytes:
            raise SourceBudgetExceeded(
                f"Source '{self.manifest.source_id}' exceeded byte budget: "
                f"requested {n} bytes, consumed {self.consumed_bytes}/{self.max_bytes}"
            )
        self.consumed_bytes += n

    def random_float(self) -> float:
        """Generate a double-precision float uniformly distributed in [0, 1)."""
        raw = self.read_bytes(8)
        val = struct.unpack(">Q", raw)[0]
        # Standard IEEE-754 53-bit resolution mantissa conversion
        return (val >> 11) * (1.0 / 9007199254740992.0)

    def read_uniform(self, low: float = 0.0, high: float = 1.0) -> float:
        """Generate float in [low, high)."""
        return low + (high - low) * self.random_float()

    def read_bits(self, n: int) -> List[int]:
        """Read n bits from byte-aligned consumption; discard final padding bits.

        Repeated short calls are not a contiguous bit stream: e.g. read_bits(1)
        advances8 physical bits per returned bit. Buffer bytes for contiguous
        bit diagnostics; Markov rho refers to contiguous emitted bits.
        """
        if isinstance(n, bool) or not isinstance(n, int) or n < 0:
            raise ValueError("bit count must be a nonnegative integer")
        num_bytes = (n + 7) // 8
        raw = self.read_bytes(num_bytes)
        bits = []
        for b in raw:
            for shift in range(7, -1, -1):
                if len(bits) < n:
                    bits.append((b >> shift) & 1)
        return bits

    def read_int(self, low: int, high: int) -> int:
        """Generate an integer uniformly in [low, high) without modulo bias."""
        if low >= high:
            raise ValueError(f"low ({low}) must be strictly less than high ({high})")
        range_size = high - low
        if (range_size & (range_size - 1)) == 0:
            # Power of two range
            bits_needed = range_size.bit_length() - 1
            if bits_needed == 0:
                return low
            num_bytes = (bits_needed + 7) // 8
            raw = self.read_bytes(num_bytes)
            val = int.from_bytes(raw, byteorder="big") & (range_size - 1)
            return low + val

        # Rejection sampling to prevent modulo bias
        mask = (1 << range_size.bit_length()) - 1
        num_bytes = (range_size.bit_length() + 7) // 8
        while True:
            raw = self.read_bytes(num_bytes)
            candidate = int.from_bytes(raw, byteorder="big") & mask
            if candidate < range_size:
                return low + candidate

    @property
    def stream_digest(self) -> str:
        """SHA-256 hash of all consumed bytes for complete stream provenance."""
        return self._sha256.hexdigest()


@dataclass
class StochasticChannels:
    """Names stochastic channels; distinct references do not prove independence.
    
    Crucial for identifying whether randomness defects compromise proposal generation,
    acceptance filtering, initial configurations, or resampling procedures.
    """
    proposal: RandomSource
    acceptance: RandomSource
    initialization: RandomSource
    resampling: Optional[RandomSource] = None

    def to_manifest_dict(self) -> Dict[str, Any]:
        return {
            "proposal": self.proposal.manifest.to_dict(),
            "acceptance": self.acceptance.manifest.to_dict(),
            "initialization": self.initialization.manifest.to_dict(),
            "resampling": self.resampling.manifest.to_dict() if self.resampling else None,
        }
