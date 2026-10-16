"""Operating-system cryptographic randomness and periodically reseeded CSPRNGs."""
from __future__ import annotations

import hashlib
import hmac
import secrets
import struct
from typing import Any

from .base import RandomSource, SourceKind


class OSCryptoSource(RandomSource):
    """Host operating-system cryptographic randomness backed by /dev/urandom or getrandom()."""

    def __init__(self, source_id: str = "crypto-os"):
        super().__init__(source_id=source_id, kind=SourceKind.OS_CRYPTOGRAPHIC)

    def _read_bytes(self, n: int) -> bytes:
        return secrets.token_bytes(n)


class ReseedingCSPRNGSource(RandomSource):
    """CSPRNG based on HMAC-SHA256 counter mode with periodic automatic reseeding from an entropy source."""

    def __init__(
        self,
        entropy_source: RandomSource,
        reseed_interval_bytes: int = 65536,
        source_id: str = "crypto-reseeding-csprng",
    ):
        super().__init__(source_id=source_id, kind=SourceKind.RESEEDING_CSPRNG)
        self.entropy_source = entropy_source
        self.reseed_interval_bytes = reseed_interval_bytes
        self.bytes_since_reseed = 0
        self._counter = 0

        # Initial seed
        initial_entropy = self.entropy_source.read_bytes(32)
        self._key = hashlib.sha256(initial_entropy + b"initial_seed").digest()
        self.metadata.extra["reseed_interval_bytes"] = reseed_interval_bytes
        self.metadata.extra["entropy_source_id"] = entropy_source.source_id

    def reseed(self, entropy: bytes | None = None) -> None:
        """Reseed the internal HMAC state with fresh entropy."""
        fresh_entropy = entropy if entropy is not None else self.entropy_source.read_bytes(32)
        self._key = hmac.new(self._key, fresh_entropy + struct.pack(">Q", self._counter), hashlib.sha256).digest()
        self.bytes_since_reseed = 0
        self.metadata.reseeding_events += 1

    def _read_bytes(self, n: int) -> bytes:
        chunks = []
        needed = n
        while needed > 0:
            if self.bytes_since_reseed >= self.reseed_interval_bytes:
                self.reseed()

            # Generate block: HMAC(K, counter)
            block = hmac.new(self._key, struct.pack(">Q", self._counter), hashlib.sha256).digest()
            self._counter += 1

            take = min(needed, len(block))
            chunks.append(block[:take])
            needed -= take
            self.bytes_since_reseed += take

        return b"".join(chunks)
