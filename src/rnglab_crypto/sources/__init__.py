"""Byte-source adapters sharing RNGLab14's audited consumption boundary."""
from __future__ import annotations

import hashlib
import json
import os
from typing import Protocol

from rnglab11.drbg import HMACDRBG
from rnglab14.sources.base import RandomSource, SourceManifest


class ByteSource(Protocol):
    def read_bytes(self, n: int) -> bytes: ...


def seed_for(namespace: str, replicate: int, role: str) -> bytes:
    """Unambiguous, public, domain-separated replay material; NOT entropy."""
    if not namespace or type(replicate) is not int or replicate < 0 or not role:
        raise ValueError("namespace, nonnegative replicate and role required")
    return hashlib.sha256(json.dumps(["BLUE-v1", namespace, replicate, role],
                                    separators=(",", ":")).encode()).digest()


class ReplaySource(RandomSource):
    """Synthetic secret-unexposed reference model; publicly replayable in reality.

    Each read is one HMAC-DRBG Generate call, so the request schedule is part
    of replay provenance. This adapter confers zero validated entropy credit.
    """
    def __init__(self, seed: bytes, *, max_bytes: int = 1_048_576):
        if len(seed) != 32:
            raise ValueError("32 byte public replay seed required")
        self.seed = bytes(seed)
        self.drbg = HMACDRBG(seed, b"BLUE replay v1!!", max_bytes=max_bytes,
                             security_strength=128)
        manifest = SourceManifest(hashlib.sha256(seed).hexdigest(), "synthetic_replay",
                                  True, False, created_at="deterministic-replay",
                                  metadata={"seed_hex": seed.hex(), "entropy_credit_bits": 0,
                                            "request_schedule_matters": True})
        super().__init__(manifest, max_bytes=max_bytes)

    def _read_bytes(self, n):
        return self.drbg.read(n)

    def restart(self):
        return ReplaySource(self.seed, max_bytes=self.max_bytes)


class OSReference(RandomSource):
    """Read-only OS CSPRNG statistical reference; no interception or replacement."""
    def __init__(self, max_bytes=1_048_576):
        super().__init__(SourceManifest("host-os", "os_csprng", False, False,
                         metadata={"replayable": False, "physical_source_validated": False}),
                         max_bytes=max_bytes)

    def _read_bytes(self, n):
        return os.urandom(n)


class CapturedSource(RandomSource):
    """Finite read-only capture. Public physical recordings are replay, not fresh entropy."""
    def __init__(self, data: bytes, *, provenance: str, physical: bool = False):
        if not provenance:
            raise ValueError("capture provenance required")
        self.data = bytes(data)
        super().__init__(SourceManifest(hashlib.sha256(data).hexdigest(), "capture_replay",
                         True, physical, created_at="recorded-in-provenance",
                         metadata={"provenance": provenance, "entropy_credit_bits": 0}),
                         max_bytes=len(data))

    def _read_bytes(self, n):
        return self.data[self.consumed_bytes:self.consumed_bytes+n]
