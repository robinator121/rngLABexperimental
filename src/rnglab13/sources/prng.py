"""Deterministic PRNG implementations for RNGLab v0.13+."""
from __future__ import annotations

import random
import struct
from typing import Any

from .base import RandomSource, SourceKind


class MersenneTwisterSource(RandomSource):
    """Deterministic PRNG backed by Mersenne Twister (MT19937)."""

    def __init__(self, seed: int, source_id: str = "prng-mt19937"):
        super().__init__(source_id=source_id, kind=SourceKind.DETERMINISTIC_PRNG)
        self.seed = seed
        self._rng = random.Random(seed)
        self.metadata.extra["seed"] = seed

    def _read_bytes(self, n: int) -> bytes:
        return self._rng.randbytes(n)

    def reseed(self, entropy: bytes) -> None:
        new_seed = int.from_bytes(entropy, "big")
        self._rng.seed(new_seed)
        self.metadata.reseeding_events += 1


class Xoshiro256PlusPlusSource(RandomSource):
    """Deterministic PRNG backed by Xoshiro256++ 64-bit generator."""

    def __init__(self, seed: int, source_id: str = "prng-xoshiro256pp"):
        super().__init__(source_id=source_id, kind=SourceKind.DETERMINISTIC_PRNG)
        self.seed = seed
        # Initialize 4x 64-bit state words using SplitMix64 on seed
        self._s = self._seed_splitmix64(seed)
        self.metadata.extra["seed"] = seed

    @staticmethod
    def _seed_splitmix64(x: int) -> list[int]:
        s = []
        cur = x & 0xFFFFFFFFFFFFFFFF
        for _ in range(4):
            cur = (cur + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
            z = cur
            z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
            z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
            z = (z ^ (z >> 31)) & 0xFFFFFFFFFFFFFFFF
            s.append(z)
        if all(word == 0 for word in s):
            s[0] = 1
        return s

    @staticmethod
    def _rotl(x: int, k: int) -> int:
        return ((x << k) & 0xFFFFFFFFFFFFFFFF) | (x >> (64 - k))

    def _next_u64(self) -> int:
        s = self._s
        result = (self._rotl((s[0] + s[3]) & 0xFFFFFFFFFFFFFFFF, 23) + s[0]) & 0xFFFFFFFFFFFFFFFF
        t = (s[1] << 17) & 0xFFFFFFFFFFFFFFFF

        s[2] ^= s[0]
        s[3] ^= s[1]
        s[1] ^= s[2]
        s[0] ^= s[3]

        s[2] ^= t
        s[3] = self._rotl(s[3], 45)

        return result

    def _read_bytes(self, n: int) -> bytes:
        chunks = []
        needed = n
        while needed > 0:
            val = self._next_u64()
            chunk = struct.pack(">Q", val)
            if needed >= 8:
                chunks.append(chunk)
                needed -= 8
            else:
                chunks.append(chunk[:needed])
                needed = 0
        return b"".join(chunks)
