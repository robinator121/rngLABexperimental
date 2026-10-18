"""RNGLab v0.14: Reference PRNG and CSPRNG Generators.

Implements standard statistical, cryptographic, and round-reduced stream generators:
- MT19937 (Mersenne Twister)
- PCG64 (Permuted Congruential Generator)
- Xoshiro256++ (Modern high-speed PRNG)
- ChaCha20 / ChaCha8 (Cryptographic stream generator with selectable round degradation)
- System CSPRNG (OS entropy pool)
"""
from __future__ import annotations

import hashlib
import os
import random
import struct
from typing import Optional

from .base import RandomSource, SourceManifest


class MersenneTwisterSource(RandomSource):
    """Deterministic reference generator using 64-bit Mersenne Twister (MT19937)."""

    def __init__(self, seed: int = 42, max_bytes: Optional[int] = None):
        manifest = SourceManifest(
            source_id=f"mt19937_s{seed}",
            family="reference_prng",
            is_deterministic=True,
            is_physical=False,
            seed=seed,
            metadata={"generator": "MT19937", "period": "2^19937-1"},
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self._rng = random.Random(seed)

    def _read_bytes(self, n: int) -> bytes:
        return self._rng.randbytes(n)


class PCG64Source(RandomSource):
    """Permuted Congruential Generator (PCG-XSL-RR 128/64) for statistical benchmarking."""

    def __init__(self, seed: int = 42, stream_id: int = 1, max_bytes: Optional[int] = None):
        manifest = SourceManifest(
            source_id=f"pcg64_s{seed}_str{stream_id}",
            family="reference_prng",
            is_deterministic=True,
            is_physical=False,
            seed=seed,
            metadata={"generator": "PCG64", "stream_id": stream_id},
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self.state = (seed & 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF)
        self.inc = ((stream_id << 1) | 1) & 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        self._buffer = bytearray()
        # Initial step
        self._step()

    def _step(self) -> int:
        # LCG 128-bit step: state = state * multiplier + increment
        mult = 0x2360ED051FC65DA44385DF649FCCF645
        self.state = (self.state * mult + self.inc) & 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        # XSL-RR 128/64 output permutation
        hi = self.state >> 64
        lo = self.state & 0xFFFFFFFFFFFFFFFF
        xored = hi ^ lo
        rot = (self.state >> 122) & 63
        # Rotate right 64-bit
        out = ((xored >> rot) | (xored << ((64 - rot) & 63))) & 0xFFFFFFFFFFFFFFFF
        return out

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            val = self._step()
            self._buffer.extend(struct.pack(">Q", val))
        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out


class Xoshiro256PlusPlusSource(RandomSource):
    """Xoshiro256++ generator by Blackman & Vigna (2019)."""

    def __init__(self, seed: int = 42, max_bytes: Optional[int] = None):
        manifest = SourceManifest(
            source_id=f"xoshiro256pp_s{seed}",
            family="reference_prng",
            is_deterministic=True,
            is_physical=False,
            seed=seed,
            metadata={"generator": "Xoshiro256++"},
        )
        super().__init__(manifest, max_bytes=max_bytes)
        # Initialize 4 x 64-bit state using SplitMix64 on seed
        s = seed & 0xFFFFFFFFFFFFFFFF
        state = []
        for _ in range(4):
            s = (s + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
            z = s
            z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
            z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
            z = (z ^ (z >> 31)) & 0xFFFFFFFFFFFFFFFF
            state.append(z)
        if all(x == 0 for x in state):
            state[0] = 1
        self.s = state
        self._buffer = bytearray()

    def _step(self) -> int:
        s0, s1, s2, s3 = self.s
        # Output function: rotl(s0 + s3, 23) + s0
        res = ((((s0 + s3) & 0xFFFFFFFFFFFFFFFF) << 23) | (((s0 + s3) & 0xFFFFFFFFFFFFFFFF) >> 41)) & 0xFFFFFFFFFFFFFFFF
        res = (res + s0) & 0xFFFFFFFFFFFFFFFF

        t = (s1 << 17) & 0xFFFFFFFFFFFFFFFF
        s2 ^= s0
        s3 ^= s1
        s1 ^= s2
        s0 ^= s3
        s2 ^= t
        # rotl(s3, 45)
        self.s[3] = ((s3 << 45) | (s3 >> 19)) & 0xFFFFFFFFFFFFFFFF
        self.s[0] = s0
        self.s[1] = s1
        self.s[2] = s2
        return res

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            val = self._step()
            self._buffer.extend(struct.pack(">Q", val))
        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out


class ChaCha20Source(RandomSource):
    """ChaCha stream cipher generator with configurable round count (e.g. 8, 12, or 20).
    
    Allows controlled degradation studies of cryptographic diffusion.
    """

    def __init__(self, key: Optional[bytes] = None, rounds: int = 20, max_bytes: Optional[int] = None):
        if rounds % 2 != 0 or rounds < 2:
            raise ValueError("ChaCha rounds must be an even integer >= 2")
        if key is None:
            key = b"\x01" * 32
        if len(key) != 32:
            raise ValueError("ChaCha20 key must be exactly 32 bytes")

        manifest = SourceManifest(
            source_id=f"chacha_{rounds}rounds",
            family="csprng",
            is_deterministic=True,
            is_physical=False,
            metadata={"rounds": rounds, "key_sha256": hashlib.sha256(key).hexdigest()[:16]},
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self.key = key
        self.rounds = rounds
        self.counter = 0
        self.nonce = b"\x00" * 12
        self._buffer = bytearray()

    def _quarter_round(self, x: list[int], a: int, b: int, c: int, d: int):
        mask = 0xFFFFFFFF
        x[a] = (x[a] + x[b]) & mask
        x[d] = ((x[d] ^ x[a]) << 16 | (x[d] ^ x[a]) >> 16) & mask
        x[c] = (x[c] + x[d]) & mask
        x[b] = ((x[b] ^ x[c]) << 12 | (x[b] ^ x[c]) >> 20) & mask
        x[a] = (x[a] + x[b]) & mask
        x[d] = ((x[d] ^ x[a]) << 8 | (x[d] ^ x[a]) >> 24) & mask
        x[c] = (x[c] + x[d]) & mask
        x[b] = ((x[b] ^ x[c]) << 7 | (x[b] ^ x[c]) >> 25) & mask

    def _block(self) -> bytes:
        constants = [0x61707865, 0x3320646E, 0x79622D32, 0x6B206574]
        key_words = list(struct.unpack("<8I", self.key))
        counter_words = [self.counter & 0xFFFFFFFF]
        nonce_words = list(struct.unpack("<3I", self.nonce))
        state = constants + key_words + counter_words + nonce_words
        initial = list(state)

        for _ in range(0, self.rounds, 2):
            # Column round
            self._quarter_round(state, 0, 4, 8, 12)
            self._quarter_round(state, 1, 5, 9, 13)
            self._quarter_round(state, 2, 6, 10, 14)
            self._quarter_round(state, 3, 7, 11, 15)
            # Diagonal round
            self._quarter_round(state, 0, 5, 10, 15)
            self._quarter_round(state, 1, 6, 11, 12)
            self._quarter_round(state, 2, 7, 8, 13)
            self._quarter_round(state, 3, 4, 9, 14)

        out_words = [(state[i] + initial[i]) & 0xFFFFFFFF for i in range(16)]
        self.counter += 1
        return struct.pack("<16I", *out_words)

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            self._buffer.extend(self._block())
        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out


class OSEntropySource(RandomSource):
    """Non-deterministic OS system entropy pool (os.urandom)."""

    def __init__(self, max_bytes: Optional[int] = None):
        manifest = SourceManifest(
            source_id="os_urandom",
            family="csprng",
            is_deterministic=False,
            is_physical=False,
            metadata={"device": "/dev/urandom"},
        )
        super().__init__(manifest, max_bytes=max_bytes)

    def _read_bytes(self, n: int) -> bytes:
        return os.urandom(n)
