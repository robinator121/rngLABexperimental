"""Conditioning and multi-source mixing implementations for RNGLab v0.13+."""
from __future__ import annotations

import collections
import hashlib
from typing import Sequence

from .base import RandomSource, SourceExhausted, SourceKind


class VonNeumannConditioner(RandomSource):
    """Classic John von Neumann unbiasing conditioner.

    Processes bits in pairs:
    01 -> 0
    10 -> 1
    00 or 11 -> discarded
    Eliminates constant stationary bit bias at the cost of variable output rate.
    """

    def __init__(self, raw_source: RandomSource, source_id: str = "conditioned-von-neumann"):
        super().__init__(
            source_id=source_id,
            kind=SourceKind.CONDITIONED,
            is_physical=raw_source.metadata.is_physical,
            is_synthetic=raw_source.metadata.is_synthetic,
        )
        self.raw_source = raw_source
        self._bit_queue = collections.deque()
        self._raw_bit_queue = collections.deque()
        self.metadata.extra["underlying_source"] = raw_source.source_id
        self.metadata.extra["conditioning_method"] = "von_neumann"

    def _fill_raw_bits(self, count: int = 128) -> None:
        data = self.raw_source.read_bytes(count // 8 + 1)
        for b in data:
            for shift in range(7, -1, -1):
                self._raw_bit_queue.append((b >> shift) & 1)

    def _extract_conditioned_bit(self) -> int:
        while True:
            while len(self._raw_bit_queue) < 2:
                self._fill_raw_bits()
            b1 = self._raw_bit_queue.popleft()
            b2 = self._raw_bit_queue.popleft()
            if b1 == 0 and b2 == 1:
                return 0
            elif b1 == 1 and b2 == 0:
                return 1
            # 00 and 11 are discarded

    def _read_bytes(self, n: int) -> bytes:
        needed_bits = n * 8
        while len(self._bit_queue) < needed_bits:
            self._bit_queue.append(self._extract_conditioned_bit())

        out = bytearray(n)
        for i in range(n):
            byte_val = 0
            for _ in range(8):
                byte_val = (byte_val << 1) | self._bit_queue.popleft()
            out[i] = byte_val
        return bytes(out)


class HashConditioner(RandomSource):
    """Cryptographic hash-based entropy conditioning (SHA-256 block compression)."""

    def __init__(
        self,
        raw_source: RandomSource,
        input_block_bytes: int = 64,
        source_id: str = "conditioned-sha256",
    ):
        super().__init__(
            source_id=source_id,
            kind=SourceKind.CONDITIONED,
            is_physical=raw_source.metadata.is_physical,
            is_synthetic=raw_source.metadata.is_synthetic,
        )
        self.raw_source = raw_source
        self.input_block_bytes = input_block_bytes
        self._buffer = bytearray()
        self.metadata.extra["input_block_bytes"] = input_block_bytes
        self.metadata.extra["compression_ratio"] = 32 / input_block_bytes

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            raw_block = self.raw_source.read_bytes(self.input_block_bytes)
            hashed = hashlib.sha256(raw_block).digest()
            self._buffer.extend(hashed)

        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out


class XORCombiner(RandomSource):
    """Combines two or more independent entropy sources via bitwise XOR.

    Under the Vazirani-Vazirani parity lemma, XORing independent biased sources
    exponentially reduces the bias: epsilon_combined = 2^(k-1) * prod(epsilon_i).
    """

    def __init__(self, sources: Sequence[RandomSource], source_id: str = "mixed-xor-combiner"):
        if len(sources) < 2:
            raise ValueError("XORCombiner requires at least two independent sources")
        is_physical = any(s.metadata.is_physical for s in sources)
        is_synthetic = any(s.metadata.is_synthetic for s in sources)
        super().__init__(
            source_id=source_id,
            kind=SourceKind.MIXED,
            is_physical=is_physical,
            is_synthetic=is_synthetic,
        )
        self.sources = list(sources)
        self.metadata.extra["combined_source_ids"] = [s.source_id for s in sources]

    def _read_bytes(self, n: int) -> bytes:
        first = bytearray(self.sources[0].read_bytes(n))
        for s in self.sources[1:]:
            other = s.read_bytes(n)
            for i in range(n):
                first[i] ^= other[i]
        return bytes(first)
