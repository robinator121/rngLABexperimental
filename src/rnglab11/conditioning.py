"""Auditable streaming treatments; none creates physical entropy."""
from __future__ import annotations

import hashlib
import time
from typing import Iterable

from .sources import ByteSource, BytesSource, SourceExhausted, SourceBudgetExceeded


class VonNeumannConditioner:
    """01 -> 0, 10 -> 1; discard 00/11 (bits are MSB first).

    Unbiased output requires suitable independence/stationarity assumptions;
    arbitrary correlated physical inputs need not satisfy them. Pending pairs
    and output bits survive feed boundaries. finish drops incomplete bytes.
    """

    name = "von_neumann_pairs_msb"

    def __init__(self):
        self.input_bits = self.output_bits = self.output_bytes = 0
        self.discarded_pairs = 0
        self._pending_bit = None
        self._output_value = self._output_count = 0
        self.finished = False
        self.discarded_final_bits = 0

    def feed_bits(self, bits: Iterable[int]) -> bytes:
        if self.finished:
            raise ValueError("conditioner already finished")
        output = bytearray()
        for bit in bits:
            if bit not in (0, 1):
                raise ValueError("input bits must be 0 or 1")
            self.input_bits += 1
            if self._pending_bit is None:
                self._pending_bit = bit
                continue
            first, self._pending_bit = self._pending_bit, None
            if first == bit:
                self.discarded_pairs += 1
                continue
            self.output_bits += 1
            self._output_value = (self._output_value << 1) | first
            self._output_count += 1
            if self._output_count == 8:
                output.append(self._output_value)
                self._output_value = self._output_count = 0
                self.output_bytes += 1
        return bytes(output)

    def feed(self, data: bytes) -> bytes:
        return self.feed_bits((byte >> shift) & 1 for byte in data for shift in range(7, -1, -1))

    def finish(self) -> bytes:
        if not self.finished:
            self.discarded_final_bits = self._output_count
            self.finished = True
        return b""

    def report(self):
        return {"treatment": self.name, "input_bits": self.input_bits,
                "input_bytes": self.input_bits // 8, "output_bits_before_byte_truncation": self.output_bits,
                "output_bytes": self.output_bytes, "discarded_equal_pairs": self.discarded_pairs,
                "pending_input_bit": int(self._pending_bit is not None),
                "discarded_final_output_bits": self.discarded_final_bits,
                "entropy_claim": "none; VN assumptions require separate assessment"}


class HashConditioner:
    """SHA256 of each fixed-size input block; incomplete last block is dropped.

    This treatment is not an SP 800-90B vetted conditioning component merely
    because SHA256 is used. Entropy accounting and an assessed source are absent.
    """

    name = "sha256_fixed_blocks"

    def __init__(self, block_bytes: int = 64):
        if not isinstance(block_bytes, int) or not 32 <= block_bytes <= 1_048_576:
            raise ValueError("block_bytes must be in [32,1048576]")
        self.block_bytes = block_bytes
        self.input_bytes = self.output_bytes = 0
        self.discarded_tail_bytes = 0
        self._pending = bytearray()
        self.finished = False

    def feed(self, data: bytes) -> bytes:
        if self.finished:
            raise ValueError("conditioner already finished")
        self.input_bytes += len(data)
        output = bytearray()
        # Only one partial input block is retained between calls.
        cursor = 0
        if self._pending:
            take = min(self.block_bytes - len(self._pending), len(data))
            self._pending.extend(data[:take])
            cursor = take
            if len(self._pending) == self.block_bytes:
                output.extend(hashlib.sha256(self._pending).digest())
                self._pending.clear()
        while cursor + self.block_bytes <= len(data):
            output.extend(hashlib.sha256(data[cursor:cursor + self.block_bytes]).digest())
            cursor += self.block_bytes
        self._pending.extend(data[cursor:])
        self.output_bytes += len(output)
        return bytes(output)

    def finish(self) -> bytes:
        if not self.finished:
            self.discarded_tail_bytes = len(self._pending)
            self._pending.clear()
            self.finished = True
        return b""

    def report(self):
        return {"treatment": self.name, "block_bytes": self.block_bytes,
                "input_bytes": self.input_bytes, "output_bytes": self.output_bytes,
                "discarded_tail_bytes": self.discarded_tail_bytes,
                "pending_input_bytes": len(self._pending), "entropy_claim": "none"}


class ConditionedSource(ByteSource):
    """Bounded producer/consumer adapter; raw and emitted-byte counts differ."""

    def __init__(self, source: ByteSource, conditioner, *, max_input_bytes: int = 1_048_576,
                 max_output_bytes: int = 1_048_576, chunk_bytes: int = 4096):
        if max_input_bytes < 0 or chunk_bytes < 1 or chunk_bytes > 1_048_576:
            raise ValueError("invalid conditioning limits")
        self.source, self.conditioner = source, conditioner
        self.max_input_bytes, self.chunk_bytes = max_input_bytes, chunk_bytes
        self.input_bytes_consumed = 0
        self._buffer = bytearray()
        self._done = False
        self.processing_seconds = 0.0
        super().__init__({"adapter": "conditioned", "origin_source_id": source.source_id,
                          "treatment": conditioner.name, "entropy_claim": "none"},
                         max_bytes=max_output_bytes)

    def _read(self, n):
        while len(self._buffer) < n and not self._done:
            remaining = self.max_input_bytes - self.input_bytes_consumed
            if self.source.max_bytes is not None:
                remaining = min(remaining, self.source.max_bytes - self.source.bytes_consumed)
            if remaining <= 0:
                self._done = True
                self._buffer.extend(self.conditioner.finish())
                break
            try:
                data = self.source.read(min(remaining, self.chunk_bytes))
            except SourceExhausted as exc:
                data = exc.partial
                self._done = True
            self.input_bytes_consumed += len(data)
            start = time.perf_counter()
            self._buffer.extend(self.conditioner.feed(data))
            if self._done:
                self._buffer.extend(self.conditioner.finish())
            self.processing_seconds += time.perf_counter() - start
        data = bytes(self._buffer[:n])
        del self._buffer[:n]
        return data

    def provenance(self):
        return {**super().provenance(), "input_bytes_consumed": self.input_bytes_consumed,
                "buffered_output_bytes": len(self._buffer), "processing_seconds": self.processing_seconds,
                "transformation": self.conditioner.report(), "origin": self.source.provenance()}

    def close(self):
        self.source.close()


def compare_treatments(source: ByteSource, *, input_bytes: int = 4096,
                       hash_block_bytes: int = 64, max_input_bytes: int = 1_048_576):
    """Construct matched finite treatments from the SAME exact source segment.

    The bounded development helper returns replayable output sources and
    transformation accounting. It does not compare security or estimate entropy.
    For large inputs use independent clones with ConditionedSource instead.
    """
    if not 0 <= input_bytes <= max_input_bytes <= 16_777_216:
        raise SourceBudgetExceeded("matched comparison exceeds explicit input budget")
    origin_offset = source.offset
    data = source.read(input_bytes)
    origin = {"origin_source_id": source.source_id, "origin_offset": origin_offset,
              "origin_length": len(data), "origin_segment_sha256": hashlib.sha256(data).hexdigest()}
    treatments = {"raw": (None, data)}
    for name, conditioner in (("debiased", VonNeumannConditioner()),
                              ("conditioned", HashConditioner(hash_block_bytes))):
        start = time.perf_counter()
        output = conditioner.feed(data) + conditioner.finish()
        treatments[name] = ({**conditioner.report(), "processing_seconds": time.perf_counter() - start}, output)
    return {name: {"source": BytesSource(output, metadata={**origin, "treatment": name}),
                   "report": report or {"treatment": "raw", "input_bytes": len(data),
                                         "output_bytes": len(data), "entropy_claim": "none"},
                   "origin": dict(origin)}
            for name, (report, output) in treatments.items()}
