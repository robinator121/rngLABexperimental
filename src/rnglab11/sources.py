"""Bounded byte streams and explicit transformations for research inputs.

Accounting is for bytes actually returned/consumed, including partial reads on
exhaustion. No class infers physical entropy from a file name or from its output.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import random
from typing import Any, Mapping

import numpy as np


class SourceExhausted(EOFError):
    """The partial bytes have already been consumed and must not be retried."""

    def __init__(self, requested: int, partial: bytes = b""):
        self.requested, self.partial = requested, partial
        super().__init__(f"source exhausted: requested {requested}, obtained {len(partial)} bytes")


class SourceBudgetExceeded(ValueError):
    """A request would exceed its hard budget; nothing was consumed."""


class ByteSource:
    """Exact-read interface; instances are single-consumer, not thread safe."""

    def __init__(self, metadata: Mapping[str, Any], *, max_bytes: int | None = None):
        if max_bytes is not None and (not isinstance(max_bytes, int) or max_bytes < 0):
            raise ValueError("max_bytes must be a nonnegative integer")
        # Freeze caller-owned nested metadata by storing its canonical encoding.
        self._metadata_json = json.dumps(dict(metadata), sort_keys=True, allow_nan=False)
        self.source_id = hashlib.sha256(self._metadata_json.encode()).hexdigest()
        self.max_bytes = max_bytes
        self.bytes_consumed = 0
        self.offset = 0
        self._consumed_digest = hashlib.sha256()

    def _read(self, n: int) -> bytes:
        raise NotImplementedError

    def read(self, n: int) -> bytes:
        if not isinstance(n, int) or isinstance(n, bool) or n < 0:
            raise ValueError("read size must be a nonnegative integer")
        if self.max_bytes is not None and self.bytes_consumed + n > self.max_bytes:
            raise SourceBudgetExceeded("request exceeds source byte budget")
        if n == 0:
            return b""
        data = self._read(n)
        if len(data) > n:
            raise RuntimeError("source returned more bytes than requested")
        self.bytes_consumed += len(data)
        self.offset += len(data)
        self._consumed_digest.update(data)
        if len(data) != n:
            raise SourceExhausted(n, data)
        return data

    def provenance(self) -> dict[str, Any]:
        return {"source_id": self.source_id, "metadata": json.loads(self._metadata_json),
                "bytes_consumed": self.bytes_consumed, "offset": self.offset,
                "consumed_bytes_sha256": self._consumed_digest.hexdigest(),
                "max_bytes": self.max_bytes}

    def clone(self) -> "ByteSource":
        """Return an independent stream positioned at its beginning."""
        raise NotImplementedError("this source cannot be replayed without a capture")

    def close(self) -> None:
        pass

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


class BytesSource(ByteSource):
    def __init__(self, data: bytes, *, metadata: Mapping[str, Any] | None = None,
                 max_bytes: int | None = None):
        self._data = bytes(data)
        self._origin_metadata = json.loads(json.dumps(dict(metadata or {}), allow_nan=False))
        super().__init__({**self._origin_metadata, "adapter": "bytes", "size_bytes": len(data),
                          "content_sha256": hashlib.sha256(data).hexdigest()}, max_bytes=max_bytes)

    def _read(self, n):
        return self._data[self.offset:self.offset + n]

    def seek(self, offset: int):
        if not isinstance(offset, int) or not 0 <= offset <= len(self._data):
            raise ValueError("invalid offset")
        self.offset = offset

    def clone(self):
        return BytesSource(self._data, metadata=self._origin_metadata, max_bytes=self.max_bytes)


class FileSource(ByteSource):
    """Binary input with a stable file-identity manifest and consumed-prefix hash.

    The identity digest hashes metadata, not unseen file contents. A declared
    content SHA256 is a caller assertion unless independently verified. Replays
    reject changed size/mtime; callers should use immutable datasets.
    """

    def __init__(self, path: str | Path, *, kind: str = "binary",
                 metadata: Mapping[str, Any] | None = None, max_bytes: int | None = None):
        self.path = Path(path).resolve()
        self.kind = kind
        self._origin_metadata = json.loads(json.dumps(dict(metadata or {}), allow_nan=False))
        self._file = self.path.open("rb")
        stat = os.fstat(self._file.fileno())
        self._identity = (stat.st_size, stat.st_mtime_ns, stat.st_dev, stat.st_ino)
        super().__init__({**self._origin_metadata, "adapter": "binary_file", "kind": kind,
                          "path": str(self.path), "size_bytes": stat.st_size,
                          "mtime_ns": stat.st_mtime_ns, "file_identity": list(self._identity),
                          "entropy_claim": "none"}, max_bytes=max_bytes)

    def _check_unchanged(self):
        stat = os.fstat(self._file.fileno())
        if (stat.st_size, stat.st_mtime_ns, stat.st_dev, stat.st_ino) != self._identity:
            raise RuntimeError("input file changed during use")

    def _read(self, n):
        self._check_unchanged()
        return self._file.read(n)

    def seek(self, offset: int):
        self._check_unchanged()
        if not isinstance(offset, int) or not 0 <= offset <= self._identity[0]:
            raise ValueError("invalid offset")
        self._file.seek(offset)
        self.offset = offset

    def clone(self):
        self._check_unchanged()
        other = FileSource(self.path, kind=self.kind, metadata=self._origin_metadata,
                           max_bytes=self.max_bytes)
        if other._identity != self._identity:
            other.close()
            raise RuntimeError("input file changed before replay")
        return other

    def close(self):
        self._file.close()


class SyntheticSource(ByteSource):
    """Chunk-independent seeded Python MT19937 stream, expressly NOT a CSPRNG."""

    def __init__(self, seed: int | str, *, max_bytes: int = 1_048_576):
        if not isinstance(seed, (int, str)):
            raise ValueError("seed must be an integer or string")
        self.seed = seed
        self._rng = random.Random(seed)
        self._buffer = bytearray()
        super().__init__({"adapter": "synthetic_mt19937", "seed": seed,
                          "cryptographic_security": False, "entropy_claim": "none"},
                         max_bytes=max_bytes)

    def _read(self, n):
        while len(self._buffer) < n:
            self._buffer.extend(self._rng.getrandbits(512).to_bytes(64, "little"))
        data = bytes(self._buffer[:n])
        del self._buffer[:n]
        return data

    def clone(self):
        return SyntheticSource(self.seed, max_bytes=self.max_bytes)


class SystemCSPRNGSource(ByteSource):
    """OS randomness, optionally captured to a NEW private file for replay."""

    def __init__(self, *, max_bytes: int = 1_048_576, capture_path: str | Path | None = None):
        self.capture_path = Path(capture_path).resolve() if capture_path else None
        self._capture = None
        if self.capture_path:
            fd = os.open(self.capture_path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
            self._capture = os.fdopen(fd, "wb")
        super().__init__({"adapter": "os_urandom", "capture_path": str(self.capture_path)
                          if self.capture_path else None, "physical_entropy_assessed": False},
                         max_bytes=max_bytes)

    def _read(self, n):
        data = os.urandom(n)
        if self._capture:
            self._capture.write(data)
            self._capture.flush()
        return data

    def close(self):
        if self._capture:
            self._capture.close()
            self._capture = None


class TimestampCSVSource(ByteSource):
    """CSV integer event ticks -> nonnegative uint64 interarrival ticks (LE).

    This is a lossless interval encoding relative to the first timestamp, NOT
    an entropy extractor. The first timestamp is recorded, not emitted. Input
    bytes/records and adapter output bytes have separate counters.
    """

    def __init__(self, path: str | Path, *, timestamp_column: str = "timestamp",
                 delimiter: str = ",", tick_unit: str = "unspecified",
                 max_bytes: int = 1_048_576, max_input_bytes: int = 16_777_216,
                 max_line_bytes: int = 4096):
        if len(delimiter) != 1 or max_line_bytes < 2 or max_input_bytes < 1:
            raise ValueError("invalid CSV adapter limits")
        self.path = Path(path).resolve()
        self._file = self.path.open("rb")
        self.timestamp_column, self.delimiter = timestamp_column, delimiter
        self.tick_unit = tick_unit
        self.max_input_bytes, self.max_line_bytes = max_input_bytes, max_line_bytes
        self.input_bytes_consumed = self.records_consumed = 0
        self._input_digest = hashlib.sha256()
        self.first_timestamp = self._previous = None
        self._buffer = bytearray()
        self._columns = None
        stat = os.fstat(self._file.fileno())
        self._identity = (stat.st_size, stat.st_mtime_ns, stat.st_ino)
        super().__init__({"adapter": "csv_integer_tick_intervals_u64le", "path": str(self.path),
                          "file_identity": list(self._identity), "timestamp_column": timestamp_column,
                          "delimiter": delimiter, "tick_unit": tick_unit,
                          "entropy_claim": "none"}, max_bytes=max_bytes)

    def _line(self):
        remaining = self.max_input_bytes - self.input_bytes_consumed
        if remaining <= 0:
            raise SourceBudgetExceeded("timestamp raw input budget exhausted")
        line = self._file.readline(min(self.max_line_bytes + 1, remaining))
        self.input_bytes_consumed += len(line)
        self._input_digest.update(line)
        if len(line) > self.max_line_bytes:
            raise ValueError("timestamp CSV line exceeds limit")
        if line and not line.endswith(b"\n") and self.input_bytes_consumed >= self.max_input_bytes:
            raise SourceBudgetExceeded("raw input budget ended inside a CSV line")
        return line

    def _read(self, n):
        stat = os.fstat(self._file.fileno())
        if (stat.st_size, stat.st_mtime_ns, stat.st_ino) != self._identity:
            raise RuntimeError("timestamp input changed during use")
        while len(self._buffer) < n:
            line = self._line()
            if not line:
                break
            row = next(csv.reader([line.decode("utf-8-sig")], delimiter=self.delimiter))
            if self._columns is None:
                self._columns = row
                if row.count(self.timestamp_column) != 1:
                    raise ValueError("timestamp column must occur exactly once")
                continue
            if len(row) != len(self._columns):
                raise ValueError("timestamp CSV column count changed")
            try:
                tick = int(row[self._columns.index(self.timestamp_column)])
            except ValueError as exc:
                raise ValueError("timestamps must be integer ticks; convert fractional units explicitly") from exc
            self.records_consumed += 1
            if self._previous is None:
                self.first_timestamp = self._previous = tick
                continue
            delta = tick - self._previous
            if not 0 <= delta < 2**64:
                raise ValueError("timestamps must be nondecreasing with uint64 intervals")
            self._previous = tick
            self._buffer.extend(delta.to_bytes(8, "little"))
        data = bytes(self._buffer[:n])
        del self._buffer[:n]
        return data

    def provenance(self):
        return {**super().provenance(), "input_bytes_consumed": self.input_bytes_consumed,
                "input_prefix_sha256": self._input_digest.hexdigest(),
                "records_consumed": self.records_consumed, "first_timestamp": self.first_timestamp}

    def clone(self):
        other = TimestampCSVSource(self.path, timestamp_column=self.timestamp_column,
                                   delimiter=self.delimiter, tick_unit=self.tick_unit,
                                   max_bytes=self.max_bytes, max_input_bytes=self.max_input_bytes,
                                   max_line_bytes=self.max_line_bytes)
        if other._identity != self._identity:
            other.close()
            raise RuntimeError("timestamp input changed before replay")
        return other

    def close(self):
        self._file.close()


class UniformStream:
    """One uniform consumes exactly eight bytes; top 53 bits map to [0,1).

    Multiple channels are consecutive draws reshaped in row-major order. They
    are not asserted to be independent. Integer conversion is little-endian.
    """

    def __init__(self, source: ByteSource, *, max_draws: int = 1_000_000):
        self.source, self.max_draws = source, max_draws
        self.uniforms_generated = 0

    def random(self, count: int, channels: int = 1) -> np.ndarray:
        if not isinstance(count, int) or count < 0 or not isinstance(channels, int) or channels < 1:
            raise ValueError("count must be nonnegative and channels positive integers")
        draws = count * channels
        if self.uniforms_generated + draws > self.max_draws:
            raise SourceBudgetExceeded("uniform draw budget exceeded")
        data = self.source.read(8 * draws)
        values = (np.frombuffer(data, dtype="<u8") >> np.uint64(11)).astype(np.float64)
        self.uniforms_generated += draws
        return (values * (1.0 / 2**53)).reshape(count, channels)


class SourceFaultSource(ByteSource):
    """Deterministic byte-position fault trajectory BEFORE any conditioning.

    These are synthetic interventions with known control seeds, not fitted
    physical degradation models. No physical entropy is claimed or created.
    Finite start/end positions represent bursts; linearly varying severity
    represents a bounded drift. Repetition uses the preceding emitted byte.
    """

    def __init__(self, source: ByteSource, *, mode: str = "stuck_zero",
                 start_offset: int = 0, end_offset: int | None = None,
                 severity: float = 1.0, final_severity: float | None = None,
                 fault_seed: int = 0, xor_mask: int = 1, max_bytes: int = 1_048_576):
        if mode not in ("stuck_zero", "bit_bias_zero", "repeat_previous", "bit_flip"):
            raise ValueError("unsupported source fault")
        if (not isinstance(start_offset, int) or start_offset < 0 or
                end_offset is not None and (not isinstance(end_offset, int) or end_offset <= start_offset)):
            raise ValueError("invalid fault interval")
        final_severity = severity if final_severity is None else final_severity
        if not 0 <= severity <= 1 or not 0 <= final_severity <= 1:
            raise ValueError("fault severity must be in [0,1]")
        if end_offset is None and final_severity != severity:
            raise ValueError("drift requires a finite end_offset")
        if not isinstance(fault_seed, int) or not isinstance(xor_mask, int) or not 0 <= xor_mask <= 255:
            raise ValueError("invalid deterministic fault controls")
        self.source, self.mode = source, mode
        self.start_offset, self.end_offset = start_offset, end_offset
        self.severity, self.final_severity = float(severity), float(final_severity)
        self.xor_mask = xor_mask
        self._rng = random.Random(fault_seed)
        self._previous = None
        self.exposed_bytes = self.modified_bytes = self.modified_bits = 0
        super().__init__({"adapter": "synthetic_source_fault_v1", "mode": mode,
                          "start_offset": start_offset, "end_offset_exclusive": end_offset,
                          "severity": self.severity, "final_severity": self.final_severity,
                          "fault_seed": fault_seed, "xor_mask": xor_mask,
                          "origin_source_id": source.source_id,
                          "entropy_claim": "none; deterministic synthetic fault intervention"}, max_bytes=max_bytes)

    def severity_at(self, offset: int) -> float:
        if offset < self.start_offset or self.end_offset is not None and offset >= self.end_offset:
            return 0.0
        if self.end_offset is None:
            return self.severity
        fraction = (offset - self.start_offset) / max(1, self.end_offset - self.start_offset - 1)
        return self.severity + fraction * (self.final_severity - self.severity)

    def _read(self, n):
        try:
            raw = self.source.read(n)
        except SourceExhausted as exc:
            raw = exc.partial
        output = bytearray(raw)
        for index, original in enumerate(raw):
            offset = self.offset + index
            active = offset >= self.start_offset and (self.end_offset is None or offset < self.end_offset)
            changed = original
            if active:
                self.exposed_bytes += 1
                probability = self.severity_at(offset)
                if self.mode == "bit_bias_zero":
                    for bit in range(8):
                        if self._rng.random() < probability:
                            changed &= ~(1 << bit)
                elif self._rng.random() < probability:
                    if self.mode == "stuck_zero":
                        changed = 0
                    elif self.mode == "bit_flip":
                        changed ^= self.xor_mask
                    elif self._previous is not None:
                        changed = self._previous
            self.modified_bytes += int(changed != original)
            self.modified_bits += (changed ^ original).bit_count()
            output[index] = changed
            self._previous = changed
        return bytes(output)

    def provenance(self):
        return {**super().provenance(), "exposed_bytes": self.exposed_bytes,
                "modified_bytes": self.modified_bytes, "modified_bits": self.modified_bits,
                "origin": self.source.provenance()}

    def close(self):
        self.source.close()


def iter_chunks(source: ByteSource, *, chunk_bytes: int = 65_536,
                total_bytes: int = 1_048_576):
    """Finite chunk iterator; yields a final partial chunk on physical EOF."""
    if chunk_bytes < 1 or total_bytes < 0:
        raise ValueError("invalid chunk limits")
    remaining = total_bytes
    while remaining:
        try:
            data = source.read(min(chunk_bytes, remaining))
        except SourceExhausted as exc:
            if exc.partial:
                yield exc.partial
            break
        yield data
        remaining -= len(data)
