"""Hardware-independent acquisition interfaces and adapters for future instrumentation.

Supports:
- Serial framing, packet sequence counters, and dropped-packet detection.
- Pulse-event interval acquisition (e.g. Geiger-Müller detector, photodiode pulses).
- Checksum validation and acquisition metadata.
"""
from __future__ import annotations

import collections
import io
import struct
import time
from typing import Any, BinaryIO

from .base import RandomSource, SourceExhausted, SourceHealth, SourceKind


SERIAL_START_BYTE = 0xAA
SERIAL_END_BYTE = 0x55


class SerialTRNGAdapter(RandomSource):
    """Adapter for framed serial/USB hardware entropy instruments.

    Packet frame format:
    [START (1B)] [SEQ_NUM (2B)] [LEN (1B)] [PAYLOAD (N B)] [CRC16 (2B)] [END (1B)]
    """

    def __init__(self, stream: BinaryIO, source_id: str = "hardware-serial-trng"):
        super().__init__(
            source_id=source_id,
            kind=SourceKind.HARDWARE_SERIAL,
            is_physical=True,
            is_synthetic=False,
        )
        self.stream = stream
        self.last_seq: int | None = None
        self.dropped_packets: int = 0
        self.crc_errors: int = 0
        self._buffer = bytearray()
        self.metadata.extra.update({
            "dropped_packets": 0,
            "crc_errors": 0,
            "framing": "START(0xAA)-SEQ16-LEN8-PAYLOAD-CRC16-END(0x55)",
        })

    @staticmethod
    def _crc16(data: bytes) -> int:
        """Standard CCITT CRC-16."""
        crc = 0xFFFF
        for b in data:
            crc ^= (b << 8)
            for _ in range(8):
                if crc & 0x8000:
                    crc = ((crc << 1) ^ 0x1021) & 0xFFFF
                else:
                    crc = (crc << 1) & 0xFFFF
        return crc

    @classmethod
    def create_frame(cls, seq: int, payload: bytes) -> bytes:
        """Utility to construct a valid framed serial packet."""
        if len(payload) > 255:
            raise ValueError("Payload length cannot exceed 255 bytes per frame")
        header = struct.pack(">BBH", SERIAL_START_BYTE, len(payload), seq & 0xFFFF)
        crc = cls._crc16(header[1:] + payload)
        footer = struct.pack(">HB", crc, SERIAL_END_BYTE)
        return header + payload + footer

    def _read_next_frame(self) -> bytes:
        """Parse one valid frame from the serial stream."""
        while True:
            # Sync to start byte
            start = self.stream.read(1)
            if not start:
                raise SourceExhausted("Serial stream reached EOF")
            if start[0] != SERIAL_START_BYTE:
                continue

            hdr_bytes = self.stream.read(3)
            if len(hdr_bytes) < 3:
                raise SourceExhausted("Premature EOF during serial header")
            length, seq = struct.unpack(">BH", hdr_bytes)

            payload = self.stream.read(length)
            if len(payload) < length:
                raise SourceExhausted("Premature EOF during serial payload")

            tail = self.stream.read(3)
            if len(tail) < 3:
                raise SourceExhausted("Premature EOF during serial footer")
            crc_received, end_byte = struct.unpack(">HB", tail)

            if end_byte != SERIAL_END_BYTE:
                continue

            expected_crc = self._crc16(hdr_bytes + payload)
            if crc_received != expected_crc:
                self.crc_errors += 1
                self.metadata.extra["crc_errors"] = self.crc_errors
                continue

            # Check sequence numbering
            if self.last_seq is not None:
                expected_seq = (self.last_seq + 1) & 0xFFFF
                if seq != expected_seq:
                    dropped = (seq - expected_seq) & 0xFFFF
                    self.dropped_packets += dropped
                    self.metadata.extra["dropped_packets"] = self.dropped_packets
            self.last_seq = seq
            return payload

    def _read_bytes(self, n: int) -> bytes:
        while len(self._buffer) < n:
            frame_data = self._read_next_frame()
            self._buffer.extend(frame_data)
        out = bytes(self._buffer[:n])
        del self._buffer[:n]
        return out


class PulseEventAdapter(RandomSource):
    """Adapter converting physical pulse arrival intervals to random bits.

    Takes a stream of arrival timestamps (e.g. from Geiger counter or photodiode),
    and applies Von Neumann interval comparison:
    t2 - t1 > t3 - t2  => bit 1
    t2 - t1 < t3 - t2  => bit 0
    t2 - t1 == t3 - t2 => discarded
    """

    def __init__(self, timestamps_iterable, source_id: str = "hardware-pulse-events"):
        super().__init__(
            source_id=source_id,
            kind=SourceKind.HARDWARE_PULSE,
            is_physical=True,
            is_synthetic=False,
        )
        self._timestamps = iter(timestamps_iterable)
        self._bit_buffer = collections.deque()
        self.total_pulses = 0
        self.metadata.extra["total_pulses"] = 0

    def _extract_next_bit(self) -> int:
        while True:
            try:
                t1 = next(self._timestamps)
                t2 = next(self._timestamps)
                t3 = next(self._timestamps)
                self.total_pulses += 3
                self.metadata.extra["total_pulses"] = self.total_pulses
            except StopIteration:
                raise SourceExhausted("Pulse event stream exhausted")

            dt1 = t2 - t1
            dt2 = t3 - t2
            if dt1 > dt2:
                return 1
            elif dt1 < dt2:
                return 0
            # Tie: discard and continue

    def _read_bytes(self, n: int) -> bytes:
        needed_bits = n * 8
        while len(self._bit_buffer) < needed_bits:
            bit = self._extract_next_bit()
            self._bit_buffer.append(bit)

        out = bytearray(n)
        for i in range(n):
            byte_val = 0
            for b in range(8):
                byte_val = (byte_val << 1) | self._bit_buffer.popleft()
            out[i] = byte_val
        return bytes(out)
