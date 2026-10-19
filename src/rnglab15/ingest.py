"""Explicit schemas only. Acquisition order is never sorted or repaired."""
from __future__ import annotations
import csv
from dataclasses import dataclass
from decimal import Decimal
import re
from .integrity import lines

INTEGER = re.compile(r"[+-]?[0-9]+\Z")


def integer(text):
    if not isinstance(text, str) or not INTEGER.fullmatch(text):
        raise ValueError("expected decimal integer")
    return int(text)


@dataclass(frozen=True)
class Waveform:
    row: int
    channel: int
    timestamp: int
    samples: tuple[int, ...]


def waveforms(path, *, samples, metadata_fields=2, min_code=0, max_code=65535):
    if metadata_fields != 2 or samples < 1:
        raise ValueError("adapter requires explicit two-field schema and sample count")
    for row, text in enumerate(lines(path), 1):
        if not text.strip():
            continue
        fields = text.split()
        if len(fields) != metadata_fields + samples:
            raise ValueError(f"row {row}: expected {metadata_fields + samples} fields, got {len(fields)}")
        try:
            values = tuple(integer(x) for x in fields)
        except ValueError as e:
            raise ValueError(f"row {row}: {e}") from e
        if any(x < min_code or x > max_code for x in values[2:]):
            raise ValueError(f"row {row}: waveform outside declared code range")
        yield Waveform(row, values[0], values[1], values[2:])


def csv_rows(path):
    reader = csv.DictReader(lines(path), strict=True)
    if not reader.fieldnames or len(reader.fieldnames) != len(set(reader.fieldnames)):
        raise ValueError("empty or duplicate CSV columns")
    for row, record in enumerate(reader, 2):
        if None in record or None in record.values():
            raise ValueError(f"row {row}: malformed CSV shape")
        yield row, record


def sensors(path, columns):
    required = {"time_us", *columns}
    for row, record in csv_rows(path):
        if not required <= record.keys():
            raise ValueError(f"row {row}: missing declared sensor columns")
        values = {}
        try:
            timestamp = integer(record["time_us"])
            for key in columns:
                values[key] = None if record[key] == "NA" else integer(record[key])
        except ValueError as e:
            raise ValueError(f"row {row}: {e}") from e
        yield row, timestamp, values, record


def events(path, timestamp_column="timestamp", channel_column="channel"):
    """Integer device ticks only. Count-window data must use a different adapter."""
    for row, record in csv_rows(path):
        if timestamp_column not in record or channel_column not in record:
            raise ValueError("explicit timestamp and channel columns required")
        yield row, record[channel_column], integer(record[timestamp_column])


def seconds_per_tick(metadata):
    value = metadata.get("seconds_per_tick")
    if value is None:
        return None
    if not metadata.get("clock_evidence"):
        raise ValueError("clock scale requires a documentation/calibration reference")
    value = Decimal(str(value))
    if not value.is_finite() or value <= 0:
        raise ValueError("tick scale must be positive and finite")
    return value
