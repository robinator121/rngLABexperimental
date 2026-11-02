"""Read-only physical export audit. Does not turn measurements into real keys."""
from __future__ import annotations

import hashlib
import math
from pathlib import Path

import numpy as np

from rnglab15.ingest import waveforms, seconds_per_tick
from rnglab15.assessment import Timing


def audit_waveforms(path, metadata):
    if not isinstance(metadata.get("samples_per_record"), int) or metadata["samples_per_record"] < 2:
        raise ValueError("explicit samples_per_record >=2 required")
    if not metadata.get("provenance"):
        raise ValueError("source provenance required, even if incomplete")
    scale = seconds_per_tick(metadata)
    timing, seen_codes = Timing(), set()
    total = ceiling = floor = clipped_records = count = pairs = 0
    sx = sy = xx = yy = xy = 0.
    first_tick = last_tick = None
    for record in waveforms(path, samples=metadata["samples_per_record"]):
        count += 1
        timing.add(record.timestamp)
        if first_tick is None:
            first_tick = record.timestamp
        last_tick = record.timestamp
        data = np.asarray(record.samples, dtype=float)
        total += len(data)
        ceiling += int(np.sum(data == 65535))
        floor += int(np.sum(data == 0))
        clipped_records += bool(np.any(data == 65535))
        seen_codes.update(record.samples)
        x, y = data[:-1], data[1:]
        pairs += len(x)
        sx += float(x.sum()); sy += float(y.sum())
        xx += float(x@x); yy += float(y@y); xy += float(x@y)
    if not count:
        raise ValueError("empty physical export")
    denominator = (pairs*xx-sx*sx)*(pairs*yy-sy*sy)
    source_hash = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            source_hash.update(chunk)
    return {"source_sha256": source_hash.hexdigest(), "source_bytes": Path(path).stat().st_size,
            "metadata": metadata, "records": count, "samples": total,
            "ceiling_fraction": ceiling/total, "floor_fraction": floor/total,
            "records_with_ceiling_fraction": clipped_records/count,
            "distinct_codes": len(seen_codes), "native_adc_resolution_bits": metadata.get("native_adc_bits"),
            "effective_resolution_bits": None, "effective_resolution_reason": "requires calibrated ADC/noise characterization",
            "within_record_raw_lag1_correlation": (pairs*xy-sx*sy)/math.sqrt(denominator) if denominator > 0 else None,
            "dependence_qualification": "raw waveform shape contributes; boundaries are never joined",
            "timing": timing.report(scale), "first_tick": first_tick, "last_tick": last_tick,
            "acquisition_start": metadata.get("acquisition_start"), "acquisition_end": metadata.get("acquisition_end"),
            "conditioning": "none in this audit; source processing requires provenance",
            "validated_entropy_credit_bits": 0,
            "crypto_use": "public replay research only; not fresh entropy or real key material"}
