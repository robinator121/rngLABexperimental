"""Synthetic interventions with frozen position units and explicit placement.

Byte intervals are absolute, zero-based, half-open byte positions. Uniform
intervals are absolute, zero-based, half-open application rows. Interventions
never imply that physical entropy has been measured or a fault is detectable.
"""
from __future__ import annotations

import hashlib
import math
from typing import Mapping

import numpy as np

from rnglab11.degradation import DegradationSource
from rnglab11.sources import ByteSource, SourceExhausted

APPLICATION_CHANNELS = {
    "integration": {"proposal": [0]},
    "pi": {"proposal": [0, 1]},
    "mcmc": {"proposal": [0], "acceptance": [1]},
    "vmc": {"proposal": [0], "acceptance": [1]},
    "ising": {"selection": [0], "acceptance": [1]},
    "harmonic_pimc": {"selection": [0], "proposal": [1], "acceptance": [2]},
    "quartic_pimc": {"selection": [0], "proposal": [1], "acceptance": [2]},
    "double_well_pimc": {"selection": [0], "proposal": [1], "acceptance": [2]},
    "tilted_double_well_pimc": {"selection": [0], "proposal": [1], "acceptance": [2]},
}
_ALIASES = {"harmonic": "harmonic_pimc", "gaussian_mcmc": "mcmc", "pimc": "harmonic_pimc",
            "polynomial_pimc": "quartic_pimc", "quartic": "quartic_pimc", "double_well": "double_well_pimc"}
FAULT_MODES = {"bias", "autocorrelation", "copy_bit", "stuck_bit", "periodicity", "zero_burst",
               "drift", "entropy_reduction", "stuck_zero", "bit_flip"}


def _integer(value, name, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or maximum is not None and value > maximum:
        raise ValueError(f"invalid {name}")
    return value


def _probability(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f"{name} must be finite in [0,1]")
    return float(value)


def _interval(spec, unit):
    singular = "row" if unit == "rows" else "byte"
    other = "byte" if unit == "rows" else "row"
    foreign_keys = {f"start_{other}", f"start_{other}s", f"end_{other}", f"end_{other}s", f"duration_{other}s"}
    if foreign_keys.intersection(spec):
        raise ValueError(f"fault position fields must use {unit}, not {other}s")
    for field in ("start", "end"):
        if f"{field}_{singular}" in spec and f"{field}_{unit}" in spec and spec[f"{field}_{singular}"] != spec[f"{field}_{unit}"]:
            raise ValueError("contradictory fault interval aliases")
    start = _integer(spec.get(f"start_{singular}", spec.get(f"start_{unit}", spec.get("start", 0))), f"start_{unit}")
    duration = spec.get(f"duration_{unit}", spec.get("duration"))
    end = spec.get(f"end_{singular}", spec.get(f"end_{unit}"))
    if end is not None and duration is not None:
        raise ValueError("provide duration or end, not both")
    if duration is not None:
        end = start + _integer(duration, f"duration_{unit}", 1)
    if end is not None:
        _integer(end, f"end_{unit}", start + 1)
    declared = spec.get("unit")
    if declared is not None and declared != unit:
        raise ValueError(f"fault position unit must be {unit}")
    return start, end


def _faults(spec):
    if not isinstance(spec, Mapping):
        raise ValueError("fault specification must be a mapping")
    if "faults" in spec:
        if "mode" in spec or not isinstance(spec["faults"], list) or not 1 <= len(spec["faults"]) <= 16:
            raise ValueError("fault combinations require a nonempty list of at most 16 interventions")
        for child in spec["faults"]:
            if not isinstance(child, Mapping) or "faults" in child:
                raise ValueError("nested fault combinations are not supported")
        return spec["faults"]
    return [spec]


class _BitFaultSource(ByteSource):
    def __init__(self, source, spec, *, max_bytes):
        self.source = source
        self.mode = spec["mode"]
        self.start, self.end = _interval(spec, "bytes")
        self.bit = _integer(spec.get("bit", 0), "bit", 0, 7)
        self.from_bit = _integer(spec.get("from_bit", 1), "from_bit", 0, 7)
        self.stuck_value = _integer(spec.get("value", 0), "value", 0, 1)
        self.keep_bits = _integer(spec.get("keep_bits", 4), "keep_bits", 0, 8)
        self.severity = _probability(spec.get("severity", 1.), "severity")
        self.final = _probability(spec.get("final_severity", self.severity), "final_severity")
        self.seed = _integer(spec.get("seed", 0), "seed")
        if self.final != self.severity and self.end is None:
            raise ValueError("drifting intervention probability needs a finite interval")
        self.modified_bytes = self.exposed_bytes = 0
        super().__init__({"adapter": "synthetic_bit_fault_v012", "mode": self.mode,
                          "start_byte": self.start, "end_byte_exclusive": self.end,
                          "bit_lsb_zero_based": self.bit, "from_bit": self.from_bit,
                          "stuck_value": self.stuck_value, "keep_high_bits": self.keep_bits,
                          "severity": self.severity, "final_severity": self.final,
                          "control_seed": self.seed,
                          "origin_source_id": source.source_id,
                          "physical_entropy_assessed": False, "entropy_claim": "none"}, max_bytes=max_bytes)

    def _read(self, n):
        try:
            raw = self.source.read(n)
        except SourceExhausted as exc:
            raw = exc.partial
        output = bytearray(raw)
        for j, value in enumerate(raw):
            absolute = self.offset + j
            if absolute < self.start or self.end is not None and absolute >= self.end:
                continue
            self.exposed_bytes += 1
            fraction = 0 if self.end is None else (absolute - self.start) / max(1, self.end - self.start - 1)
            probability = self.severity + fraction * (self.final - self.severity)
            control = int.from_bytes(hashlib.sha256(f"byte-gate:{self.seed}:{absolute}".encode()).digest()[:8], "little") / 2**64
            if control >= probability:
                continue
            if self.mode == "copy_bit":
                changed = (value & ~(1 << self.bit)) | (((value >> self.from_bit) & 1) << self.bit)
            elif self.mode == "stuck_bit":
                changed = (value & ~(1 << self.bit)) | (self.stuck_value << self.bit)
            else:
                changed = value & (255 << (8 - self.keep_bits) & 255)
            output[j] = changed
            self.modified_bytes += changed != value
        return bytes(output)

    def provenance(self):
        return {**super().provenance(), "exposed_bytes": self.exposed_bytes,
                "modified_bytes": self.modified_bytes, "origin": self.source.provenance()}

    def close(self):
        self.source.close()


class _PlacedFaultSource(ByteSource):
    """Record the caller's explicit wrapper position and ordered interventions."""

    def __init__(self, source, spec, max_bytes):
        self.source = source
        super().__init__({"adapter": "explicit_byte_fault_placement_v012",
                          "ordered_faults": [dict(item) for item in _faults(spec)],
                          "position_unit": "absolute_zero_based_bytes",
                          "physical_entropy_assessed": False,
                          "entropy_reduction_control_is_not_an_entropy_measurement": True}, max_bytes=max_bytes)

    def _read(self, n):
        try:
            return self.source.read(n)
        except SourceExhausted as exc:
            return exc.partial

    def provenance(self):
        return {**super().provenance(), "origin": self.source.provenance()}

    def close(self):
        self.source.close()


def make_byte_fault(source: ByteSource, spec: Mapping, *, max_bytes=1_048_576):
    """Wrap byte faults in listed order, preserving the consumed native stream.

    DegradationSource retains its cycle/previous-byte/random-control state
    across reads. copy_bit copies a bit WITHIN the same delivered byte; stream
    autocorrelation uses previous-byte copying. Quantization is a synthetic
    predictability control, not an estimate of source min-entropy.
    """
    for fault in _faults(spec):
        mode = fault.get("mode")
        if mode not in FAULT_MODES:
            raise ValueError("unsupported byte fault mode")
        if fault.get("scope", "all") != "all":
            raise ValueError("byte faults precede application channels and require scope=all")
        if fault.get("placement", "source_bytes") not in ("source_bytes", "generator_output_bytes", "conditioner_input_bytes", "conditioner_output_bytes", "reseed_input_bytes"):
            raise ValueError("invalid byte fault placement")
        start, end = _interval(fault, "bytes")
        if mode in ("copy_bit", "stuck_bit", "entropy_reduction"):
            source = _BitFaultSource(source, fault, max_bytes=max_bytes)
        else:
            mode = {"zero_burst": "stuck_zero"}.get(mode, mode)
            source = DegradationSource(source, mode=mode, start=start, end=end,
                                       severity=_probability(fault.get("severity", 1.), "severity"),
                                       final_severity=_probability(fault.get("final_severity", fault.get("severity", 1.)), "final_severity"),
                                       seed=fault.get("seed", 0), period=fault.get("period", 32),
                                       max_bytes=max_bytes)
    return _PlacedFaultSource(source, spec, max_bytes)


def _control_uniforms(rows, channels, seed, purpose="pattern"):
    """Counter-addressed SHA256 controls avoid hidden mutable PRNG state."""
    _integer(seed, "seed")
    return np.array([[int.from_bytes(hashlib.sha256(f"rnglab12:{purpose}:{seed}:{int(row)}:{int(c)}".encode()).digest()[:8], "little") >> 11
                      for c in channels] for row in rows], dtype=np.float64).reshape(-1, len(channels)) / 2**53


def _uniform_pattern(fault, rows, channels, mode, start):
    period = _integer(fault.get("period", 4), "period", 1, 4096)
    pattern = fault.get("pattern")
    relative = rows - start
    if mode == "autocorrelation":
        return _control_uniforms(relative // period, channels, fault.get("seed", 0))
    if pattern is None:
        return _control_uniforms(relative % period, channels, fault.get("seed", 0))
    pattern = np.asarray(pattern, dtype=float)
    if pattern.ndim == 1:
        pattern = np.repeat(pattern[:, None], len(channels), axis=1)
    if pattern.shape != (period, len(channels)) or not np.isfinite(pattern).all() or np.any((pattern < 0) | (pattern >= 1)):
        raise ValueError("pattern must contain period x selected_channels finite values in [0,1)")
    return pattern[relative % period]


def apply_uniform_fault(u, start_row: int, spec: Mapping | None, application: str):
    """Return a copy with explicit post-uniform channel interventions.

    Periodicity repeats a supplied frozen pattern or a deterministic synthetic
    control pattern. It does not silently read unseen earlier input chunks.
    Autocorrelation replaces values by deterministic repeated control blocks.
    To study repetition of native variates, provide the per-seed first-period
    native values as pattern. This contract makes arbitrary partitions exact.
    """
    _integer(start_row, "start_row")
    out = np.array(u, dtype=np.float64, copy=True)
    if out.ndim != 2 or not np.isfinite(out).all() or np.any((out < 0) | (out >= 1)):
        raise ValueError("uniforms must be a finite two-dimensional array in [0,1)")
    application = _ALIASES.get(application, application)
    mapping = APPLICATION_CHANNELS.get(application)
    if mapping is None:
        raise ValueError("application channel semantics are not registered")
    expected = max(c for values in mapping.values() for c in values) + 1
    if out.shape[1] != expected:
        raise ValueError("uniform channel count does not match application semantics")
    if spec is None:
        return out
    absolute_rows = np.arange(start_row, start_row + len(out), dtype=np.int64)
    for fault in _faults(spec):
        mode = fault.get("mode")
        if mode not in FAULT_MODES:
            raise ValueError("unsupported uniform fault mode")
        if fault.get("placement", "post_uniform") not in ("post_uniform", "uniform"):
            raise ValueError("uniform faults require post_uniform placement")
        scope = fault.get("scope", "all")
        if scope == "all":
            channels = list(range(expected))
        elif scope in mapping:
            channels = mapping[scope]
        else:
            raise ValueError(f"scope {scope!r} is not defined for {application}")
        start, end = _interval(fault, "rows")
        active = absolute_rows >= start
        if end is not None:
            active &= absolute_rows < end
        indices = np.flatnonzero(active)
        rows = absolute_rows[active]
        # Validate even empty intervals so invalid settings cannot become valid
        # solely because the requested chunk does not intersect the fault.
        severity = _probability(fault.get("severity", 1.), "severity")
        final = _probability(fault.get("final_severity", severity), "final_severity")
        if final != severity and end is None:
            raise ValueError("drift requires a finite duration/end")
        if mode == "drift" and end is None:
            raise ValueError("drift requires a finite duration/end")
        fraction = np.zeros(len(rows)) if end is None else (rows - start) / max(1, end - start - 1)
        strength = severity + fraction * (final - severity)
        values = out[np.ix_(indices, channels)]
        if mode in ("periodicity", "autocorrelation"):
            changed = _uniform_pattern(fault, rows, channels, mode, start)
        elif mode in ("zero_burst", "stuck_zero"):
            changed = np.zeros_like(values)
        elif mode in ("bias", "drift"):
            changed = values * (1. - strength[:, None])
        else:
            bits = np.floor(values * 2**53).astype(np.uint64)
            if mode == "entropy_reduction":
                keep = _integer(fault.get("keep_bits", 8), "keep_bits", 0, 53)
                shift = np.uint64(53 - keep)
                bits = (bits >> shift) << shift
            elif mode == "bit_flip":
                bit = _integer(fault.get("bit", 0), "bit", 0, 52)
                bits ^= np.uint64(1 << bit)
            else:
                bit = _integer(fault.get("bit", 0), "bit", 0, 52)
                mask = np.uint64(1 << bit)
                replacement = (_integer(fault.get("value", 0), "value", 0, 1)
                               if mode == "stuck_bit" else None)
                if mode == "copy_bit":
                    source_bit = _integer(fault.get("from_bit", 1), "from_bit", 0, 52)
                    replacement = (bits >> np.uint64(source_bit)) & np.uint64(1)
                bits = (bits & ~mask) | (np.asarray(replacement, dtype=np.uint64) << np.uint64(bit))
            changed = bits.astype(float) / 2**53
        # Bias/drift severity defines the shrink transformation. For remaining
        # interventions it is the deterministic addressable substitution rate.
        if mode not in ("bias", "drift") and (severity < 1 or final < 1):
            control = _control_uniforms(rows, channels, fault.get("seed", 0), purpose="substitution_gate")
            changed = np.where(control < strength[:, None], changed, values)
        out[np.ix_(indices, channels)] = changed
    return out
