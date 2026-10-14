"""Explicit research-source families; declarations do not authenticate physics.

The unmodified v0.11 byte sources implement bounded I/O, conditioning, and the
research HMAC_DRBG. This module adds a strict acquisition-manifest boundary and
fixed-size internal Generate requests, so application read partitioning does
not change a deterministic comparison stream or its reseed schedule.
"""
from __future__ import annotations

from datetime import datetime
import hashlib
import math
from pathlib import Path
import re
from typing import Any, Mapping

from rnglab11.conditioning import ConditionedSource, HashConditioner, VonNeumannConditioner
from rnglab11.drbg import ReseedingDRBG
from rnglab11.sources import (ByteSource, FileSource, SourceExhausted,
                             SyntheticSource, SystemCSPRNGSource, UniformStream)

SOURCE_FAMILIES = (
    "raw_physical", "conditioned_physical", "standalone_csprng",
    "physical_seeded_csprng", "physical_reseeded_csprng",
    "deterministic_seeded_or_reseeded_csprng", "synthetic_fault",
)
SOURCE_TYPES = SOURCE_FAMILIES
_PLACEHOLDERS = {"", "todo", "unknown", "unspecified", "none", "null", "n/a", "na",
                 "placeholder", "test", "dummy", "example", "synthetic", "fixture"}


def _integer(value, name, minimum=0, maximum=None):
    if type(value) is not int or value < minimum or (maximum is not None and value > maximum):
        raise ValueError(f"{name} must be an integer in the permitted range")
    return value


def source_type(spec: Mapping[str, Any]) -> str:
    """Return the declared architecture, rejecting ambiguous classifications."""
    if not isinstance(spec, Mapping):
        raise ValueError("source specification must be a mapping")
    family = spec.get("family", spec.get("source_type"))
    if family not in SOURCE_FAMILIES:
        raise ValueError(f"family must be one of {SOURCE_FAMILIES}")
    if "family" in spec and "source_type" in spec and spec["family"] != spec["source_type"]:
        raise ValueError("contradictory source family labels")
    return family


def verify_content_file(path: str | Path, content_sha256: str, *, max_input_bytes=67_108_864):
    """Check an entire finite file, not merely a consumed prefix or its metadata.

    This generic integrity check is also testable on explicitly synthetic test
    fixtures. It establishes neither acquisition authenticity nor entropy.
    """
    _integer(max_input_bytes, "max_input_bytes", 1)
    if not isinstance(content_sha256, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", content_sha256):
        raise ValueError("a declared full-file SHA256 is required")
    if len(set(content_sha256.lower())) == 1:
        raise ValueError("placeholder digest rejected")
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise ValueError("physical input must be a regular file")
    before = path.stat()
    if not 0 < before.st_size <= max_input_bytes:
        raise ValueError("input is empty or exceeds its declared verification budget")
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while block := handle.read(1_048_576):
            digest.update(block)
    after = path.stat()
    identity = lambda s: (s.st_size, s.st_mtime_ns, s.st_dev, s.st_ino)
    if identity(before) != identity(after):
        raise ValueError("input changed during full-file verification")
    if digest.hexdigest() != content_sha256.lower():
        raise ValueError("physical input content SHA256 does not match declaration")
    return path, identity(after)


def validate_physical_manifest(spec: Mapping[str, Any]):
    """Validate a caller's acquisition declaration and independently verify bytes.

    Even a complete truthful-looking declaration is not independent evidence of
    physical origin. Explicit synthetic fixtures can never enter this family.
    """
    if spec.get("test_fixture") or spec.get("synthetic"):
        raise ValueError("synthetic/test fixtures cannot be classified as physical observations")
    acquisition = spec.get("acquisition")
    if not isinstance(acquisition, Mapping):
        raise ValueError("acquisition instrument, time, units and origin are required")
    for key in ("instrument", "time", "units", "origin"):
        value = acquisition.get(key)
        if not isinstance(value, str) or value.strip().lower() in _PLACEHOLDERS:
            raise ValueError(f"nonplaceholder acquisition {key} is required")
        lowered = value.lower()
        if any(marker in lowered for marker in ("todo", "placeholder", "synthetic fixture", "test fixture")):
            raise ValueError(f"placeholder acquisition {key} rejected")
    try:
        timestamp = datetime.fromisoformat(acquisition["time"].replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError("acquisition time must be ISO8601 with a timezone") from exc
    if timestamp.tzinfo is None:
        raise ValueError("acquisition time requires a timezone")
    if not spec.get("path"):
        raise ValueError("physical data path is required")
    path, identity = verify_content_file(spec["path"], spec.get("content_sha256"),
                                         max_input_bytes=spec.get("max_input_bytes", 67_108_864))
    return {"path": str(path), "content_sha256": spec["content_sha256"].lower(),
            "acquisition": dict(acquisition), "verified_file_identity": list(identity),
            "integrity_status": "full_file_sha256_verified",
            "physical_origin_status": "declared_by_operator_not_independently_authenticated",
            "physical_entropy_assessed": False, "entropy_claim": "none"}


class ClassifiedSource(ByteSource):
    """Keep classification separate from the underlying mechanism's metadata."""

    def __init__(self, source: ByteSource, family: str, *, details=None, max_bytes=1_048_576):
        self.source = source
        super().__init__({"family": family, "implementation": "rnglab12.sources.v1",
                          "origin_source_id": source.source_id, "entropy_claim": "none",
                          **(details or {})}, max_bytes=max_bytes)

    def _read(self, n):
        try:
            return self.source.read(n)
        except SourceExhausted as exc:
            return exc.partial

    def provenance(self):
        origin = self.source.provenance()
        result = {**super().provenance(), "family": self.provenance_family,
                  "origin": origin}
        if hasattr(self.source, "request_bytes"):
            result["internal_generate_request_bytes"] = self.source.request_bytes
            result["generator"] = origin
        return result

    @property
    def provenance_family(self):
        import json
        return json.loads(self._metadata_json)["family"]

    def close(self):
        self.source.close()


class FixedRequestDRBG(ByteSource):
    """One Generate unit is request_bytes internally, never an application step.

    Ordinary read sizes only drain this fixed request buffer. Reseeding occurs
    after N completed internal Generate requests and before request N+1, using
    the v0.11 safe HMAC_DRBG update. Pending bytes and key state persist across
    reads. No secret state is put in routine provenance.
    """

    def __init__(self, seed_source, *, reseed_source=None, interval=None,
                 request_bytes=4096, max_bytes=1_048_576,
                 initial_seed_label="deterministic_test_seed", reseed_label="none"):
        _integer(request_bytes, "request_bytes", 1, 65_536)
        if interval is not None:
            _integer(interval, "reseed_interval_generate_calls", 1)
        self.request_bytes = request_bytes
        self._buffer = bytearray()
        self.external_read_calls = 0
        self.seed_source = seed_source
        self.reseed_source = reseed_source if reseed_source is not None else seed_source
        internal_budget = math.ceil(max_bytes / request_bytes) * request_bytes
        self.drbg = ReseedingDRBG(seed_source, reseed_source=self.reseed_source,
                                 interval=interval, configuration="R0", mixing="safe_mix",
                                 max_bytes=internal_budget,
                                 max_cycles=max(1, math.ceil(max_bytes / request_bytes)),
                                 initial_seed_label=initial_seed_label, reseed_label=reseed_label)
        super().__init__({"adapter": "fixed_request_research_hmac_drbg_sha256",
                          "internal_generate_request_bytes": request_bytes,
                          "reseed_interval_generate_calls": interval,
                          "reseed_unit": "internal_nonempty_Generate_requests",
                          "application_read_partition_invariant": True,
                          "initial_seed_label": initial_seed_label,
                          "reseed_label": reseed_label, "safe_state_mix": True,
                          "validated_cryptographic_product": False,
                          "physical_entropy_assessed": False}, max_bytes=max_bytes)

    def _read(self, n):
        while len(self._buffer) < n:
            self._buffer.extend(self.drbg.read(self.request_bytes))
        result = bytes(self._buffer[:n])
        del self._buffer[:n]
        self.external_read_calls += 1
        return result

    def provenance(self):
        return {**super().provenance(), "external_read_calls": self.external_read_calls,
                "buffered_generated_bytes": len(self._buffer),
                "internal_generate_calls": self.drbg.cycles,
                "state_continuity": "same_DRBG_instance_and_pending_buffer_across_reads",
                "generator": self.drbg.provenance()}

    def close(self):
        self.seed_source.close()
        if self.reseed_source is not self.seed_source:
            self.reseed_source.close()


def _physical_source(spec, max_bytes):
    manifest = validate_physical_manifest(spec)
    source = FileSource(manifest["path"], kind="declared_physical_acquisition",
                        metadata=manifest, max_bytes=max_bytes)
    if list(source._identity) != manifest["verified_file_identity"]:
        source.close()
        raise ValueError("physical file changed between verification and opening")
    return source, manifest


def _physical_seed(spec, max_bytes):
    if source_type(spec) not in ("raw_physical", "conditioned_physical"):
        raise ValueError("physical seed material requires a validated physical source family")
    return make_source(spec, max_bytes=max_bytes)


def make_source(spec: Mapping[str, Any], *, max_bytes=1_048_576) -> ByteSource:
    """Build a bounded, explicitly classified source without hidden fallbacks.

    No physical file is supplied by this release. Runtime comparisons with a
    deterministic seed use the deterministic family, even though the output
    mechanism is a cryptographic generator. MT19937 itself is never a CSPRNG.
    """
    _integer(max_bytes, "max_bytes", 0)
    family = source_type(spec)
    details = {}
    if family in ("raw_physical", "conditioned_physical"):
        input_budget = spec.get("max_input_bytes", 67_108_864)
        source, details = _physical_source(spec, input_budget if family == "conditioned_physical" else max_bytes)
        if family == "conditioned_physical":
            treatment = spec.get("conditioner", "sha256")
            if treatment == "sha256":
                conditioner = HashConditioner(spec.get("block_bytes", 64))
            elif treatment == "von_neumann":
                conditioner = VonNeumannConditioner()
            else:
                source.close()
                raise ValueError("conditioner must be sha256 or von_neumann")
            source = ConditionedSource(source, conditioner, max_input_bytes=input_budget,
                                       max_output_bytes=max_bytes,
                                       chunk_bytes=spec.get("chunk_bytes", 4096))
            details["conditioning_is_not_entropy_certification"] = True
    elif family == "standalone_csprng":
        if "seed" in spec:
            raise ValueError("a known deterministic seed requires the deterministic source family")
        source = SystemCSPRNGSource(max_bytes=max_bytes)
        details["mechanism"] = "operating_system_csprng"
    elif family in ("physical_seeded_csprng", "physical_reseeded_csprng",
                    "deterministic_seeded_or_reseeded_csprng"):
        interval = spec.get("reseed_interval_generate_calls")
        if family == "physical_seeded_csprng" and interval is not None:
            raise ValueError("periodic reseeding requires physical_reseeded_csprng")
        if family == "physical_reseeded_csprng" and interval is None:
            raise ValueError("physical_reseeded_csprng requires an explicit positive interval")
        if interval is not None:
            _integer(interval, "reseed_interval_generate_calls", 1)
        seed_budget = spec.get("max_seed_source_bytes", 1_048_576)
        _integer(seed_budget, "max_seed_source_bytes", 48)
        if family == "deterministic_seeded_or_reseeded_csprng":
            if "seed" not in spec or type(spec["seed"]) not in (int, str):
                raise ValueError("an explicit deterministic int/string research seed is required")
            seed_source = SyntheticSource(spec["seed"], max_bytes=seed_budget)
            reseed_source = (SyntheticSource(spec["reseed_seed"], max_bytes=seed_budget)
                             if "reseed_seed" in spec else seed_source)
            initial_label = "deterministic_MT19937_test_material_not_assessed_entropy"
            reseed_label = initial_label if interval is not None else "none"
        else:
            if not isinstance(spec.get("seed_source"), Mapping):
                raise ValueError("an explicit validated seed_source manifest is required")
            seed_source = _physical_seed(spec["seed_source"], seed_budget)
            try:
                reseed_source = (_physical_seed(spec["reseed_source"], seed_budget)
                                 if "reseed_source" in spec else seed_source)
            except Exception:
                seed_source.close()
                raise
            initial_label = "declared_physical_material_integrity_verified_entropy_unassessed"
            reseed_label = initial_label if interval is not None else "none"
        try:
            source = FixedRequestDRBG(seed_source, reseed_source=reseed_source, interval=interval,
                                      request_bytes=spec.get("request_bytes", 4096), max_bytes=max_bytes,
                                      initial_seed_label=initial_label, reseed_label=reseed_label)
        except Exception:
            seed_source.close()
            if reseed_source is not seed_source:
                reseed_source.close()
            raise
    else:
        from .faults import make_byte_fault
        base = spec.get("base", {"family": "deterministic_seeded_or_reseeded_csprng", "seed": spec.get("seed", 0)})
        if source_type(base) == "synthetic_fault":
            raise ValueError("nest combinations in fault.faults, not recursive source families")
        source = make_source(base, max_bytes=max_bytes)
        try:
            source = make_byte_fault(source, spec.get("fault", {"mode": "zero_burst"}), max_bytes=max_bytes)
        except Exception:
            source.close()
            raise
        details["synthetic_intervention_not_measured_entropy_degradation"] = True
    return ClassifiedSource(source, family, details=details, max_bytes=max_bytes)
