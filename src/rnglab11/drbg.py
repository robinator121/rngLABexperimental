"""Research HMAC_DRBG (SP 800-90A SHA256 mechanism) and reseeding controls.

This module is not a validated cryptographic product or complete 90C RBG. Seed
length checks cannot establish entropy. Use the OS/library RNG in production.
"""
from __future__ import annotations

import hashlib
import hmac
from typing import Any

from .sources import ByteSource, SourceBudgetExceeded, SourceExhausted


RESEED_CONFIGURATIONS = {"R0": None, "R1": 100, "R2": 1_000, "R3": 10_000,
                         "R4": 100_000, "R5": 1_000_000}


class InitializationUnavailable(RuntimeError):
    """Initialization consumed an incomplete source segment; no seed is exposed."""
    def __init__(self, consumed_bytes):
        self.consumed_bytes=consumed_bytes
        super().__init__(f"Initial seed source exhausted after {consumed_bytes} bytes; no generator created")


class ReseedUnavailable(RuntimeError):
    """Reseed input was incomplete; no generator output or seed is exposed."""

    def __init__(self, consumed_bytes: int):
        self.consumed_bytes = consumed_bytes
        super().__init__(f"reseed source exhausted after consuming {consumed_bytes} input bytes; no output generated")


class HMACDRBG(ByteSource):
    """SP 800-90A §10.1.2 instantiate/update/reseed/generate, SHA256.

    Each successful nonempty read is one Generate call. Splitting a request
    changes the stream because Generate updates state after each call.
    Prediction resistance and a production entropy-source assessment are absent.
    """

    max_request_bytes = 65_536  # 2**19 bits, 90A per-request limit

    def __init__(self, entropy_input: bytes, nonce: bytes, *, personalization: bytes = b"",
                 security_strength: int = 256, max_bytes: int = 1_048_576):
        if security_strength not in (128, 192, 256):
            raise ValueError("security_strength must be 128, 192, or 256")
        if len(entropy_input) * 8 < security_strength or len(nonce) * 8 < security_strength // 2:
            raise ValueError("seed/nonce length is below selected strength's length requirement")
        self._key, self._value = b"\0" * 32, b"\1" * 32
        self.reseed_counter = 1
        self.generate_calls = self.reseed_events = 0
        self.security_strength = security_strength
        self._update(bytes(entropy_input) + bytes(nonce) + bytes(personalization))
        # Never place seed hashes or state in routine metadata: degraded seeds
        # can have a small dictionary and hashes would disclose them.
        super().__init__({"adapter": "research_hmac_drbg_sha256", "mechanism": "SP800-90A-10.1.2",
                          "selected_strength_bits": security_strength, "entropy_assessed": False,
                          "validated_implementation": False}, max_bytes=max_bytes)

    def _mac(self, data: bytes):
        return hmac.new(self._key, data, hashlib.sha256).digest()

    def _update(self, provided_data: bytes = b""):
        self._key = self._mac(self._value + b"\0" + provided_data)
        self._value = self._mac(self._value)
        if provided_data:
            self._key = self._mac(self._value + b"\1" + provided_data)
            self._value = self._mac(self._value)

    def reseed(self, entropy_input: bytes, *, additional_input: bytes = b""):
        if len(entropy_input) * 8 < self.security_strength:
            raise ValueError("reseed material too short for selected strength")
        self._update(bytes(entropy_input) + bytes(additional_input))
        self.reseed_counter = 1
        self.reseed_events += 1

    def read(self, n: int, *, additional_input: bytes = b"") -> bytes:
        if not isinstance(n, int) or isinstance(n, bool) or not 0 <= n <= self.max_request_bytes:
            raise ValueError("Generate request must be between 0 and 65536 bytes")
        if n == 0:
            return b""
        if self.reseed_counter > 2**48:
            raise RuntimeError("90A mechanism reseed interval exhausted")
        if len(additional_input) > 65_536:
            raise ValueError("research additional-input limit exceeded")
        if self.max_bytes is not None and self.bytes_consumed + n > self.max_bytes:
            raise SourceBudgetExceeded("DRBG output budget exceeded")
        # Keep additional input local to this call; update before and after.
        if additional_input:
            self._update(additional_input)
        output = bytearray()
        while len(output) < n:
            self._value = self._mac(self._value)
            output.extend(self._value)
        self._update(additional_input)
        self.reseed_counter += 1
        self.generate_calls += 1
        data = bytes(output[:n])
        self.bytes_consumed += n
        self.offset += n
        self._consumed_digest.update(data)
        return data

    generate = read

    def provenance(self):
        return {**super().provenance(), "generate_calls": self.generate_calls,
                "reseed_events": self.reseed_events, "reseed_counter": self.reseed_counter,
                "consumed_digest_scope": "since_state_restore" if
                getattr(self, "_digest_resumed_without_prefix", False) else "instance_lifetime"}

    def export_sensitive_state(self, *, allow_sensitive_research_state: bool = False) -> dict[str, Any]:
        if not allow_sensitive_research_state:
            raise PermissionError("plaintext state export requires explicit research-only opt-in")
        return {"mechanism": "HMAC_DRBG_SHA256", "key_hex": self._key.hex(),
                "value_hex": self._value.hex(), "reseed_counter": self.reseed_counter,
                "generate_calls": self.generate_calls, "reseed_events": self.reseed_events,
                "security_strength": self.security_strength,
                "generated_bytes": self.bytes_consumed, "max_bytes": self.max_bytes,
                "warning": "SENSITIVE deterministic research state; never publish real keys"}

    @classmethod
    def from_sensitive_state(cls, state, *, allow_sensitive_research_state: bool = False):
        if not allow_sensitive_research_state:
            raise PermissionError("plaintext state import requires explicit research-only opt-in")
        if state.get("mechanism") != "HMAC_DRBG_SHA256":
            raise ValueError("unsupported state mechanism")
        key, value = bytes.fromhex(state["key_hex"]), bytes.fromhex(state["value_hex"])
        if len(key) != 32 or len(value) != 32 or not 1 <= state["reseed_counter"] <= 2**48 + 1:
            raise ValueError("invalid DRBG state")
        if (not isinstance(state["generated_bytes"], int) or state["generated_bytes"] < 0
                or state["generated_bytes"] > state["max_bytes"]
                or state["generate_calls"] < 0 or state["reseed_events"] < 0):
            raise ValueError("invalid DRBG state counters")
        obj = cls(b"\0" * 32, b"\0" * 16, security_strength=state["security_strength"],
                  max_bytes=state["max_bytes"])
        obj._key, obj._value = key, value
        obj.reseed_counter = state["reseed_counter"]
        obj.generate_calls, obj.reseed_events = state["generate_calls"], state["reseed_events"]
        obj.bytes_consumed = obj.offset = state["generated_bytes"]
        obj._digest_resumed_without_prefix = True
        return obj


class ReseedingDRBG(ByteSource):
    """HMAC_DRBG architecture with explicit finite source/output/cycle budgets.

    A cycle is one successful nonempty Generate call, independent of byte count
    or application iterations. After N completed cycles reseed occurs BEFORE
    call N+1. Source exhaustion blocks the call; it does not silently skip reseed.
    """

    def __init__(self, seed_source: ByteSource, *, reseed_source: ByteSource | None = None,
                 configuration: str = "R0", interval: int | None = None,
                 mixing: str = "safe_mix", allow_unsafe_control: bool = False,
                 max_bytes: int = 1_048_576, max_cycles: int = 100_000,
                 seed_bytes: int = 32, nonce_bytes: int = 16, reseed_bytes: int = 32,
                 initial_seed_label: str = "unassessed", reseed_label: str = "unassessed"):
        for label, value, minimum in (("max_bytes", max_bytes, 0), ("max_cycles", max_cycles, 1),
                                      ("seed_bytes", seed_bytes, 32), ("nonce_bytes", nonce_bytes, 16),
                                      ("reseed_bytes", reseed_bytes, 32)):
            if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
                raise ValueError(f"{label} must be an integer at least {minimum}")
        if configuration not in RESEED_CONFIGURATIONS:
            raise ValueError("unknown reseeding configuration")
        selected_interval = interval if interval is not None else RESEED_CONFIGURATIONS[configuration]
        if selected_interval is not None and (not isinstance(selected_interval, int) or selected_interval < 1):
            raise ValueError("reseed interval must be a positive integer or R0")
        if mixing not in ("safe_mix", "unsafe_replace"):
            raise ValueError("unknown reseed mixing mode")
        if mixing == "unsafe_replace" and not allow_unsafe_control:
            raise ValueError("unsafe replacement requires explicit negative-control opt-in")
        if seed_bytes + nonce_bytes > 65_536 or reseed_bytes > 65_536:
            raise ValueError("research seed request limit exceeded")
        self.seed_source = seed_source
        self.reseed_source = reseed_source if reseed_source is not None else seed_source
        self._seed_source_before = seed_source.bytes_consumed
        self._reseed_source_before = self.reseed_source.bytes_consumed
        self.interval, self.mixing = selected_interval, mixing
        self.reseed_bytes, self.max_cycles = reseed_bytes, max_cycles
        self.cycles = self.application_iterations = self.reseed_events = 0
        self._last_reseed_cycle = 0
        self.events = []
        self._reseed_blocked = False
        # One exact read prevents a half-created initial DRBG from being used.
        start = seed_source.offset
        try:
            material = seed_source.read(seed_bytes + nonce_bytes)
        except SourceExhausted as exc:
            raise InitializationUnavailable(len(exc.partial)) from None
        self.drbg = HMACDRBG(material[:seed_bytes], material[seed_bytes:], max_bytes=max_bytes)
        self.initial_source_segment = {"source_id": seed_source.source_id, "start": start,
                                       "end": start + len(material), "bytes": len(material)}
        super().__init__({"adapter": "reseeding_research_hmac_drbg", "configuration": configuration,
                          "interval_generate_calls": selected_interval, "mixing": mixing,
                          "unsafe_negative_control": mixing == "unsafe_replace",
                          "initial_seed_label": initial_seed_label, "reseed_label": reseed_label,
                          "seed_source_id": seed_source.source_id,
                          "reseed_source_id": self.reseed_source.source_id,
                          "entropy_assessed": False, "validated_implementation": False}, max_bytes=max_bytes)

    @property
    def source_bytes_consumed(self):
        if self.seed_source is self.reseed_source:
            return self.seed_source.bytes_consumed - self._seed_source_before
        return (self.seed_source.bytes_consumed - self._seed_source_before
                + self.reseed_source.bytes_consumed - self._reseed_source_before)

    def next_reseed(self):
        if self.interval is None:
            return {"status": "disabled", "generate_calls_until_reseed": None,
                    "source_segment": None}
        next_boundary = self._last_reseed_cycle + self.interval
        return {"status": "blocked_source_exhausted" if self._reseed_blocked else "scheduled",
                "after_completed_cycle": next_boundary,
                "before_generate_call": next_boundary + 1,
                "generate_calls_until_reseed": max(0, next_boundary - self.cycles),
                "source_segment": {"source_id": self.reseed_source.source_id,
                                   "start": self.reseed_source.offset,
                                   "end": self.reseed_source.offset + self.reseed_bytes,
                                   "bytes": self.reseed_bytes},
                "exposure_is_not_application_failure": True}

    def _perform_reseed(self):
        if self._reseed_blocked:
            raise RuntimeError("reseed source previously exhausted; explicitly reconstruct/recover before retry")
        start = self.reseed_source.offset
        try:
            incoming = self.reseed_source.read(self.reseed_bytes)
        except SourceExhausted as exc:
            self._reseed_blocked = True
            self.events.append({"event": "reseed_blocked", "cycle": self.cycles,
                                "source_start": start, "source_end": self.reseed_source.offset,
                                "partial_bytes_consumed": len(exc.partial)})
            # The underlying partial bytes are SEED material, not partial DRBG
            # output; never expose them through SourceExhausted.partial.
            raise ReseedUnavailable(len(exc.partial)) from None
        if self.mixing == "safe_mix":
            self.drbg.reseed(incoming)
        else:
            # Deliberately discard healthy previous state. This is a negative
            # experimental control, not the standardized reseed operation.
            remaining_budget = self.max_bytes - self.bytes_consumed
            self.drbg = HMACDRBG(incoming, b"\0" * 16, max_bytes=remaining_budget)
        self.reseed_events += 1
        self._last_reseed_cycle = self.cycles
        self.events.append({"event": "reseed", "cycle": self.cycles, "mixing": self.mixing,
                            "source_id": self.reseed_source.source_id,
                            "source_start": start, "source_end": self.reseed_source.offset,
                            "source_bytes": len(incoming)})

    def read(self, n):
        if not isinstance(n, int) or isinstance(n, bool) or not 0 <= n <= HMACDRBG.max_request_bytes:
            raise ValueError("Generate request must be between 0 and 65536 bytes")
        if n == 0:
            return b""
        if self.cycles >= self.max_cycles:
            raise SourceBudgetExceeded("DRBG cycle budget exceeded")
        return super().read(n)

    generate = read

    def _read(self, n):
        if self.interval is not None and self.cycles - self._last_reseed_cycle >= self.interval:
            self._perform_reseed()
        data = self.drbg.read(n)
        self.cycles += 1
        return data

    def record_application_iterations(self, count: int):
        if not isinstance(count, int) or isinstance(count, bool) or count < 0:
            raise ValueError("application iterations must be nonnegative")
        self.application_iterations += count

    def provenance(self):
        return {**super().provenance(), "generated_bytes": self.bytes_consumed,
                "source_bytes_consumed": self.source_bytes_consumed,
                "cycles": self.cycles, "application_iterations": self.application_iterations,
                "reseed_events": self.reseed_events, "events": list(self.events),
                "initial_source_segment": self.initial_source_segment,
                "next_reseed": self.next_reseed()}


class OutputFaultSource(ByteSource):
    """Post-generator faults leave internal generator state untouched."""

    def __init__(self, source: ByteSource, *, mode: str = "zero", start_offset: int = 0,
                 xor_mask: int = 1, max_bytes: int = 1_048_576):
        if mode not in ("zero", "xor") or not 0 <= xor_mask <= 255 or start_offset < 0:
            raise ValueError("invalid output fault")
        self.source, self.mode = source, mode
        self.start_offset, self.xor_mask = start_offset, xor_mask
        self.faulted_bytes = 0
        super().__init__({"adapter": "output_side_fault", "mode": mode,
                          "start_offset": start_offset, "xor_mask": xor_mask,
                          "origin_source_id": source.source_id, "research_control": True}, max_bytes=max_bytes)

    def _transform(self, raw):
        out = bytearray(raw)
        first = max(0, self.start_offset - self.offset)
        for i in range(first, len(out)):
            out[i] = 0 if self.mode == "zero" else out[i] ^ self.xor_mask
        self.faulted_bytes += max(0, len(out) - first)
        return bytes(out)

    def _read(self, n):
        try:
            raw = self.source.read(n)
        except SourceExhausted as exc:
            return self._transform(exc.partial)
        return self._transform(raw)

    def provenance(self):
        return {**super().provenance(), "faulted_bytes": self.faulted_bytes,
                "origin": self.source.provenance()}
