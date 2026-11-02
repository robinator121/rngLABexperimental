"""Explicit opt-in, operation-indexed synthetic faults. Never patches an OS RNG."""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import math

import numpy as np

from rnglab14.sources.base import RandomSource, SourceManifest
from rnglab_crypto.sources import ReplaySource


PLACEMENTS = {
    "reduced_entropy": "drbg_initialization",
    "seed_reuse": "drbg_initialization",
    "fixed_nonce": "downstream_nonce_api",
    "bias": "raw_source",
    "serial": "raw_source",
    "periodic": "raw_source",
    "rollback": "drbg_state",
    "starvation": "raw_source",
    "drift": "raw_source",
    "partial_exposure": "bounded_seed_state_model",
}


@dataclass(frozen=True)
class FaultSpec:
    kind: str
    onset: int = 0
    duration: int = 256
    intensity: float = 1.0
    entropy_bits: int = 8
    exposed_bits: int = 4
    period_bytes: int = 2
    recovery: str = "resume_shadow"
    target_layer: str | None = None

    def __post_init__(self):
        if self.kind not in PLACEMENTS:
            raise ValueError("unknown fault")
        for name, lo, hi in (("onset", 0, 100_000), ("duration", 1, 100_000),
                             ("entropy_bits", 0, 256), ("exposed_bits", 0, 256),
                             ("period_bytes", 1, 4096)):
            v = getattr(self, name)
            if type(v) is not int or not lo <= v <= hi:
                raise ValueError(f"invalid {name}")
        if type(self.intensity) not in (float, int) or not math.isfinite(self.intensity) or not 0 <= self.intensity <= 1:
            raise ValueError("intensity must be finite in [0,1]")
        if self.exposed_bits > self.entropy_bits or self.recovery != "resume_shadow":
            raise ValueError("invalid exposure/recovery")
        if self.target_layer is not None and self.target_layer not in set(PLACEMENTS.values()):
            raise ValueError("unknown target layer")

    @property
    def placement(self):
        return self.target_layer or PLACEMENTS[self.kind]


class EntropyUnavailable(RuntimeError):
    pass


class FaultAdapter(RandomSource):
    """One read = one operation; exact onset/duration including failed reads.

    The healthy shadow advances throughout the fault. Recovery returns that
    shadow's next block; this is a simulation intervention, not a secure repair
    recipe. Bias/serial/drift parameters describe bits, not byte shifts.
    Rollback restores the INITIAL synthetic DRBG state at onset. Seed reuse
    repeats one seed block. Partial exposure leaves output unchanged and records
    the known prefix of a deliberately bounded seed, not native DRBG internals.
    """
    def __init__(self, base: ReplaySource, fault: FaultSpec, *, lab_only: bool = False):
        if lab_only is not True or not isinstance(base, ReplaySource):
            raise ValueError("faults require an explicitly opted-in synthetic ReplaySource")
        self.base, self.fault = base, fault
        self.operation = 0
        self.last_bit = None
        self.reused = None
        self.rollback = None
        self.leak = None
        self.rng = np.random.Generator(np.random.PCG64(int.from_bytes(base.seed, "big")))
        self.pattern = hashlib.shake_256(base.seed).digest(fault.period_bytes)
        self.pattern_offset = 0
        super().__init__(SourceManifest("fault-"+base.manifest.source_id, "synthetic_fault",
                         True, False, created_at="deterministic-replay",
                         metadata={**asdict(fault), "placement": fault.placement,
                                   "operation_units": "nonempty read attempts", "entropy_credit_bits": 0}),
                         max_bytes=base.max_bytes)

    def _read_bytes(self, n):
        f = self.fault
        # Validate shape before advancing any state.
        if f.kind in ("reduced_entropy", "partial_exposure") and f.entropy_bits > 8*n:
            raise ValueError("seed width smaller than bounded entropy support")
        if f.kind == "seed_reuse" and self.reused is not None and len(self.reused) != n:
            raise ValueError("seed-reuse fault requires fixed-size reads")
        raw = self.base.read_bytes(n)
        op = self.operation
        self.operation += 1
        if not f.onset <= op < f.onset + f.duration or f.intensity == 0:
            return raw
        if f.kind == "starvation":
            raise EntropyUnavailable("synthetic source unavailable; no substitute output")
        if f.kind in ("reduced_entropy", "partial_exposure"):
            value = int.from_bytes(raw, "big") & ((1 << f.entropy_bits)-1)
            if f.kind == "partial_exposure":
                self.leak = {"known_prefix": value >> (f.entropy_bits-f.exposed_bits),
                             "unknown_suffix_bits": f.entropy_bits-f.exposed_bits}
            return value.to_bytes(n, "big")
        if f.kind == "seed_reuse":
            if self.reused is None:
                self.reused = raw
            return self.reused
        if f.kind == "fixed_nonce":
            return bytes(n)
        if f.kind == "rollback":
            if self.rollback is None:
                self.rollback = self.base.restart()
            return self.rollback.read_bytes(n)
        if f.kind == "periodic":
            out = bytes(self.pattern[(self.pattern_offset+i) % len(self.pattern)] for i in range(n))
            self.pattern_offset += n
            return out
        uniforms = self.rng.random(8*n)
        strength = f.intensity
        if f.kind == "drift":
            strength *= (op-f.onset+1)/f.duration
        if f.kind in ("bias", "drift"):
            bits = uniforms < .5 + .5*strength
        else:
            bits = np.empty(8*n, dtype=np.uint8)
            for j, u in enumerate(uniforms):
                p = .5 if self.last_bit is None else (.5 + strength*(self.last_bit-.5))
                self.last_bit = int(u < p)
                bits[j] = self.last_bit
        return np.packbits(bits).tobytes()
