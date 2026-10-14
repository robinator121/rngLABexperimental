"""Stateful, chunk-invariant faults on source bits; no floating-point faults."""
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class FaultSpec:
    kind: str
    strength: float
    onset_bit: int = 0
    ramp_bits: int = 0
    pattern: str = "0011"
    frozen_bit: int = 0

    def __post_init__(self):
        if self.kind not in {"bias", "freeze", "periodic", "autocorrelation"}:
            raise ValueError(f"Unknown fault: {self.kind}")
        if not 0 <= self.strength <= 1:
            raise ValueError("strength must be in [0, 1]")
        if not isinstance(self.onset_bit, int) or not isinstance(self.ramp_bits, int) or min(self.onset_bit, self.ramp_bits) < 0:
            raise ValueError("onset_bit and ramp_bits must be nonnegative integers")
        if not self.pattern or set(self.pattern) - {"0", "1"} or self.frozen_bit not in (0, 1):
            raise ValueError("pattern and frozen_bit must be binary")


class Fault:
    def __init__(self, spec: FaultSpec, seed):
        self.spec = spec
        self.rng = np.random.Generator(np.random.PCG64(seed))
        self.previous = None

    def apply(self, bits, offset):
        s = self.spec
        positions = np.arange(offset, offset + len(bits), dtype=np.int64)
        if s.ramp_bits:
            envelope = np.clip((positions - s.onset_bit) / s.ramp_bits, 0, 1)
        else:
            envelope = positions >= s.onset_bit
        mask = self.rng.random(len(bits)) < s.strength * envelope
        out = bits.copy()
        if s.kind == "bias":
            out[mask] = 1  # isolated IID input: P(1) = (1 + strength)/2
        elif s.kind == "freeze":
            out[mask] = s.frozen_bit
        elif s.kind == "periodic":
            pattern = np.array(list(s.pattern), dtype=np.uint8)
            out[mask] = pattern[positions[mask] % len(pattern)]
        elif len(out):
            # Each masked bit copies the PREVIOUS OUTPUT bit, including across calls.
            if self.previous is None:
                mask[0] = False
            refresh = np.where(~mask, np.arange(1, len(out) + 1), 0)
            origin = np.maximum.accumulate(refresh)
            out = np.concatenate(([0 if self.previous is None else self.previous], out))[origin].astype(np.uint8)
        if len(out):
            self.previous = int(out[-1])
        return out
