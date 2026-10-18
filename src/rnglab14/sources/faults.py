"""RNGLab v0.14: Controlled Fault Injection Sources.

Defines standardized adversarial, physical, and degradation fault models:
- Constant Bias (skewed P(1) != 0.5)
- First-order and Higher-order Autocorrelation (Markov coupling)
- Linear and Non-linear Drift
- Stuck Bit / Zero Burst hardware faults
- Periodic Coupling / Power-supply hum
- Intermittent / Recovering faults
"""
from __future__ import annotations

import math
import random
from typing import Optional

from .base import RandomSource, SourceManifest


class DegradedFaultSource(RandomSource):
    """Wraps a base random source and injects mathematically controlled defects."""

    def __init__(
        self,
        base_source: RandomSource,
        bias_delta: float = 0.0,
        autocorrelation_rho: float = 0.0,
        stuck_bit_index: Optional[int] = None,
        stuck_bit_value: int = 1,
        periodic_coupling_freq: float = 0.0,
        periodic_coupling_amp: float = 0.0,
        max_bytes: Optional[int] = None,
    ):
        if not all(math.isfinite(v) for v in [bias_delta, autocorrelation_rho, periodic_coupling_freq, periodic_coupling_amp]):
            raise ValueError("fault parameters must be finite")
        if periodic_coupling_freq < 0 or periodic_coupling_amp < 0:
            raise ValueError("periodic frequency/amplitude must be nonnegative")
        p_one = .5 + bias_delta
        if not 0 < p_one < 1 or not -1 <= autocorrelation_rho <= 1:
            raise ValueError("invalid bias or rho")
        p11 = p_one + autocorrelation_rho * (1 - p_one)
        p01 = p_one * (1 - autocorrelation_rho)
        if autocorrelation_rho and not (0 <= p11 <= 1 and 0 <= p01 <= 1):
            raise ValueError("rho and marginal bit probability are incompatible")
        if stuck_bit_index is not None and stuck_bit_index not in range(8):
            raise ValueError("stuck bit index must lie in [0,7]")
        if stuck_bit_value not in (0, 1):
            raise ValueError("stuck bit value must be binary")
        manifest = SourceManifest(
            source_id=f"fault_{base_source.manifest.source_id}",
            family="synthetic_fault",
            is_deterministic=base_source.manifest.is_deterministic,
            is_physical=False,
            metadata={
                "base_source": base_source.manifest.source_id,
                "fault_version": 2,
                "placement": "emitted_bytes_all_consumption_APIs",
                "rho_mechanism": "stationary_binary_markov" if autocorrelation_rho else "none",
                "rho_scope": "lag1 of contiguous emitted bits; short read_bits calls discard padding",
                "bias_mechanism": "bit_marginal" if autocorrelation_rho else "byte_shift_clipping",
                "bit_marginal_when_markov": p_one,
                "bias_delta": bias_delta,
                "autocorrelation_rho": autocorrelation_rho,
                "stuck_bit_index": stuck_bit_index,
                "periodic_coupling_freq": periodic_coupling_freq,
                "periodic_coupling_amp": periodic_coupling_amp,
            },
        )
        super().__init__(manifest, max_bytes=max_bytes)
        self.base = base_source
        self.bias = bias_delta
        self.rho = autocorrelation_rho
        self.stuck_bit = stuck_bit_index
        self.stuck_val = stuck_bit_value
        self.freq = periodic_coupling_freq
        self.amp = periodic_coupling_amp
        self.step_idx = 0
        self.last_bit = None
        self.p_one = p_one
        self.p11, self.p01 = p11, p01

    def _read_bytes(self, n: int) -> bytes:
        if self.rho:
            # One persistent Markov bit sequence, packed into bytes; floats and
            # direct bytes now see the same defect as read_bits. Channel effects
            # after packing need separate study; rho is a bit-lag parameter.
            raw = bytearray()
            for _ in range(n):
                value = 0
                for _ in range(8):
                    probability = self.p_one if self.last_bit is None else (self.p11 if self.last_bit else self.p01)
                    self.last_bit = int(self.base.random_float() < probability)
                    value = (value << 1) | self.last_bit
                raw.append(value)
        else:
            raw = bytearray(self.base.read_bytes(n))

        # 1. Stuck bit fault injection
        if self.stuck_bit is not None and 0 <= self.stuck_bit < 8:
            mask = 1 << self.stuck_bit
            for i in range(n):
                if self.stuck_val:
                    raw[i] |= mask
                else:
                    raw[i] &= ~mask

        # 2. Bias and Periodic Coupling
        if (not self.rho and abs(self.bias) > 1e-6) or self.amp > 1e-6:
            for i in range(n):
                harmonic = self.amp * math.sin(2.0 * math.pi * self.freq * self.step_idx) if self.freq > 0 else 0.0
                effective_shift = int(round(((self.bias if not self.rho else 0.0) + harmonic) * 255))
                val = raw[i] + effective_shift
                raw[i] = max(0, min(255, val))
                self.step_idx += 1

        return bytes(raw)

