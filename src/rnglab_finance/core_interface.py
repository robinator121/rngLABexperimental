"""
Core Interface Adapter for RNGLab-Finance.

This module provides a strictly isolated boundary between RNGLab-Finance
and the underlying scientific core (RNGLab v0.14). It enables future
extraction of RNGLab-Finance into a proprietary commercial codebase without
entangling application-layer financial models with core scientific primitives.
"""

from typing import Dict, Any, List, Optional, Tuple, Protocol
import numpy as np

# Core source abstractions
from rnglab14.sources.base import RandomSource, StochasticChannels
from rnglab14.sources.reference import (
    PCG64Source,
    ChaCha20Source,
    Xoshiro256PlusPlusSource,
    MersenneTwisterSource,
    OSEntropySource,
)
from rnglab14.sources.conditioning import (
    HashBlockConditioner,
    VonNeumannConditioner,
    ReseedingDRBGPipeline,
)
from rnglab14.sources.faults import DegradedFaultSource
from rnglab14.monitoring.dual_monitor import DualIntegrityMonitor, MonitoringQuadrant
from rnglab14.standards.sp800_22 import frequency_monobit_test, runs_test, SP800_22_Result
from rnglab14.standards.health_90b import RepetitionCountTest, AdaptiveProportionTest
from rnglab14.statistics.multitest import (
    bonferroni_correction,
    holm_bonferroni_correction,
    benjamini_hochberg_fdr,
)


class FinancialRNGProvider:
    """
    Standardized provider wrapping RNGLab sources for financial Monte Carlo engines.
    Supports single-stream, multi-stream channel isolation, and antithetic pair generation.
    """

    def __init__(self, source: RandomSource, name: str = "FinancialRNG"):
        self.source = source
        self.name = name
        self._total_drawn = 0

    def draw_uniforms(self, n: int) -> np.ndarray:
        """Draw n i.i.d. uniforms in [0, 1) using IEEE-754 53-bit resolution."""
        if n <= 0:
            return np.empty(0, dtype=np.float64)
        self._total_drawn += n
        raw = self.source.read_bytes(8 * n)
        u64 = np.frombuffer(raw, dtype=">u8")
        return (u64 >> 11) * (1.0 / 9007199254740992.0)

    def draw_normals(self, n: int, method: str = "box_muller") -> np.ndarray:
        """
        Draw n standard normal variates Z ~ N(0, 1) using a verifiable transform.
        Defaults to Box-Muller transformation of uniform pairs.
        """
        if method == "box_muller":
            # Box-Muller requires even number of uniforms
            m = (n + 1) // 2
            u1 = np.maximum(self.draw_uniforms(m), 1e-15)  # avoid log(0)
            u2 = self.draw_uniforms(m)
            
            radius = np.sqrt(-2.0 * np.log(u1))
            theta = 2.0 * np.pi * u2
            
            z0 = radius * np.cos(theta)
            z1 = radius * np.sin(theta)
            z = np.empty(2 * m, dtype=np.float64)
            z[0::2] = z0
            z[1::2] = z1
            return z[:n]
        elif method == "scipy_ppf":
            # Inverse CDF method: Z = Phi^{-1}(U)
            from scipy.stats import norm
            u = np.clip(self.draw_uniforms(n), 1e-15, 1.0 - 1e-15)
            return norm.ppf(u)
        else:
            raise ValueError(f"Unsupported normal transform method: {method}")

    def draw_antithetic_normals(self, n: int) -> Tuple[np.ndarray, np.ndarray]:
        """
        Draw n pairs of antithetic standard normal variates (Z, -Z).
        Guarantees exact zero sample mean for paired innovations.
        """
        z = self.draw_normals(n)
        return z, -z

    @property
    def total_drawn(self) -> int:
        return self._total_drawn

    def get_source_digest(self) -> str:
        return self.source.digest()
