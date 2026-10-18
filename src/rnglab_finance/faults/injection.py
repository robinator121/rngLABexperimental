"""
Financial Fault Injection Engine.

Provides a documented, reproducible catalogue of computational failure modes
encountered in production quantitative finance Monte Carlo pipelines:
1. Reused Seeds across distributed/threaded workers (sample size collapse).
2. Duplicated or overlapping simulation paths.
3. Biased normal innovations (Martingale violation).
4. Unintended cross-stream correlation (diversification falsification).
5. Periodic corruption / clock harmonics.
6. Gradual degradation / parameter drift.
"""

import math
from typing import Dict, Any, List, Optional, Tuple
import numpy as np

from rnglab_finance.core_interface import FinancialRNGProvider
from rnglab14.sources.reference import PCG64Source


class ReusedSeedProvider(FinancialRNGProvider):
    """
    Simulates distributed workers initialized with identical random seeds.
    Divides total paths into M workers, each producing identical random sequences.
    Consequence: Effective sample size collapses to N / M; standard error is
    falsely deflated by sqrt(M), destroying confidence interval coverage.
    """

    def __init__(self, base_seed: int = 42, num_workers: int = 4, name: str = "ReusedSeed"):
        super().__init__(source=PCG64Source(seed=base_seed), name=name)
        self.base_seed = base_seed
        self.num_workers = num_workers

    def draw_normals(self, n: int, method: str = "box_muller") -> np.ndarray:
        # Each worker draws n // num_workers paths
        paths_per_worker = max(n // self.num_workers, 1)
        
        # Single worker generates paths
        worker_source = PCG64Source(seed=self.base_seed)
        worker_provider = FinancialRNGProvider(worker_source)
        batch = worker_provider.draw_normals(paths_per_worker, method=method)
        
        # Duplicate across workers
        duplicated = np.tile(batch, self.num_workers)[:n]
        self._total_drawn += len(duplicated)
        return duplicated


class DuplicatedPathProvider(FinancialRNGProvider):
    """
    Simulates worker buffer caching where a chunk of K simulation paths
    is repeatedly replayed cyclically.
    """

    def __init__(self, base_seed: int = 2026, chunk_size: int = 500, name: str = "DuplicatedPath"):
        super().__init__(source=PCG64Source(seed=base_seed), name=name)
        self.chunk_size = chunk_size
        self._cached_chunk: Optional[np.ndarray] = None

    def draw_normals(self, n: int, method: str = "box_muller") -> np.ndarray:
        if self._cached_chunk is None:
            self._cached_chunk = super().draw_normals(self.chunk_size, method=method)
        
        repeats = (n + self.chunk_size - 1) // self.chunk_size
        stream = np.tile(self._cached_chunk, repeats)[:n]
        self._total_drawn += n
        return stream


class BiasedInnovationProvider(FinancialRNGProvider):
    """
    Injects a mean offset delta or variance scaling into standard normal innovations.
    Z_biased = (Z * scale) + delta.
    Consequence: Directly breaks risk-neutral Martingale invariance E[e^{-rT} S_T] = S0.
    """

    def __init__(self, delta: float = 0.05, scale: float = 1.0, base_seed: int = 2026, name: str = "BiasedInnovation"):
        super().__init__(source=PCG64Source(seed=base_seed), name=name)
        self.delta = delta
        self.scale = scale

    def draw_normals(self, n: int, method: str = "box_muller") -> np.ndarray:
        raw_z = super().draw_normals(n, method=method)
        return raw_z * self.scale + self.delta


class CorrelatedStreamProvider(FinancialRNGProvider):
    """
    Injects unintended cross-correlation between supposedly independent asset streams.
    Given K assets, couples each asset innovation with a common latent noise factor:
        Z_k = sqrt(1 - rho) * Z_k_indep + sqrt(rho) * Z_common
    Consequence: Subverts portfolio diversification, underestimating VaR/ES.
    """

    def __init__(self, coupling_rho: float = 0.30, base_seed: int = 2026, name: str = "CorrelatedStreams"):
        super().__init__(source=PCG64Source(seed=base_seed), name=name)
        self.coupling_rho = float(coupling_rho)

    def draw_correlated_asset_normals(self, n_scenarios: int, n_assets: int) -> np.ndarray:
        """
        Draw (n_scenarios, n_assets) normal innovations where columns share
        unintended cross-correlation coupling_rho.
        """
        # Common factor per scenario
        common_z = super().draw_normals(n_scenarios).reshape(n_scenarios, 1)
        # Independent factors
        indep_z = super().draw_normals(n_scenarios * n_assets).reshape(n_scenarios, n_assets)
        
        # Coupled innovations
        sqrt_1_rho = math.sqrt(max(1.0 - self.coupling_rho, 0.0))
        sqrt_rho = math.sqrt(max(self.coupling_rho, 0.0))
        coupled_z = sqrt_1_rho * indep_z + sqrt_rho * common_z
        return coupled_z.flatten()

    def draw_normals(self, n: int, method: str = "box_muller") -> np.ndarray:
        return super().draw_normals(n, method=method)


class PeriodicCorruptionProvider(FinancialRNGProvider):
    """
    Injects periodic glitches / clock harmonics into the RNG stream at interval P.
    """

    def __init__(self, period: int = 64, glitch_magnitude: float = 2.5, base_seed: int = 2026, name: str = "PeriodicCorruption"):
        super().__init__(source=PCG64Source(seed=base_seed), name=name)
        self.period = period
        self.glitch_magnitude = glitch_magnitude

    def draw_normals(self, n: int, method: str = "box_muller") -> np.ndarray:
        z = super().draw_normals(n, method=method)
        # Apply periodic spike
        indices = np.arange(0, n, self.period)
        z[indices] += self.glitch_magnitude
        return z


class GradualDegradationProvider(FinancialRNGProvider):
    """
    Injects a non-stationary linear drift into innovation variance or mean over batches:
        delta(t) = drift_rate * (t / total_samples)
    """

    def __init__(self, max_drift: float = 0.15, base_seed: int = 2026, name: str = "GradualDegradation"):
        super().__init__(source=PCG64Source(seed=base_seed), name=name)
        self.max_drift = max_drift

    def draw_normals(self, n: int, method: str = "box_muller") -> np.ndarray:
        z = super().draw_normals(n, method=method)
        t = np.linspace(0.0, self.max_drift, n)
        return z + t
