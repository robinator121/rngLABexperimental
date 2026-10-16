"""2D Ising Model Monte Carlo Simulation on square lattice L x L using Metropolis algorithm."""
from __future__ import annotations

from dataclasses import dataclass
import math
import random
from typing import Any
import numpy as np

from rnglab13.sources.base import RandomSource


@dataclass(frozen=True)
class IsingObservableResult:
    L: int
    temperature: float
    total_sweeps: int
    burn_in: int
    energy_mean: float
    energy_var: float
    magnetization_mean: float
    abs_magnetization_mean: float
    magnetization_var: float
    specific_heat: float
    susceptibility: float
    autocorrelation_time_energy: float
    effective_sample_size_energy: float
    healthy: bool


class Ising2D:
    """2D Ising Model on L x L lattice with periodic boundary conditions."""

    def __init__(self, L: int = 16, temperature: float = 2.269185, seed: int = 42):
        self.L = L
        self.N = L * L
        self.temperature = temperature
        self.beta = 1.0 / temperature
        self.seed = seed
        self.rng = random.Random(seed)
        # Random initial spin configuration in {-1, +1}
        self.spins = np.array(
            [[1 if self.rng.random() < 0.5 else -1 for _ in range(L)] for _ in range(L)],
            dtype=np.int8,
        )

    def total_energy(self) -> float:
        """Compute exact total Hamiltonian energy: H = - sum_<i,j> s_i s_j."""
        L = self.L
        s = self.spins
        # Periodic neighbor interactions
        right = np.roll(s, -1, axis=1)
        down = np.roll(s, -1, axis=0)
        return float(-np.sum(s * (right + down)))

    def total_magnetization(self) -> float:
        """Total magnetization M = sum s_i."""
        return float(np.sum(self.spins))

    def run_simulation(
        self,
        source: RandomSource,
        sweeps: int = 2000,
        burn_in: int = 500,
    ) -> IsingObservableResult:
        """Execute Metropolis Monte Carlo simulation consuming random numbers from source."""
        L = self.L
        N = self.N
        beta = self.beta

        # Precompute Boltzmann factors for delta E in {-8, -4, 0, 4, 8}
        boltzmann = {dE: math.exp(-beta * dE) for dE in (-8, -4, 0, 4, 8)}

        energies = []
        magnetizations = []
        abs_mags = []

        cur_energy = self.total_energy()
        cur_mag = self.total_magnetization()

        for sweep in range(sweeps):
            # One sweep = N spin flip attempts
            for _ in range(N):
                # Pick random site
                i = int(source.random_float() * L) % L
                j = int(source.random_float() * L) % L

                s = self.spins[i, j]
                # Neighbor sum with periodic boundary conditions
                nb_sum = (
                    self.spins[(i + 1) % L, j]
                    + self.spins[(i - 1) % L, j]
                    + self.spins[i, (j + 1) % L]
                    + self.spins[i, (j - 1) % L]
                )
                delta_E = 2 * s * nb_sum

                # Metropolis acceptance criterion
                if delta_E <= 0 or source.random_float() < boltzmann[delta_E]:
                    self.spins[i, j] = -s
                    cur_energy += delta_E
                    cur_mag += -2 * s

            if sweep >= burn_in:
                # Per-spin observables
                energies.append(cur_energy / N)
                magnetizations.append(cur_mag / N)
                abs_mags.append(abs(cur_mag) / N)

        # Statistical analysis
        e_arr = np.array(energies)
        m_arr = np.array(magnetizations)
        abs_m_arr = np.array(abs_mags)

        e_mean = float(np.mean(e_arr))
        e_var = float(np.var(e_arr, ddof=1))
        m_mean = float(np.mean(m_arr))
        abs_m_mean = float(np.mean(abs_m_arr))
        m_var = float(np.var(m_arr, ddof=1))

        # Specific heat C_v = beta^2 * N * Var(E)
        specific_heat = (beta ** 2) * N * e_var
        # Susceptibility chi = beta * N * Var(|M|)
        susceptibility = beta * N * float(np.var(abs_m_arr, ddof=1))

        # Integrated autocorrelation time of energy
        tau_int = self._estimate_autocorrelation_time(e_arr)
        ess = len(e_arr) / max(1.0, 2.0 * tau_int)

        return IsingObservableResult(
            L=L,
            temperature=self.temperature,
            total_sweeps=sweeps,
            burn_in=burn_in,
            energy_mean=round(e_mean, 6),
            energy_var=round(e_var, 8),
            magnetization_mean=round(m_mean, 6),
            abs_magnetization_mean=round(abs_m_mean, 6),
            magnetization_var=round(m_var, 8),
            specific_heat=round(specific_heat, 6),
            susceptibility=round(susceptibility, 6),
            autocorrelation_time_energy=round(tau_int, 3),
            effective_sample_size_energy=round(ess, 1),
            healthy=(source.health.value == "healthy"),
        )

    @staticmethod
    def _estimate_autocorrelation_time(series: np.ndarray, max_lag: int = 50) -> float:
        """Estimate integrated autocorrelation time tau_int using windowed ACF."""
        n = len(series)
        if n < 20:
            return 1.0
        s = series - np.mean(series)
        var0 = np.dot(s, s) / n
        if var0 < 1e-15:
            return 1.0

        max_lag = min(max_lag, n // 4)
        acf = []
        for lag in range(1, max_lag):
            c = np.dot(s[:-lag], s[lag:]) / (n - lag)
            rho = c / var0
            if rho <= 0:
                break
            acf.append(rho)

        return 0.5 + sum(acf)
