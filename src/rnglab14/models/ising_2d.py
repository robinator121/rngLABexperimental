"""RNGLab v0.14: 2D Ising Model Simulation with Decoupled Stochastic Channels.

Isolates:
- Site proposal channel: random lattice coordinate selection (i, j)
- Acceptance channel: uniform random float u in [0, 1) tested against Boltzmann factor e^(-beta * dE)
- Initialization channel: random initial spin orientation (+1 or -1)
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Optional, Union

import numpy as np

from rnglab14.sources.base import RandomSource, StochasticChannels


@dataclass(frozen=True)
class IsingObservableResult:
    L: int
    temperature: float
    total_sweeps: int
    burn_in_sweeps: int
    energy_mean: float
    energy_variance: float
    specific_heat: float
    magnetization_mean: float
    magnetization_variance: float
    susceptibility: float
    acceptance_rate: float
    effective_sample_size_energy: float
    healthy: bool


class Ising2D:
    """2D square lattice Ising model with periodic boundary conditions."""

    def __init__(self, L: int = 8, temperature: float = 2.269, seed: Optional[int] = None):
        self.L = L
        self.N = L * L
        self.temperature = temperature
        self.beta = 1.0 / temperature
        self._seed = seed
        self.spins = np.ones((L, L), dtype=np.int8)

    def total_energy(self) -> float:
        s = self.spins
        right = np.roll(s, -1, axis=1)
        down = np.roll(s, -1, axis=0)
        return float(-np.sum(s * (right + down)))

    def total_magnetization(self) -> float:
        return float(np.sum(self.spins))

    def run_simulation(
        self,
        source_or_channels: Union[RandomSource, StochasticChannels],
        sweeps: int = 1000,
        burn_in: int = 200,
    ) -> IsingObservableResult:
        if isinstance(source_or_channels, StochasticChannels):
            prop_src = source_or_channels.proposal
            acc_src = source_or_channels.acceptance
            init_src = source_or_channels.initialization
        else:
            prop_src = source_or_channels
            acc_src = source_or_channels
            init_src = source_or_channels

        # Initialize random spin configuration using init_src
        for i in range(self.L):
            for j in range(self.L):
                self.spins[i, j] = 1 if init_src.random_float() < 0.5 else -1

        L = self.L
        N = self.N
        beta = self.beta
        boltzmann = {dE: math.exp(-beta * dE) for dE in (-8, -4, 0, 4, 8)}

        energies = []
        magnetizations = []
        accepted_flips = 0
        total_proposals = 0

        current_E = self.total_energy()
        current_M = self.total_magnetization()

        for sweep in range(sweeps):
            for _ in range(N):
                # 1. Site proposal
                i = prop_src.read_int(0, L)
                j = prop_src.read_int(0, L)
                s_ij = self.spins[i, j]

                # Nearest neighbor sum with periodic boundary conditions
                nb_sum = (
                    self.spins[(i + 1) % L, j]
                    + self.spins[(i - 1) % L, j]
                    + self.spins[i, (j + 1) % L]
                    + self.spins[i, (j - 1) % L]
                )
                dE = int(2 * s_ij * nb_sum)
                total_proposals += 1

                # 2. Acceptance decision
                if dE <= 0:
                    self.spins[i, j] = -s_ij
                    current_E += dE
                    current_M -= 2 * s_ij
                    accepted_flips += 1
                else:
                    u = acc_src.random_float()
                    if u < boltzmann.get(dE, math.exp(-beta * dE)):
                        self.spins[i, j] = -s_ij
                        current_E += dE
                        current_M -= 2 * s_ij
                        accepted_flips += 1

            if sweep >= burn_in:
                energies.append(current_E)
                magnetizations.append(abs(current_M) / N)

        e_arr = np.array(energies, dtype=float)
        m_arr = np.array(magnetizations, dtype=float)

        e_mean = float(np.mean(e_arr)) / N
        e_var = float(np.var(e_arr)) / N
        m_mean = float(np.mean(m_arr))
        m_var = float(np.var(m_arr))

        # Fluctuation-dissipation observables
        c_v = float(np.var(e_arr)) / (N * (self.temperature ** 2))
        chi = (N / self.temperature) * m_var
        acc_rate = accepted_flips / max(1, total_proposals)

        # Autocorrelation time estimation for energy
        n_meas = len(e_arr)
        if n_meas > 20:
            e_centered = e_arr - np.mean(e_arr)
            c0 = float(np.dot(e_centered, e_centered) / n_meas)
            tau_int = 0.5
            if c0 > 1e-12:
                for lag in range(1, min(50, n_meas // 4)):
                    c_lag = float(np.dot(e_centered[:-lag], e_centered[lag:]) / (n_meas - lag))
                    rho_lag = c_lag / c0
                    if rho_lag <= 0.05:
                        break
                    tau_int += rho_lag
            ess_e = max(1.0, float(n_meas / (2.0 * tau_int)))
        else:
            ess_e = float(n_meas)

        # Sanity bounds near critical point Tc ~ 2.269
        healthy = (abs(e_mean) < 2.0 and 0.0 <= m_mean <= 1.0 and ess_e > 5.0)

        return IsingObservableResult(
            L=L,
            temperature=self.temperature,
            total_sweeps=sweeps,
            burn_in_sweeps=burn_in,
            energy_mean=round(e_mean, 6),
            energy_variance=round(e_var, 6),
            specific_heat=round(c_v, 6),
            magnetization_mean=round(m_mean, 6),
            magnetization_variance=round(m_var, 6),
            susceptibility=round(chi, 6),
            acceptance_rate=round(acc_rate, 4),
            effective_sample_size_energy=round(ess_e, 2),
            healthy=healthy,
        )
