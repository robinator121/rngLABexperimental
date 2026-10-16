"""Variational Monte Carlo (VMC) for Quantum Harmonic Oscillator.

Hamiltonian: H = -1/2 d^2/dx^2 + 1/2 x^2 (hbar = m = omega = 1)
Exact ground state: E_0 = 0.5, psi_0(x) = exp(-0.5 * x^2)

Trial wavefunction: psi_alpha(x) = exp(-alpha * x^2)
Local energy: E_L(x) = alpha + x^2 * (1/2 - 2 * alpha^2)

CRITICAL PHENOMENON (Zero-Variance Principle vs Observable Masking):
When alpha = 0.5, E_L(x) = 0.5 identically for ALL x!
Therefore, the energy estimate <E> remains exactly 0.5 EVEN IF the random
sampling is completely corrupted, while physical moments like <x^2> and <x^4>
become severely biased!
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from rnglab13.sources.base import RandomSource


@dataclass(frozen=True)
class VMCResult:
    alpha: float
    total_steps: int
    burn_in: int
    acceptance_rate: float
    energy_mean: float
    energy_std: float
    x2_mean: float
    x4_mean: float
    exact_energy: float
    exact_x2: float
    energy_bias: float
    x2_bias: float
    observable_masking_detected: bool


class HarmonicOscillatorVMC:
    """VMC simulation of 1D Harmonic Oscillator using Metropolis sampling."""

    def __init__(self, alpha: float = 0.5, delta: float = 1.0):
        self.alpha = alpha
        self.delta = delta  # Metropolis proposal step size [-delta, delta]
        self.exact_energy = alpha + (0.5 - 2.0 * alpha ** 2) / (4.0 * alpha)
        self.exact_x2 = 1.0 / (4.0 * alpha)

    def run_simulation(
        self,
        source: RandomSource,
        steps: int = 10000,
        burn_in: int = 2000,
    ) -> VMCResult:
        """Run VMC Metropolis sampling driven by source."""
        x = 0.0
        accepted = 0

        energies = []
        x2_vals = []
        x4_vals = []

        alpha = self.alpha
        delta = self.delta

        for step in range(steps):
            # Propose move: uniform in [-delta, delta] using two floats from source
            u1 = source.random_float()
            dx = (u1 * 2.0 - 1.0) * delta
            x_prop = x + dx

            # Metropolis acceptance probability: |psi(x_prop)|^2 / |psi(x)|^2
            # log ratio = -2 * alpha * (x_prop^2 - x^2)
            log_prob = -2.0 * alpha * (x_prop ** 2 - x ** 2)

            u2 = source.random_float()
            if log_prob >= 0.0 or u2 < math.exp(log_prob):
                x = x_prop
                accepted += 1

            if step >= burn_in:
                # Local energy: E_L(x) = alpha + x^2 * (1/2 - 2 * alpha^2)
                e_loc = alpha + (x ** 2) * (0.5 - 2.0 * (alpha ** 2))
                energies.append(e_loc)
                x2_vals.append(x ** 2)
                x4_vals.append(x ** 4)

        e_arr = np.array(energies)
        x2_arr = np.array(x2_vals)
        x4_arr = np.array(x4_vals)

        e_mean = float(np.mean(e_arr))
        e_std = float(np.std(e_arr, ddof=1)) if len(e_arr) > 1 else 0.0
        x2_mean = float(np.mean(x2_arr))
        x4_mean = float(np.mean(x4_arr))

        e_bias = e_mean - self.exact_energy
        x2_bias = x2_mean - self.exact_x2

        # Observable masking: energy bias is small (< 0.01) but x^2 bias is large (> 0.05)
        masking = (abs(e_bias) < 0.02 and abs(x2_bias) > 0.05)

        return VMCResult(
            alpha=alpha,
            total_steps=steps,
            burn_in=burn_in,
            acceptance_rate=round(accepted / steps, 4),
            energy_mean=round(e_mean, 6),
            energy_std=round(e_std, 6),
            x2_mean=round(x2_mean, 6),
            x4_mean=round(x4_mean, 6),
            exact_energy=round(self.exact_energy, 6),
            exact_x2=round(self.exact_x2, 6),
            energy_bias=round(e_bias, 6),
            x2_bias=round(x2_bias, 6),
            observable_masking_detected=masking,
        )
