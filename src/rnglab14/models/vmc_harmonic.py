"""RNGLab v0.14: Variational Monte Carlo Harmonic Oscillator with Observable Masking.

Channel Decoupling:
- Proposal Channel: Coordinate step delta x in [-delta, delta]
- Acceptance Channel: Metropolis acceptance ratio comparison u in [0, 1)

Observable Masking Analysis (MATH-02):
- At optimal alpha = 0.5: Local energy E_L(x) = 0.5 identically, with zero variance.
- When proposal stream is corrupted, spatial moments <x^2>, <x^4> deviate drastically,
  while energy expectation remains pristine (Observable Masking).
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Optional, Union

import numpy as np

from rnglab14.sources.base import RandomSource, StochasticChannels
from rnglab12.uncertainty import summarize_series
from rnglab14.monitoring.calibration import FrozenRankCalibration


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
    exploratory_masking_flag: bool
    masking_p_value: Optional[float]
    masking_calibration_hash: Optional[str]
    monitoring_status: str
    energy_uncertainty: dict
    x2_uncertainty: dict


class HarmonicOscillatorVMC:
    """Variational Monte Carlo for 1D Harmonic Oscillator: H = -1/2 d^2/dx^2 + 1/2 x^2.
    
    Trial wavefunction: Psi_alpha(x) = exp(-alpha * x^2).
    Exact ground state corresponds to alpha = 0.5, with E_0 = 0.5 and <x^2>_0 = 0.5.
    """

    def __init__(self, alpha: float = 0.5, delta: float = 1.0, seed: Optional[int] = None):
        if not math.isfinite(alpha) or alpha <= 0 or not math.isfinite(delta) or delta <= 0:
            raise ValueError("finite positive alpha and delta required")
        self.alpha = alpha
        self.delta = delta
        # Exact values under uncorrupted Psi_alpha
        self.exact_energy = alpha + (0.5 - 2.0 * alpha ** 2) / (4.0 * alpha)
        self.exact_x2 = 1.0 / (4.0 * alpha)

    def monitoring_protocol(self, steps: int, burn_in: int) -> dict:
        """Whole-chain statistic and sampler settings; independent seeds vary."""
        return {"model": "harmonic_vmc", "alpha": float(self.alpha), "delta": float(self.delta),
                "steps": steps, "burn_in": burn_in,
                "score": "abs(post_burn_in_mean_x_squared - 1/(4*alpha))",
                "unit": "one independently seeded complete chain",
                "null": "same healthy proposal/acceptance law and initialization protocol"}

    def run_simulation(
        self,
        source_or_channels: Union[RandomSource, StochasticChannels],
        steps: int = 10000,
        burn_in: int = 2000,
        masking_calibration: Optional[FrozenRankCalibration] = None,
    ) -> VMCResult:
        if (type(steps) is not int or type(burn_in) is not int or burn_in < 0
                or steps - burn_in < 32):
            raise ValueError("integer steps/burn_in and at least 32 retained steps required")
        protocol = self.monitoring_protocol(steps, burn_in)
        if masking_calibration is not None:
            # Fail before spending compute if the frozen protocol is incompatible.
            masking_calibration.p_value(0.0, protocol)
        if isinstance(source_or_channels, StochasticChannels):
            prop_src = source_or_channels.proposal
            acc_src = source_or_channels.acceptance
        else:
            prop_src = source_or_channels
            acc_src = source_or_channels

        x = 0.0
        accepted = 0
        energies = []
        x2_vals = []
        x4_vals = []

        alpha = self.alpha
        delta = self.delta

        for step in range(steps):
            # 1. Propose coordinate displacement
            u_prop = prop_src.random_float()
            dx = (u_prop * 2.0 - 1.0) * delta
            x_prop = x + dx

            # 2. Acceptance decision
            log_prob = -2.0 * alpha * (x_prop ** 2 - x ** 2)
            u_acc = acc_src.random_float()

            if log_prob >= 0.0 or u_acc < math.exp(log_prob):
                x = x_prop
                accepted += 1

            if step >= burn_in:
                # Local energy: E_L(x) = alpha + x^2 * (1/2 - 2 * alpha^2)
                # At alpha = 0.5: E_L(x) = 0.5 + x^2 * 0 = 0.5 identically!
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

        # Retain the historical threshold as an explicitly exploratory diagnostic.
        exploratory_masking = abs(e_bias) < 0.02 and abs(x2_bias) > 0.05
        p_value = None if masking_calibration is None else masking_calibration.p_value(abs(x2_bias), protocol)
        # The exact energy blindness identity applies only at the exact eigenstate.
        masking = (self.alpha == 0.5 and masking_calibration is not None
                   and p_value <= masking_calibration.alpha)
        status = ("UNCALIBRATED_EXPLORATORY_ONLY" if masking_calibration is None else
                  "CALIBRATED_EXCHANGEABLE_CHAIN_NULL" if self.alpha == 0.5 else
                  "CALIBRATED_POSITION_SCORE_NOT_EXACT_ENERGY_MASKING_TARGET")

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
            exploratory_masking_flag=exploratory_masking,
            masking_p_value=p_value,
            masking_calibration_hash=None if masking_calibration is None else masking_calibration.digest,
            monitoring_status=status,
            energy_uncertainty=summarize_series(e_arr, reference=self.exact_energy),
            x2_uncertainty=summarize_series(x2_arr, reference=self.exact_x2),
        )
