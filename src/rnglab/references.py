"""RNGLab v0.18 namespace: exact references.

Re-exports, unchanged:

- exact_finite_bead_energy: primitive-action harmonic PIMC energy at finite M,
  m = hbar = omega = 1 (rnglab14.models.pimc_harmonic).
- ising_exact_enumeration: exact finite-lattice Ising observables by full
  enumeration, periodic L in {2, 3, 4} only (rnglab16.core.enumeration).
  Susceptibility there is the signed-M variance / (N T).
- Two-state chain stationary law and finite-N mean variance, and the
  augmented driver kernel (rnglab16.core).
- ISING_TC: the infinite-lattice critical temperature, a regime marker only.
"""
from rnglab14.models.pimc_harmonic import exact_finite_bead_energy
from rnglab16.core import TC as ISING_TC
from rnglab16.core import augmented as augmented_driver_kernel
from rnglab16.core import enumeration as ising_exact_enumeration
from rnglab16.core import stationary, two_state
from rnglab16.core import variance_mean as two_state_mean_variance

__all__ = [
    "exact_finite_bead_energy",
    "ising_exact_enumeration",
    "ISING_TC",
    "stationary",
    "two_state",
    "two_state_mean_variance",
    "augmented_driver_kernel",
]
