"""RNGLab v0.18 namespace: Monte Carlo models.

Re-exports, unchanged:

- Harmonic PIMC (rnglab14.models.pimc_harmonic): the package sampler with the
  optional centroid heat bath, the corrected thermodynamic and centroid-virial
  estimators, and the exact finite-bead Gaussian path sampler.
- Harmonic VMC (rnglab14.models.vmc_harmonic).
- 2D Ising: the source-API sampler Ising2D (rnglab14.models.ising_2d) and the
  vectorized NumPy Metropolis/Wolff kernels used by v0.16/v0.17
  (rnglab16.core.ising, rnglab16.core.wolff). The vectorized kernels draw
  from NumPy SeedSequence children, not from RandomSource objects.
- Monte Carlo decision support (rnglab14.models.decision_support).

Coverage of PIMC intervals is validated only for the Phase II standalone
implementation (research/phase2-20261030); validation through this API is the
v0.18 Phase 2 milestone.
"""
from rnglab14.models.decision_support import DecisionEvaluationResult, DecisionOption, MonteCarloDecisionSupport
from rnglab14.models.ising_2d import Ising2D, IsingObservableResult
from rnglab14.models.pimc_harmonic import (
    PathIntegralHarmonicOscillator,
    PIMCResult,
    harmonic_centroid_heatbath,
    harmonic_energy_estimators,
    sample_exact_finite_bead_paths,
)
from rnglab14.models.vmc_harmonic import HarmonicOscillatorVMC, VMCResult
from rnglab16.core import ising as ising_vectorized
from rnglab16.core import wolff as wolff_vectorized

__all__ = [
    "PathIntegralHarmonicOscillator",
    "PIMCResult",
    "harmonic_centroid_heatbath",
    "harmonic_energy_estimators",
    "sample_exact_finite_bead_paths",
    "HarmonicOscillatorVMC",
    "VMCResult",
    "Ising2D",
    "IsingObservableResult",
    "ising_vectorized",
    "wolff_vectorized",
    "MonteCarloDecisionSupport",
    "DecisionOption",
    "DecisionEvaluationResult",
]
