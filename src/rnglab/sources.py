"""RNGLab v0.18 namespace: random sources and stochastic channels.

Re-exports the audited source hierarchy from rnglab14 unchanged; see
docs/v018/ARCHITECTURE.md. Every read passes through RandomSource.read_bytes,
which counts and SHA-256 hashes the emitted bytes. A digest identifies the
consumed stream; it is not an entropy or independence certificate.
StochasticChannels names the proposal/acceptance/initialization/resampling
roles; distinct objects do not prove independent origin.

Not re-exported: the physical emulators (GeigerMullerPulseSource,
PhotodiodeNoiseSource, MicrocontrollerSerialSource) and the conditioners, which
are synthetic research treatments (docs/archive/V0141_CORRECTIONS.md). Import
them from rnglab14.sources if needed.
"""
from rnglab14.sources.base import RandomSource, SourceBudgetExceeded, SourceManifest, StochasticChannels
from rnglab14.sources.reference import (
    ChaCha20Source,
    MersenneTwisterSource,
    OSEntropySource,
    PCG64Source,
    Xoshiro256PlusPlusSource,
)

__all__ = [
    "RandomSource",
    "SourceManifest",
    "SourceBudgetExceeded",
    "StochasticChannels",
    "PCG64Source",
    "Xoshiro256PlusPlusSource",
    "MersenneTwisterSource",
    "ChaCha20Source",
    "OSEntropySource",
]
