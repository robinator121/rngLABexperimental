"""RNGLab v0.20 nuclear module: synthetic stochastic nuclear-physics workloads.

Purpose: measure how defects in random inputs propagate into the outputs of
stochastic nuclear-physics simulations. This is research into computational
reliability. Every nuclide, rate and configuration is synthetic and
non-operational; nothing here is reactor-physics, neutron-transport or
reactor-design software, and passing tests do not make it validated physics
software.

Status (2026-11-16): Layer A (event-level radioactive decay) is implemented
and unit-tested; its verification against analytic laws is the frozen study
research/v020-nuclear-a-verify. Layer B (branching processes and reduced
kinetics), Layer C (detector response), fault-injection campaigns and HPC
execution are not implemented yet. Architecture and status:
docs/v020/nuclear/NUCLEAR_ARCHITECTURE.md.
"""
from .bateman import (
    bateman_chain,
    decay_probabilities,
    expected_populations,
    generator_matrix,
    path_state_probabilities,
    state_probabilities,
    survival_probability,
    transient_equilibrium_ratio,
    truncated_exponential_cdf,
)
from .decay import (
    LAYOUT,
    MODE_SPECS,
    PARTICLE_CODES,
    Channel,
    DecayRun,
    DecayScheme,
    Mode,
    Nuclide,
    simulate_decays,
    simulate_decays_reference,
)
from .streams import NumpyPCG64Source, WordStream, labelled_source, words_to_uniforms

__all__ = [
    "LAYOUT", "MODE_SPECS", "PARTICLE_CODES", "Channel", "DecayRun", "DecayScheme", "Mode", "Nuclide",
    "simulate_decays", "simulate_decays_reference",
    "bateman_chain", "decay_probabilities", "expected_populations", "generator_matrix",
    "path_state_probabilities", "state_probabilities", "survival_probability",
    "transient_equilibrium_ratio", "truncated_exponential_cdf",
    "NumpyPCG64Source", "WordStream", "labelled_source", "words_to_uniforms",
]
