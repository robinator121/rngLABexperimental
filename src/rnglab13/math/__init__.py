"""RNGLab v0.13+ Mathematical Research Modules (MATH-01 to MATH-06)."""
from .entropy_equivalence import (
    MarkovEquivalenceExperiment,
    binary_entropy,
    effective_sample_size,
    exact_variance_closed_form,
    exact_variance_summation,
    generate_markov_bits,
    markov_entropy_rate,
    run_entropy_equivalence_experiment,
)
from .impossibility import (
    IndistinguishablePairExperiment,
    compute_lipschitz_drift_bound,
    run_indistinguishable_histories_counterexample,
)
from .information_functional import (
    FunctionalStabilityReport,
    fundamental_equation_of_information_defect,
    perturbed_information_function,
    shannon_information_function,
    test_hyers_ulam_stability_of_fei,
    xor_composition_cauchy_identity,
)
from .perturbation_stability import (
    PerturbationAnalysisResult,
    analyze_markov_perturbation,
    kl_divergence,
    markov_fundamental_matrix,
    markov_stationary_distribution,
    total_variation_distance,
)
from .recovery_memory import (
    RecoveryTrajectoryResult,
    predict_recovery_horizon,
    simulate_recovery_trajectory,
)

__all__ = [
    "binary_entropy",
    "markov_entropy_rate",
    "exact_variance_summation",
    "exact_variance_closed_form",
    "effective_sample_size",
    "generate_markov_bits",
    "run_entropy_equivalence_experiment",
    "MarkovEquivalenceExperiment",
    "total_variation_distance",
    "kl_divergence",
    "markov_stationary_distribution",
    "markov_fundamental_matrix",
    "analyze_markov_perturbation",
    "PerturbationAnalysisResult",
    "run_indistinguishable_histories_counterexample",
    "IndistinguishablePairExperiment",
    "compute_lipschitz_drift_bound",
    "fundamental_equation_of_information_defect",
    "shannon_information_function",
    "perturbed_information_function",
    "test_hyers_ulam_stability_of_fei",
    "FunctionalStabilityReport",
    "xor_composition_cauchy_identity",
    "simulate_recovery_trajectory",
    "predict_recovery_horizon",
    "RecoveryTrajectoryResult",
]
