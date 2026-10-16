"""MATH-01: Entropy Equivalence and Application Reliability.

Investigates stationary symmetric binary Markov sources with opposite autocorrelation (+rho vs -rho).

Key Mathematical Theorems & Identities:
1. Entropy Rate Symmetry:
   H_rate(rho) = h_2((1 + rho) / 2) = H_rate(-rho)
   where h_2(p) = -p*log2(p) - (1-p)*log2(1-p) is binary entropy.

2. Finite-Sample Variance Identity:
   For X_i in {0, 1} with E[X_i] = 1/2 and Corr(X_i, X_{i+k}) = rho^k:
   Var(X_mean) = (1 / (4n)) * [ 1 + 2 * sum_{k=1}^{n-1} (1 - k/n) * rho^k ]
               = (1 / (4n)) * [ (1 + rho)/(1 - rho) - (2*rho*(1 - rho^n)) / (n*(1 - rho)^2) ]

3. Asymptotic Variance & Effective Sample Size (ESS):
   lim_{n -> inf} n * Var(X_mean) = (1/4) * (1 + rho)/(1 - rho)
   tau_int = (1 + rho) / (1 - rho)
   ESS(rho) = n / tau_int = n * (1 - rho) / (1 + rho)

Crucial Scientific Finding:
   While +rho and -rho have IDENTICAL entropy rates, their Monte Carlo variances
   differ by a factor of ((1+rho)/(1-rho))^2. Thus, entropy rate is fundamentally
   insufficient to characterize or predict application-level stochastic reliability.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import random
from typing import Any, Sequence


def binary_entropy(p: float) -> float:
    """Binary Shannon entropy function h_2(p) in bits."""
    if p <= 0.0 or p >= 1.0:
        return 0.0
    return -p * math.log2(p) - (1.0 - p) * math.log2(1.0 - p)


def markov_entropy_rate(rho: float) -> float:
    """Entropy rate H_rate(rho) = h_2((1 + rho)/2) in bits per sample."""
    if not -1.0 < rho < 1.0:
        raise ValueError("Autocorrelation rho must be in (-1, 1)")
    p_trans = (1.0 + rho) / 2.0
    return binary_entropy(p_trans)


def exact_variance_summation(n: int, rho: float) -> float:
    """Compute finite-sample variance via exact summation:

    Var(X_mean) = (1 / (4n)) * [ 1 + 2 * sum_{k=1}^{n-1} (1 - k/n) * rho^k ]
    """
    if n <= 0:
        raise ValueError("Sample count n must be positive")
    if abs(rho) >= 1.0:
        raise ValueError("rho must be in (-1, 1)")
    if n == 1:
        return 0.25

    sum_term = sum((1.0 - k / n) * (rho ** k) for k in range(1, n))
    return (0.25 / n) * (1.0 + 2.0 * sum_term)


def exact_variance_closed_form(n: int, rho: float) -> float:
    """Compute finite-sample variance via closed-form algebraic formula:

    Var(X_mean) = (1 / (4n)) * [ (1 + rho)/(1 - rho) - (2*rho*(1 - rho^n)) / (n*(1 - rho)^2) ]
    """
    if n <= 0:
        raise ValueError("Sample count n must be positive")
    if abs(rho) >= 1.0:
        raise ValueError("rho must be in (-1, 1)")
    if n == 1 or abs(rho) < 1e-15:
        return 0.25 / n

    term1 = (1.0 + rho) / (1.0 - rho)
    term2 = (2.0 * rho * (1.0 - rho ** n)) / (n * ((1.0 - rho) ** 2))
    return (0.25 / n) * (term1 - term2)


def effective_sample_size(n: int, rho: float) -> float:
    """Asymptotic effective sample size ESS = n / tau_int."""
    tau_int = (1.0 + rho) / (1.0 - rho)
    return n / tau_int


def generate_markov_bits(n: int, rho: float, seed: int) -> list[int]:
    """Generate n bits from stationary symmetric binary Markov chain.

    Stationary distribution: P(X=0) = P(X=1) = 0.5.
    Transition probabilities:
    P(0|0) = (1 + rho)/2, P(1|0) = (1 - rho)/2
    P(1|1) = (1 + rho)/2, P(0|1) = (1 - rho)/2
    """
    rng = random.Random(seed)
    bits = [0] * n
    # Initial state from stationary distribution
    state = 1 if rng.random() < 0.5 else 0
    bits[0] = state

    p_same = (1.0 + rho) / 2.0
    for i in range(1, n):
        if rng.random() < p_same:
            # Stay in same state
            pass
        else:
            state = 1 - state
        bits[i] = state
    return bits


@dataclass(frozen=True)
class MarkovEquivalenceExperiment:
    n: int
    rho: float
    entropy_rate: float
    theoretical_variance: float
    empirical_mean: float
    empirical_variance: float
    mean_bias: float
    relative_variance_error: float
    effective_sample_size: float
    coverage_95: float
    trials: int


def run_entropy_equivalence_experiment(
    n: int,
    rho: float,
    trials: int = 1000,
    base_seed: int = 42,
) -> MarkovEquivalenceExperiment:
    """Run empirical Monte Carlo trials to verify analytical variance and coverage."""
    h_rate = markov_entropy_rate(rho)
    var_theory = exact_variance_closed_form(n, rho)
    ess = effective_sample_size(n, rho)
    std_theory = math.sqrt(var_theory)

    sample_means = []
    covered_count = 0
    true_mean = 0.5

    for t in range(trials):
        seed = base_seed + t * 7919
        bits = generate_markov_bits(n, rho, seed)
        m = sum(bits) / n
        sample_means.append(m)

        # 95% Confidence interval check using theoretical standard deviation
        ci_lower = m - 1.95996 * std_theory
        ci_upper = m + 1.95996 * std_theory
        if ci_lower <= true_mean <= ci_upper:
            covered_count += 1

    emp_mean = sum(sample_means) / trials
    emp_var = sum((m - emp_mean) ** 2 for m in sample_means) / (trials - 1)
    bias = emp_mean - true_mean
    rel_err = abs(emp_var - var_theory) / var_theory
    coverage = covered_count / trials

    return MarkovEquivalenceExperiment(
        n=n,
        rho=rho,
        entropy_rate=round(h_rate, 6),
        theoretical_variance=round(var_theory, 8),
        empirical_mean=round(emp_mean, 6),
        empirical_variance=round(emp_var, 8),
        mean_bias=round(bias, 6),
        relative_variance_error=round(rel_err, 4),
        effective_sample_size=round(ess, 2),
        coverage_95=round(coverage, 4),
        trials=trials,
    )
