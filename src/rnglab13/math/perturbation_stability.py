"""MATH-02: Perturbation Sensitivity and Stability.

Investigates the relationship between random-source perturbations (theta)
and application-level observable errors:
S(theta) = |E_theta[f] - E_0[f]| <= K * |theta|

Mathematical Rigor:
1. Total Variation Distance Bounds:
   For any bounded observable f with ||f||_inf <= M:
   |E_theta[f] - E_0[f]| <= 2 * M * d_TV(P_theta, P_0)

2. Relative Entropy (KL Divergence) & Pinsker's Inequality:
   d_TV(P_theta, P_0) <= sqrt( (1/2) * D_KL(P_0 || P_theta) )

3. Markov Chain Transition Kernel Perturbations:
   Let P_0 be ergodic with stationary distribution pi_0 and fundamental matrix:
   Z = (I - P_0 + 1 * pi_0)^(-1)
   Then for perturbed kernel P_theta = P_0 + theta * Delta_P:
   pi_theta - pi_0 = theta * pi_0 * Delta_P * Z + O(theta^2)
   The sensitivity constant K is bounded by ||Z||_1 * ||Delta_P||_1.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable, Sequence
import numpy as np


def total_variation_distance(p: np.ndarray, q: np.ndarray) -> float:
    """Compute total variation distance d_TV(P, Q) = (1/2) * sum |p_i - q_i|."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    return 0.5 * float(np.sum(np.abs(p - q)))


def kl_divergence(p: np.ndarray, q: np.ndarray) -> float:
    """Kullback-Leibler divergence D_KL(P || Q) in nats."""
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    mask = (p > 0)
    if np.any((q <= 0) & mask):
        return float("inf")
    return float(np.sum(p[mask] * np.log(p[mask] / q[mask])))


def markov_stationary_distribution(P: np.ndarray) -> np.ndarray:
    """Compute stationary distribution pi = pi * P for ergodic transition matrix P."""
    P = np.asarray(P, dtype=float)
    dim = P.shape[0]
    # Solve (P^T - I) pi^T = 0 with sum(pi) = 1
    A = np.vstack([P.T - np.eye(dim), np.ones(dim)])
    b = np.append(np.zeros(dim), [1.0])
    pi, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    return pi


def markov_fundamental_matrix(P: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    """Compute stationary distribution pi, fundamental matrix Z, and spectral gap gamma."""
    P = np.asarray(P, dtype=float)
    dim = P.shape[0]
    pi = markov_stationary_distribution(P)
    one_pi = np.outer(np.ones(dim), pi)
    # Z = (I - P + 1*pi)^(-1)
    Z = np.linalg.inv(np.eye(dim) - P + one_pi)

    # Eigenvalues for spectral gap
    eigenvals = np.sort(np.abs(np.linalg.eigvals(P)))[::-1]
    spectral_gap = 1.0 - float(eigenvals[1]) if len(eigenvals) > 1 else 1.0

    return pi, Z, spectral_gap


@dataclass(frozen=True)
class PerturbationAnalysisResult:
    theta: float
    observable_baseline: float
    observable_perturbed: float
    observed_sensitivity: float
    theoretical_tv_bound: float
    lipschitz_ratio: float
    linear_regime_valid: bool


def analyze_markov_perturbation(
    P0: np.ndarray,
    Delta_P: np.ndarray,
    observable: np.ndarray,
    theta_values: Sequence[float],
) -> list[PerturbationAnalysisResult]:
    """Analyze observable shift under kernel perturbation P(theta) = P0 + theta * Delta_P."""
    P0 = np.asarray(P0, dtype=float)
    Delta_P = np.asarray(Delta_P, dtype=float)
    observable = np.asarray(observable, dtype=float)

    pi0, Z, gamma = markov_fundamental_matrix(P0)
    e0 = float(np.dot(pi0, observable))
    m_norm = float(np.max(np.abs(observable)))

    results = []
    # Theoretical first-order sensitivity slope: dE/dtheta = pi0 * Delta_P * Z * observable
    first_order_slope = float(np.dot(np.dot(pi0, Delta_P), np.dot(Z, observable)))

    for theta in theta_values:
        P_theta = P0 + theta * Delta_P
        # Ensure row stochasticity
        P_theta = np.clip(P_theta, 0.0, 1.0)
        P_theta = P_theta / np.sum(P_theta, axis=1, keepdims=True)

        pi_theta = markov_stationary_distribution(P_theta)
        e_theta = float(np.dot(pi_theta, observable))
        s_theta = abs(e_theta - e0)

        # TV distance between stationary distributions
        d_tv = total_variation_distance(pi_theta, pi0)
        tv_bound = 2.0 * m_norm * d_tv

        ratio = s_theta / abs(theta) if abs(theta) > 1e-12 else abs(first_order_slope)
        is_linear = (s_theta <= tv_bound + 1e-12)

        results.append(PerturbationAnalysisResult(
            theta=round(theta, 6),
            observable_baseline=round(e0, 6),
            observable_perturbed=round(e_theta, 6),
            observed_sensitivity=round(s_theta, 6),
            theoretical_tv_bound=round(tv_bound, 6),
            lipschitz_ratio=round(ratio, 4),
            linear_regime_valid=is_linear,
        ))

    return results
