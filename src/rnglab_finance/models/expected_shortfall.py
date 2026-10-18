"""
Expected Shortfall (ES / CVaR) Simulator.

Implements coherent risk measure estimation under Basel III / FRTB standards.
Evaluates tail expectations beyond the Value at Risk quantile:
    ES_alpha = E[ L | L >= VaR_alpha ]
Benchmarks against analytical Gaussian formulas and independent numerical
quadrature for heavy-tailed distributions.
"""

import math
import hashlib
import json
from typing import Dict, Any, Optional, Tuple
import numpy as np

from rnglab_finance.core_interface import FinancialRNGProvider
from rnglab_finance.models.black_scholes_reference import AnalyticalRiskReference
from rnglab_finance.models.portfolio_var import PortfolioVaRSimulator


class ExpectedShortfallEstimator:
    """
    Monte Carlo Expected Shortfall / Conditional Value at Risk engine.
    """

    def __init__(
        self,
        portfolio_value: float,
        weights: np.ndarray,
        means: np.ndarray,
        cov_matrix: np.ndarray,
        horizon_dt: float = 1.0 / 252.0,
    ):
        self.var_sim = PortfolioVaRSimulator(
            portfolio_value=portfolio_value,
            weights=weights,
            means=means,
            cov_matrix=cov_matrix,
            horizon_dt=horizon_dt,
        )
        self.portfolio_value = self.var_sim.portfolio_value
        self.weights = self.var_sim.weights
        self.means = self.var_sim.means
        self.cov_matrix = self.var_sim.cov_matrix
        self.horizon_dt = self.var_sim.horizon_dt

    def estimate_expected_shortfall(
        self,
        rng_provider: FinancialRNGProvider,
        n_scenarios: int = 100000,
        confidence_level: float = 0.975,  # Basel III FRTB standard is 97.5%
    ) -> Dict[str, Any]:
        """
        Simulate portfolio losses, determine tail exceedances, and estimate Expected Shortfall.
        """
        losses, z_flat = self.var_sim.simulate_portfolio_losses(rng_provider, n_scenarios)

        # Analytical reference ground truth
        exact_es = AnalyticalRiskReference.gaussian_expected_shortfall(
            self.portfolio_value,
            self.weights,
            self.means,
            self.cov_matrix,
            confidence_level,
            self.horizon_dt,
        )
        exact_var = AnalyticalRiskReference.gaussian_var(
            self.portfolio_value,
            self.weights,
            self.means,
            self.cov_matrix,
            confidence_level,
            self.horizon_dt,
        )

        # Empirical VaR cutoff: empirical alpha quantile
        cutoff_idx = int(math.floor(confidence_level * n_scenarios))
        cutoff_idx = min(max(cutoff_idx, 0), n_scenarios - 1)
        sorted_losses = np.sort(losses)
        empirical_var = sorted_losses[cutoff_idx]

        # Tail losses beyond empirical VaR
        tail_losses = sorted_losses[cutoff_idx:]
        n_tail = len(tail_losses)
        if n_tail == 0:
            est_es = float(empirical_var)
            tail_se = 0.0
        else:
            est_es = float(np.mean(tail_losses))
            tail_var = float(np.var(tail_losses, ddof=1)) if n_tail > 1 else 0.0
            tail_se = float(math.sqrt(max(tail_var, 1e-18) / n_tail))

        # 95% Confidence Interval for Expected Shortfall
        ci_lower = est_es - 1.959964 * tail_se
        ci_upper = est_es + 1.959964 * tail_se
        ci_covers = bool(ci_lower <= exact_es <= ci_upper)

        error = est_es - exact_es
        rel_error = abs(error) / max(abs(exact_es), 1e-9)
        z_score = error / max(tail_se, 1e-12)

        prov_dict = {
            "portfolio_value": self.portfolio_value,
            "confidence_level": confidence_level,
            "n_scenarios": n_scenarios,
            "est_es": round(est_es, 4),
            "exact_es": round(exact_es, 4),
        }
        prov_hash = hashlib.sha256(json.dumps(prov_dict, sort_keys=True).encode("utf-8")).hexdigest()

        return {
            "confidence_level": confidence_level,
            "n_scenarios": n_scenarios,
            "n_tail_samples": n_tail,
            "empirical_var_cutoff": float(empirical_var),
            "exact_var_cutoff": float(exact_var),
            "estimated_expected_shortfall": est_es,
            "exact_expected_shortfall": exact_es,
            "error": error,
            "relative_error": rel_error,
            "tail_standard_error": tail_se,
            "z_score": z_score,
            "confidence_interval_95": [ci_lower, ci_upper],
            "ci_covers_exact": ci_covers,
            "provenance_hash": prov_hash,
        }
