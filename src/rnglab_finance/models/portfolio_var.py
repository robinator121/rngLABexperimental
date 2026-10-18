"""
Multi-Asset Portfolio Value at Risk (VaR) Simulator.

Implements Monte Carlo portfolio loss simulation with correlated asset returns
via Cholesky decomposition. Benchmarks against analytical Gaussian VaR ground truth,
evaluates empirical quantile errors, order-statistic confidence intervals,
and executes the Kupiec Proportion of Failures (POF) test.
"""

import math
import hashlib
import json
from typing import Dict, Any, List, Optional, Tuple
import numpy as np
from scipy.stats import norm, chi2

from rnglab_finance.core_interface import FinancialRNGProvider
from rnglab_finance.models.black_scholes_reference import AnalyticalRiskReference


class PortfolioVaRSimulator:
    """
    Monte Carlo Value at Risk (VaR) engine for linear multi-asset portfolios.
    """

    def __init__(
        self,
        portfolio_value: float,
        weights: np.ndarray,
        means: np.ndarray,
        cov_matrix: np.ndarray,
        horizon_dt: float = 1.0 / 252.0,  # 1-day horizon by default
    ):
        if portfolio_value <= 0.0:
            raise ValueError("portfolio_value must be > 0")

        self.portfolio_value = float(portfolio_value)
        self.weights = np.asarray(weights, dtype=np.float64)
        self.means = np.asarray(means, dtype=np.float64)
        self.cov_matrix = np.asarray(cov_matrix, dtype=np.float64)
        self.horizon_dt = float(horizon_dt)

        self.num_assets = len(self.weights)
        if len(self.means) != self.num_assets or self.cov_matrix.shape != (self.num_assets, self.num_assets):
            raise ValueError("Dimensions of weights, means, and cov_matrix must match.")

        # Ensure positive-definiteness and compute Cholesky factor
        try:
            self.chol = np.linalg.cholesky(self.cov_matrix)
        except np.linalg.LinAlgError:
            # Add tiny regularization if on boundary
            jitter = 1e-12 * np.eye(self.num_assets)
            self.chol = np.linalg.cholesky(self.cov_matrix + jitter)

        # Analytical benchmarks
        self.exact_var_95 = AnalyticalRiskReference.gaussian_var(
            self.portfolio_value, self.weights, self.means, self.cov_matrix, 0.95, self.horizon_dt
        )
        self.exact_var_99 = AnalyticalRiskReference.gaussian_var(
            self.portfolio_value, self.weights, self.means, self.cov_matrix, 0.99, self.horizon_dt
        )

    def simulate_portfolio_losses(
        self,
        rng_provider: FinancialRNGProvider,
        n_scenarios: int,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Simulate n_scenarios of correlated multi-asset returns and return portfolio losses.
        Loss L = - Portfolio_Return * Portfolio_Value.
        """
        # Draw uncorrelated standard normals: (n_scenarios, num_assets)
        total_normals = n_scenarios * self.num_assets
        z_flat = rng_provider.draw_normals(total_normals)
        z = z_flat.reshape(n_scenarios, self.num_assets)

        # Correlated innovations: (n_scenarios, num_assets)
        # R_i = mu_i * dt + sqrt(dt) * (L * Z)_i
        sqrt_dt = math.sqrt(self.horizon_dt)
        drift = self.means * self.horizon_dt
        diffusion = np.dot(z, self.chol.T) * sqrt_dt
        asset_returns = drift + diffusion

        # Portfolio return: R_p = sum(w_i * R_i)
        port_returns = np.dot(asset_returns, self.weights)
        # Loss: L = -V0 * R_p
        losses = -self.portfolio_value * port_returns
        return losses, z_flat

    def calculate_var(
        self,
        rng_provider: FinancialRNGProvider,
        n_scenarios: int = 100000,
        confidence_level: float = 0.95,
    ) -> Dict[str, Any]:
        """
        Run Monte Carlo simulation, compute empirical VaR, and perform backtesting checks.
        """
        losses, z_flat = self.simulate_portfolio_losses(rng_provider, n_scenarios)
        losses_sorted = np.sort(losses)

        # Exact analytical reference
        exact_var = AnalyticalRiskReference.gaussian_var(
            self.portfolio_value,
            self.weights,
            self.means,
            self.cov_matrix,
            confidence_level,
            self.horizon_dt,
        )

        # Empirical quantile: index = floor(alpha * N)
        index = int(math.floor(confidence_level * n_scenarios))
        index = min(max(index, 0), n_scenarios - 1)
        est_var = float(losses_sorted[index])

        # Order statistic 95% confidence interval for the alpha-quantile
        # Standard error of quantile: se_q = sqrt(alpha * (1 - alpha) / N) / f(q)
        q_se = math.sqrt(confidence_level * (1.0 - confidence_level) / n_scenarios)
        idx_lower = max(int(math.floor((confidence_level - 1.96 * q_se) * n_scenarios)), 0)
        idx_upper = min(int(math.ceil((confidence_level + 1.96 * q_se) * n_scenarios)), n_scenarios - 1)
        ci_lower = float(losses_sorted[idx_lower])
        ci_upper = float(losses_sorted[idx_upper])
        ci_covers = bool(ci_lower <= exact_var <= ci_upper)

        error = est_var - exact_var
        rel_error = abs(error) / max(abs(exact_var), 1e-9)

        # Kupiec Proportion of Failures (POF) test on theoretical VaR
        # Exceedance: loss > exact_var
        exceedances = int(np.sum(losses > exact_var))
        expected_exceedances = (1.0 - confidence_level) * n_scenarios
        actual_rate = exceedances / float(n_scenarios)
        nominal_rate = 1.0 - confidence_level

        # Likelihood ratio test statistic
        if exceedances == 0:
            lr_pof = -2.0 * n_scenarios * math.log(max(confidence_level, 1e-15))
        elif exceedances == n_scenarios:
            lr_pof = -2.0 * n_scenarios * math.log(max(nominal_rate, 1e-15))
        else:
            term1 = (n_scenarios - exceedances) * math.log(confidence_level / (1.0 - actual_rate))
            term2 = exceedances * math.log(nominal_rate / actual_rate)
            lr_pof = -2.0 * (term1 + term2)

        p_value_pof = float(1.0 - chi2.cdf(max(lr_pof, 0.0), df=1))
        kupiec_passed = bool(p_value_pof >= 0.05)

        # Innovation health checks
        z_mean = float(np.mean(z_flat))
        z_var = float(np.var(z_flat, ddof=1))

        prov_dict = {
            "portfolio_value": self.portfolio_value,
            "n_assets": self.num_assets,
            "n_scenarios": n_scenarios,
            "confidence_level": confidence_level,
            "est_var": round(est_var, 4),
            "exact_var": round(exact_var, 4),
        }
        prov_hash = hashlib.sha256(json.dumps(prov_dict, sort_keys=True).encode("utf-8")).hexdigest()

        return {
            "confidence_level": confidence_level,
            "n_scenarios": n_scenarios,
            "estimated_var": est_var,
            "exact_var": exact_var,
            "error": error,
            "relative_error": rel_error,
            "quantile_ci_95": [ci_lower, ci_upper],
            "ci_covers_exact": ci_covers,
            "exceedance_count": exceedances,
            "expected_exceedances": expected_exceedances,
            "actual_exceedance_rate": actual_rate,
            "nominal_exceedance_rate": nominal_rate,
            "kupiec_lr_statistic": float(lr_pof),
            "kupiec_p_value": p_value_pof,
            "kupiec_test_passed": kupiec_passed,
            "innovation_mean": z_mean,
            "innovation_variance": z_var,
            "provenance_hash": prov_hash,
        }
