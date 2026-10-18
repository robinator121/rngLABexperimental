"""
Analytical Reference Benchmarks for Financial Models.

Provides exact closed-form solutions for Black-Scholes option pricing,
Greeks, Martingale discounted price invariants, Put-Call Parity,
and analytical Gaussian Value at Risk (VaR) and Expected Shortfall (ES).
"""

import math
from typing import Dict, Any, Tuple
import numpy as np
from scipy.stats import norm
import scipy.integrate as integrate


class BlackScholesReference:
    """
    Exact analytical Black-Scholes formula for European options.
    Serves as the infallible mathematical ground truth for estimator validation.
    """

    @staticmethod
    def d1_d2(s0: float, k: float, r: float, sigma: float, t: float) -> Tuple[float, float]:
        """Compute d1 and d2 parameters."""
        if s0 <= 0.0 or k <= 0.0 or sigma <= 0.0 or t <= 0.0:
            raise ValueError("All parameters s0, k, sigma, t must be strictly positive.")
        sqrt_t = math.sqrt(t)
        d1 = (math.log(s0 / k) + (r + 0.5 * sigma * sigma) * t) / (sigma * sqrt_t)
        d2 = d1 - sigma * sqrt_t
        return d1, d2

    @classmethod
    def call_price(cls, s0: float, k: float, r: float, sigma: float, t: float) -> float:
        """Exact European Call price: C = S0 * N(d1) - K * exp(-r*T) * N(d2)."""
        d1, d2 = cls.d1_d2(s0, k, r, sigma, t)
        discount = math.exp(-r * t)
        return float(s0 * norm.cdf(d1) - k * discount * norm.cdf(d2))

    @classmethod
    def put_price(cls, s0: float, k: float, r: float, sigma: float, t: float) -> float:
        """Exact European Put price: P = K * exp(-r*T) * N(-d2) - S0 * N(-d1)."""
        d1, d2 = cls.d1_d2(s0, k, r, sigma, t)
        discount = math.exp(-r * t)
        return float(k * discount * norm.cdf(-d2) - s0 * norm.cdf(-d1))

    @classmethod
    def greeks(cls, s0: float, k: float, r: float, sigma: float, t: float) -> Dict[str, float]:
        """Exact Black-Scholes Greeks."""
        d1, d2 = cls.d1_d2(s0, k, r, sigma, t)
        sqrt_t = math.sqrt(t)
        discount = math.exp(-r * t)
        pdf_d1 = norm.pdf(d1)

        delta_call = float(norm.cdf(d1))
        delta_put = float(delta_call - 1.0)
        gamma = float(pdf_d1 / (s0 * sigma * sqrt_t))
        vega = float(s0 * pdf_d1 * sqrt_t)
        theta_call = float(
            -(s0 * pdf_d1 * sigma) / (2.0 * sqrt_t) - r * k * discount * norm.cdf(d2)
        )
        rho_call = float(k * t * discount * norm.cdf(d2))

        return {
            "delta_call": delta_call,
            "delta_put": delta_put,
            "gamma": gamma,
            "vega": vega,
            "theta_call": theta_call,
            "rho_call": rho_call,
        }

    @classmethod
    def put_call_parity_rhs(cls, s0: float, k: float, r: float, t: float) -> float:
        """Put-Call Parity identity: C - P = S0 - K * exp(-r*T)."""
        return float(s0 - k * math.exp(-r * t))

    @classmethod
    def martingale_asset_expectation(cls, s0: float) -> float:
        """Under risk-neutral measure Q, E[exp(-r*T) * S_T] = S0."""
        return float(s0)


class AnalyticalRiskReference:
    """
    Exact analytical formulas for Gaussian Portfolio Value at Risk and Expected Shortfall.
    """

    @staticmethod
    def portfolio_moments(
        weights: np.ndarray,
        means: np.ndarray,
        cov_matrix: np.ndarray,
        horizon_dt: float = 1.0,
    ) -> Tuple[float, float]:
        """
        Compute mean and standard deviation of linear portfolio returns.
        mu_p = w^T mu * dt
        sigma_p = sqrt(w^T Sigma w * dt)
        """
        w = np.asarray(weights, dtype=np.float64)
        mu = np.asarray(means, dtype=np.float64)
        cov = np.asarray(cov_matrix, dtype=np.float64)

        port_mean = float(np.dot(w, mu) * horizon_dt)
        port_variance = float(np.dot(w, np.dot(cov, w)) * horizon_dt)
        port_std = math.sqrt(max(port_variance, 1e-18))
        return port_mean, port_std

    @classmethod
    def gaussian_var(
        cls,
        portfolio_value: float,
        weights: np.ndarray,
        means: np.ndarray,
        cov_matrix: np.ndarray,
        confidence_level: float = 0.95,
        horizon_dt: float = 1.0,
    ) -> float:
        """
        Exact Gaussian Value at Risk (Loss quantile).
        Loss L = -Portfolio_Return * Portfolio_Value
        VaR_alpha = - (mu_p - z_alpha * sigma_p) * V0
        where z_alpha = Phi^{-1}(alpha).
        """
        mu_p, sigma_p = cls.portfolio_moments(weights, means, cov_matrix, horizon_dt)
        z_alpha = float(norm.ppf(confidence_level))
        loss_quantile = -(mu_p - z_alpha * sigma_p) * portfolio_value
        return float(loss_quantile)

    @classmethod
    def gaussian_expected_shortfall(
        cls,
        portfolio_value: float,
        weights: np.ndarray,
        means: np.ndarray,
        cov_matrix: np.ndarray,
        confidence_level: float = 0.95,
        horizon_dt: float = 1.0,
    ) -> float:
        """
        Exact Gaussian Expected Shortfall (CVaR).
        ES_alpha = E[L | L >= VaR_alpha]
        For normal loss: ES_alpha = mu_L + sigma_L * (phi(z_alpha) / (1 - alpha))
        """
        mu_p, sigma_p = cls.portfolio_moments(weights, means, cov_matrix, horizon_dt)
        z_alpha = float(norm.ppf(confidence_level))
        pdf_z = float(norm.pdf(z_alpha))
        tail_factor = pdf_z / (1.0 - confidence_level)

        mu_loss = -mu_p * portfolio_value
        sigma_loss = sigma_p * portfolio_value
        es = mu_loss + sigma_loss * tail_factor
        return float(es)

    @classmethod
    def student_t_expected_shortfall_quad(
        cls,
        portfolio_value: float,
        df: float,
        scale: float,
        confidence_level: float = 0.95,
    ) -> float:
        """
        Expected Shortfall under Student-t distribution via numerical quadrature.
        Used as independent non-Gaussian ground truth.
        """
        from scipy.stats import t as student_t
        t_alpha = student_t.ppf(confidence_level, df=df, scale=scale)
        # integrate x * f(x) from t_alpha to infinity
        integral, _ = integrate.quad(
            lambda x: x * student_t.pdf(x, df=df, scale=scale),
            t_alpha,
            np.inf,
            limit=200,
        )
        tail_es = (integral / (1.0 - confidence_level)) * portfolio_value
        return float(tail_es)
