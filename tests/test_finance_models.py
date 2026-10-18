"""
Unit and Convergence Tests for RNGLab-Finance Simulation Models.
"""

import unittest
import math
import numpy as np

from rnglab_finance.core_interface import FinancialRNGProvider
from rnglab14.sources.reference import PCG64Source
from rnglab_finance.models.black_scholes_reference import (
    BlackScholesReference,
    AnalyticalRiskReference,
)
from rnglab_finance.models.european_option import EuropeanOptionPricer
from rnglab_finance.models.portfolio_var import PortfolioVaRSimulator
from rnglab_finance.models.expected_shortfall import ExpectedShortfallEstimator


class TestFinancialModels(unittest.TestCase):
    """Verify analytical reference formulas and Monte Carlo model convergence."""

    def setUp(self):
        self.s0 = 100.0
        self.k = 100.0
        self.r = 0.05
        self.sigma = 0.20
        self.t = 1.0

    def test_black_scholes_analytical_values(self):
        """Check known analytical values for Black-Scholes European options."""
        call = BlackScholesReference.call_price(self.s0, self.k, self.r, self.sigma, self.t)
        put = BlackScholesReference.put_price(self.s0, self.k, self.r, self.sigma, self.t)
        
        # Standard textbook value: S0=100, K=100, r=0.05, sigma=0.2, T=1 -> Call ~ 10.4506
        self.assertAlmostEqual(call, 10.4506, places=3)
        # Put-Call parity: C - P = S0 - K * exp(-r*T)
        rhs = BlackScholesReference.put_call_parity_rhs(self.s0, self.k, self.r, self.t)
        self.assertAlmostEqual(call - put, rhs, places=7)

    def test_greeks_consistency(self):
        """Check Delta and Gamma properties."""
        greeks = BlackScholesReference.greeks(self.s0, self.k, self.r, self.sigma, self.t)
        self.assertGreater(greeks["delta_call"], 0.5)
        self.assertLess(greeks["delta_call"], 0.7)
        self.assertAlmostEqual(greeks["delta_call"] - greeks["delta_put"], 1.0, places=7)
        self.assertGreater(greeks["gamma"], 0.0)
        self.assertGreater(greeks["vega"], 0.0)

    def test_european_call_monte_carlo_convergence(self):
        """Monte Carlo European call pricer should converge within 95% CI to analytical value."""
        source = PCG64Source(seed=20261018)
        provider = FinancialRNGProvider(source)
        pricer = EuropeanOptionPricer(self.s0, self.k, self.r, self.sigma, self.t, option_type="call")
        
        res = pricer.price(provider, n_paths=50000, antithetic=True)
        self.assertTrue(res["ci_95_covers"])
        self.assertLess(res["relative_error"], 0.01)  # within 1% error
        self.assertTrue(res["martingale_preserved"])

    def test_european_put_monte_carlo_convergence(self):
        """Monte Carlo European put pricer should converge within 95% CI."""
        source = PCG64Source(seed=20261018)
        provider = FinancialRNGProvider(source)
        pricer = EuropeanOptionPricer(self.s0, self.k, self.r, self.sigma, self.t, option_type="put")
        
        res = pricer.price(provider, n_paths=50000, antithetic=True)
        self.assertTrue(res["ci_95_covers"])
        self.assertLess(res["relative_error"], 0.015)

    def test_portfolio_var_simulation(self):
        """Portfolio VaR Monte Carlo should match analytical Gaussian VaR."""
        weights = np.array([0.4, 0.3, 0.3])
        means = np.array([0.08, 0.12, 0.05])
        cov = np.array([
            [0.04, 0.01, 0.005],
            [0.01, 0.09, 0.02],
            [0.005, 0.02, 0.0225],
        ])
        sim = PortfolioVaRSimulator(portfolio_value=1_000_000.0, weights=weights, means=means, cov_matrix=cov)
        
        source = PCG64Source(seed=20261018)
        provider = FinancialRNGProvider(source)
        res = sim.calculate_var(provider, n_scenarios=50000, confidence_level=0.95)
        
        self.assertLess(res["relative_error"], 0.03)  # within 3% error on quantile
        self.assertTrue(res["ci_covers_exact"])
        self.assertTrue(res["kupiec_test_passed"])

    def test_expected_shortfall_simulation(self):
        """Expected shortfall estimator should match analytical Gaussian ground truth."""
        weights = np.array([0.5, 0.5])
        means = np.array([0.10, 0.10])
        cov = np.array([
            [0.04, 0.01],
            [0.01, 0.04],
        ])
        es_est = ExpectedShortfallEstimator(portfolio_value=500_000.0, weights=weights, means=means, cov_matrix=cov)
        
        source = PCG64Source(seed=20261018)
        provider = FinancialRNGProvider(source)
        res = es_est.estimate_expected_shortfall(provider, n_scenarios=50000, confidence_level=0.975)
        
        self.assertLess(res["relative_error"], 0.04)  # within 4% error on 97.5% tail
        self.assertTrue(res["ci_covers_exact"])


if __name__ == "__main__":
    unittest.main()
