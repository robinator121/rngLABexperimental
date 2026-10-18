"""
Tests for Financial Fault Injection Engine.
"""

import unittest
import numpy as np

from rnglab_finance.models.european_option import EuropeanOptionPricer
from rnglab_finance.models.portfolio_var import PortfolioVaRSimulator
from rnglab_finance.faults.injection import (
    ReusedSeedProvider,
    DuplicatedPathProvider,
    BiasedInnovationProvider,
    CorrelatedStreamProvider,
    PeriodicCorruptionProvider,
    GradualDegradationProvider,
)


class TestFinancialFaults(unittest.TestCase):
    """Verify reproducible behavior and mathematical consequences of injected defects."""

    def test_reused_seed_path_duplication(self):
        """ReusedSeedProvider should produce identical sub-sequences across workers."""
        provider = ReusedSeedProvider(base_seed=42, num_workers=4)
        z = provider.draw_normals(4000)
        self.assertEqual(len(z), 4000)
        
        # Quarter 1 and Quarter 2 should be identical
        q1 = z[0:1000]
        q2 = z[1000:2000]
        np.testing.assert_array_equal(q1, q2)
        
        # Unique elements should be exactly 1000, not 4000!
        self.assertEqual(len(np.unique(z)), 1000)

    def test_biased_innovation_breaks_martingale(self):
        """Biased innovations (delta > 0) should break Martingale preservation and bias call price upward."""
        provider = BiasedInnovationProvider(delta=0.08, scale=1.0, base_seed=42)
        pricer = EuropeanOptionPricer(s0=100.0, k=100.0, r=0.05, sigma=0.20, t=1.0, option_type="call")
        
        res = pricer.price(provider, n_paths=30000, antithetic=False)
        # Martingale error should be positive and z-score should exceed threshold
        self.assertGreater(res["martingale_error"], 0.5)
        self.assertFalse(res["martingale_preserved"])
        self.assertGreater(res["z_score"], 3.0)  # Call price significantly overestimated

    def test_duplicated_path_provider(self):
        """DuplicatedPathProvider should repeat chunks of size chunk_size."""
        provider = DuplicatedPathProvider(base_seed=123, chunk_size=200)
        z = provider.draw_normals(600)
        self.assertEqual(len(z), 600)
        np.testing.assert_array_equal(z[0:200], z[200:400])
        np.testing.assert_array_equal(z[0:200], z[400:600])

    def test_correlated_stream_provider(self):
        """CorrelatedStreamProvider should induce empirical correlation between asset columns."""
        provider = CorrelatedStreamProvider(coupling_rho=0.50, base_seed=999)
        z_flat = provider.draw_correlated_asset_normals(n_scenarios=10000, n_assets=2)
        z_mat = z_flat.reshape(10000, 2)
        
        corr = np.corrcoef(z_mat[:, 0], z_mat[:, 1])[0, 1]
        # Should be close to coupling_rho = 0.50
        self.assertAlmostEqual(corr, 0.50, delta=0.05)

    def test_periodic_corruption_provider(self):
        """PeriodicCorruptionProvider should introduce spikes at regular intervals."""
        provider = PeriodicCorruptionProvider(period=50, glitch_magnitude=4.0, base_seed=100)
        z = provider.draw_normals(200)
        # At indices 0, 50, 100, 150, z should be noticeably shifted
        for idx in [0, 50, 100, 150]:
            self.assertGreater(z[idx], 1.0)

    def test_gradual_degradation_provider(self):
        """GradualDegradationProvider should exhibit an increasing mean trend."""
        provider = GradualDegradationProvider(max_drift=0.30, base_seed=200)
        z = provider.draw_normals(10000)
        first_quarter_mean = np.mean(z[:2500])
        last_quarter_mean = np.mean(z[7500:])
        self.assertGreater(last_quarter_mean, first_quarter_mean + 0.15)


if __name__ == "__main__":
    unittest.main()
