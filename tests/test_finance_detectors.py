"""
Tests for Financial Integrity Monitoring and 4-Quadrant Classifier.
"""

import unittest
import numpy as np

from rnglab_finance.core_interface import FinancialRNGProvider
from rnglab14.sources.reference import PCG64Source
from rnglab_finance.models.european_option import EuropeanOptionPricer
from rnglab_finance.faults.injection import (
    ReusedSeedProvider,
    BiasedInnovationProvider,
)
from rnglab_finance.detectors.financial_monitor import FinancialIntegrityMonitor


class TestFinancialDetectors(unittest.TestCase):
    """Verify detection accuracy, false positive immunity, and 4-quadrant classification."""

    def setUp(self):
        self.monitor = FinancialIntegrityMonitor(
            alpha_significance=0.01,
            martingale_z_threshold=3.0,
            max_acceptable_pricing_error_pct=0.02,
        )
        self.pricer = EuropeanOptionPricer(
            s0=100.0, k=100.0, r=0.05, sigma=0.20, t=1.0, option_type="call"
        )

    def test_healthy_simulation_produces_q1_nominal(self):
        """Under a healthy reference source (PCG64), monitor should report Q1_Nominal with zero alarms."""
        source = PCG64Source(seed=20261018)
        provider = FinancialRNGProvider(source)
        st, z = self.pricer.simulate_terminal_prices(provider, n_paths=20000, antithetic=True)
        pricing_res = self.pricer.price(provider, n_paths=20000, antithetic=True)
        
        assessment = self.monitor.assess_financial_run(provider, pricing_res, innovations=z)
        self.assertEqual(assessment["operational_quadrant"], "Q1_Nominal")
        self.assertEqual(assessment["risk_level"], "LOW")
        self.assertFalse(assessment["source_alarm"])
        self.assertFalse(assessment["application_alarm"])

    def test_reused_seed_detection_exposes_baseline_blindness(self):
        """Reused seeds should be detected by path duplication while conventional SE check is blind."""
        provider = ReusedSeedProvider(base_seed=42, num_workers=4)
        st, z = self.pricer.simulate_terminal_prices(provider, n_paths=20000, antithetic=False)
        pricing_res = self.pricer.price(provider, n_paths=20000, antithetic=False)
        
        assessment = self.monitor.assess_financial_run(provider, pricing_res, innovations=z)
        self.assertTrue(assessment["duplication_detected"])
        self.assertTrue(assessment["application_alarm"])
        # Conventional SE check passed because duplicate paths falsely deflate SE
        self.assertTrue(assessment["baseline_se_check_passed"])
        self.assertTrue(assessment["baseline_blind_to_defect"])
        self.assertTrue(assessment["rnglab_superior_to_baseline"])

    def test_biased_innovation_triggers_martingale_alarm(self):
        """Biased innovations should trigger Martingale violation alarm."""
        provider = BiasedInnovationProvider(delta=0.08, scale=1.0, base_seed=42)
        st, z = self.pricer.simulate_terminal_prices(provider, n_paths=20000, antithetic=False)
        pricing_res = self.pricer.price(provider, n_paths=20000, antithetic=False)
        
        assessment = self.monitor.assess_financial_run(provider, pricing_res, innovations=z)
        self.assertTrue(assessment["martingale_violated"])
        self.assertTrue(assessment["application_alarm"])
        self.assertIn("MartingaleViolation", assessment["application_alarms_detail"][0])


if __name__ == "__main__":
    unittest.main()
