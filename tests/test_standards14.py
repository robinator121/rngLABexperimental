"""Unit tests for RNGLab v0.14 Standards Traceability and Health Tests."""
import unittest

from rnglab14.standards.drbg_90a import verify_kat_hmac_drbg
from rnglab14.standards.entropy_90b import assess_min_entropy
from rnglab14.standards.health_90b import AdaptiveProportionTest, RepetitionCountTest
from rnglab14.standards.sp800_22 import run_sp800_22_battery
from rnglab14.standards.traceability import get_standards_catalog
from rnglab14.sources import MersenneTwisterSource


class StandardsV14Tests(unittest.TestCase):
    def test_standards_catalog_coverage(self):
        cat = get_standards_catalog()
        self.assertGreater(len(cat), 8)
        standards = set(r["standard"] for r in cat)
        self.assertIn("NIST SP 800-90A Rev. 1", standards)
        self.assertIn("NIST SP 800-90B", standards)
        self.assertIn("NIST SP 800-22 Rev. 1a", standards)
        self.assertIn("BSI AIS 20/31", standards)
        self.assertIn("FIPS 140-3", standards)

    def test_kat_hmac_drbg(self):
        self.assertTrue(verify_kat_hmac_drbg())

    def test_sp800_90b_rct_and_apt(self):
        rct = RepetitionCountTest(min_entropy_per_sample=1.0, alpha=2**-10)
        # Healthy alternating stream
        passed = True
        for i in range(100):
            res = rct.step(i % 2)
            if not res.passed:
                passed = False
        self.assertTrue(passed)

        # Stuck bit stream
        stuck_triggered = False
        for _ in range(20):
            res = rct.step(1)
            if not res.passed:
                stuck_triggered = True
                break
        self.assertTrue(stuck_triggered)

    def test_sp800_22_battery(self):
        src = MersenneTwisterSource(seed=777)
        raw = src.read_bytes(2048)
        battery = run_sp800_22_battery(raw)
        self.assertTrue(all(r.passed for r in battery))


if __name__ == "__main__":
    unittest.main()
