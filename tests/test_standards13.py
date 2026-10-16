"""Unit tests for RNGLab v0.13+ Standards Integration (SP 800-90A/B/C, SP 800-22)."""
import random
import unittest

from rnglab13.sources import MersenneTwisterSource
from rnglab13.standards import (
    AdaptiveProportionTest,
    HMAC_DRBG,
    MinEntropyReport,
    RepetitionCountTest,
    StandardsRequirement,
    assess_min_entropy,
    estimate_mcv,
    get_traceability_matrix,
    run_sp800_22_battery,
    verify_kat_hmac_drbg,
)


class StandardsV13Tests(unittest.TestCase):
    """Test standards compliance tooling, health tests, KAT vectors, and min-entropy."""

    def test_repetition_count_test_healthy_and_stuck(self):
        # Min-entropy = 1.0 bit/sample, alpha = 2^-20 => cutoff = 1 + ceil(20/1) = 21
        rct = RepetitionCountTest(min_entropy_per_sample=1.0, alpha=2**-20)
        self.assertEqual(rct.cutoff, 21)

        # Alternating stream: passes
        for i in range(100):
            res = rct.step(i % 2)
            self.assertTrue(res.passed)

        # Stuck stream: 21 identical samples triggers failure
        rct_stuck = RepetitionCountTest(min_entropy_per_sample=1.0, alpha=2**-20)
        last_res = None
        for _ in range(20):
            res = rct_stuck.step(1)
            self.assertTrue(res.passed)
        # 21st identical sample: fails
        res_fail = rct_stuck.step(1)
        self.assertFalse(res_fail.passed)
        self.assertGreaterEqual(res_fail.current_value, 21)

    def test_adaptive_proportion_test_healthy_and_biased(self):
        apt = AdaptiveProportionTest(window_size=512, min_entropy_per_sample=1.0, alpha=2**-20)
        # On balanced random stream: should pass
        rng = random.Random(42)
        for _ in range(512):
            res = apt.step(rng.randint(0, 1))
            self.assertTrue(res.passed)

        # On heavily biased stream (95% ones): triggers failure within window
        apt_biased = AdaptiveProportionTest(window_size=512, min_entropy_per_sample=1.0, alpha=2**-20)
        failed = False
        for _ in range(512):
            sample = 1 if rng.random() < 0.95 else 0
            res = apt_biased.step(sample)
            if not res.passed:
                failed = True
                break
        self.assertTrue(failed, "Adaptive Proportion Test should fail on 95% biased stream")

    def test_min_entropy_estimation_uniform_vs_biased(self):
        rng = random.Random(1234)
        # 1. Near-uniform bits (high entropy)
        uniform_samples = [rng.randint(0, 1) for _ in range(4000)]
        rep_unif: MinEntropyReport = assess_min_entropy(uniform_samples)
        self.assertGreater(rep_unif.min_entropy, 0.50)
        self.assertGreater(rep_unif.mcv_estimate, 0.85)

        # 2. Biased binary stream (p = 0.9 -> theoretical min-entropy = -log2(0.9) ~ 0.152)
        biased_samples = [1 if rng.random() < 0.90 else 0 for _ in range(4000)]
        rep_biased: MinEntropyReport = assess_min_entropy(biased_samples)
        self.assertLess(rep_biased.min_entropy, 0.4)
        self.assertLess(rep_biased.mcv_estimate, 0.3)

    def test_kat_hmac_drbg_cavp_vector(self):
        """Verify HMAC_DRBG matches official NIST CAVP Known Answer Test vector."""
        kat_passed = verify_kat_hmac_drbg()
        self.assertTrue(kat_passed, "HMAC_DRBG must match NIST CAVP test vector exactly")

    def test_sp800_22_battery_passes_good_source(self):
        prng = MersenneTwisterSource(seed=777)
        data = prng.read_bytes(2048)  # 16,384 bits
        results = run_sp800_22_battery(data)
        for res in results:
            self.assertTrue(res.passed, f"{res.test_name} should pass with p={res.p_value}")
            self.assertGreaterEqual(res.p_value, 0.01)

    def test_sp800_22_battery_fails_stuck_bits(self):
        stuck_data = bytes([0xFF] * 1024)
        results = run_sp800_22_battery(stuck_data)
        # Monobit must fail
        monobit = next(r for r in results if r.test_name == "Frequency Monobit")
        self.assertFalse(monobit.passed)
        self.assertLess(monobit.p_value, 0.01)

    def test_traceability_matrix_catalog(self):
        matrix = get_traceability_matrix()
        self.assertGreaterEqual(len(matrix), 5)
        req_ids = {r["req_id"] for r in matrix}
        self.assertIn("REQ-90B-RCT", req_ids)
        self.assertIn("REQ-90B-APT", req_ids)
        self.assertIn("REQ-90B-MIN-ENTROPY", req_ids)
        self.assertIn("REQ-90A-HMAC-DRBG", req_ids)
        self.assertIn("REQ-SP800-22", req_ids)


if __name__ == "__main__":
    unittest.main()
