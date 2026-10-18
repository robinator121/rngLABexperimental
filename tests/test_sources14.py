"""Unit tests for RNGLab v0.14 Sources, Physical Adapters, and Conditioning Pipelines."""
import unittest

from rnglab14.sources import (
    ChaCha20Source,
    DegradedFaultSource,
    GeigerMullerPulseSource,
    HashBlockConditioner,
    MersenneTwisterSource,
    MicrocontrollerSerialSource,
    OSEntropySource,
    PCG64Source,
    PhotodiodeNoiseSource,
    ReseedingDRBGPipeline,
    SourceBudgetExceeded,
    StochasticChannels,
    VonNeumannConditioner,
    XORIndependenceCheckedConditioner,
    Xoshiro256PlusPlusSource,
)


class SourcesV14Tests(unittest.TestCase):
    def test_reference_prng_determinism(self):
        # PCG64
        p1 = PCG64Source(seed=123)
        p2 = PCG64Source(seed=123)
        self.assertEqual(p1.read_bytes(32), p2.read_bytes(32))

        # Xoshiro256++
        x1 = Xoshiro256PlusPlusSource(seed=456)
        x2 = Xoshiro256PlusPlusSource(seed=456)
        self.assertEqual(x1.read_bytes(32), x2.read_bytes(32))

        # ChaCha20
        c1 = ChaCha20Source(key=b"K" * 32, rounds=20)
        c2 = ChaCha20Source(key=b"K" * 32, rounds=20)
        self.assertEqual(c1.read_bytes(64), c2.read_bytes(64))

        # ChaCha8 (round-reduced)
        c_r8 = ChaCha20Source(key=b"K" * 32, rounds=8)
        self.assertNotEqual(c1.read_bytes(64), c_r8.read_bytes(64))

    def test_source_budget_exceeded(self):
        src = MersenneTwisterSource(seed=1, max_bytes=16)
        raw = src.read_bytes(10)
        self.assertEqual(len(raw), 10)
        with self.assertRaises(SourceBudgetExceeded):
            src.read_bytes(10)  # Total 20 > 16

    def test_geiger_muller_dead_time_modes(self):
        # Non-paralyzable
        gm_np = GeigerMullerPulseSource(mean_cps=1000.0, dead_time_seconds=100e-6, paralyzable=False, seed=42)
        bytes_np = gm_np.read_bytes(64)
        self.assertEqual(len(bytes_np), 64)
        self.assertGreater(len(set(bytes_np)), 10)

        # Paralyzable
        gm_p = GeigerMullerPulseSource(mean_cps=1000.0, dead_time_seconds=100e-6, paralyzable=True, seed=42)
        bytes_p = gm_p.read_bytes(64)
        self.assertEqual(len(bytes_p), 64)

    def test_photodiode_adc_sampling(self):
        pd = PhotodiodeNoiseSource(adc_bits=12, noise_std_dev=20.0, seed=77)
        data = pd.read_bytes(128)
        self.assertEqual(len(data), 128)
        self.assertGreater(len(set(data)), 15)

    def test_mcu_serial_framing_and_loss_recovery(self):
        mcu = MicrocontrollerSerialSource(packet_payload_size=32, simulated_loss_rate=0.2, seed=888)
        data = mcu.read_bytes(100)
        self.assertEqual(len(data), 100)
        self.assertGreater(mcu.valid_packets, 0)

    def test_von_neumann_debiasing_efficiency(self):
        # Feed biased bits: P(1) = 0.8 -> extraction efficiency is p*(1-p) = 0.16
        raw = DegradedFaultSource(MersenneTwisterSource(seed=55), bias_delta=0.30)
        vn = VonNeumannConditioner(raw)
        data = vn.read_bytes(32)
        self.assertEqual(len(data), 32)
        self.assertGreater(vn.extraction_efficiency, 0.05)
        self.assertLess(vn.extraction_efficiency, 0.30)

    def test_xor_independence_check_rejects_correlated_streams(self):
        # Create two identical (perfectly correlated) streams
        s1 = MersenneTwisterSource(seed=1234)
        s2 = MersenneTwisterSource(seed=1234)
        with self.assertRaises(ValueError) as ctx:
            XORIndependenceCheckedConditioner([s1, s2], check_independence=True)
        self.assertIn("diagnostic rejected", str(ctx.exception))

        # Independent streams should pass
        s3 = MersenneTwisterSource(seed=101)
        s4 = MersenneTwisterSource(seed=202)
        xor_valid = XORIndependenceCheckedConditioner([s3, s4], check_independence=True)
        out = xor_valid.read_bytes(32)
        self.assertEqual(len(out), 32)

    def test_reseeding_drbg_pipeline_and_failure(self):
        raw = MersenneTwisterSource(seed=999)
        drbg = ReseedingDRBGPipeline(raw, reseed_interval_samples=4, simulate_reseed_failure=False)
        out1 = drbg.read_bytes(256)
        self.assertEqual(len(out1), 256)
        self.assertGreater(drbg.reseed_count, 1)

    def test_stochastic_channels_container(self):
        prop = PCG64Source(seed=1)
        acc = PCG64Source(seed=2)
        init = PCG64Source(seed=3)
        channels = StochasticChannels(proposal=prop, acceptance=acc, initialization=init)
        m = channels.to_manifest_dict()
        self.assertIn("proposal", m)
        self.assertIn("acceptance", m)
        self.assertIn("initialization", m)


if __name__ == "__main__":
    unittest.main()
