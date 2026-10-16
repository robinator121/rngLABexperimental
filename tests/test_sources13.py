"""Unit tests for RNGLab v0.13+ Randomness Source Infrastructure."""
import io
from pathlib import Path
import tempfile
import unittest

from rnglab13.sources import (
    BinaryFileSource,
    HashConditioner,
    MersenneTwisterSource,
    OSCryptoSource,
    PulseEventAdapter,
    RandomSource,
    ReseedingCSPRNGSource,
    SerialTRNGAdapter,
    SourceExhausted,
    SourceHealth,
    SourceKind,
    SyntheticTRNGSource,
    VonNeumannConditioner,
    XORCombiner,
    Xoshiro256PlusPlusSource,
)


class SourcesV13Tests(unittest.TestCase):
    """Test the complete randomness-source abstraction suite."""

    def test_prng_determinism_and_manifest(self):
        s1 = MersenneTwisterSource(seed=12345, source_id="mt-test-1")
        s2 = MersenneTwisterSource(seed=12345, source_id="mt-test-2")
        data1 = s1.read_bytes(100)
        data2 = s2.read_bytes(100)
        self.assertEqual(data1, data2)
        self.assertEqual(s1.total_bytes_read, 100)
        self.assertEqual(s2.total_bytes_read, 100)

        man1 = s1.reproducibility_manifest()
        man2 = s2.reproducibility_manifest()
        self.assertEqual(man1["stream_sha256"], man2["stream_sha256"])
        self.assertEqual(man1["total_bytes_read"], 100)

    def test_xoshiro256pp_determinism(self):
        s1 = Xoshiro256PlusPlusSource(seed=999)
        s2 = Xoshiro256PlusPlusSource(seed=999)
        s3 = Xoshiro256PlusPlusSource(seed=1000)
        d1 = s1.read_bytes(256)
        d2 = s2.read_bytes(256)
        d3 = s3.read_bytes(256)
        self.assertEqual(d1, d2)
        self.assertNotEqual(d1, d3)

    def test_typed_readers(self):
        s = MersenneTwisterSource(seed=42)
        u32 = s.read_u32()
        self.assertIsInstance(u32, int)
        self.assertTrue(0 <= u32 <= 0xFFFFFFFF)

        u64 = s.read_u64()
        self.assertIsInstance(u64, int)
        self.assertTrue(0 <= u64 <= 0xFFFFFFFFFFFFFFFF)

        flt = s.random_float()
        self.assertIsInstance(flt, float)
        self.assertTrue(0.0 <= flt < 1.0)

    def test_os_crypto_source(self):
        crypto = OSCryptoSource()
        b1 = crypto.read_bytes(64)
        b2 = crypto.read_bytes(64)
        self.assertEqual(len(b1), 64)
        self.assertEqual(len(b2), 64)
        self.assertNotEqual(b1, b2)
        self.assertEqual(crypto.kind, SourceKind.OS_CRYPTOGRAPHIC)

    def test_reseeding_csprng_automatic_reseed(self):
        entropy = MersenneTwisterSource(seed=777)
        csprng = ReseedingCSPRNGSource(entropy_source=entropy, reseed_interval_bytes=64)
        self.assertEqual(csprng.metadata.reseeding_events, 0)

        # Read less than interval
        _ = csprng.read_bytes(40)
        self.assertEqual(csprng.metadata.reseeding_events, 0)

        # Read past interval -> triggers reseed
        _ = csprng.read_bytes(50)
        self.assertGreaterEqual(csprng.metadata.reseeding_events, 1)

    def test_synthetic_trng_flags_not_physical(self):
        syn = SyntheticTRNGSource(seed=101)
        data = syn.read_bytes(128)
        self.assertEqual(len(data), 128)
        self.assertTrue(syn.metadata.is_synthetic)
        self.assertFalse(syn.metadata.is_physical)
        self.assertIn("SYNTHETIC SIMULATION ONLY", syn.metadata.extra["scientific_notice"])

    def test_binary_file_source_and_exhaustion(self):
        with tempfile.NamedTemporaryFile("wb", delete=False) as f:
            test_bytes = b"RNGLAB_TEST_ENTROPY_DATA_12345678"
            f.write(test_bytes)
            f_path = f.name

        try:
            bf = BinaryFileSource(f_path)
            read_1 = bf.read_bytes(10)
            self.assertEqual(read_1, test_bytes[:10])
            read_2 = bf.read_bytes(23)
            self.assertEqual(read_2, test_bytes[10:33])

            # Exhaustion
            with self.assertRaises(SourceExhausted):
                bf.read_bytes(1)
            self.assertEqual(bf.health, SourceHealth.EXHAUSTED)
        finally:
            Path(f_path).unlink(missing_ok=True)

    def test_serial_trng_adapter_framing_and_dropped_packets(self):
        buf = io.BytesIO()
        # Create frames 0, 1, 3 (simulating dropped frame 2)
        payload0 = b"abc"
        payload1 = b"defgh"
        payload3 = b"ijkl"
        buf.write(SerialTRNGAdapter.create_frame(0, payload0))
        buf.write(SerialTRNGAdapter.create_frame(1, payload1))
        # Frame 2 is skipped: jump from seq 1 to 3
        buf.write(SerialTRNGAdapter.create_frame(3, payload3))
        buf.seek(0)

        adapter = SerialTRNGAdapter(buf)
        data = adapter.read_bytes(len(payload0) + len(payload1) + len(payload3))
        self.assertEqual(data, payload0 + payload1 + payload3)
        self.assertEqual(adapter.dropped_packets, 1)

    def test_pulse_event_adapter(self):
        # Timestamps with alternating intervals: (dt1 > dt2 -> 1), (dt1 < dt2 -> 0)
        # t: 0, 10, 12 => dt1=10, dt2=2 => bit 1
        # t: 20, 22, 35 => dt1=2, dt2=13 => bit 0
        timestamps = []
        # Generate 8 bits = 1 byte (10101010 = 0xAA)
        base = 0.0
        for i in range(8):
            bit = 1 if (i % 2 == 0) else 0
            if bit == 1:
                dt1, dt2 = 1.0, 0.2
            else:
                dt1, dt2 = 0.2, 1.0
            t1 = base
            t2 = t1 + dt1
            t3 = t2 + dt2
            timestamps.extend([t1, t2, t3])
            base = t3 + 0.1

        pulse_src = PulseEventAdapter(timestamps)
        out = pulse_src.read_bytes(1)
        self.assertEqual(out[0], 0xAA)
        self.assertEqual(pulse_src.total_pulses, 24)

    def test_von_neumann_conditioner(self):
        raw = MersenneTwisterSource(seed=123)
        cond = VonNeumannConditioner(raw)
        out = cond.read_bytes(32)
        self.assertEqual(len(out), 32)
        self.assertEqual(cond.metadata.kind, SourceKind.CONDITIONED)

    def test_hash_conditioner_compression(self):
        raw = MersenneTwisterSource(seed=555)
        h_cond = HashConditioner(raw, input_block_bytes=64)
        out = h_cond.read_bytes(48)
        self.assertEqual(len(out), 48)
        # 48 bytes requires 2x 32-byte hash blocks -> 2x 64 = 128 raw bytes read
        self.assertEqual(raw.total_bytes_read, 128)

    def test_xor_combiner_mixing(self):
        s1 = MersenneTwisterSource(seed=1)
        s2 = MersenneTwisterSource(seed=2)
        combiner = XORCombiner([s1, s2])
        out = combiner.read_bytes(16)
        self.assertEqual(len(out), 16)
        self.assertEqual(combiner.metadata.kind, SourceKind.MIXED)


if __name__ == "__main__":
    unittest.main()
