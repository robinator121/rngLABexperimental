"""Unit tests for the Stage C (H1) library, research/v018-h1-detection/h1lib.py.

SP 800-22 expectations are the worked examples of NIST SP 800-22 Rev 1a
(sections 2.4.8, 2.11.4/2.11.6, 2.12.4/2.12.8, 2.13.4/2.13.8), read from the
standard on 2026-11-06. Two of its printed values are internally inconsistent and
are handled explicitly:
- serial: step 5 prints 0.9057 and 0.8805, but recomputing from the stated psi^2
  values gives the section 2.11.6 values 0.808792 and 0.670320, which are tested;
- longest run: section 2.4.4 gives 0.180598 and section 2.4.8 gives 0.180609
  (chi^2 4.882605 vs 4.882457); the tolerance covers both.
Runs use only the smoke namespace and short runs.
"""
import importlib.util
import math
import struct
import sys
from pathlib import Path

import numpy as np
import pytest

from rnglab11.statistics import HealthTests
from rnglab14.sources import PCG64Source

DIR = Path(__file__).resolve().parents[2] / "research" / "v018-h1-detection"
spec = importlib.util.spec_from_file_location("v018_h1lib", DIR / "h1lib.py")
h1 = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = h1
spec.loader.exec_module(h1)

SMOKE = "rnglab-v018-h1-detection-smoke"
E100 = "1100100100001111110110101010001000100001011010001100001000110100110001001100011001100010100010111000"
LR128 = ("11001100000101010110110001001100111000000000001001" "00110101010001000100111101011010000000110101111100"
         "1100111001101101100010110010")


def bits(s):
    return np.array([int(c) for c in s], dtype=np.int8)


def test_serial_worked_example():
    p1, p2 = h1.serial_test(bits("0011011101"), 3)
    assert (round(p1, 6), round(p2, 6)) == (0.808792, 0.670320)


def test_approximate_entropy_worked_examples():
    assert round(h1.approximate_entropy_test(bits("0100110101"), 3), 6) == 0.261961
    assert round(h1.approximate_entropy_test(bits(E100), 2), 6) == 0.235301


def test_cumulative_sums_forward_worked_examples():
    assert abs(h1.cumulative_sums_forward(bits("1011010111")) - 0.4116588) < 5e-7
    assert round(h1.cumulative_sums_forward(bits(E100)), 6) == 0.219194


def test_lower_order_pattern_counts_are_exact():
    b = np.unpackbits(np.frombuffer(bytes(PCG64Source(5).read_bytes(4096)), dtype=np.uint8))
    nine = h1.pattern_counts(b, 9)
    for m in range(1, 9):
        assert np.array_equal(nine[m], h1.pattern_counts(b, m)[m])
    assert nine[0].tolist() == [len(b)]


def test_longest_run_worked_example_and_tables():
    p, nu = h1.longest_run_test(bits(LR128))
    assert nu == [4, 9, 3, 0]
    assert abs(p - 0.1806) < 2e-5
    assert [h1.longest_run_block_size(n) for n in (128, 6272, 750000)] == [8, 128, 10000]
    for k, edges, pi in h1.LONGEST_RUN_TABLES.values():
        assert len(edges) == len(pi) == k + 1 and abs(sum(pi) - 1) < 2e-4


def test_word_tests_detect_a_single_duplicate_and_accept_random_words():
    data = bytes(PCG64Source(7).read_bytes(8 * 4096))
    clean = h1.word_tests(data)
    assert clean["duplicate_count"] == 0 and clean["adjacent_duplicates"] == 1.0 and clean["lag1_correlation"] > 1e-4
    dup = data[:800] + data[792:800] + data[800 + 8:]
    hit = h1.word_tests(dup)
    assert hit["duplicate_count"] == 1
    assert math.isclose(hit["adjacent_duplicates"], 4095 * 2.0**-64, rel_tol=1e-6)


def test_conventional_window_returns_every_p_value():
    window = h1.conventional_window(bytes(PCG64Source(11).read_bytes(108000)))
    assert set(window) == set(h1.CONV_P_KEYS)
    assert all(v is None or 0 <= v <= 1 for v in window.values())


def test_vectorized_health_tests_match_the_repository_implementation():
    rng = np.random.default_rng(3)
    x = rng.integers(0, 256, size=20000)
    x[1000:1006] = 42                       # repetition-count event
    x[4096:4096 + 512:3] = x[4096]          # adaptive-proportion event in one window
    tests = HealthTests(min_entropy=7, alpha=2**-20, alphabet_size=256)
    events = HealthTests(min_entropy=7, alpha=2**-20, alphabet_size=256).update(x)
    rct, apt = h1.health_test_events(x, tests)
    assert rct.tolist() == [e["sample"] for e in events if e["test"] == "repetition_count"]
    assert apt.tolist() == [e["sample"] for e in events if e["test"] == "adaptive_proportion"]
    assert len(rct) >= 1 and len(apt) >= 1


def test_trinomial_irreducibility():
    assert h1.trinomial_irreducible(4, 1) and h1.trinomial_irreducible(7, 3)
    assert not h1.trinomial_irreducible(4, 2)  # (x^2 + x + 1)^2
    assert all(h1.trinomial_irreducible(p, q) for p, q in h1.GFSR_PAIRS)


def test_gfsr_recurrence_and_determinism():
    g = h1.GFSR(31, 13, seed=5)
    state, i = list(g.state), g.i
    value = g.next()
    assert value == state[i] ^ state[(i + 31 - 13) % 31] and value < 2**32
    assert h1.GFSR(250, 103, 9).word() == h1.GFSR(250, 103, 9).word()


def words(source, n):
    return [source.read_bytes(8) for _ in range(n)]


@pytest.mark.parametrize("family,severity", [("SHIFT", 0.01), ("LAG1", 0.3), ("DUP", 0.2), ("GFSR", (31, 13))])
def test_prefix_identical_and_onset_exact(family, severity):
    onset = 8 * 100
    healthy = h1.FaultyWords(PCG64Source(1), "healthy", 0, 2, 3, onset)
    faulty = h1.FaultyWords(PCG64Source(1), family, severity, 2, 3, onset)
    a, b = words(healthy, 400), words(faulty, 400)
    assert a[:100] == b[:100] and a[100:] != b[100:]
    assert bytes(faulty.emitted) == b"".join(b) and bytes(healthy.emitted) == b"".join(a)


def test_shift_map_is_exact():
    base = h1.FaultyWords(PCG64Source(4), "healthy", 0, 0, 0, 0)
    shifted = h1.FaultyWords(PCG64Source(4), "SHIFT", 0.25, 0, 0, 0)
    for raw, out in zip(words(base, 2000), words(shifted, 2000)):
        w = struct.unpack(">Q", raw)[0]
        u = min((w >> 11) * 2.0**-53 + 0.25, h1.TOP)
        assert struct.unpack(">Q", out)[0] == (int(u * 2.0**53) << 11) | (w & 0x7FF)


def test_dup_follows_the_fault_generator():
    base = words(h1.FaultyWords(PCG64Source(8), "healthy", 0, 0, 0, 0), 3000)
    out = words(h1.FaultyWords(PCG64Source(8), "DUP", 0.1, 21, 0, 0), 3000)
    decisions = np.random.default_rng(21).random(3000) < 0.1
    expected, prev = [], None
    for raw, d in zip(base, decisions):
        prev = prev if (d and prev is not None) else raw
        expected.append(prev)
    assert out == expected and decisions.sum() > 200


def test_lag1_keeps_uniform_margins_and_induces_correlation():
    out = words(h1.FaultyWords(PCG64Source(12), "LAG1", 0.5, 0, 0, 0), 20000)
    u = np.array([(struct.unpack(">Q", w)[0] >> 11) * 2.0**-53 for w in out])
    assert abs(u.mean() - 0.5) < 0.01 and abs(u.var() - 1 / 12) < 0.003
    expected = 6 / math.pi * math.asin(0.5 / 2)  # Pearson correlation of a Gaussian copula's uniforms
    assert abs(np.corrcoef(u[:-1], u[1:])[0, 1] - expected) < 0.03


def test_short_run_records_and_pairing():
    healthy, arrays = h1.run_one(SMOKE, "unit", "healthy", 0, 0, sweeps=1300, onset_bytes=2**14)
    shifted, _ = h1.run_one(SMOKE, "unit", "SHIFT", 0.01, 0, sweeps=1300, onset_bytes=2**14)
    assert healthy["per_sweep_bytes_ok"] and len(arrays["marks"]) == 1300
    assert healthy["prefix_sha256"] == shifted["prefix_sha256"] and healthy["local_sha256"] != shifted["local_sha256"]
    marks = arrays["marks"]
    k = healthy["onset_sweep"]
    assert marks[k - 1] <= 2**14 < marks[k]
    assert [c["end"] for c in healthy["conventional"]] == [512, 768, 1024, 1280]
    assert [c["end"] for c in healthy["application"]] == [1024, 1280]
    assert "harm" not in healthy  # harm needs the full 10,240 sweeps


def test_arm_statistics_threshold_and_detection():
    cal = [{"application": [{"mean_D": d, "uphill_fraction": 0.4 + d}]} for d in np.linspace(-0.1, 0.1, 21)]
    constants = h1.application_constants(cal)
    assert constants["mu_f"] == pytest.approx(0.4)
    assert h1.threshold([3, 1, 2, 5, 4], 4) == 4.0
    series = [(512, 1.0), (768, 9.0), (1024, 2.0), (1280, 9.0)]
    out = h1.detection(series, tau=5.0, onset_sweep=900, marks=np.arange(1, 1281) * 200, onset_bytes=180000)
    assert out == {"whole_run_alarm": True, "post_onset_alarm": True, "alarm_end": 1280,
                   "latency_sweeps": 380, "latency_bytes": 1280 * 200 - 180000}
    assert h1.conventional_statistic({k: 0.5 for k in h1.CONV_P_KEYS} | {"runs": None, "frequency": 1e-4}) == pytest.approx(4)
