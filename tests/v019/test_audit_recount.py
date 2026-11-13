"""Independent recount library: known answers, cross-implementation agreement, edge cases.

Reference values are from NIST SP 800-22 Rev 1a (worked examples of sections
2.1-2.4 and 2.11-2.13), from closed forms, or constructed by hand. Agreement
with the primary code is a cross-check, not the definition of correctness.
"""
import json
import math
from fractions import Fraction

import numpy as np
import pytest
from scipy import optimize, special, stats

import importlib.util as _ilu
import sys as _sys
from pathlib import Path as _Path

_spec = _ilu.spec_from_file_location("_audit_loader", _Path(__file__).with_name("_audit_loader.py"))
au = _sys.modules.get("_audit_loader") or _ilu.module_from_spec(_spec)
if "_audit_loader" not in _sys.modules:
    _sys.modules["_audit_loader"] = au
    _spec.loader.exec_module(au)
from rnglab11.statistics import evaluate_bits
from rnglab12.uncertainty import summarize_series
from rnglab14.models.pimc_harmonic import exact_finite_bead_energy
from rnglab14.sources import PCG64Source

rc = au.recount()
h1 = au.h1lib()
PI_100 = ("1100100100001111110110101010001000100001011010001100001000110100"
          "110001001100011001100010100010111000")
LONGEST_128 = ("11001100000101010110110001001100111000000000001001001101010100010001001111010110"
               "100000001101011111001100111001101101100010110010")


def bits(text):
    return np.array([int(c) for c in text], dtype=np.uint8)


# --------------------------------------------------------------------------- battery: known answers
def test_nist_frequency_block_runs_examples():
    assert rc.frequency_test(bits("1011010101"))[1] == pytest.approx(0.527089, abs=5e-7)
    b = bits(PI_100)
    assert len(b) == 100 and int(b.sum()) == 42  # S_100 = -16
    z, p = rc.frequency_test(b)
    assert z == pytest.approx(-1.6) and p == pytest.approx(0.109599, abs=5e-7)
    chi2, p = rc.block_frequency_test(bits("0110011010"), 3)
    assert chi2 == pytest.approx(1.0) and p == pytest.approx(0.801252, abs=5e-7)
    chi2, p = rc.block_frequency_test(b, 10)
    assert chi2 == pytest.approx(7.2) and p == pytest.approx(0.706438, abs=5e-7)
    assert rc.runs_test(bits("1001101011"))[3] == pytest.approx(0.147232, abs=5e-7)
    v, _, _, p = rc.runs_test(b)
    assert v == 52 and p == pytest.approx(0.500798, abs=5e-7)


def test_nist_serial_entropy_cusum_longest_run_examples():
    s = rc.serial_test(bits("0011011101"), 3)
    assert s["psi2"] == pytest.approx([0.4, 1.2, 2.8])
    assert (s["p1"], s["p2"]) == pytest.approx((0.808792, 0.670320), abs=5e-7)
    assert rc.approximate_entropy_test(bits("0100110101"), 3)["p"] == pytest.approx(0.261961, abs=5e-7)
    a = rc.approximate_entropy_test(bits(PI_100), 2)
    assert a["chi2"] == pytest.approx(5.550792, abs=5e-6) and a["p"] == pytest.approx(0.235301, abs=5e-7)
    assert rc.cumulative_sums_forward_test(bits("1011010111")) == (4, pytest.approx(0.4116588, abs=5e-7))
    z, p = rc.cumulative_sums_forward_test(bits(PI_100))
    assert z == 16 and p == pytest.approx(0.219194, abs=5e-7)
    longest = rc.longest_run_test(bits(LONGEST_128))
    assert longest["block"] == 8 and longest["nu"] == [4, 9, 3, 0]
    assert longest["chi2"] == pytest.approx(4.882457, abs=5e-4)  # the document prints 4.882457
    assert longest["p"] == pytest.approx(0.1806, abs=2e-5)       # printed as 0.180598 and as 0.180609


def test_longest_run_per_block_matches_direct_count():
    rng = np.random.default_rng(11)
    b = (rng.random(5000) < 0.7).astype(np.uint8)
    direct = [max((len(run) for run in "".join(map(str, b[i:i + 128])).split("0")), default=0)
              for i in range(0, 5000 - 5000 % 128, 128)]
    assert rc.longest_runs_per_block(b, 128).tolist() == direct
    assert rc.longest_runs_per_block(np.ones(16, np.uint8), 8).tolist() == [8, 8]  # runs do not cross blocks
    assert rc.longest_runs_per_block(np.zeros(16, np.uint8), 8).tolist() == [0, 0]


@pytest.mark.parametrize("order", [1, 3, 9])
def test_pattern_counts_direct_and_shifted_agree(order):
    b = (np.random.default_rng(order).random(777) < 0.4).astype(np.uint8)
    text = "".join(map(str, b)) + "".join(map(str, b[:order - 1]))
    direct = np.zeros(2**order, dtype=np.int64)
    for i in range(len(b)):
        direct[int(text[i:i + order], 2)] += 1
    assert np.array_equal(rc.cyclic_pattern_counts(b, order), direct)
    assert np.array_equal(rc.cyclic_pattern_counts(b, order, (9, rc.cyclic_pattern_values(b, 9))), direct)


# --------------------------------------------------------------------------- battery: agreement with Stage C
@pytest.mark.parametrize("words,bias", [(2048, 0.0), (13000, 0.0), (13000, 0.02), (16384, 0.0)])
def test_battery_agrees_with_stage_c_window(words, bias):
    rng = np.random.default_rng(words + int(1000 * bias))
    w = rng.integers(0, 2**64, size=words, dtype=np.uint64)
    if bias:
        w[rng.random(words) < bias] = np.uint64(0xFFFFFFFFFFFFF800)  # saturated uniforms, as SHIFT produces
        w[100], w[101] = w[99], w[99]                               # adjacent duplicates
    mine = rc.battery_checkpoint(w)
    data = rc.words_to_bytes(w)
    assert 8 * len(data) > rc.BIT_CAP_BITS if words == 16384 else 8 * len(data) <= rc.BIT_CAP_BITS
    reference = h1.conventional_window(data[:125000])
    reference.update({k: h1.word_tests(data)[k] for k in ("lag1_correlation", "adjacent_duplicates")})
    for key, p in zip(rc.P_KEYS, mine["p_values"]):
        if reference[key] is None:
            assert p is None
        elif reference[key] <= 1e-290:
            assert p <= 1e-290
        else:
            assert p == pytest.approx(reference[key], rel=1e-9, abs=1e-300)
    assert mine["was_cut"] == (words == 16384) and mine["bit_used_bits"] == min(64 * words, 1_000_000)
    assert mine["word_count"] == words and mine["original_bits"] == 64 * words
    assert mine["S"] == pytest.approx(-math.log10(max(min(v for v in reference.values() if v is not None), 1e-300)))


def test_bit_cap_boundary_is_exactly_one_million_bits():
    rng = np.random.default_rng(5)
    exact = rng.integers(0, 2**64, size=15625, dtype=np.uint64)  # 1,000,000 bits: reachable, 524288 + 64*7433
    assert (1_000_000 - 512 * 128 * 8) % 64 == 0
    assert not rc.battery_checkpoint(exact)["was_cut"]
    over = np.concatenate([exact, rng.integers(0, 2**64, size=1, dtype=np.uint64)])
    cut = rc.battery_checkpoint(over)
    assert cut["was_cut"] and cut["bit_used_bits"] == 1_000_000 and cut["original_bits"] == 1_000_064
    kept = rc.battery_checkpoint(exact)
    assert cut["p_values"][:8] == kept["p_values"][:8]      # bit tests see only the first 1,000,000 bits
    assert cut["pair_count"] == kept["pair_count"] + 1      # word tests see every word
    over[-1] = over[-2]
    assert rc.battery_checkpoint(over)["duplicate_count"] == 1 and rc.battery_checkpoint(over)["S"] > 15


def test_runs_prerequisite_and_short_windows_are_not_applicable():
    w = np.full(4096, np.uint64(0xFFFFFFFFFFFFFFFF))
    out = rc.battery_checkpoint(w)
    base = evaluate_bits(np.ones(64 * 4096, np.uint8), block_size=max(128, 64 * 4096 // 99 + 1))["tests"]
    assert out["runs"] is None and base["runs"]["p_value"] is None and out["prerequisite_failed"]
    assert out["runs_signed_z"] is None  # pi = 1: no finite z
    assert out["lag1_correlation"] is None and out["lag1_r"] is None and out["duplicate_count"] == 4095  # zero variance
    assert out["S"] == 300.0  # floored
    assert rc.bit_statistics(np.zeros(0, np.uint8))["frequency"] is None
    assert rc.bit_statistics(bits(PI_100)[:99])["frequency"] is None  # n < 100


def test_alarm_statistic_rules():
    assert rc.alarm_statistic([None, float("nan"), 2.0, -0.1, float("inf")]) == (0.0, 0)
    assert rc.alarm_statistic([0.5, None, 1e-4]) == (pytest.approx(4.0), 2)
    assert rc.alarm_statistic([0.0, 0.3]) == (300.0, 2)
    assert rc.alarm_statistic([1.0]) == (0.0, 1)


# --------------------------------------------------------------------------- streams
def test_numpy_pcg64_reproduces_repository_source():
    for seed in (0, 1, 2**63 + 12345, 2**64 - 1):
        words = rc.pcg64_words(seed, 3000)
        assert rc.words_to_bytes(words) == PCG64Source(seed).read_bytes(24000)
        assert np.array_equal(rc.bytes_to_words(rc.words_to_bytes(words)), words)
    with pytest.raises(rc.RecountError):
        rc.bytes_to_words(b"1234567")


def test_seed_recipe_matches_frozen_definition():
    from rnglab.seeds import namespace_seed
    for cond, role, i in (("calibration", "local", 0), ("fault/SHIFT", "profile", 999), ("healthy", "start", 1999)):
        assert rc.namespace_seed("rnglab-v019-d-pilot", cond, role, i) == namespace_seed(
            "rnglab-v019-d-pilot", cond, role, f"{i:05d}")
    assert rc.analysis_seed("ns", "fault/DUP", "bootstrap-L20480") == namespace_seed(
        "ns", "fault/DUP", "analysis", "bootstrap-L20480")


def test_sweep_of_word_and_ledger():
    marks = np.array([128, 264, 520, 648], dtype=np.int64)
    assert rc.sweep_of_word(marks).tolist() == [0] * 16 + [1] * 17 + [2] * 32 + [3] * 16
    assert rc.check_ledger(marks, 4) == []
    assert rc.window_words(marks, 1, 3) == (16, 65) and rc.window_words(marks, 0, 1) == (0, 16)
    assert any("expected 5" in p for p in rc.check_ledger(marks, 5))
    assert rc.check_ledger(marks.astype(np.int32))
    for bad in ([128, 128], [128, 257], [120], [128, 392], [128, 261]):
        assert rc.check_ledger(np.array(bad, dtype=np.int64))
    with pytest.raises(rc.RecountError):
        rc.sweep_of_word(np.array([128, 129], dtype=np.int64))


def test_severity_and_closed_form_mean():
    caps = rc.protocol_caps()
    assert caps == {"SHIFT": 0.046415888336127774, "LAG1": 0.5, "DUP": 0.5}
    assert not rc.severity_path(None, None, None, 64).any()
    for t0, t_cap, s_max in ((3, 7.5, 0.5), (2048, 2048.0, 0.5), (8192, 163840.0, caps["SHIFT"]), (5, 4.0, 0.25)):
        path = rc.severity_path(t0, t_cap, s_max, 20480)
        assert path[t0] == 0 and path[t0 + 1] > 0 and not path[:t0].any()
        assert path.max() <= s_max and np.all(np.diff(path) >= 0)
        first_cap = t0 + math.ceil(t_cap)
        if first_cap < len(path):
            assert path[first_cap] == s_max and path[first_cap - 1] < s_max or t_cap == int(t_cap)
        for length in (t0, t0 + 1, t0 + 2, 10240, 20480):
            assert rc.s_bar_closed_form(length, t0, t_cap, s_max) == pytest.approx(
                math.fsum(path[:length]) / length, rel=1e-13, abs=1e-18)


def test_dup_forward_fill_equals_sequential_definition():
    rng = np.random.default_rng(3)
    base = rng.integers(0, 2**64, size=400, dtype=np.uint64)
    sweep = np.repeat(np.arange(20), 20)
    severity = rc.severity_path(0, 4.0, 0.9, 20)
    out = rc.apply_fault("DUP", base, sweep, severity, 0, 77)
    q = np.random.default_rng(77)
    expected, previous = [], None
    for i, w in enumerate(base.tolist()):
        draw = q.random()
        word = previous if draw < severity[sweep[i]] and previous is not None else w
        expected.append(word)
        previous = word
    assert out.tolist() == expected and len(set(expected)) < 400
    assert np.array_equal(rc.apply_fault("healthy", base, sweep, severity, None, 0), base)
    with pytest.raises(rc.RecountError):
        rc.apply_fault("GFSR", base, sweep, severity, 0, 0)


# --------------------------------------------------------------------------- endpoint
def test_reference_energy_normal_modes_equal_closed_form():
    for beads, beta in ((16, 2.0), (2, 0.5), (64, 4.0), (16, 0.1)):
        assert rc.reference_energy(beads, beta) == pytest.approx(exact_finite_bead_energy(beads, beta), rel=2e-15)
    assert rc.reference_energy(4096, 2.0) == pytest.approx(0.5 / math.tanh(1.0), rel=1e-7)  # continuum limit
    assert rc.T975_15 == pytest.approx(2.131450, abs=5e-7)  # tabulated t quantile
    assert rc.T975_15 == pytest.approx(stats.t.ppf(0.975, 15), rel=1e-14)


def ar1(seed, n, phi=0.9, mean=0.0):
    rng = np.random.default_rng(seed)
    e = rng.normal(size=n)
    x = np.empty(n)
    x[0] = e[0] / math.sqrt(1 - phi * phi)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x + mean


@pytest.mark.parametrize("n", [512, 1024, 10240, 20480])
def test_endpoint_agrees_with_primary_summary(n):
    for seed in range(12):
        x = ar1(seed, n, mean=0.3)
        mine = rc.endpoint(x, n, 0.3)
        ref = summarize_series(x[:n], reference=0.3, batch_sizes=(n // 16,))["batch_means"][str(n // 16)]
        assert ref["status"] == "ESTIMATED"
        assert mine["valid"] is ref["covered"] and mine["covered"] is ref["covered"]
        assert mine["estimate"] == pytest.approx(ref["prefix_mean"], rel=1e-13)
        assert mine["mcse"] == pytest.approx(ref["mcse"], rel=1e-12)
        assert mine["interval_95"] == pytest.approx(ref["interval_95"], rel=1e-12)
        assert mine["margin"] == pytest.approx(mine["half_width"] - abs(mine["signed_error"]))
        assert (mine["margin"] >= 0) is mine["covered"]


def test_endpoint_edge_cases():
    e = rc.endpoint(np.full(64, 1.5), 64, 1.5)
    assert not e["valid"] and e["reason"] == "DEGENERATE_BATCH_VARIANCE" and e["estimate"] == 1.5
    x = ar1(1, 64)
    assert rc.endpoint(np.r_[x[:63], np.nan], 64, 0.0)["reason"] == "NONFINITE_INPUT"
    assert rc.endpoint(np.r_[x[:63], np.inf], 64, 0.0)["reason"] == "NONFINITE_INPUT"
    assert rc.endpoint(x[:48], 64, 0.0)["reason"] == "INCOMPLETE_TRACE"
    assert rc.endpoint(x, 16, 0.0)["reason"] == "EXCEPTION"  # 16 samples: the prescribed call refuses fewer than 32
    with pytest.raises(ValueError):
        summarize_series(x[:16], reference=0.0, batch_sizes=(1,))
    assert rc.endpoint(x, None, 0.0) == {**rc.endpoint(x, None, 0.0), "valid": None, "reason": "NO_ALARM"}
    for bad in (0, -16, 17, 40):
        with pytest.raises(rc.RecountError):
            rc.endpoint(x, bad, 0.0)
    assert rc.issuance_length(None) is None and rc.issuance_length(1023) == 1008 and rc.issuance_length(768) == 768


def test_endpoint_coverage_boundary_is_inclusive_and_exact():
    # Batch means 0,1,...,15 exactly: estimate 7.5, variance 22.666..., so the upper bound is known.
    x = np.repeat(np.arange(16.0), 4)
    e = rc.endpoint(x, 64, 0.0)
    half = rc.T975_15 * math.sqrt((sum((k - 7.5)**2 for k in range(16)) / 15) / 16)
    assert e["estimate"] == 7.5 and e["half_width"] == pytest.approx(half, rel=1e-15)
    upper, lower = e["interval_95"][1], e["interval_95"][0]
    assert rc.endpoint(x, 64, upper)["covered"] and rc.endpoint(x, 64, lower)["covered"]
    assert not rc.endpoint(x, 64, np.nextafter(upper, np.inf))["covered"]
    assert not rc.endpoint(x, 64, np.nextafter(lower, -np.inf))["covered"]
    assert rc.endpoint(x, 64, upper)["margin"] == pytest.approx(0.0, abs=1e-15)
    primary = summarize_series(x, reference=upper, batch_sizes=(4,))["batch_means"]["4"]
    assert primary["covered"] is True and primary["interval_95"][1] == pytest.approx(upper, rel=1e-15)


# --------------------------------------------------------------------------- calibration, issuance, rates
def test_clopper_pearson_root_form_matches_beta_form():
    for k, n in ((0, 20), (20, 20), (1, 2000), (100, 2000), (1864, 2000), (1946, 2000), (3867, 4096), (199, 4096)):
        lo, hi = rc.clopper_pearson(k, n)
        ref_lo = 0.0 if k == 0 else stats.beta.ppf(0.025, k, n - k + 1)
        ref_hi = 1.0 if k == n else stats.beta.ppf(0.975, k + 1, n - k)
        assert (lo, hi) == pytest.approx((ref_lo, ref_hi), rel=1e-11, abs=1e-13)
    assert rc.clopper_pearson(0, 0) is None
    assert rc.clopper_pearson(3867, 4096) == pytest.approx([0.9366110953196926, 0.9509308125864511], abs=1e-12)  # Stage C V1
    with pytest.raises(rc.RecountError):
        rc.clopper_pearson(5, 4)
    lo, hi = rc.clopper_pearson(7, 50, 0.9)
    assert stats.binom.sf(6, 50, lo) == pytest.approx(0.05, rel=1e-9) and stats.binom.cdf(7, 50, hi) == pytest.approx(0.05, rel=1e-9)


def test_threshold_is_the_950th_order_statistic_and_alarm_is_strict():
    values = np.random.default_rng(0).permutation(np.arange(999.0))
    assert rc.order_statistic_threshold(values) == 949.0
    tied = np.r_[np.zeros(940), np.full(59, 7.0)]
    assert rc.order_statistic_threshold(tied) == 7.0  # ties retained; strict exceedance then gives no alarm at 7.0
    for bad in (values[:998], np.r_[values[:998], np.nan]):
        with pytest.raises(rc.RecountError):
            rc.order_statistic_threshold(bad)
    s = np.zeros(159)
    assert rc.first_alarm(s, 0.0) is None
    s[10] = 7.0
    assert rc.first_alarm(s, 7.0) is None and rc.first_alarm(s, np.nextafter(7.0, 0)) == 512 + 256 * 10
    s[3], s[40] = 9.0, 99.0
    assert rc.first_alarm(s, 7.0) == 512 + 256 * 3  # first exceedance, not the largest
    with pytest.raises(rc.RecountError):
        rc.first_alarm(np.r_[s[:158], np.nan], 7.0)
    with pytest.raises(rc.RecountError):
        rc.first_alarm(s[:158], 7.0)


def test_onset_class_and_eligibility_boundaries():
    assert rc.onset_class(None, 3000) == "no_alarm" and rc.onset_class(2048, None) is None
    assert rc.onset_class(3072, 3072) == "pre_onset"   # window ends at sweep 3071: no onset word observed
    assert rc.onset_class(3328, 3072) == "post_onset" and rc.onset_class(512, 8192) == "pre_onset"
    assert rc.eligible(10239, 10240) and not rc.eligible(10240, 10240) and not rc.eligible(None, 10240)
    assert rc.gate_v2d(100)["pass"] and not rc.gate_v2d(174)["pass"]  # passes iff at most 173 of 2000 alarm
    assert rc.gate_v2d(173)["cp95"][1] <= 0.10 < rc.gate_v2d(174)["cp95"][1]


# --------------------------------------------------------------------------- features and models
def test_trailing_slope_rules():
    ends = rc.BATTERY_ENDS.tolist()
    v = [0.5 * e / 256 + (0.3 if i % 2 else -0.3) for i, e in enumerate(ends)]
    for j in (3, 10, 158):
        slope, count = rc.trailing_slope(ends, v, j)
        assert count == 4 and slope == pytest.approx(np.polyfit(ends[j - 3:j + 1], v[j - 3:j + 1], 1)[0], rel=1e-9)
    assert rc.trailing_slope(ends, v, 0) == (0.0, 1)
    assert rc.trailing_slope(ends, v, 1) == (pytest.approx((v[1] - v[0]) / 256), 2)
    assert rc.trailing_slope(ends, v, 2)[1] == 3
    assert rc.trailing_slope(ends, [1.0, None, 2.0, 3.0, 4.0], 4) == (None, 4)
    assert rc.trailing_slope(ends, [None, 1.0, 2.0, 3.0, 4.0], 4)[0] == pytest.approx(1 / 256)  # only the last four count
    x = [Fraction(e) for e in ends[:4]]
    assert rc.trailing_slope(ends, [1e-300, 2e-300, 3e-300, 4e-300], 3)[0] == pytest.approx(1e-300 / 256)
    assert len(x) == 4


def test_application_alignment():
    assert rc.application_index(512) is None and rc.application_index(768) is None
    assert rc.application_index(1024) == 0 and rc.application_index(1279) == 0 and rc.application_index(1280) == 1
    assert rc.application_index(40960) == 156 == len(rc.APPLICATION_ENDS) - 1
    assert len(rc.BATTERY_ENDS) == 159 and rc.BATTERY_ENDS[[0, -1]].tolist() == [512, 40960]


def frozen_fit(z, y):
    """The fit exactly as protocol.json forecasters.logistic.call specifies it."""
    d = z.shape[1]
    return optimize.minimize(lambda t: rc.logistic_objective(t, z, y), np.zeros(d + 1), jac=True, method="L-BFGS-B",
                             bounds=None, options={"gtol": 1e-4, "ftol": 64 * np.finfo(np.float64).eps,
                                                   "maxiter": 100, "maxfun": 15000, "maxls": 50, "maxcor": 10})


@pytest.mark.parametrize("seed,n,d", [(0, 400, 1), (1, 300, 3), (2, 600, 5), (3, 240, 3)])
def test_newton_optimum_and_frozen_lbfgsb_agree(seed, n, d):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(n, d)) * rng.uniform(0.1, 30, size=d) + rng.normal(size=d)
    y = (rng.random(n) < special.expit(x @ rng.normal(size=d) / 10 + 0.3)).astype(float)
    mu, scale = rc.standardise(x)
    z = (x - mu) / scale
    assert z.mean(axis=0) == pytest.approx(0, abs=1e-12) and z.std(axis=0) == pytest.approx(1)
    newton = rc.fit_logistic_newton(z, y)
    assert newton["converged"] and newton["gradient_inf"] <= 1e-10
    numeric = optimize.approx_fprime(newton["theta"], lambda t: rc.logistic_objective(t, z, y)[0], 1e-6)
    assert numeric == pytest.approx(rc.logistic_objective(newton["theta"], z, y)[1], abs=1e-3)
    result = frozen_fit(z, y)
    assert result.success and result.status == 0 and result.nit <= 100
    assert result.message.startswith("CONVERGENCE:")  # the two accepted terminations in this SciPy
    assert result.fun >= newton["objective"] - 1e-9   # Newton reaches the minimum of the convex objective
    q_frozen = rc.logistic_predict(result.x, mu, scale, x)
    q_newton = rc.logistic_predict(newton["theta"], mu, scale, x)
    assert np.max(np.abs(q_frozen - q_newton)) < 2e-4  # gtol 1e-4 on the summed loss leaves this much room


def test_standardise_constant_feature_and_fit_guards():
    x = np.column_stack([np.arange(6.0), np.full(6, 3.0)])
    mu, scale = rc.standardise(x)
    assert scale[1] == 1.0 and mu[1] == 3.0 and scale[0] == pytest.approx(np.std(np.arange(6.0)))
    with pytest.raises(rc.RecountError):
        rc.fit_logistic_newton(x, np.ones(6))  # one class: the protocol falls back to F0
    with pytest.raises(rc.RecountError):
        rc.fit_logistic_newton(x[:0], np.zeros(0))


# --------------------------------------------------------------------------- scores
def test_brier_skill_and_bins():
    y = [1, 0, 1, 1]
    assert rc.brier([1, 0, 1, 1], y) == 0.0 and rc.brier([0.5] * 4, y) == 0.25
    assert rc.brier([], []) is None and rc.brier([float("nan"), 0, 0, 0], y) is None
    assert rc.skill(0.1, 0.2) == pytest.approx(0.5) and rc.skill(0.1, 0.0) is None and rc.skill(None, 0.2) is None
    rows = rc.calibration_bins([0.0, 0.19999, 0.2, 0.4, 0.6, 0.8, 1.0, 0.7999999], [0, 0, 1, 1, 0, 1, 1, 0])
    assert [r["count"] for r in rows] == [2, 1, 1, 2, 2]
    assert rows[4]["observed_fraction"] == 1.0 and rows[3]["observed_fraction"] == 0.0
    assert rc.calibration_bins([0.5], [1])[0] == {"bin": 0, "count": 0, "mean_predicted": None,
                                                  "observed_fraction": None, "cp95": None}
    protocol_index = lambda q: min(np.searchsorted([0.0, 0.2, 0.4, 0.6, 0.8, 1.0], q, side="right") - 1, 4)
    for q in np.r_[np.linspace(0, 1, 101), 0.2, 0.4, 0.6, 0.8, np.nextafter(0.2, 0), np.nextafter(0.8, 1)]:
        mine = [r["count"] for r in rc.calibration_bins([q], [1])].index(1)
        assert mine == protocol_index(q)
    with pytest.raises(rc.RecountError):
        rc.calibration_bins([1.2], [1])


def test_bootstrap_block_draw_equals_sequential_draws():
    for n in (1, 7, 240, 400):
        rng = np.random.default_rng(99)
        sequential = np.stack([rng.integers(0, n, size=n) for _ in range(200)])
        assert np.array_equal(rc.bootstrap_indices(99, n, 200), sequential)


def test_paired_bootstrap_against_manual_resampling():
    rng = np.random.default_rng(4)
    y = (rng.random(300) < 0.4).astype(int)
    q = {"F0": np.full(300, y.mean()), "F0+": np.clip(0.4 + 0.1 * rng.normal(size=300), 0, 1),
         "F-B-lite": np.clip(0.4 + 0.5 * (y - 0.4) + 0.1 * rng.normal(size=300), 0, 1)}
    out = rc.paired_bootstrap(y, q, seed=123, resamples=500)
    draws = np.random.default_rng(123)
    manual = []
    for _ in range(500):
        idx = draws.integers(0, 300, size=300)
        manual.append(1 - np.mean((q["F-B-lite"][idx] - y[idx])**2) / np.mean((q["F0+"][idx] - y[idx])**2))
    level = 1 - 0.05 / 3
    expected = np.quantile(manual, [(1 - level) / 2, 1 - (1 - level) / 2], method="linear")
    assert out["models"]["F-B-lite"]["skill_ci_criterion"] == pytest.approx(expected.tolist(), rel=1e-12)
    assert out["models"]["F0+"]["skill"] == 0.0 and out["models"]["F0+"]["skill_ci95"] == [0.0, 0.0]
    assert out["models"]["F-B-lite"]["skill"] == pytest.approx(
        1 - np.mean((q["F-B-lite"] - y)**2) / np.mean((q["F0+"] - y)**2))
    assert out["undefined_skill_resamples"] == 0 and out["n"] == 300


def test_paired_bootstrap_undefined_skill_is_not_discarded():
    y = np.array([1] * 5 + [0])
    q = {"F0+": np.array([1.0] * 5 + [0.5]), "F-B-lite": np.array([0.9] * 5 + [0.1])}
    out = rc.paired_bootstrap(y, q, seed=1, resamples=400)
    assert out["undefined_skill_resamples"] > 0  # resamples without the sixth run have a zero baseline score
    assert out["models"]["F-B-lite"]["skill_ci_criterion"] is None and out["models"]["F-B-lite"]["brier_ci95"] is not None
    assert rc.paired_bootstrap([], {"F0+": []}, seed=1) == {"n": 0, "models": {}}


# --------------------------------------------------------------------------- oracle
def test_oracle_pieces():
    assert rc.slope_through_origin([1, 2, 3], [2, 4, 6]) == pytest.approx(2.0)
    assert rc.slope_through_origin([0, 0], [1, 2]) is None
    assert rc.sigma_eff(10240, [0.01, 0.03]) == pytest.approx(math.sqrt(10240) * 0.02)
    assert rc.validity_probability_stable(0.0) == pytest.approx(0.95, abs=1e-12)
    assert rc.validity_probability_stable(-3.0) == pytest.approx(rc.validity_probability_stable(3.0), rel=1e-12)
    values = [rc.validity_probability_stable(d) for d in np.arange(0, 60, 0.5)]
    assert all(0 <= v <= 1 for v in values) and all(a >= b for a, b in zip(values, values[1:]))


def test_stable_validity_probability_matches_simulation():
    rng = np.random.default_rng(8)
    for delta in (0.0, 1.5, 4.0):
        z = rng.normal(size=(200000, 16)) + delta / 4.0  # sixteen unit-variance batch means, mean shift delta/4
        t = z.mean(axis=1) / (z.std(axis=1, ddof=1) / 4.0)
        observed = np.mean(np.abs(t) <= rc.T975_15)
        p = rc.validity_probability_stable(delta)
        assert abs(observed - p) < 5 * math.sqrt(p * (1 - p) / 200000) + 1e-5


def test_frozen_oracle_call_is_nan_or_correct():
    """protocol.json oracle.probability fixes scipy.stats.nct. In SciPy 1.18.1 its lower-tail term
    is NaN for scattered noncentralities from about 8.75; where it is finite it is accurate."""
    nan_points = []
    for delta in np.arange(0.0, 45.01, 0.25):
        frozen = rc.validity_probability_frozen(float(delta))
        if math.isnan(frozen):
            nan_points.append(float(delta))
        else:
            assert frozen == pytest.approx(rc.validity_probability_stable(float(delta)), abs=1e-12)
    assert all(d >= 8.0 for d in nan_points)  # never in the range where the probability is not negligible
    assert all(rc.validity_probability_stable(d) < 1e-9 for d in nan_points)


# --------------------------------------------------------------------------- criteria
def test_criterion_1_counts_and_design_probabilities():
    accepted = [k for k in range(1800, 2001) if rc.criterion_1({10240: k})["result"] == "PASS"]
    assert (accepted[0], accepted[-1]) == (1864, 1946) and accepted == list(range(1864, 1947))
    design = json.loads((au.REPO / "research/v019-d-pilot/design.json").read_text())
    option = next(o for o in design["size_options"] if o["healthy"] == 2000 and o["fault_per_family"] == 1000)
    for row in option["criterion_1"]["probabilities"]:
        p = row["true_coverage"]
        independent = stats.binom.cdf(1946, 2000, p) - stats.binom.cdf(1863, 2000, p)
        assert row["p_pass"] == pytest.approx(independent, rel=1e-9)
    c = rc.criterion_1({10240: 1900, 20480: 1863, 40960: 1900})
    assert c["result"] == "INSUFFICIENT_INFORMATION" and not c["stops_D"]
    assert c["per_length"][20480]["label"] == "STRADDLES_BAND"
    c = rc.criterion_1({10240: 1900, 20480: 1800, 40960: 1900})
    assert c["result"] == "FAIL" and c["stops_D"] and c["per_length"][20480]["label"] == "OUTSIDE_BAND"
    assert rc.criterion_1({10240: 1990})["per_length"][10240]["label"] == "OUTSIDE_BAND"
    assert rc.criterion_1({1: 1900, 2: 1900, 3: 1900})["result"] == "PASS"


def family_rows(n_alarm, valid_fraction, invalid_fraction, alarm=12000, total=1000):
    alarms = [alarm] * n_alarm + [None] * (total - n_alarm)
    valid = [i < valid_fraction * n_alarm for i in range(total)]
    invalid = [i < invalid_fraction * n_alarm for i in range(total)]
    return {"alarm": alarms, "valid": valid, "already_invalid": invalid}


def test_criteria_2_and_3_votes():
    good, decided = family_rows(900, 0.5, 0.1), family_rows(900, 0.95, 0.1)
    out = rc.criteria_2_3({"SHIFT": good, "LAG1": good, "DUP": decided})
    assert out["criterion_2"] == "PASS" and out["criterion_3"] == "PASS"
    assert out["per_family"]["DUP"]["criterion_2"] is False and out["per_family"]["SHIFT"]["median_remaining"] == 8480
    out = rc.criteria_2_3({"SHIFT": good, "LAG1": decided, "DUP": decided})
    assert out["criterion_2"] == "FAIL" and out["criterion_3"] == "PASS"
    few = family_rows(299, 0.5, 0.1)
    out = rc.criteria_2_3({"SHIFT": good, "LAG1": few, "DUP": decided})
    assert out["criterion_2"] == "INSUFFICIENT_INFORMATION" and out["per_family"]["LAG1"]["criterion_2"] is None
    late = family_rows(900, 0.5, 0.1, alarm=15361)
    assert rc.criteria_2_3({"SHIFT": late, "LAG1": late, "DUP": good})["criterion_3"] == "FAIL"
    edge = family_rows(900, 0.5, 0.1, alarm=15360)
    assert rc.criteria_2_3({"SHIFT": edge, "LAG1": edge, "DUP": edge})["criterion_3"] == "PASS"  # median exactly 5120
    half_invalid = family_rows(900, 0.5, 0.5)
    assert rc.criteria_2_3({"SHIFT": half_invalid, "LAG1": half_invalid, "DUP": good})["criterion_2"] == "FAIL"
    unknown = {**good, "alarm": ["UNKNOWN"] + good["alarm"][1:]}
    assert rc.criteria_2_3({"SHIFT": unknown, "LAG1": good, "DUP": decided})["criterion_2"] == "INSUFFICIENT_INFORMATION"
    at_length = family_rows(900, 0.5, 0.1, alarm=20480)
    assert rc.criteria_2_3({"SHIFT": at_length, "LAG1": good, "DUP": good})["per_family"]["SHIFT"]["eligible"] == 0


def test_criterion_4_votes():
    def result(n, skill, lower, baseline=0.2):
        return {"n": n, "models": {"F0+": {"brier": baseline}, "F-B-lite": {
            "skill": skill, "skill_ci_criterion": None if lower is None else [lower, 0.5]}}}
    assert rc.criterion_4({"SHIFT": result(400, 0.12, 0.01), "LAG1": result(400, 0.0, -0.1)})["result"] == "PASS"
    assert rc.criterion_4({"SHIFT": result(400, 0.0999, 0.05), "LAG1": result(400, 0.3, 0.0)})["result"] == "FAIL"
    assert rc.criterion_4({"SHIFT": result(239, 0.5, 0.2), "LAG1": result(400, 0.0, -0.1)})["result"] == "INSUFFICIENT_INFORMATION"
    assert rc.criterion_4({"SHIFT": result(400, 0.5, None)})["result"] == "INSUFFICIENT_INFORMATION"
    assert rc.criterion_4({"SHIFT": result(400, 0.5, 0.2, baseline=0.0)})["result"] == "INSUFFICIENT_INFORMATION"


def test_compare_is_exact_for_flags_and_null_patterns():
    assert rc.compare(1.0, 1.0 + 5e-13) and not rc.compare(1.0, 1.0 + 5e-11)
    assert rc.compare(None, None) and not rc.compare(None, 0.0) and not rc.compare(0.0, None)
    assert rc.compare(True, True) and not rc.compare(True, 1) and not rc.compare(1.0, True)
    assert rc.compare([1.0, None, "a"], (1.0, None, "a")) and not rc.compare([1.0], [1.0, 2.0])
    assert rc.compare({"a": [1, 2]}, {"a": [1, 2]}) and not rc.compare({"a": 1}, {"b": 1})
    assert rc.compare(3, 3) and not rc.compare(3, 4) and not rc.compare(float("nan"), float("nan"))


def test_criterion_5_label_follows_the_two_deadlines():
    assert rc.criterion_5("2027-04-30", "2026-12-22", "2027-01-30")["result"] == "PASS"
    assert rc.criterion_5("2027-04-30", "2026-12-23", "2027-01-30")["result"] == "FAIL"      # engineering one day late
    assert rc.criterion_5("2027-04-30", "2026-12-22", "2027-01-31")["result"] == "FAIL"      # results one day late
    for missing in range(3):
        for absent in (None, ""):
            dates = ["2027-04-30", "2026-12-01", "2027-01-10"]
            dates[missing] = absent
            out = rc.criterion_5(*dates)
            assert out["result"] == "INSUFFICIENT_INFORMATION" and out["date_problems"] == []
    out = rc.criterion_5("2027-04-30", "2026-12-01", "2027-01-10")
    assert out["engineering_deadline"] == "2026-12-22" and out["timezone"] == "Europe/Budapest"
    # the frozen text states the same date: approval 2026-11-10 plus six weeks
    assert "gives 2026-12-22" in rc.load_protocol()["stop_criteria"]["5"]["calculation"]


@pytest.mark.parametrize("deadline, expected", [
    ("2027-04-30", "2027-01-30"), ("2027-05-31", "2027-02-28"), ("2028-05-31", "2028-02-29"),
    ("2027-05-29", "2027-02-28"), ("2027-02-15", "2026-11-15"), ("2027-03-31", "2026-12-31"),
    ("2027-01-31", "2026-10-31"), ("2027-07-31", "2027-04-30"), ("2027-12-31", "2027-09-30"),
    ("2028-03-01", "2027-12-01")])
def test_criterion_5_three_calendar_months_with_clamping(deadline, expected):
    assert rc.criterion_5(deadline, "2026-12-01", "2026-12-02")["results_deadline"] == expected


@pytest.mark.parametrize("position, value", [
    (0, "2027-4-30"), (0, "30/04/2027"), (0, 20270430), (0, "2027-02-30"), (1, "2026-12-22T10:00:00"),
    (2, "soon"), (2, 5.0), (1, True), (0, "20270430"), (0, "2027-W17-5"), (2, "2027-01-10\n")])
def test_criterion_5_gives_no_label_for_a_text_that_is_not_a_date(position, value):
    dates = ["2027-04-30", "2026-12-01", "2027-01-10"]
    dates[position] = value
    out = rc.criterion_5(*dates)
    assert out["result"] is None and len(out["date_problems"]) == 1 and "not a calendar date" in out["date_problems"][0]


def test_criterion_5_dates_against_the_manifest_date():
    stamp = "2026-12-05T10:00:00+0100"
    clean = rc.criterion_5("2027-04-30", "2026-12-05", "2026-12-05", stamp)       # the same day is consistent
    assert clean["result"] == "PASS" and clean["date_problems"] == [] and clean["manifest_date"] == "2026-12-05"
    late = rc.criterion_5("2027-04-30", "2026-12-06", "2026-12-20", stamp)        # engineering after the manifest
    assert late["result"] == "PASS" and len(late["date_problems"]) == 1 and "freeze steps 4 and 6" in late["date_problems"][0]
    early = rc.criterion_5("2027-04-30", "2026-12-01", "2026-12-04", stamp)       # results before the manifest
    assert early["result"] == "PASS" and len(early["date_problems"]) == 1 and "freeze steps 6 and 7" in early["date_problems"][0]
    # the date is taken in the protocol's time zone: 23:30 UTC is the next day in Budapest
    assert rc.criterion_5("2027-04-30", "2026-12-01", "2026-12-09", "2026-12-05T23:30:00+0000")["manifest_date"] == "2026-12-06"
    assert rc.manifest_date("2026-12-05T10:00:00", "Europe/Budapest").isoformat() == "2026-12-05"   # no offset: as written
    broken = rc.criterion_5("2027-04-30", "2026-12-01", "2026-12-09", "yesterday")
    assert broken["manifest_date"] is None and "written_at" in broken["date_problems"][0]
