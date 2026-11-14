"""Per-run Stage D engineering: smoke namespace, fixed random bytes, no campaign."""
import hashlib
import importlib.util
import inspect
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from scipy.special import gammaincc
from scipy.stats import norm, poisson, t

from rnglab12.uncertainty import summarize_series

REPO = Path(__file__).resolve().parents[2]
_spec = importlib.util.spec_from_file_location("v019_measure", REPO / "research/v019-d-pilot/measure.py")
m = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(m)
h1 = m.h1lib
SMOKE = "rnglab-v019-d-pilot-smoke"
SPEC = json.loads((REPO / "research/v019-d-pilot/protocol.json").read_text())


def marks_for(size):
    # Evenly distribute complete words over 512 sweeps, including an exact cap.
    return (np.arange(1, 513, dtype=np.int64) * (size // 8) // 512) * 8


@pytest.mark.parametrize("size", [65536, 125000, 131072])
def test_cap_exact_calls_and_metadata(size):
    window = np.random.default_rng(831).bytes(size)
    c = m.battery_checkpoint(window, marks_for(size), 512)
    p = h1.conventional_window(window[:125000])
    if size > 125000:
        p.update({k: h1.word_tests(window)[k] for k in h1.CONV_P_KEYS[8:]})
    assert c["p_values"] == [p[k] for k in h1.CONV_P_KEYS]
    # Compare binary representations too, including underflow zeros.
    for actual, expected in zip(c["p_values"], p.values()):
        assert actual is None if expected is None else np.float64(actual).tobytes() == np.float64(expected).tobytes()
    assert c["was_cut"] == (size > 125000)
    assert c["original_bytes"] == c["word_used_bytes"] == size
    assert c["original_bits"] == size * 8
    assert c["bit_used_bytes"] == min(size, 125000)
    assert c["bit_used_bits"] == min(size * 8, 1000000)
    assert c["word_count"] == size // 8 and c["pair_count"] == size // 8 - 1
    assert c["byte_start"] == 0 and c["byte_end"] == size
    assert c["complete"] and not c["failures"]
    base = m.evaluate_bits(np.unpackbits(np.frombuffer(window[:125000], dtype=np.uint8)),
                           block_size=max(128, c["bit_used_bits"] // 99 + 1))["tests"]
    assert c["p_statuses"][:3] == [base[k]["status"] for k in h1.CONV_P_KEYS[:3]]
    assert c["p_reasons"][:3] == [base[k].get("reason") for k in h1.CONV_P_KEYS[:3]]
    for entry in SPEC["signal_schema"]["checkpoint_metadata"]:
        assert set(entry["name"].split(",")) <= c.keys()


def fault_window(family):
    raw = np.random.default_rng(619).bytes(131072)
    if family == "random":
        return raw
    cap = m.dlib.CAPS[family]
    p = {"t0": 0, "T_cap": 1.0, "r": cap, "s_max": cap}
    seeds = m.dlib.seeds_for(SMOKE, "fault/" + family, 0)
    return m.dlib.reapply(family, p, seeds, marks_for(len(raw)), len(raw), step=True, base_stream=raw)


def cumsum_p(z, n):
    # Independent truncation toward zero, deliberately not Python floor division.
    div = lambda a, b: int(a / b)
    hi = div(div(n, z) - 1, 4)
    k1 = np.arange(div(div(-n, z) + 1, 4), hi + 1)
    k2 = np.arange(div(div(-n, z) - 3, 4), hi + 1)
    sqn = math.sqrt(n)
    return float(1 - np.sum(norm.cdf((4*k1+1)*z/sqn) - norm.cdf((4*k1-1)*z/sqn))
                 + np.sum(norm.cdf((4*k2+3)*z/sqn) - norm.cdf((4*k2+1)*z/sqn)))


@pytest.mark.parametrize("family", ["random", "SHIFT", "LAG1", "DUP"])
def test_all_statistics_and_p_formulas(family):
    data = fault_window(family)
    c = m.battery_checkpoint(data, marks_for(len(data)), 512)
    bits = np.unpackbits(np.frombuffer(data[:125000], dtype=np.uint8))
    n = len(bits)
    expected = h1.conventional_window(data[:125000])
    words = h1.word_tests(data)
    expected.update({k: words[k] for k in h1.CONV_P_KEYS[8:]})
    assert c["p_values"] == [expected[k] for k in h1.CONV_P_KEYS]
    formulas = [math.erfc(abs(c["frequency_signed_z"]) / math.sqrt(2)),
                gammaincc(c["blocks"] / 2, c["block_frequency_chi2"] / 2),
                math.erfc(abs(c["runs_signed_z"]) / math.sqrt(2)),
                gammaincc(2**6, c["serial_delta1"] / 2), gammaincc(2**5, c["serial_delta2"] / 2),
                gammaincc(2**7, c["approximate_entropy_chi2"] / 2),
                cumsum_p(c["cumulative_sums_max_abs"], n),
                gammaincc(h1.LONGEST_RUN_TABLES[c["longest_run_table_key"]][0] / 2, c["longest_run_chi2"] / 2),
                math.erfc(abs(c["lag1_test_z"]) / math.sqrt(2)),
                1.0 if c["duplicate_count"] == 0 else max(float(poisson.sf(c["duplicate_count"] - 1,
                                                                                 c["duplicate_lambda"])), h1.P_FLOOR)]
    for p, reconstructed in zip(c["p_values"], formulas):
        if p is not None:
            assert reconstructed == pytest.approx(p, rel=2e-12, abs=1e-300)
    assert c["ones"] == int(bits.sum())
    assert c["frequency_signed_z"] == (2 * c["ones"] - n) / math.sqrt(n)
    assert c["V"] == 1 + np.count_nonzero(np.diff(bits))
    assert c["pi"] == c["ones"] / n
    assert c["runs_residual"] == c["V"] - 2 * n * c["pi"] * (1 - c["pi"])
    assert c["runs_signed_z"] == c["runs_residual"] / (2 * math.sqrt(n) * c["pi"] * (1 - c["pi"]))
    assert c["prerequisite_failed"] == (abs(c["pi"] - .5) >= 2 / math.sqrt(n))
    counts = h1.pattern_counts(bits, 9)
    for k in (6, 7, 8, 9):
        assert c[f"pattern_counts_{k}"] == counts[k].tolist()
        assert sum(c[f"pattern_counts_{k}"]) == n
    assert c["serial_psi2"] == [2**k / n * np.sum(counts[k].astype(float)**2) - n for k in (6, 7, 8)]
    for k in (8, 9):
        fractions = counts[k][counts[k] > 0] / n
        assert c[f"approximate_entropy_phi{k}"] == np.sum(fractions * np.log(fractions))
    assert c["approximate_entropy"] == c["approximate_entropy_phi8"] - c["approximate_entropy_phi9"]
    assert c["cumulative_sums_max_abs"] == max(abs(np.cumsum(2 * bits.astype(np.int64) - 1)))
    M = c["longest_run_block_size"]
    K, edges, probs = h1.LONGEST_RUN_TABLES[M]
    longest = []
    for row in bits[:c["longest_run_blocks"] * M].reshape(-1, M):
        # Independent byte-string longest ones recount (no Stage C helper).
        longest.append(max(map(len, "".join(map(str, row)).split("0"))))
    categories = np.clip(longest, edges[0], edges[-1])
    assert c["longest_run_categories"] == [int(np.sum(categories == e)) for e in edges]
    assert c["longest_run_edges"] == list(edges)
    assert c["longest_run_unused_bits"] == n % M
    assert c["unused_bits"] == n % c["block_size"]
    assert c["lag1_r"] == words["lag1_r"]
    assert c["lag1_test_z"] == math.sqrt(c["word_count"]) * c["lag1_r"]
    assert c["lag1_feature"] == math.sqrt(c["pair_count"]) * c["lag1_r"]
    assert c["duplicate_count"] == words["duplicate_count"]
    assert c["duplicate_rate"] == c["duplicate_count"] / c["pair_count"]
    assert c["duplicate_lambda"] == c["pair_count"] * 2.0**-64
    for entry in SPEC["signal_schema"]["battery_statistics"]:
        for key in entry["name"].split(","):
            assert key in c and key in c["statistic_statuses"] and key in c["statistic_reasons"]
    if family == "SHIFT":
        assert c["frequency_signed_z"] > 0 and c["prerequisite_failed"]
        assert c["p_values"][2] is None and c["p_statuses"][2] == "NOT_APPLICABLE"
    if family == "LAG1":
        assert c["lag1_r"] > 0
    if family == "DUP":
        assert c["duplicate_rate"] > .4


def test_statistic_guard_and_subset():
    assert m.statistic([]) == (0.0, 0)
    assert m.statistic([None, np.nan, np.inf, -1, 1.1]) == (0.0, 0)
    assert m.statistic([0, .5, None]) == (300.0, 2)
    assert m.statistic([1e-310]) == (300.0, 1)
    data = fault_window("DUP")
    c = m.battery_checkpoint(data, marks_for(len(data)), 512)
    assert (c["S"], c["available_p_count"]) == m.statistic(c["p_values"])
    assert (c["S_prime"], c["available_p_count_prime"]) == m.statistic(c["p_values"][:8])
    # Explicit synthetic column vector separates subset and whole battery.
    ps = [.5] * 8 + [1e-100, 0]
    assert m.statistic(ps)[0] == 300 and m.statistic(ps[:8])[0] == -math.log10(.5)


def test_unavailable_statistics_and_missing_checkpoint():
    c = m.battery_checkpoint(bytes(65536), marks_for(65536), 512)
    assert c["runs_signed_z"] is None and c["runs_residual"] == 1
    assert c["lag1_r"] is None and c["p_values"][8] is None
    assert c["p_statuses"][8] == "NOT_APPLICABLE" and c["complete"]
    c = m.battery_checkpoint(b"", np.array([], dtype=np.int64), 512)
    assert c["p_values"] == [None] * 10 and c["S"] == 0
    assert c["byte_start"] is None and not c["complete"]
    assert all(v == "ERROR" for v in c["p_statuses"])
    json.dumps(c, allow_nan=False)
    for byte, sign in ((b"\x0f", -1), (b"\x55", 1)):
        c = m.battery_checkpoint(byte * 65536, marks_for(65536), 512)
        assert sign * c["runs_residual"] > 0 and sign * c["runs_signed_z"] > 0


def test_test_failures_and_original_unavailable_p(monkeypatch):
    data = np.random.default_rng(182).bytes(65536)
    original = m.evaluate_bits

    def broken_values(*args, **kwargs):
        result = original(*args, **kwargs)
        for key, p in zip(h1.CONV_P_KEYS[:3], [np.nan, np.inf, 2.0]):
            result["tests"][key]["p_value"] = p
        return result

    monkeypatch.setattr(m, "evaluate_bits", broken_values)
    c = m.battery_checkpoint(data, marks_for(len(data)), 512)
    assert c["p_values"][:3] == [None] * 3 and c["p_statuses"][:3] == ["ERROR"] * 3
    assert c["p_original_values"][:3] == ["nan", "inf", 2.0]
    assert not c["complete"] and len(c["failures"]) == 3
    json.dumps(c, allow_nan=False)
    monkeypatch.setattr(m, "evaluate_bits", original)

    def failed_word_test(data):
        raise ArithmeticError("injected engineering failure")

    # Replace only the pilot's module reference, never h1lib itself.
    monkeypatch.setattr(m, "h1lib", SimpleNamespace(**{**vars(h1), "word_tests": failed_word_test}))
    c = m.battery_checkpoint(data, marks_for(len(data)), 512)
    assert not c["complete"] and c["complete_prime"]
    assert c["p_statuses"][8:] == ["ERROR", "ERROR"]
    assert c["failures"][0]["exception"]["class"] == "ArithmeticError"
    cs = checkpoints()
    cs[0] = c
    assert m.issuance(cs, -1)["status"] == "UNKNOWN"
    assert m.issuance(cs, -1, prime=True)["alarm_sweep"] == 512


@pytest.mark.parametrize("L", m.LENGTHS)
def test_endpoint_selected_batch_and_all_fields(L):
    cv = np.random.default_rng(23).normal(m.E_M, .1, m.H)
    e = m.endpoint(cv, L, m.E_M)
    batch = summarize_series(cv[:L], reference=m.E_M, batch_sizes=(L // 16,))["batch_means"][str(L // 16)]
    assert e["source_status"] == batch["status"]
    assert e["valid"] == (batch["status"] == "ESTIMATED" and batch["covered"] is True)
    assert e["y"] == int(e["valid"])
    assert e["covered"] == batch["covered"]
    assert e["estimate"] == batch["prefix_mean"] and e["mcse"] == batch["mcse"]
    assert e["interval_95"] == batch["interval_95"]
    assert e["half_width"] == t.ppf(.975, 15) * batch["mcse"]
    assert e["batch_means"] == cv[:L].reshape(16, L // 16).mean(axis=1).tolist()
    assert e["signed_error"] == e["estimate"] - m.E_M
    assert e["input_count"] == L and e["severity_path"] == {"array": "severity", "start": 0, "stop": L}
    assert set(e) == {"L", "input_count", "source_status", "valid", "y", "covered", "invalid_reason", "exception",
                      "estimate", "signed_error", "mcse", "half_width", "interval_95", "batch_means", "severity_path"}


def test_endpoint_invalid_and_precedence():
    cases = [(np.arange(512, dtype=float), 512, "UNCOVERED"),
             (np.full(512, m.E_M), 512, "DEGENERATE_BATCH_VARIANCE"),
             (np.r_[np.nan, np.zeros(511)], 512, "NONFINITE_INPUT"),
             (np.zeros(511), 512, "INCOMPLETE_TRACE"),
             (np.r_[np.nan, np.zeros(15)], 512, "NONFINITE_INPUT"),
             (np.zeros((512, 1)), 512, "INCOMPLETE_TRACE"),
             (np.full(512, 1e308), 512, "EXCEPTION"),
             (np.arange(16, dtype=float), 16, "EXCEPTION")]
    for cv, n, reason in cases:
        e = m.endpoint(cv, n, m.E_M)
        assert not e["valid"] and e["invalid_reason"] == reason
        if reason == "EXCEPTION":
            assert e["exception"]["class"] == "ValueError" and e["interval_95"] is None
        json.dumps(e, allow_nan=False)
    # A later nonfinite sample must not invalidate a complete earlier prefix.
    cv = np.r_[np.random.default_rng(5).normal(m.E_M, .1, 512), np.nan]
    assert m.endpoint(cv, 512, m.E_M)["invalid_reason"] != "NONFINITE_INPUT"
    for boundary in m.endpoint(cv, 512, m.E_M)["interval_95"]:
        assert m.endpoint(cv, 512, boundary)["valid"]  # inclusive coverage


def test_endpoint_unreachable_status_and_malformed_contract(monkeypatch):
    # Only replace the new module's binding; never mutate frozen modules.
    cv = np.arange(512, dtype=float)
    result = {"status": "INSUFFICIENT_BATCHES", "prefix_mean": None, "mcse": None,
              "interval_95": None, "covered": None}
    monkeypatch.setattr(m, "summarize_series", lambda *a, **k: {"batch_means": {"32": result}})
    assert m.endpoint(cv, 512, 0)["invalid_reason"] == "INSUFFICIENT_BATCHES"
    result["status"] = "FUTURE_STATUS"
    e = m.endpoint(cv, 512, 0)
    assert e["invalid_reason"] == "OTHER_STATUS" and e["source_status"] == "FUTURE_STATUS"
    for change in ({"mcse": np.inf}, {"covered": 1}, {"status": "ESTIMATED"}, {"interval_95": [0]}):
        original = dict(result)
        result.update(change)
        assert m.endpoint(cv, 512, 0)["invalid_reason"] == "MALFORMED_RESULT"
        result.clear()
        result.update(original)
    result.pop("covered")
    assert m.endpoint(cv, 512, 0)["invalid_reason"] == "MALFORMED_RESULT"


def checkpoints():
    return [{"end": end, "S": 2.0, "S_prime": 1.0, "complete": True} for end in m.BATTERY_ENDS]


def test_issuance_strict_first_onset_blind_and_eligibility():
    assert "t0" not in inspect.signature(m.issuance).parameters
    cs = checkpoints()
    assert m.issuance(cs, 2)["status"] == "NO_ALARM"
    cs[1]["S"], cs[3]["S"] = 2.1, 9.0
    result = m.issuance(list(reversed(cs)), 2)
    assert result["alarm_sweep"] == 768 and result["S_at_alarm"] == 2.1
    assert result["whole_horizon_max"] == 9 and result["whole_horizon_alarm"]
    assert m.issuance(cs, 2, prime=True)["status"] == "NO_ALARM"
    cs[0]["S_prime"] = 3
    assert m.issuance(cs, 2, prime=True)["alarm_sweep"] == 512
    for onset, label, latency in ((769, "pre_onset", -1), (768, "pre_onset", 0), (767, "post_onset", 1)):
        assert m.onset_classification(768, onset) == {"onset_classification": label, "alarm_minus_t0": latency}
    for L in (512, 768, 1024):
        assert m.eligibility(768, L)["eligible"] == (768 < L)
    assert m.eligibility(None, 1024)["reason"] == "no_alarm"
    assert m.eligibility(1024, 1024)["reason"] == "alarm_at_or_after_L"
    assert m.onset_classification(None, 768)["onset_classification"] == "no_alarm"
    assert m.issuance_endpoint([], None) == {"n": None, "endpoint": None, "already_invalid": None}
    cs = checkpoints()
    result = m.issuance(cs, 2)
    assert result["alarm_sweep"] is result["S_at_alarm"] is None
    assert not result["whole_horizon_alarm"]
    assert m.issuance(cs[1:], 2)["status"] == "UNKNOWN"
    assert m.eligibility(None, 1024, known=False)["eligible"] is None
    assert m.onset_classification(None, 768, known=False)["onset_classification"] == "unknown"
    cs[0]["S"] = 3
    assert m.issuance(cs[:-1], 2)["alarm_sweep"] == 512
    assert m.issuance(cs[:-1], 2)["whole_horizon_max"] is None
    cs[0]["complete"] = False
    assert m.issuance(cs, 2)["status"] == "UNKNOWN"
    cs[1]["S"] = 4
    assert m.issuance(cs, 2)["whole_horizon_alarm"] is True


def test_issuance_endpoint_and_healthy_short_intervals():
    cv = np.random.default_rng(111).normal(m.E_M, .1, m.H)
    result = m.issuance_endpoint(cv, 777)
    assert result["n"] == 768 and result["endpoint"] == m.endpoint(cv, 768, m.E_M)
    assert result["already_invalid"] == (not result["endpoint"]["valid"])
    short = m.short_interval_check(cv)
    assert [e["L"] for e in short] == list(m.BATTERY_ENDS)
    assert short[0] == m.endpoint(cv, 512, m.E_M)
    assert short[-1] == m.endpoint(cv, m.H, m.E_M)


def test_application_ledger_sign_and_constants():
    k = np.arange(2048, dtype=np.int64) % 17
    marks = np.cumsum(128 + 8 * k)
    cv = np.ones(2048)
    td = cv + np.linspace(-2, 4, 2048)
    c = m.application_checkpoint(cv, td, marks, 1536)
    assert c["mean_D"] == np.mean(td[512:1536] - cv[512:1536]) > 0
    assert c["count_D"] == c["finite_input_count"] == 1024
    assert c["uphill_count"] == int(k[512:1536].sum())
    assert c["uphill_count"] == sum((marks[t] - marks[t - 1] - 128) / 8 for t in range(512, 1536))
    assert c["proposal_count"] == 16384
    assert c["uphill_fraction"] == c["uphill_count"] / 16384
    assert m.application_checkpoint(td, cv, marks, 1536)["mean_D"] < 0
    calibration = [{"application": [m.application_checkpoint(cv, td, marks, e) for e in (1024, 1280, 1536)]}]
    constants = h1.application_constants(calibration)
    z = m.standardise_application(c, constants)
    assert z["A1_over_sigma_D"] == c["mean_D"] / constants["sigma_D"]
    assert z["A2_standardised"] == (c["uphill_fraction"] - constants["mu_f"]) / constants["sigma_f"]
    for constants in ({}, {"sigma_D": 0, "sigma_f": -1, "mu_f": 0},
                      {"sigma_D": np.inf, "sigma_f": 1, "mu_f": np.nan}):
        z = m.standardise_application(c, constants)
        assert z["A1_over_sigma_D"] is z["A2_standardised"] is None


def test_application_missing_and_bad_ledger():
    cv, td = np.ones(1024), np.zeros(1024)
    marks = np.arange(1, 1025, dtype=np.int64) * 128
    for bad in (marks[:-1], marks + 1, marks * 3):
        c = m.application_checkpoint(cv, td, bad, 1024)
        assert c["uphill_count"] is c["proposal_count"] is c["uphill_fraction"] is None
        assert c["statistic_statuses"]["uphill_count"] == "ERROR"
    td[0] = np.nan
    c = m.application_checkpoint(cv, td, marks, 1024)
    assert c["mean_D"] is None and c["finite_input_count"] == 1023
    c = m.application_checkpoint(cv[:-1], np.zeros(1024), marks, 1024)
    assert c["mean_D"] is None and c["finite_input_count"] == 1023


@pytest.fixture(scope="module")
def full_run():
    # Supplied engineering stress profile, NOT a draw from the pilot design.
    # Saturating SHIFT forces the high-consumption path and exercises the cap
    # end to end; sampled pilot-cap profiles need not ever produce a cut.
    p = {"t0": 0, "T_cap": 1.0, "r": 1.0, "s_max": 1.0}
    run = m.dlib.generate_run(SMOKE, "fault/SHIFT", "SHIFT", 0, profile=p)
    record, arrays = m.measure_run(run)
    return run, record, arrays


def test_full_horizon_smoke(full_run):
    run, record, arrays = full_run
    assert record["namespace"] == SMOKE
    assert record["profile_source"] == "supplied"
    assert record["status"] == "SUCCESS", record["failures"]
    assert len(record["conventional"]) == 159 and len(record["application"]) == 157
    assert [c["end"] for c in record["conventional"]] == list(m.BATTERY_ENDS)
    assert [c["end"] for c in record["application"]] == list(m.APPLICATION_ENDS)
    assert set(record["endpoints"]) == set(map(str, m.LENGTHS))
    assert any(c["was_cut"] for c in record["conventional"])
    assert record["measurement_seconds"] > 0 and "generation_seconds" not in record
    assert record["completed_battery_checkpoints"] == 159 and record["completed_application_checkpoints"] == 157
    assert record["split"] == "fitting" and record["source_manifests"]
    assert record["validity_pattern"] == "".join(str(int(e["valid"])) for e in record["endpoints"].values())
    for key, dtype in m.ARRAY_DTYPES.items():
        assert arrays[key].shape == (m.H,) and arrays[key].dtype.str == dtype
        assert np.array_equal(arrays[key], run[key])
    assert all(k not in record for k in (*m.ARRAY_DTYPES, "emitted", "raw_base"))


def test_storage_round_trip_canonical_hashes(tmp_path, full_run):
    _, record, arrays = full_run
    entry = m.write_run(record, arrays, tmp_path)
    loaded, restored = m.load_run(tmp_path, record["condition"], record["replica"], manifest=entry)
    trace = tmp_path / "traces/fault/SHIFT/00000.npz"
    path = tmp_path / "results/runs/fault/SHIFT/00000.json"
    encoded = path.read_bytes()
    assert encoded == (json.dumps(loaded, sort_keys=True, indent=2, ensure_ascii=False, allow_nan=False) + "\n").encode()
    assert b"NaN" not in encoded and b"Infinity" not in encoded
    assert loaded["trace_sha256"] == hashlib.sha256(trace.read_bytes()).hexdigest()
    assert entry["files"][str(path.relative_to(tmp_path))] == hashlib.sha256(encoded).hexdigest()
    for key, dtype in m.ARRAY_DTYPES.items():
        assert restored[key].tobytes() == arrays[key].tobytes()
        assert loaded["array_hashes"][key] == {
            "sha256": hashlib.sha256(np.ascontiguousarray(arrays[key], dtype=dtype).tobytes(order="C")).hexdigest(),
            "dtype": dtype, "shape": [m.H]}
    assert loaded["compression"]["numpy"] == np.__version__
    trace.write_bytes(trace.read_bytes() + b"corruption")
    with pytest.raises(ValueError, match="hash mismatch"):
        m.load_run(tmp_path, record["condition"], record["replica"], manifest=entry)


def test_incomplete_run_keeps_earlier_endpoints(full_run):
    run, _, _ = full_run
    partial = {**run, **{k: run[k][:10240] for k in m.ARRAY_DTYPES}, "current_sweep": 10240}
    partial["local_bytes"] = int(partial["marks"][-1])
    partial["emitted"] = run["emitted"][:partial["local_bytes"]]
    record, _ = m.measure_run(partial)
    assert record["status"] == "FAILED" and record["failure"] is not None
    assert record["endpoints"]["10240"] == m.endpoint(run["CV"], 10240, m.E_M)
    assert record["endpoints"]["20480"]["invalid_reason"] == "INCOMPLETE_TRACE"
    assert m.issuance(record["conventional"], 300)["status"] == "UNKNOWN"


def test_failed_storage_safe_null_and_marks_hash(tmp_path):
    arrays = {k: np.zeros(32, dtype=dtype) for k, dtype in m.ARRAY_DTYPES.items()}
    arrays["marks"][0] = 2**53 + 1  # must never pass through float64
    record = {"namespace": SMOKE, "condition": "healthy", "replica": 2, "status": "FAILED",
              "missing": np.nan, "infinity": np.inf, "failure": {"phase": "test"}}
    entry = m.write_run(record, arrays, tmp_path)
    loaded, restored = m.load_run(tmp_path, "healthy", 2, manifest=entry)
    assert loaded["missing"] is loaded["infinity"] is None
    assert restored["marks"][0] == 2**53 + 1
    with pytest.raises(ValueError, match="dtype"):
        m.write_run(record, {**arrays, "marks": arrays["marks"].astype(float)}, tmp_path)
