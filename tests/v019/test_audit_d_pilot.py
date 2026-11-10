"""Adversarial checks of the v0.19 Stage D profile and fault adapter (research/v019-d-pilot/dlib.py).

Each test states the frozen requirement it attacks (protocol.json key). The
reference for every comparison is the independent code in research/v019-audit/,
which shares no adapter, clock, source class or sampler class with dlib.py.
Runs use the audit-only namespace; nothing here draws from a pilot or smoke stream.
"""
import hashlib
import importlib.util as _ilu
import json
import math
import multiprocessing
import os
import struct
import sys as _sys
from pathlib import Path as _Path

import numpy as np
import pytest

_spec = _ilu.spec_from_file_location("_audit_loader", _Path(__file__).with_name("_audit_loader.py"))
au = _sys.modules.get("_audit_loader") or _ilu.module_from_spec(_spec)
if "_audit_loader" not in _sys.modules:
    _sys.modules["_audit_loader"] = au
    _spec.loader.exec_module(au)

rc, rp, d = au.recount(), au.replay(), au.dlib()
NS = au.AUDIT_NAMESPACE
FAMILIES = ("SHIFT", "LAG1", "DUP")


def hand_profile(t0=3, t_cap=7.5, s_max=0.5):
    return {"t0": t0, "T_cap": np.float64(t_cap), "r": np.float64(s_max / np.float64(t_cap)), "s_max": np.float64(s_max)}


def condition(family):
    return "healthy" if family == "healthy" else "fault/" + family


def assert_same_run(run, mine):
    for key in ("marks", "CV", "TD", "severity"):
        assert run[key].dtype == mine[key].dtype and np.array_equal(run[key], mine[key]), key
    for key in ("emitted_sha256", "base_sha256", "prefix_sha256", "start_sha256", "B0", "local_bytes"):
        assert run[key] == mine[key], key
    assert rc.words_to_bytes(mine["emitted_words"]) == run["emitted"]


# --- application, sweep_clock, fault_transforms: an independent regeneration reproduces every stored quantity
@pytest.mark.parametrize("family", ("healthy",) + FAMILIES)
def test_independent_replay_reproduces_short_runs(family):
    profile = None if family == "healthy" else hand_profile(s_max=d.CAPS[family])
    run = d.generate_run(NS, condition(family), family, 2, sweeps=96, profile=profile)
    assert_same_run(run, rp.replay(NS, condition(family), family, 2, sweeps=96, profile=profile))
    assert rc.check_ledger(run["marks"], 96) == [] and run["I-1"]["passed"]
    assert len(run["emitted"]) == len(run["raw_base"]) == run["local_bytes"] == run["marks"][-1]


@pytest.mark.parametrize("family", FAMILIES)
def test_independent_replay_reproduces_seeded_profile_past_onset(family):
    seeds = rc.run_seeds(NS, condition(family), 0)
    assert seeds == d.seeds_for(NS, condition(family), 0)
    mine_profile = rc.draw_profile(seeds["profile"], rc.protocol_caps()[family])
    theirs = d.draw_profile(family, seeds)
    assert mine_profile["t0"] == theirs["t0"]
    for field in ("T_cap", "r", "s_max"):
        assert struct.pack(">d", mine_profile[field]) == struct.pack(">d", theirs[field])
    sweeps = mine_profile["t0"] + 48
    run = d.generate_run(NS, condition(family), family, 0, sweeps=sweeps)
    assert_same_run(run, rp.replay(NS, condition(family), family, 0, sweeps=sweeps))
    assert run["B0"] == run["marks"][mine_profile["t0"] - 1]


# --- sweep_clock.callback and profile.severity: one severity per sweep, switching exactly at the marks
def test_shift_severity_switches_at_sweep_boundaries():
    profile = hand_profile(t0=4, t_cap=6.0, s_max=0.25)
    run = d.generate_run(NS, "fault/SHIFT", "SHIFT", 3, sweeps=64, profile=profile)
    raw, out = rc.bytes_to_words(run["raw_base"]), rc.bytes_to_words(run["emitted"])
    sweep = rc.sweep_of_word(run["marks"])
    u = (raw >> np.uint64(11)).astype(float) * 2.0**-53
    v = (out >> np.uint64(11)).astype(float) * 2.0**-53
    implied = v - u
    severity = rc.severity_path(4, 6.0, 0.25, 64)
    assert severity[4] == 0 and severity[5] > 0 and severity[10] == 0.25 and severity[9] < 0.25
    unsaturated = u + severity[sweep] < rc.TOP
    assert np.all(np.abs(implied[unsaturated] - severity[sweep][unsaturated]) <= 2.0**-52)
    assert np.all(v[~unsaturated] == rc.TOP) and (~unsaturated).sum() > 0
    assert np.array_equal(out[sweep <= 4], raw[sweep <= 4])   # sweep t0 has severity 0: SHIFT leaves it unchanged
    first_of_5 = int(np.flatnonzero(sweep == 5)[0])
    assert first_of_5 == run["marks"][4] // 8 and out[first_of_5] != raw[first_of_5]
    assert np.all((out & np.uint64(0x7FF)) == (raw & np.uint64(0x7FF)))  # the low 11 bits are never touched
    assert run["B0"] == run["marks"][3]


def test_onset_sweep_consumes_fault_draws_and_keeps_lag1_state():
    profile = hand_profile(t0=4, t_cap=6.0, s_max=0.5)
    run = d.generate_run(NS, "fault/DUP", "DUP", 3, sweeps=40, profile=profile)
    sweep = rc.sweep_of_word(run["marks"])
    raw, out = rc.bytes_to_words(run["raw_base"]), rc.bytes_to_words(run["emitted"])
    assert np.array_equal(out[sweep <= 4], raw[sweep <= 4])  # zero severity cannot duplicate
    post = int((sweep >= 4).sum())
    draws = np.random.default_rng(run["seeds"]["fault"]).random(post)
    duplicate = draws < rc.severity_path(4, 6.0, 0.5, 40)[sweep[sweep >= 4]]
    assert np.array_equal(out[sweep >= 4] != raw[sweep >= 4], duplicate)  # one draw per word from B0, sweep t0 included
    assert duplicate.sum() > 20
    run = d.generate_run(NS, "fault/LAG1", "LAG1", 3, sweeps=40, profile=profile)
    sweep = rc.sweep_of_word(run["marks"])
    raw, out = rc.bytes_to_words(run["raw_base"]), rc.bytes_to_words(run["emitted"])
    assert np.array_equal(out[sweep < 4], raw[sweep < 4])
    at_onset = sweep == 4
    du = ((out[at_onset] >> np.uint64(11)).astype(float) - (raw[at_onset] >> np.uint64(11)).astype(float)) * 2.0**-53
    assert np.all(np.abs(du) <= 2.0**-50) and np.any(du != 0)  # normal-score round trip at zero severity
    vector = rc.apply_fault("LAG1", raw, sweep, rc.severity_path(4, 6.0, 0.5, 40), 4, run["seeds"]["fault"])
    assert np.array_equal(vector, out)


@pytest.mark.parametrize("family", FAMILIES)
def test_vectorised_transform_reproduces_emitted_stream(family):
    profile = hand_profile(t0=2, t_cap=20.0, s_max=d.CAPS[family])
    run = d.generate_run(NS, condition(family), family, 5, sweeps=128, profile=profile)
    base = rc.pcg64_words(run["seeds"]["local"], run["local_bytes"] // 8)
    assert rc.words_to_bytes(base) == run["raw_base"]
    out = rc.apply_fault(family, base, rc.sweep_of_word(run["marks"]),
                         rc.severity_path(2, 20.0, d.CAPS[family], 128), 2, run["seeds"]["fault"])
    assert rc.words_to_bytes(out) == run["emitted"]
    # Moving one word across a sweep boundary must change the stream. For DUP the onset boundary is
    # used (the draw sequence shifts by one); for SHIFT and LAG1 a later boundary (the severity differs).
    moved_marks = run["marks"].copy()
    moved_marks[1 if family == "DUP" else 6] += 8
    moved = rc.apply_fault(family, base, rc.sweep_of_word(moved_marks),
                           rc.severity_path(2, 20.0, d.CAPS[family], 128), 2, run["seeds"]["fault"])
    assert rc.words_to_bytes(moved) != run["emitted"]


# --- fault_transforms.LAG1: the "midpoint" u + 2^-54 is exact only below 1/2
def test_lag1_midpoint_rounding_and_domain_failure_set():
    for k in (1, 2**51 + 1, 2**52 - 1):
        u = k * 2.0**-53
        assert Fraction_exact(u + 2.0**-54) == 2 * k + 1  # an exact midpoint, in units of 2^-54
    for k in (2**52, 2**52 + 2, 2**53 - 2):               # even mantissa at or above 1/2: the offset is rounded away
        assert k * 2.0**-53 + 2.0**-54 == k * 2.0**-53
    for k in (2**52 + 1, 2**53 - 3):                      # odd mantissa: rounded up to the next uniform
        assert k * 2.0**-53 + 2.0**-54 == (k + 1) * 2.0**-53
    assert (2**53 - 1) * 2.0**-53 + 2.0**-54 == 1.0       # the single word value outside the inverse-CDF domain
    with pytest.raises(ValueError):
        d.STD.inv_cdf((2**53 - 1) * 2.0**-53 + 2.0**-54)
    assert d.STD.inv_cdf((2**53 - 2) * 2.0**-53 + 2.0**-54) < 8.3


def Fraction_exact(value):
    from fractions import Fraction
    scaled = Fraction(value) * 2**54
    assert scaled.denominator == 1
    return scaled.numerator


# --- storage.parallelism: results do not depend on process scheduling
def _worker(job):
    family, replica = job
    module = au.dlib()
    profile = None if family == "healthy" else hand_profile(s_max=module.CAPS[family])
    run = module.generate_run(NS, condition(family), family, replica, sweeps=64, profile=profile)
    digest = hashlib.sha256()
    for key in ("CV", "TD", "marks", "severity"):
        digest.update(np.ascontiguousarray(run[key]).tobytes())
    return family, replica, digest.hexdigest(), run["emitted_sha256"], run["start_sha256"]


def test_parallel_execution_is_deterministic():
    jobs = [(family, replica) for family in ("healthy",) + FAMILIES for replica in (10, 11)]
    serial = [_worker(job) for job in jobs]
    with multiprocessing.get_context("spawn").Pool(4) as pool:
        first = pool.map(_worker, jobs, chunksize=1)
        second = pool.map(_worker, list(reversed(jobs)), chunksize=3)
    assert first == serial and list(reversed(second)) == serial
    assert len({row[2] for row in serial}) == len(serial)  # distinct runs give distinct traces


# --- data.seeds: no stream is shared inside the pilot, with the smoke, with the audit, or with Stage C
def test_seed_namespaces_are_disjoint():
    plan = {"calibration": 999, "healthy": 2000, "fault/SHIFT": 1000, "fault/LAG1": 1000, "fault/DUP": 1000}
    roles = ("local", "centroid", "start", "fault", "profile")
    seen = {}
    for namespace, sizes in (("rnglab-v019-d-pilot", plan), ("rnglab-v019-d-pilot-smoke", {k: 2 for k in plan}),
                             (NS, {k: 16 for k in plan})):
        for cond, count in sizes.items():
            for replica in range(count):
                for role in roles:
                    seed = rc.namespace_seed(namespace, cond, role, replica)
                    assert seed not in seen, (namespace, cond, role, replica, seen.get(seed))
                    seen[seed] = (namespace, cond, role, replica)
    assert len(seen) == 5 * (5999 + 10 + 80)
    for cond in ("calibration", "healthy", "SHIFT-5", "LAG1-5", "DUP-5"):
        for replica in range(0, 4096, 37):
            for role in roles + ("gfsr",):
                assert rc.namespace_seed("rnglab-v018-h1-detection", cond, role, replica) not in seen
    for purpose in [f"bootstrap-L{L:05d}" for L in rc.LENGTHS] + ["I4-replay-sample"]:
        for cond in plan:
            assert rc.analysis_seed("rnglab-v019-d-pilot", cond, purpose) not in seen


def test_i4_replay_sample_is_reproducible_and_inside_the_plan():
    """implementation_checks.I-4.sample: five replicas per condition, drawn from the planned IDs only."""
    plan = {"calibration": 999, "healthy": 2000, "fault/SHIFT": 1000, "fault/LAG1": 1000, "fault/DUP": 1000}
    total = 0
    for cond, count in plan.items():
        seed = rc.analysis_seed("rnglab-v019-d-pilot", cond, "I4-replay-sample")
        chosen = np.sort(np.random.default_rng(seed).choice(np.arange(count, dtype=np.int64), size=5, replace=False))
        again = np.sort(np.random.default_rng(seed).choice(np.arange(count, dtype=np.int64), size=5, replace=False))
        assert np.array_equal(chosen, again) and len(set(chosen.tolist())) == 5 and chosen.max() < count
        total += len(chosen)
    assert total == 25


# --- design.json: numbers the frozen protocol quotes, recomputed without the design script
def test_design_ramp_shares_agree_with_the_exact_distribution():
    design = json.loads((au.REPO / "research/v019-d-pilot/design.json").read_text())
    onsets = np.arange(2048, 8193)
    for row in design["ramp_profiles"]["lengths"]:
        # cap reached before L  <=>  t0 + ceil(T_cap) < L  <=>  T_cap <= L - t0 - 1, T_cap log-uniform on [2048, 163840)
        room = np.clip(row["L"] - onsets - 1, 2048, 163840).astype(float)
        exact = float(np.mean(np.log(room / 2048) / math.log(163840 / 2048)))
        se = math.sqrt(exact * (1 - exact) / design["ramp_profiles"]["draws"])
        assert abs(row["share_reaching_cap_before_L"] - exact) < 4 * se
    assert design["size_options"][2]["total_runs"] == 5999 and design["storage"]["bytes_per_run"] == 3 * 8 * 40960


# --- operational_alarm.bit_cap: the cap exists but real consumption stays far below it
def test_full_length_windows_stay_far_below_the_bit_cap():
    worst = 0
    for family in ("healthy", "SHIFT"):
        run = rp.replay(NS, condition(family), family, 0)
        assert rc.check_ledger(run["marks"], rc.H) == []
        bits = np.array([64 * (hi - lo) for lo, hi in (rc.window_words(run["marks"], e - 512, e) for e in rc.BATTERY_ENDS)])
        worst = max(worst, int(bits.max()))
        assert bits.min() >= 750_000  # so the longest-run test always uses 10,000-bit blocks
    assert worst < 900_000  # a cut needs more than 1,000,000 bits, i.e. more than 14.5 of 16 uphill moves per sweep
    assert 8 * 512 * (128 + 8 * 14) < rc.BIT_CAP_BITS < 8 * 512 * (128 + 8 * 15)


# --- measurement layer, when present: cross-checks against the independent code
def synthetic_marks(per_sweep_bytes, sweeps=rc.H):
    return np.cumsum(np.full(sweeps, per_sweep_bytes, dtype=np.int64))


def test_measure_battery_matches_independent_battery_including_cut_windows():
    m = au.measure()
    rng = np.random.default_rng(21)
    for per_sweep in (128, 208, 248, 256):  # 256 bytes per sweep gives 1,048,576-bit windows
        marks = synthetic_marks(per_sweep, 1024)
        words = rng.integers(0, 2**64, size=int(marks[-1]) // 8, dtype=np.uint64)
        emitted = rc.words_to_bytes(words)
        for end in (512, 768, 1024):
            theirs = m.battery_checkpoint(emitted, marks, end)
            lo, hi = rc.window_words(marks, end - 512, end)
            mine = rc.battery_checkpoint(words[lo:hi])
            assert theirs["was_cut"] == mine["was_cut"] == (per_sweep * 512 * 8 > 1_000_000)
            for key in ("original_bits", "bit_used_bits", "word_count", "pair_count", "duplicate_count"):
                assert theirs[key] == mine[key], key
            for a, b in zip(theirs["p_values"], mine["p_values"]):
                assert (a is None) == (b is None) and (a is None or a == pytest.approx(b, rel=1e-9, abs=1e-300))
            for key in ("S", "S_prime", "frequency_signed_z", "runs_signed_z", "serial_delta1", "serial_delta2",
                        "approximate_entropy_chi2", "longest_run_chi2", "lag1_feature", "duplicate_rate"):
                assert rc.compare(theirs[key], mine[key], rtol=1e-9, atol=1e-9), key


@pytest.mark.xfail(strict=True, reason="zero-variance window: np.corrcoef returns r = 1 from a rounding residue, "
                                       "so lag1 is logged as available with p = 0 (signal_schema says unavailable)")
def test_measure_lag1_is_unavailable_for_a_constant_window():
    m = au.measure()
    logged_as_available = []
    for per_sweep in (128, 208):
        marks = synthetic_marks(per_sweep, 512)
        for value in (0.1, 0.3, 0.7, rc.TOP):
            constant = np.full(int(marks[-1]) // 8, np.uint64(int(value * 2**53)) << np.uint64(11))
            assert rc.word_statistics(constant)["lag1_r"] is None
            theirs = m.battery_checkpoint(rc.words_to_bytes(constant), marks, 512)
            if theirs["lag1_r"] is not None or theirs["p_values"][8] is not None:
                logged_as_available.append((per_sweep, value, theirs["lag1_r"], theirs["p_values"][8]))
    assert not logged_as_available, logged_as_available


def test_measure_endpoint_and_issuance_match_independent_recount():
    m = au.measure()
    rng = np.random.default_rng(31)
    reference = rc.reference_energy()
    assert m.E_M == pytest.approx(reference, rel=2e-15)
    for trial in range(6):
        cv = reference + 0.02 * trial + rng.normal(size=rc.H) * 0.5
        for n in (512, 1024, 10240, 20480, 40960):
            theirs, mine = m.endpoint(cv, n, m.E_M), rc.endpoint(cv, n, m.E_M)
            assert theirs["valid"] is mine["valid"] and theirs["covered"] is mine["covered"]
            for key in ("estimate", "signed_error", "mcse", "half_width", "interval_95", "batch_means"):
                assert rc.compare(theirs[key], mine[key]), key
    assert m.endpoint(np.full(64, 2.0), 64, 2.0)["invalid_reason"] == rc.endpoint(np.full(64, 2.0), 64, 2.0)["reason"]
    assert m.endpoint(np.r_[np.ones(63), np.nan], 64, 1.0)["invalid_reason"] == "NONFINITE_INPUT"
    assert m.endpoint(np.ones(48), 64, 1.0)["invalid_reason"] == "INCOMPLETE_TRACE"
    statistic = rng.exponential(size=159)
    checkpoints = [{"end": int(e), "S": float(s), "S_prime": float(s) / 2, "complete": True, "complete_prime": True}
                   for e, s in zip(rc.BATTERY_ENDS, statistic)]
    for tau in (0.5, 2.0, float(statistic.max()), 50.0):
        assert m.issuance(checkpoints, tau)["alarm_sweep"] == rc.first_alarm(statistic, tau)
        assert m.issuance(checkpoints, tau, prime=True)["alarm_sweep"] == rc.first_alarm(statistic / 2, tau)
    for alarm, t0 in ((None, 3000), (3072, 3072), (3328, 3072), (512, 8192)):
        assert m.onset_classification(alarm, t0)["onset_classification"] == rc.onset_class(alarm, t0)
    for alarm in (None, 10239, 10240, 20479, 20480):
        assert m.eligibility(alarm, 20480)["eligible"] == rc.eligible(alarm, 20480)
    assert m.issuance_endpoint(cv, 1023)["n"] == rc.issuance_length(1023) == 1008


def test_measure_storage_round_trip_is_exact(tmp_path):
    m = au.measure()
    run = m.dlib.generate_run(NS, "fault/LAG1", "LAG1", 7, sweeps=rc.H) if os.environ.get("RNGLAB_AUDIT_FULL") else None
    if run is None:
        pytest.skip("full-length storage round trip runs only with RNGLAB_AUDIT_FULL=1 (about 30 s)")
    record, arrays = m.measure_run(run)
    entry = m.write_run(record, arrays, tmp_path)
    loaded, loaded_arrays = m.load_run(tmp_path, "fault/LAG1", 7, manifest=entry)
    for key, value in arrays.items():
        assert loaded_arrays[key].dtype == value.dtype and loaded_arrays[key].tobytes() == value.tobytes()
    assert json.loads(json.dumps(record)) == {k: loaded[k] for k in json.loads(json.dumps(record))}
    assert_same_run(run, rp.replay(NS, "fault/LAG1", "LAG1", 7))


# --- Stage C regression: the independent functions reproduce the stored v0.18 gates
def test_independent_recount_reproduces_stage_c_gates():
    root = au.REPO / "research/v018-h1-detection/results"
    gates = json.loads((root / "gates.json").read_text())

    def maxima(name):
        p = np.load(root / "cells" / f"{name}.npz")["conv_p"]
        return np.array([max(rc.alarm_statistic(row)[0] for row in run) for run in p])

    tau = rc.order_statistic_threshold(maxima("calibration"))
    assert tau == gates["thresholds"]["C"]
    alarms = int(np.sum(maxima("healthy") > tau))
    assert alarms == gates["validity"]["V2"]["rate_C"]["alarms"] == 199
    assert rc.clopper_pearson(alarms, 4096) == pytest.approx(gates["validity"]["V2"]["rate_C"]["cp95"], abs=1e-12)
    rows = [json.loads(line) for line in (root / "cells/healthy.jsonl").open()]
    covered = sum(r["harm"]["covered"] for r in rows)
    assert covered == gates["validity"]["V1"]["covered"] == 3867
    assert rc.clopper_pearson(covered, 4096) == pytest.approx(gates["validity"]["V1"]["cp95"], abs=1e-12)
    for cell in ("SHIFT-5", "LAG1-5", "DUP-5"):
        assert int(np.sum(maxima(cell) > tau)) == gates["cells"][cell]["C_alarm"]["k"]
        assert gates["cells"][cell]["severity"] == rc.protocol_caps()[cell.split("-")[0]]


def test_independent_endpoint_reproduces_stage_c_harm_intervals_from_traces():
    override = os.environ.get("RNGLAB_STAGEC_TRACES")
    traces = _Path(override) if override else au.REPO / "research/v018-h1-detection/traces"
    block = traces / "healthy" / "block-0000.npz"
    if not block.is_file():
        pytest.skip("Stage C traces are not in this tree (set RNGLAB_STAGEC_TRACES to their directory)")
    data = np.load(block, allow_pickle=False)
    rows = {r["replica"]: r for r in (json.loads(line) for line in
                                      (au.REPO / "research/v018-h1-detection/results/cells/healthy.jsonl").open())}
    reference = rc.reference_energy()
    for replica, cv in zip(data["replica"].tolist(), data["CV"]):
        mine = rc.endpoint(cv[2560:10240], 7680, reference)  # sixteen batches of 480, the Stage C harm window
        stored = rows[replica]["harm"]
        assert mine["covered"] is stored["covered"]
        assert mine["estimate"] == pytest.approx(stored["mean"], rel=1e-12)
        assert mine["mcse"] == pytest.approx(stored["mcse"], rel=1e-12)
        assert mine["interval_95"] == pytest.approx(stored["interval"], rel=1e-12)
