"""C2 engineering checks only: audit-only seeds and hand-built profiles, no campaign output.

PROTOCOL.md sections 2.2 and 7; protocol.json implementation_checks I-1–I-3.
Stage C modules are imported unmodified, as in tests/v018/test_h1_lib.py.
"""
import hashlib
import importlib.util
import json
import math
import struct
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from rnglab.seeds import namespace_seed
from rnglab14.sources import PCG64Source
from rnglab14.sources.base import RandomSource, SourceManifest

REPO = Path(__file__).resolve().parents[2]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


d = load("v019_dlib", REPO / "research/v019-d-pilot/dlib.py")
h1 = load("v019_test_h1", REPO / "research/v018-h1-detection/h1lib.py")
c = load("v019_test_campaign", REPO / "research/v018-h1-detection/run_campaign.py")
NAMESPACE = "rnglab-v019-audit-smoke"   # audit engineering runs only; never a pilot or smoke namespace
FAMILIES = ("SHIFT", "LAG1", "DUP")


def profile(t0=3, T_cap=7.5, s_max=0.5):
    return {"t0": t0, "T_cap": np.float64(T_cap),
            "r": np.float64(s_max / np.float64(T_cap)), "s_max": np.float64(s_max)}


def drive(family, p, seeds, marks, *, step=False):
    local = d.FaultyWords(PCG64Source(seeds["local"]), family, p, seeds["fault"], step=step)
    clock = d.SweepClock(PCG64Source(seeds["centroid"]), local)
    for boundary in marks:
        while local.consumed_bytes < boundary:
            local.read_bytes(8)
        clock.read_bytes(16)
    assert clock.marks == list(marks)
    assert local.current_sweep == len(marks)
    return local


def record(local, seeds, marks):
    return {"family": local.family, "profile": local.profile, "seeds": seeds,
            "marks": marks, "sweeps": len(marks), "severity": d.severity_path(local.profile, len(marks)),
            "local_bytes": local.consumed_bytes, "B0": local.B0,
            "emitted": bytes(local.emitted), "raw_base": bytes(local.raw_base),
            "emitted_sha256": local.stream_digest, "base_sha256": local.base.stream_digest,
            "prefix_sha256": local.prefix_sha256}


@pytest.mark.parametrize("family", FAMILIES)
def test_profile_draws_exact_and_seed_roles(family):
    frozen = json.loads((REPO / "research/v019-d-pilot/protocol.json").read_text())
    gates = json.loads((REPO / "research/v018-h1-detection/results/gates.json").read_text())
    assert d.CAPS[family] == frozen["profile"]["caps"][family]["s_max"] == gates["cells"][family + "-5"]["severity"]
    for i in range(12):
        seeds = d.seeds_for(NAMESPACE, "fault/" + family, i)
        assert set(seeds) == {"local", "centroid", "start", "fault", "profile"}
        assert all(seeds[k] == namespace_seed(NAMESPACE, "fault/" + family, k, f"{i:05d}") for k in seeds)
        assert all(seeds[k] == h1.seeds_for(NAMESPACE, "fault/" + family, i)[k] for k in seeds if k != "profile")
        rng = np.random.default_rng(seeds["profile"])
        t0 = int(rng.integers(2048, 8193))
        T_cap = np.float64(np.exp(rng.uniform(np.log(2048.0), np.log(163840.0))))
        p = d.draw_profile(family, seeds)
        assert p == d.draw_profile(family, seeds)
        assert p["t0"] == t0 and 2048 <= t0 <= 8192
        assert 2048 <= p["T_cap"] <= 163840
        for field, expected in (("T_cap", T_cap), ("r", np.float64(d.CAPS[family] / T_cap)),
                                ("s_max", np.float64(d.CAPS[family]))):
            assert isinstance(p[field], np.float64)
            assert p[field].tobytes() == expected.tobytes()


def test_healthy_profile_draws_nothing(monkeypatch):
    def forbidden(*args, **kwargs):
        pytest.fail("healthy requested profile randomness")
    # Replace only the pilot module's binding, never NumPy or Stage C objects.
    monkeypatch.setattr(d, "np", SimpleNamespace(random=SimpleNamespace(default_rng=forbidden)))
    assert d.draw_profile("healthy", {}) == dict.fromkeys(d.PROFILE_FIELDS)


@pytest.mark.parametrize("p", [profile(0, 1), profile(3, 7.5), profile(2, 9),
                               profile(5, 10000, d.CAPS["SHIFT"]), profile(80, 2)])
def test_severity_formula_and_oracle_sum(p):
    path = d.severity_path(p, 128)
    expected = np.array([np.float64(0) if t < p["t0"] else
                         np.minimum(np.float64(p["s_max"]), np.float64(p["r"] * np.float64(t - p["t0"])))
                         for t in range(128)], dtype=np.float64)
    assert path.dtype == np.float64 and path.tobytes() == expected.tobytes()
    assert path[p["t0"]] == 0 and path[p["t0"] + 1] > 0
    assert np.flatnonzero(path > 0)[0] == p["t0"] + 1
    cap_sweep = p["t0"] + math.ceil(p["T_cap"])
    if cap_sweep < len(path):
        assert np.all(path[cap_sweep:] == p["s_max"])
    for L in (1, 3, 10, 64, 128):
        m = max(L - 1 - p["t0"], 0)
        q = min(m, math.floor(p["T_cap"]))
        closed = (p["r"] * q * (q + 1) / 2 + p["s_max"] * (m - q)) / L
        assert np.sum(path[:L], dtype=np.float64) / L == pytest.approx(closed, rel=1e-14, abs=1e-16)


def test_healthy_identity_and_read_contract():
    seeds = d.seeds_for(NAMESPACE, "healthy", 0)
    local = d.FaultyWords(PCG64Source(seeds["local"]), "healthy", d.draw_profile("healthy", seeds), seeds["fault"])
    clock = d.SweepClock(PCG64Source(seeds["centroid"]), local)
    base = PCG64Source(seeds["local"])
    for _ in range(64):
        assert local.read_bytes(8) == base.read_bytes(8)
    assert local.raw_base == local.emitted
    assert local.B0 is None and local.prefix_sha256 is None and local.current_sweep == 0
    assert local.fault_rng.random() == np.random.default_rng(seeds["fault"]).random()
    for n in (0, 1, 7, 16, -1, True, 8.0):
        with pytest.raises(ValueError, match="8-byte"):
            local.read_bytes(n)
    for n in (0, 8, 32, True, 16.0):
        with pytest.raises(ValueError, match="16 bytes"):
            clock.read_bytes(n)
    assert local.consumed_bytes == 512 and clock.marks == []
    clock.read_bytes(16)
    assert clock.marks == [512] and local.current_sweep == 1


@pytest.fixture(scope="module")
def healthy_run():
    return d.generate_run(NAMESPACE, "healthy", "healthy", 0, sweeps=64)


def test_real_healthy_run(healthy_run):
    run = healthy_run
    assert run["I-1"]["passed"] and run["I-1"]["counter_checked"]
    assert run["current_sweep"] == 64 and run["marks"].dtype == np.int64
    k = (np.diff(np.r_[0, run["marks"]]) - 128) / 8
    assert np.all((k >= 0) & (k <= 16) & (k == np.floor(k)))
    assert run["emitted"] == run["raw_base"] and run["base_sha256"] == run["emitted_sha256"]
    assert run["profile"] == dict.fromkeys(d.PROFILE_FIELDS)
    assert not run["severity"].any() and run["severity"].dtype == np.float64
    assert run["CV"].shape == run["TD"].shape == (64,)
    assert np.isfinite(run["CV"]).all() and np.isfinite(run["TD"]).all()
    assert d.check_i2(run)["passed"] and d.check_i3(run)["passed"]
    # Same exact start and sampler contract as frozen Stage C, without its driver.
    seeds = run["seeds"]
    start = h1.sample_exact_finite_bead_paths(PCG64Source(seeds["start"]), 16, 2.0, 2)[0]
    reference = h1.PathIntegralHarmonicOscillator(beads=16, temperature=0.5, displacement_delta=0.4).run_simulation(
        PCG64Source(seeds["local"]), sweeps=64, burn_in=0, initial_path=start,
        centroid_source=PCG64Source(seeds["centroid"]),
        initial_path_provenance={"namespace": NAMESPACE, "condition": "healthy", "role": "start",
                                 "replica": 0, "seed": seeds["start"]},
        uncertainty_batch_sizes=(512,), record_trace=True)
    assert run["start_sha256"] == hashlib.sha256(start.astype("<f8").tobytes()).hexdigest()
    for key in ("CV", "TD"):
        assert np.array_equal(run[key], reference.trace[key])


@pytest.mark.parametrize("family", FAMILIES)
def test_step_equivalence_all_four_paths(family, healthy_run):
    seeds = d.seeds_for(NAMESPACE, "fault/" + family, 0)
    p, marks = profile(s_max=d.CAPS[family]), healthy_run["marks"]
    local = drive(family, p, seeds, marks, step=True)
    B0, n = int(marks[p["t0"] - 1]), int(marks[-1])
    stage_c = h1.FaultyWords(PCG64Source(seeds["local"]), family, p["s_max"], seeds["fault"], 0, onset_bytes=B0)
    reference = b"".join(stage_c.read_bytes(8) for _ in range(n // 8))
    assert bytes(local.emitted) == reference
    assert d.reapply(family, p, seeds, marks, n, step=True) == reference
    assert c.expected_stream(family, p["s_max"], {**seeds, "gfsr": 0}, n, onset=B0) == reference
    assert d.check_i3(record(local, seeds, marks), step=True)["passed"]


@pytest.fixture(scope="module", params=FAMILIES)
def ramp_run(request):
    family = request.param
    p = profile(s_max=d.CAPS[family])
    return d.generate_run(NAMESPACE, "fault/" + family, family, 0, sweeps=64, profile=p)


def test_real_hand_built_ramp(ramp_run):
    run = ramp_run
    assert run["I-1"]["passed"]
    checked = d.check_i2(run, expected_profile=run["profile"])
    assert checked["passed"] and checked["profile_check"] == "supplied"
    assert not d.check_i2(run)["passed"]  # hand profiles cannot pass a seeded redraw
    assert d.check_i3(run)["passed"]
    assert run["emitted"] == d.reapply(run["family"], run["profile"], run["seeds"], run["marks"], run["local_bytes"])
    B0 = run["B0"]
    assert B0 == run["marks"][2]
    assert run["emitted"][:B0] == run["raw_base"][:B0]
    assert run["prefix_sha256"] == hashlib.sha256(run["raw_base"][:B0]).hexdigest()
    assert run["emitted"][B0:] != run["raw_base"][B0:]


def test_real_draw_run_past_onset():
    seeds = d.seeds_for(NAMESPACE, "fault/SHIFT", 0)
    p = d.draw_profile("SHIFT", seeds)
    run = d.generate_run(NAMESPACE, "fault/SHIFT", "SHIFT", 0, sweeps=p["t0"] + 16)
    assert run["profile_source"] == "seeded" and run["I-1"]["passed"]
    assert d.check_i2(run)["passed"] and d.check_i3(run)["passed"]
    assert run["severity"][p["t0"]] == 0 and run["severity"][p["t0"] + 1] > 0


def test_dup_draw_schedule_and_emitted_previous():
    seeds = d.seeds_for(NAMESPACE, "fault/DUP", 0)
    p = profile(t0=1, T_cap=2, s_max=0.9)
    marks = np.arange(1, 6, dtype=np.int64) * 128
    local = drive("DUP", p, seeds, marks)
    rng = np.random.default_rng(seeds["fault"])
    out, prev, repeated = bytearray(), None, 0
    for b in range(0, int(marks[-1]), 8):
        raw = bytes(local.raw_base[b:b + 8])
        t = int(np.searchsorted(marks, b, side="right"))
        if t >= p["t0"]:
            q = rng.random()  # includes all sixteen zero-severity onset words
            if q < min(p["s_max"], p["r"] * (t - p["t0"])) and prev is not None:
                repeated += prev != bytes(local.raw_base[b - 8:b])
                raw = prev
        out += raw
        prev = raw
    assert repeated > 0 and out == local.emitted
    assert local.fault_rng.random() == rng.random()  # exact number of draws
    assert local.base.consumed_bytes == local.consumed_bytes == marks[-1]


def test_lag1_midpoint_and_state_survive_severity_change():
    seeds = d.seeds_for(NAMESPACE, "fault/LAG1", 0)
    p = profile(t0=0, T_cap=2)
    marks = np.arange(1, 4, dtype=np.int64) * 128
    local = drive("LAG1", p, seeds, marks)
    z, expected = None, bytearray()
    for b in range(0, int(marks[-1]), 8):
        w = struct.unpack(">Q", local.raw_base[b:b + 8])[0]
        e = h1.STD.inv_cdf((w >> 11) * 2.0**-53 + 2.0**-54)
        rho = min(p["s_max"], p["r"] * (b // 128))
        old = z
        z = e if z is None else rho * z + math.sqrt(1 - rho * rho) * e
        if b == 128:
            assert old is not None and z != e  # neither reset nor skipped zero-severity state
        expected += struct.pack(">Q", (int(min(h1.STD.cdf(z), h1.TOP) * 2.0**53) << 11) | (w & 0x7FF))
    assert local.emitted == expected and local.z == z and local.B0 == 0
    assert local.prefix_sha256 == hashlib.sha256(b"").hexdigest()
    assert local.fault_rng.random() == np.random.default_rng(seeds["fault"]).random()


class ConstantWord(RandomSource):
    def __init__(self, word):
        super().__init__(SourceManifest("constant", "test", True, False))
        self.word = word

    def _read_bytes(self, n):
        assert n == 8
        return struct.pack(">Q", self.word)


def test_lag1_zero_word_and_domain_failure_are_stage_c_rules():
    for w in (0, 2**64 - 1):
        local = d.FaultyWords(ConstantWord(w), "LAG1", profile(t0=0), 0)
        stage_c = h1.FaultyWords(ConstantWord(w), "LAG1", 0.0, 0, 0, onset_bytes=0)
        if w == 0:
            assert local.read_bytes(8) == stage_c.read_bytes(8)
            assert local.z == h1.STD.inv_cdf(2.0**-54)
        else:
            with pytest.raises(ValueError):
                local.read_bytes(8)
            with pytest.raises(ValueError):
                stage_c.read_bytes(8)


def test_i1_failed_runs_and_mismatch_details():
    marks = np.array([128, 264, 520], dtype=np.int64)
    assert d.check_i1(marks, 520, 3, False, current_sweep=3)["passed"]
    assert d.check_i1(marks, 528, 10, True, current_sweep=3)["passed"]
    assert d.check_i1(np.array([], dtype=np.int64), 8, 10, True)["passed"]
    for bad_marks, size, horizon, counter in ((marks, 520, 4, 3), (marks, 528, 3, 3),
                                              (marks.astype(float), 520, 3, 3), (marks, 520, 3, 2),
                                              (np.array([128, 128, 520]), 520, 3, 3),
                                              (np.array([128, 265, 520]), 520, 3, 3)):
        checked = d.check_i1(bad_marks, size, horizon, False, current_sweep=counter)
        assert not checked["passed"] and checked["mismatches"]


def test_i2_i3_detect_corruption(ramp_run):
    run = ramp_run
    for field in ("T_cap", "r", "s_max"):
        p = {**run["profile"], field: np.nextafter(run["profile"][field], np.inf)}
        checked = d.check_i2({**run, "profile": p}, expected_profile=run["profile"])
        assert not checked["passed"] and checked["mismatches"]
    for field, value in (("B0", run["B0"] + 8), ("prefix_sha256", "wrong"),
                         ("severity", run["severity"] + 1e-15)):
        assert not d.check_i2({**run, field: value}, expected_profile=run["profile"])["passed"]
    for field in ("emitted", "raw_base"):
        bad = bytearray(run[field])
        bad[0] ^= 1
        checked = d.check_i3({**run, field: bytes(bad)})
        assert not checked["passed"] and checked["mismatches"][0]["first_byte"] == 0
    for field in ("base_sha256", "emitted_sha256", "prefix_sha256"):
        assert not d.check_i3({**run, field: "wrong"})["passed"]


def test_shortened_run_before_onset():
    p = profile(t0=64)
    run = d.generate_run(NAMESPACE, "fault/DUP", "DUP", 1, sweeps=64, profile=p)
    assert run["B0"] is None and run["prefix_sha256"] is None
    assert run["emitted"] == run["raw_base"]
    checked = d.check_i2(run, expected_profile=p)
    assert checked["passed"] and not checked["onset_observed"]
    assert d.check_i3(run)["passed"]
