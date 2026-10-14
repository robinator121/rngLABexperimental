"""Bounded mechanics tests; temporary bytes are explicitly synthetic fixtures."""
import hashlib

import numpy as np
import pytest

from rnglab11.sources import BytesSource, SourceBudgetExceeded
from rnglab12.faults import apply_uniform_fault, make_byte_fault
from rnglab12.sources import (SOURCE_FAMILIES, make_source, source_type,
                             validate_physical_manifest, verify_content_file)


def test_source_inventory_and_classification():
    assert len(SOURCE_FAMILIES) == len(set(SOURCE_FAMILIES)) == 7
    with pytest.raises(ValueError):
        source_type({"family": "pcg64_csprng"})
    with pytest.raises(ValueError):
        source_type({"family": "raw_physical", "source_type": "standalone_csprng"})
    with pytest.raises(ValueError):
        make_source({"family": "standalone_csprng", "seed": 12})


def test_full_content_integrity_and_tamper(tmp_path):
    # These bytes are a synthetic unit-test fixture, never a physical dataset.
    path = tmp_path / "SYNTHETIC_FIXTURE.bin"
    data = bytes(range(256))
    path.write_bytes(data)
    digest = hashlib.sha256(data).hexdigest()
    verified, identity = verify_content_file(path, digest)
    assert verified == path and identity[0] == 256
    path.write_bytes(data[:-1] + b"x")
    with pytest.raises(ValueError, match="does not match"):
        verify_content_file(path, digest)
    with pytest.raises(ValueError, match="budget"):
        verify_content_file(path, hashlib.sha256(path.read_bytes()).hexdigest(), max_input_bytes=255)
    with pytest.raises(ValueError):
        verify_content_file(path, "0" * 64)


def test_physical_placeholders_and_synthetic_fixtures_rejected(tmp_path):
    path = tmp_path / "SYNTHETIC_FIXTURE.bin"
    path.write_bytes(bytes(range(64)))
    spec = {"family": "raw_physical", "path": str(path),
            "content_sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
    with pytest.raises(ValueError, match="acquisition"):
        make_source(spec)
    with pytest.raises(ValueError, match="fixtures"):
        make_source({**spec, "test_fixture": True})
    with pytest.raises(ValueError, match="nonplaceholder"):
        validate_physical_manifest({**spec, "acquisition": {
            "instrument": "TODO", "time": "2026-10-16T00:00:00Z", "units": "bytes", "origin": "unknown"}})
    # No positive physical ingestion test invents a genuine acquisition record.


@pytest.mark.parametrize("interval", [None, 1, 3])
def test_deterministic_hmac_replay_chunk_independence_and_reseed(interval):
    spec = {"family": "deterministic_seeded_or_reseeded_csprng", "seed": 37,
            "request_bytes": 64, "reseed_interval_generate_calls": interval}
    whole = make_source(spec, max_bytes=512)
    split = make_source(spec, max_bytes=512)
    expected = whole.read(512)
    actual = b"".join(split.read(n) for n in (1, 17, 80, 127, 287))
    assert actual == expected
    prov = split.provenance()["generator"]
    assert prov["internal_generate_calls"] == 8
    assert prov["buffered_generated_bytes"] == 0
    events = prov["generator"]["reseed_events"]
    assert events == (0 if interval is None else 7 // interval)
    assert split.provenance()["family"] == spec["family"]
    assert "not_assessed_entropy" in prov["metadata"]["initial_seed_label"]
    with pytest.raises(SourceBudgetExceeded):
        split.read(1)
    if interval is not None:
        no_reseed = make_source({**spec, "reseed_interval_generate_calls": None}, max_bytes=512)
        assert expected != no_reseed.read(512)


def test_reseed_boundary_preserves_pending_output():
    spec = {"family": "deterministic_seeded_or_reseeded_csprng", "seed": 91,
            "request_bytes": 64, "reseed_interval_generate_calls": 1}
    source = make_source(spec, max_bytes=129)
    source.read(63)
    assert source.provenance()["generator"]["generator"]["reseed_events"] == 0
    source.read(1)
    assert source.provenance()["generator"]["generator"]["reseed_events"] == 0
    source.read(1)
    assert source.provenance()["generator"]["generator"]["reseed_events"] == 1


@pytest.mark.parametrize("mode", ["bias", "autocorrelation", "copy_bit", "stuck_bit", "periodicity",
                                 "zero_burst", "drift", "entropy_reduction", "bit_flip"])
def test_byte_fault_chunk_invariance_and_finite_interval(mode):
    native = bytes(range(256)) * 2
    spec = {"mode": mode, "start_byte": 31, "duration_bytes": 290,
            "severity": .4, "final_severity": .8, "seed": 4}
    first = make_byte_fault(BytesSource(native), spec, max_bytes=512)
    second = make_byte_fault(BytesSource(native), spec, max_bytes=512)
    whole = first.read(512)
    split = b"".join(second.read(n) for n in (3, 17, 171, 321))
    assert whole == split
    assert whole[:31] == native[:31] and whole[321:] == native[321:]


@pytest.mark.parametrize("mode", ["bias", "autocorrelation", "copy_bit", "stuck_bit", "periodicity",
                                 "zero_burst", "drift", "entropy_reduction", "bit_flip"])
def test_uniform_fault_chunk_invariance_scope_and_boundaries(mode):
    native = np.random.default_rng(128).random((96, 3))
    spec = {"mode": mode, "scope": "acceptance", "start_row": 11, "duration_rows": 53,
            "severity": .5, "final_severity": .8, "seed": 33}
    whole = apply_uniform_fault(native, 0, spec, "harmonic_pimc")
    split = np.vstack([apply_uniform_fault(native[a:b], a, spec, "harmonic_pimc")
                       for a, b in ((0, 7), (7, 30), (30, 53), (53, 96))])
    np.testing.assert_array_equal(whole, split)
    np.testing.assert_array_equal(whole[:, :2], native[:, :2])
    np.testing.assert_array_equal(whole[:11], native[:11])
    np.testing.assert_array_equal(whole[64:], native[64:])
    assert not np.shares_memory(whole, native)


def test_empty_uniform_chunk_and_combination_order():
    empty = np.zeros((0, 2))
    assert apply_uniform_fault(empty, 0, {"mode": "periodicity"}, "mcmc").shape == (0, 2)
    data = np.full((9, 2), .7)
    faults = [{"mode": "zero_burst", "scope": "acceptance", "start_row": 2, "duration_rows": 3},
              {"mode": "entropy_reduction", "keep_bits": 2, "scope": "proposal"}]
    expected = apply_uniform_fault(apply_uniform_fault(data, 0, faults[0], "vmc"), 0, faults[1], "vmc")
    np.testing.assert_array_equal(apply_uniform_fault(data, 0, {"faults": faults}, "vmc"), expected)


def test_uniform_invalid_scope_unit_and_range_rejected():
    native = np.full((4, 2), .5)
    for spec in ({"mode": "zero_burst", "scope": "selection"},
                 {"mode": "zero_burst", "unit": "bytes"},
                 {"mode": "zero_burst", "start_byte": 2},
                 {"mode": "drift"}, {"mode": "entropy_reduction", "keep_bits": 54}):
        with pytest.raises(ValueError):
            apply_uniform_fault(native, 0, spec, "vmc")
    with pytest.raises(ValueError):
        apply_uniform_fault(np.ones((4, 2)), 0, {"mode": "zero_burst"}, "vmc")
    with pytest.raises(ValueError):
        apply_uniform_fault(native, 0, {"mode": "zero_burst"}, "unregistered")


def test_native_period_pattern_and_zero_acceptance_placement():
    native = np.random.default_rng(41).random((12, 2))
    pattern = native[3:5, 1].tolist()
    output = apply_uniform_fault(native, 0, {"mode": "periodicity", "scope": "acceptance",
        "start_row": 3, "duration_rows": 6, "period": 2, "pattern": pattern}, "mcmc")
    np.testing.assert_array_equal(output[3:9, 1], np.tile(pattern, 3))
    zero = apply_uniform_fault(native, 0, {"mode": "zero_burst", "scope": "acceptance"}, "mcmc")
    np.testing.assert_array_equal(zero[:, 0], native[:, 0])
    assert np.all(zero[:, 1] == 0)


def test_synthetic_fault_family_never_claims_entropy():
    source = make_source({"family": "synthetic_fault", "seed": 6,
        "fault": {"mode": "entropy_reduction", "keep_bits": 0}}, max_bytes=128)
    assert source.read(128) == bytes(128)
    p = source.provenance()
    assert p["family"] == "synthetic_fault"
    assert p["metadata"]["synthetic_intervention_not_measured_entropy_degradation"] is True


def test_byte_combination_and_zero_severity():
    native = bytes(range(64))
    no_op = make_byte_fault(BytesSource(native), {"mode": "copy_bit", "severity": 0.}, max_bytes=64)
    assert no_op.read(64) == native
    spec = {"faults": [{"mode": "stuck_bit", "bit": 0, "value": 1},
                        {"mode": "entropy_reduction", "keep_bits": 0,
                         "start_byte": 20, "duration_bytes": 11}]}
    combined = make_byte_fault(BytesSource(native), spec, max_bytes=64)
    output = combined.read(64)
    assert output[20:31] == bytes(11)
    assert all(value & 1 for value in output[:20] + output[31:])
    assert combined.provenance()["metadata"]["position_unit"] == "absolute_zero_based_bytes"
