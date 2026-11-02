import hashlib
import pytest

from rnglab14.sources.base import SourceBudgetExceeded
from rnglab_crypto.sources import CapturedSource, ReplaySource, OSReference, seed_for
from rnglab_crypto.sources.physical import audit_waveforms


def test_replay_and_audited_consumption():
    seed = seed_for("test", 0, "reference")
    a, b = ReplaySource(seed), ReplaySource(seed)
    stream = b""
    for n in (0, 3, 32, 19):
        block = a.read_bytes(n)
        assert block == b.read_bytes(n)
        stream += block
    assert a.consumed_bytes == 54
    assert a.stream_digest == hashlib.sha256(stream).hexdigest()
    assert a.manifest.metadata["entropy_credit_bits"] == 0
    assert a.manifest.to_dict() == b.manifest.to_dict()
    assert seed != seed_for("test/calibration", 0, "reference")


def test_capture_fails_closed_and_never_loops():
    src = CapturedSource(b"1234", provenance="synthetic fixture")
    assert src.read_bytes(4) == b"1234"
    with pytest.raises(SourceBudgetExceeded):
        src.read_bytes(1)
    assert src.consumed_bytes == 4


@pytest.mark.parametrize("n", [-1, True, 1.5])
def test_bad_read(n):
    with pytest.raises(ValueError):
        ReplaySource(bytes(32)).read_bytes(n)


def test_os_adapter_calls_os_only_when_requested(monkeypatch):
    monkeypatch.setattr("rnglab_crypto.sources.os.urandom", lambda n: b"z"*n)
    assert OSReference().read_bytes(10) == b"z"*10


def test_physical_audit_retains_quantization_and_saturation(tmp_path):
    p = tmp_path/"synthetic.txt"
    p.write_text("0 80 1 2 3 65535\n0 160 0 2 3 4\n0 320 1 2 3 4\n")
    before = p.read_bytes()
    result = audit_waveforms(p, {"samples_per_record": 4, "provenance": "synthetic fixture"})
    assert result["timing"]["positive_gap_gcd_raw_units"] == 80
    assert result["records_with_ceiling_fraction"] == 1/3
    assert result["ceiling_fraction"] == 1/12
    assert result["effective_resolution_bits"] is None
    assert result["validated_entropy_credit_bits"] == 0
    assert p.read_bytes() == before
    with pytest.raises(ValueError):
        audit_waveforms(p, {"samples_per_record": 3, "provenance": "bad schema"})
