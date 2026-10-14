from hashlib import sha256
import pytest
from rnglab12.review import assess_liveness, evidence_receipt, verify_receipt


def test_saved_running_status_never_establishes_liveness():
    assert assess_liveness(process_observed=False, heartbeat_at=99, checked_at=100, saved_status="RUNNING").status == "STOPPED"
    assert assess_liveness(process_observed=None, heartbeat_at=99, checked_at=100, saved_status="ACTIVE").status == "UNKNOWN"
    assert assess_liveness(process_observed=True, heartbeat_at=1, checked_at=100).status == "PROCESS_OBSERVED_NO_FRESH_HEARTBEAT"
    assert assess_liveness(process_observed=True, heartbeat_at=101, checked_at=100).status != "LIVE"
    assert assess_liveness(process_observed=True, heartbeat_at=99, checked_at=100).status == "LIVE"


def test_completed_review_requires_exact_files_and_inspection(tmp_path):
    path = tmp_path / "evidence.json"
    path.write_bytes(b'{"result": 2}')
    requested = {str(path): sha256(path.read_bytes()).hexdigest()}
    seen = []
    def inspect(path, raw):
        seen.append(raw)
        assert b'"result": 2' in raw
        return {"status": "PASS", "checks": ["Recorded result equals independently derived value 2"]}
    receipt = evidence_receipt(requested, reviewer="deterministic-check", inspector=inspect, checked_at=100)
    assert seen == [path.read_bytes()]
    assert verify_receipt(receipt, requested)
    assert not verify_receipt(receipt, {})
    receipt["files"][0]["inspection"]["checks"] = ["changed"]
    assert not verify_receipt(receipt, requested)
    path.write_bytes(b'{"result": 3}')
    with pytest.raises(ValueError, match="changed"):
        evidence_receipt(requested, reviewer="test", inspector=inspect, checked_at=100)


def test_skipped_or_unchecked_evidence_cannot_complete(tmp_path):
    path = tmp_path / "evidence.txt"
    path.write_text("bytes")
    requested = {str(path): sha256(path.read_bytes()).hexdigest()}
    for result in ({}, {"status": "PASS", "checks": []}, {"status": "PASS", "checks": [""]}):
        with pytest.raises(ValueError):
            evidence_receipt(requested, reviewer="test", inspector=lambda p, b: result, checked_at=100)
    with pytest.raises(ValueError):
        evidence_receipt({}, reviewer="test", inspector=lambda p, b: {}, checked_at=100)
