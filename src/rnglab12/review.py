"""Honest observed liveness and byte-specific deterministic review receipts.

This module never invokes a provider, starts services, or equates a prior
successful call with a currently running service. A completed receipt records
inspection, not scientific endorsement, model consensus, or confirmation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Callable, Mapping


@dataclass(frozen=True)
class Liveness:
    status: str
    process_observed: bool | None
    heartbeat_at: float | None
    checked_at: float
    heartbeat_age_seconds: float | None
    saved_status: str | None
    evidence: str

    def to_dict(self) -> dict:
        return asdict(self)


def assess_liveness(*, process_observed: bool | None, heartbeat_at: float | None,
                    checked_at: float, saved_status: str | None = None,
                    max_age_seconds: float = 60.) -> Liveness:
    """LIVE requires both an actually observed process and a fresh heartbeat.

    The caller must observe the process at checked_at; a saved PID alone is not
    evidence. None means no process check occurred. Future-dated heartbeats do
    not establish freshness, even if a saved status says ACTIVE/RUNNING.
    """
    if not math.isfinite(checked_at) or max_age_seconds <= 0 or not math.isfinite(max_age_seconds):
        raise ValueError("Invalid freshness policy")
    if process_observed not in (True, False, None):
        raise ValueError("process_observed must be true, false, or unknown")
    age = None
    if heartbeat_at is not None and math.isfinite(heartbeat_at):
        age = checked_at - heartbeat_at
    fresh = age is not None and 0 <= age <= max_age_seconds
    if process_observed is False:
        status, evidence = "STOPPED", "Process check found no matching process; saved status is historical."
    elif process_observed is True and fresh:
        status, evidence = "LIVE", "Matching process observed with a heartbeat within the freshness limit."
    elif process_observed is True:
        status, evidence = "PROCESS_OBSERVED_NO_FRESH_HEARTBEAT", "A process exists, but responsiveness is not established."
    else:
        status, evidence = "UNKNOWN", "No current process observation; heartbeat/saved status alone cannot prove liveness."
    return Liveness(status, process_observed, heartbeat_at, checked_at, age, saved_status, evidence)


def evidence_receipt(requested_hashes: Mapping[str, str], *, reviewer: str,
                     inspector: Callable[[Path, bytes], dict], checked_at: float) -> dict:
    """Inspect exactly requested evidence, binding each result to actual bytes.

    The inspector is a deterministic check or an external review adapter supplied
    by the caller. It receives the complete bytes and must return a nonempty
    list of substantive checks and status PASS/FAIL/INCONCLUSIVE. The receipt
    validates byte identity and inspection presence, not the truth of free-text
    scientific reasoning. Missing, changed, skipped, or unchecked files prevent
    completion. No external calls occur unless the caller's inspector makes one.
    """
    if not reviewer or not requested_hashes or not math.isfinite(checked_at):
        raise ValueError("A named reviewer, nonempty evidence request and timestamp are required")
    files = []
    for name, expected in sorted(requested_hashes.items()):
        if len(expected) != 64 or any(c not in "0123456789abcdef" for c in expected):
            raise ValueError("Expected lowercase SHA256 for every requested evidence file")
        path = Path(name)
        raw = path.read_bytes()
        actual = sha256(raw).hexdigest()
        if actual != expected:
            raise ValueError(f"Evidence changed before inspection: {path}")
        result = inspector(path, raw)
        if not isinstance(result, dict) or result.get("status") not in {"PASS", "FAIL", "INCONCLUSIVE"}:
            raise ValueError("Each file needs an explicit inspection verdict")
        checks = result.get("checks")
        if not isinstance(checks, list) or not checks or any(not isinstance(check, str) or not check.strip() for check in checks):
            raise ValueError("Each file needs nonempty substantive inspection checks")
        # Require a JSON-serializable immutable snapshot of the inspector output.
        result = json.loads(json.dumps(result, allow_nan=False))
        files.append({"path": str(path), "sha256": actual, "bytes": len(raw), "inspection": result})
    receipt = {"reviewer": reviewer, "checked_at": checked_at, "status": "COMPLETED",
               "requested_files": len(requested_hashes), "inspected_files": len(files),
               "files": files, "scope": "Exact requested bytes inspected; completion is not scientific approval."}
    receipt["receipt_id"] = sha256(json.dumps(receipt, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
    return receipt


def verify_receipt(receipt: dict, requested_hashes: Mapping[str, str]) -> bool:
    """Check receipt integrity and exact request coverage without rerunning review."""
    try:
        payload = {k: v for k, v in receipt.items() if k != "receipt_id"}
        identity = sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
        files = receipt["files"]
        if identity != receipt["receipt_id"] or receipt.get("status") != "COMPLETED":
            return False
        if len(files) != len(requested_hashes) or receipt["requested_files"] != len(files) or receipt["inspected_files"] != len(files):
            return False
        inspected = {f["path"]: f["sha256"] for f in files}
        if len(inspected) != len(files) or inspected != dict(requested_hashes):
            return False
        return all(f["inspection"]["status"] in {"PASS", "FAIL", "INCONCLUSIVE"} and f["inspection"]["checks"] for f in files)
    except (KeyError, TypeError, ValueError):
        return False
