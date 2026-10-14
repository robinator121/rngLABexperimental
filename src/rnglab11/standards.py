"""Standards coverage and provenance-preserving external report integration.

External batteries are not bundled, installed, launched, or certified by this
module. Their original reports can be attached without turning an absent or
unfinished test suite into a passing result.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
from typing import Any


SUITES = {
    "nist_sts": {"name": "NIST SP 800-22 STS", "executables": ["assess"],
                 "reference": "https://csrc.nist.gov/projects/random-bit-generation/documentation-and-software"},
    "nist_90b": {"name": "NIST SP 800-90B EntropyAssessment", "executables": ["ea_iid", "ea_non_iid"],
                 "reference": "https://github.com/usnistgov/SP800-90B_EntropyAssessment"},
    "testu01": {"name": "TestU01", "executables": [],
                "reference": "https://github.com/umontreal-simul/TestU01-2009"},
    "practrand": {"name": "PractRand", "executables": ["RNG_test"],
                  "reference": "https://pracrand.sourceforge.net/"},
    "dieharder": {"name": "Dieharder", "executables": ["dieharder"],
                  "reference": "https://webhome.phy.duke.edu/~rgb/General/dieharder.php"},
    "ais31": {"name": "AIS 20/31", "executables": [],
              "reference": "https://csrc.nist.gov/pubs/ir/8446/final"},
}


def standards_inventory() -> dict[str, Any]:
    external = {}
    for key, suite in SUITES.items():
        paths = {name: shutil.which(name) for name in suite["executables"]}
        external[key] = {**suite, "discovered_executables": paths,
                         "status": "NOT_RUN", "support": "REPORT_IMPORT_ONLY",
                         "certified": False,
                         "installation_verified": False,
                         "note": "PATH detection does not verify tool identity, version, configuration, or completeness"}
    return {"builtin": {"sp800_22": "THREE_TEST_SUBSET: frequency/block-frequency/runs",
                        "sp800_90b": "RCT/APT-informed checks; full entropy assessment NOT_IMPLEMENTED",
                        "sp800_90a_90c": "Architecture documentation only; no conformance certification"},
            "external": external, "certified": False,
            "interpretation": "Statistical acceptance is not physical-source validation or cryptographic assurance"}


class ExternalValidationAdapter:
    """Attach an operator-created external report to its exact input digest.

    This deliberately does not expose shell commands, install tools, upload
    datasets, or launch substantial batteries. ``completed`` is the operator's
    declaration; RNGLab preserves it without independently asserting validity.
    Unrecognized output remains UNPARSED rather than being called a pass.
    """
    def __init__(self, suite: str, maximum_report_bytes: int = 2_000_000):
        if suite not in SUITES:
            raise ValueError(f"unknown external suite: {suite}")
        if (isinstance(maximum_report_bytes, bool) or not isinstance(maximum_report_bytes, int)
                or not 1 <= maximum_report_bytes <= 10_000_000):
            raise ValueError("bounded report allowance from 1 to 10000000 bytes required")
        self.suite = suite
        self.maximum_report_bytes = maximum_report_bytes

    def import_report(self, report_path: str | Path, *, dataset_sha256: str,
                      tool_version: str, configuration: dict[str, Any], completed: bool,
                      operator_notes: str = "") -> dict[str, Any]:
        if not isinstance(dataset_sha256, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", dataset_sha256):
            raise ValueError("exact input dataset SHA256 is required")
        if not isinstance(tool_version, str) or not tool_version.strip() or len(tool_version) > 200:
            raise ValueError("explicit bounded external tool version required")
        if not isinstance(configuration, dict) or not configuration:
            raise ValueError("external battery/test/bit-order settings must be recorded")
        try:
            configuration_json = json.dumps(configuration, sort_keys=True, allow_nan=False)
        except (TypeError, ValueError) as exc:
            raise ValueError("configuration must be finite JSON-compatible metadata") from exc
        if len(configuration_json.encode()) > 64_000:
            raise ValueError("configuration exceeds the 64000-byte metadata allowance")
        if not isinstance(completed, bool) or not isinstance(operator_notes, str) or len(operator_notes)>10000:
            raise ValueError("explicit completion declaration and bounded notes required")
        path = Path(report_path)
        if not path.is_file() or path.stat().st_size > self.maximum_report_bytes:
            raise ValueError("regular report file within configured byte allowance required")
        with path.open("rb") as stream:
            data = stream.read(self.maximum_report_bytes+1)
        if len(data)>self.maximum_report_bytes:
            raise ValueError("report exceeds byte allowance")
        report = data.decode("utf-8", errors="replace")
        evidence = []
        # Do not infer a complete battery verdict from isolated p-values or
        # absence of the word FAIL. These markers are evidence snippets only.
        for number, line in enumerate(report.splitlines(), start=1):
            marker = re.search(r"\b(PASSED|WEAK|FAILED|FAIL)\b", line)
            if marker:
                evidence.append({"line": number, "marker": marker.group(1), "text": line[:1000]})
                if len(evidence)>=1000:
                    break
        return {"version": 1, "suite": self.suite, "suite_name": SUITES[self.suite]["name"],
                "status": "IMPORTED_EXTERNAL_REPORT" if completed else "IMPORTED_INCOMPLETE_REPORT",
                "assessment": "UNVERIFIED_EXTERNAL_EVIDENCE", "certified": False,
                "dataset_sha256": dataset_sha256.lower(), "report_sha256": hashlib.sha256(data).hexdigest(),
                "report_bytes": len(data), "report_path": str(path.resolve()),
                "tool_version_declared": tool_version, "configuration": json.loads(configuration_json),
                "completion_declared_by_operator": completed,
                "parsed_markers": evidence, "parse_status": "MARKERS_ONLY" if evidence else "UNPARSED",
                "operator_notes": operator_notes,
                "imported_at_utc": datetime.now(timezone.utc).isoformat(),
                "limitations": "Report provenance is recorded; tool authenticity, full execution, scientific interpretation and certification require independent review"}
