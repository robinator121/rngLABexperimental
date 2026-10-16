"""Standards Traceability Matrix linking requirements to source code and tests."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class StandardsRequirement:
    req_id: str
    standard: str
    section: str
    title: str
    implementation_module: str
    test_module: str
    status: str
    limitations: str


TRACEABILITY_CATALOG: list[StandardsRequirement] = [
    StandardsRequirement(
        req_id="REQ-90B-RCT",
        standard="NIST SP 800-90B",
        section="4.4.1",
        title="Repetition Count Test (RCT)",
        implementation_module="rnglab13.standards.health_90b.RepetitionCountTest",
        test_module="tests.test_standards13.StandardsV13Tests.test_repetition_count_test",
        status="VERIFIED",
        limitations="Detects immediate identical repetitions; does not detect periodic or Markov pattern repetitions.",
    ),
    StandardsRequirement(
        req_id="REQ-90B-APT",
        standard="NIST SP 800-90B",
        section="4.4.2",
        title="Adaptive Proportion Test (APT)",
        implementation_module="rnglab13.standards.health_90b.AdaptiveProportionTest",
        test_module="tests.test_standards13.StandardsV13Tests.test_adaptive_proportion_test",
        status="VERIFIED",
        limitations="Windowed assessment; severe faults occurring late in the window may take up to W samples to trigger.",
    ),
    StandardsRequirement(
        req_id="REQ-90B-MIN-ENTROPY",
        standard="NIST SP 800-90B",
        section="6.3",
        title="Min-Entropy Estimation Suite (MCV, Collision, Markov, Compression)",
        implementation_module="rnglab13.standards.entropy_90b.assess_min_entropy",
        test_module="tests.test_standards13.StandardsV13Tests.test_min_entropy_estimation",
        status="VERIFIED",
        limitations="Min-entropy is a lower bound estimate under i.i.d. or 1st-order Markov assumptions; does not guarantee arbitrary higher-order security.",
    ),
    StandardsRequirement(
        req_id="REQ-90A-HMAC-DRBG",
        standard="NIST SP 800-90A Rev.1",
        section="10.1.2",
        title="HMAC_DRBG with SHA-256 and Reseeding",
        implementation_module="rnglab13.standards.drbg_90a.HMAC_DRBG",
        test_module="tests.test_standards13.StandardsV13Tests.test_kat_hmac_drbg",
        status="VERIFIED",
        limitations="Reference implementation verified against official NIST CAVP test vectors. Not FIPS/CMVP certified.",
    ),
    StandardsRequirement(
        req_id="REQ-SP800-22",
        standard="NIST SP 800-22 Rev.1a",
        section="2.1-2.4",
        title="Statistical Randomness Battery (Monobit, Block Frequency, Runs)",
        implementation_module="rnglab13.standards.sp800_22.run_sp800_22_battery",
        test_module="tests.test_standards13.StandardsV13Tests.test_sp800_22_battery",
        status="VERIFIED",
        limitations="Passing statistical tests is a necessary but never sufficient condition for cryptographic security.",
    ),
]


def get_traceability_matrix() -> list[dict[str, Any]]:
    """Return all requirements as structured dictionaries."""
    return [asdict(r) for r in TRACEABILITY_CATALOG]
