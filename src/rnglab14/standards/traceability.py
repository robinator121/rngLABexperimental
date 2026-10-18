"""RNGLab v0.14: Comprehensive Standards Traceability Catalog.

Maps formal international and national cryptographic/statistical standards to implementation
modules, validation tests, and explicit capability gaps:
- NIST SP 800-90A Rev. 1 (DRBG mechanisms)
- NIST SP 800-90B (Entropy source validation & health tests)
- NIST SP 800-90C (RBG constructions: RBG1, RBG2, RBG3)
- NIST SP 800-22 Rev. 1a (Statistical test suite)
- NIST IR 8446 (Physical entropy hardware specifications)
- BSI AIS 20/31 (German Federal Office for Information Security criteria for PTG.2 / PTG.3)
- ISO/IEC 18031:2011 (Random bit generation architecture)
- FIPS 140-3 (Security requirements for cryptographic modules)

CRITICAL DISCLAIMER:
Statistical test-suite success does NOT establish cryptographic security, entropy-source
validation, NIST CAVP certification, or FIPS 140-3 module validation.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any, List


@dataclass(frozen=True)
class StandardsRequirement:
    req_id: str
    standard: str
    section: str
    title: str
    status: str  # 'VERIFIED', 'IMPLEMENTED_DIAGNOSTIC', 'PARTIAL', 'GAP'
    implementation_module: str
    test_module: str
    limitations_and_gaps: str


STANDARDS_CATALOG: List[StandardsRequirement] = [
    # --- NIST SP 800-90B ---
    StandardsRequirement(
        req_id="REQ-90B-RCT",
        standard="NIST SP 800-90B",
        section="4.4.1",
        title="Repetition Count Test (RCT)",
        status="IMPLEMENTED_DIAGNOSTIC",
        implementation_module="rnglab14.standards.health_90b.RepetitionCountTest",
        test_module="tests.test_standards14.StandardsV14Tests.test_sp800_90b_rct_and_apt",
        limitations_and_gaps="Detects immediate identical repetitions (stuck bits); does not detect pattern repetitions.",
    ),
    StandardsRequirement(
        req_id="REQ-90B-APT",
        standard="NIST SP 800-90B",
        section="4.4.2",
        title="Adaptive Proportion Test (APT)",
        status="IMPLEMENTED_DIAGNOSTIC",
        implementation_module="rnglab14.standards.health_90b.AdaptiveProportionTest",
        test_module="tests.test_standards14.StandardsV14Tests.test_sp800_90b_rct_and_apt",
        limitations_and_gaps="Windowed assessment (W=512 or 1024); detects severe proportion shifts, blind to subtle drift.",
    ),
    StandardsRequirement(
        req_id="REQ-90B-MIN-ENTROPY",
        standard="NIST SP 800-90B",
        section="6.3",
        title="Min-Entropy Estimation Suite",
        status="IMPLEMENTED_DIAGNOSTIC",
        implementation_module="rnglab14.standards.entropy_90b.assess_min_entropy",
        test_module="None",
        limitations_and_gaps="Research heuristics, including Shannon-based compression proxy; not validated90B min-entropy estimators or conservative entropy guarantees.",
    ),
    StandardsRequirement(
        req_id="REQ-90B-RESTART",
        standard="NIST SP 800-90B",
        section="3.1.4",
        title="Restart Tests (Row/Column Sanity)",
        status="GAP",
        implementation_module="None",
        test_module="tests.test_standards14.not_implemented",
        limitations_and_gaps="Not implemented. No restart matrix or physical cold-boot validation.",
    ),

    # --- NIST SP 800-90A ---
    StandardsRequirement(
        req_id="REQ-90A-HMAC-DRBG",
        standard="NIST SP 800-90A Rev. 1",
        section="10.1.2",
        title="HMAC_DRBG with SHA-256",
        status="IMPLEMENTED_DIAGNOSTIC",
        implementation_module="rnglab14.standards.drbg_90a.HMAC_DRBG",
        test_module="tests.test_standards14.StandardsV14Tests.test_kat_hmac_drbg",
        limitations_and_gaps="Verified against official NIST CAVP Known Answer Test vectors. Not CMVP certified.",
    ),
    StandardsRequirement(
        req_id="REQ-90A-CTR-DRBG",
        standard="NIST SP 800-90A Rev. 1",
        section="10.2.1",
        title="CTR_DRBG with AES-256",
        status="GAP",
        implementation_module="None",
        test_module="None",
        limitations_and_gaps="Identified gap. Currently HMAC_DRBG is the primary reference DRBG.",
    ),

    # --- NIST SP 800-90C ---
    StandardsRequirement(
        req_id="REQ-90C-CONSTRUCTION",
        standard="NIST SP 800-90C (2025 final; construction not validated)",
        section="4.1",
        title="RBG Construction Classes (RBG1, RBG2, RBG3)",
        status="IMPLEMENTED_DIAGNOSTIC",
        implementation_module="rnglab14.sources.conditioning.ReseedingDRBGPipeline",
        test_module="tests.test_sources14.test_hybrid_reseeding_drbg",
        limitations_and_gaps="Implements RBG2-type architecture (DRBG instantiated and reseeded from physical/entropy source).",
    ),

    # --- NIST SP 800-22 ---
    StandardsRequirement(
        req_id="REQ-SP800-22-BATTERY",
        standard="NIST SP 800-22 Rev. 1a",
        section="2.1 - 2.15",
        title="Statistical Test Battery for Randomness",
        status="IMPLEMENTED_DIAGNOSTIC",
        implementation_module="rnglab14.standards.sp800_22.run_sp800_22_battery",
        test_module="tests.test_standards14.StandardsV14Tests.test_sp800_22_battery",
        limitations_and_gaps="Implements core 5 tests (Frequency, Block Frequency, Runs, Longest Run, Approximate Entropy). Full 15-test battery is optional diagnostic.",
    ),

    # --- BSI AIS 20/31 ---
    StandardsRequirement(
        req_id="REQ-AIS-31-PTG2",
        standard="BSI AIS 20/31",
        section="3.3",
        title="Physical True Random Number Generator Class PTG.2 Criteria",
        status="GAP",
        implementation_module="None",
        test_module="tests.test_standards14.not_implemented",
        limitations_and_gaps="No implemented PTG.2 validation; no physical entropy assurance.",
    ),

    # --- ISO/IEC 18031 & FIPS 140-3 ---
    StandardsRequirement(
        req_id="REQ-ISO-18031-ARCH",
        standard="ISO/IEC 18031:2011",
        section="7",
        title="Non-Deterministic and Deterministic Random Bit Generator Architecture",
        status="IMPLEMENTED_DIAGNOSTIC",
        implementation_module="rnglab14.sources.base.RandomSource",
        test_module="tests.test_sources14.test_source_architecture",
        limitations_and_gaps="Decouples entropy source, conditioning, and state expansion per ISO/IEC 18031 architectural specifications.",
    ),
    StandardsRequirement(
        req_id="REQ-FIPS-140-3-CAVEAT",
        standard="FIPS 140-3",
        section="9.2",
        title="Entropy Source and Health Testing Caveats",
        status="IMPLEMENTED_DIAGNOSTIC",
        implementation_module="rnglab14.standards.traceability",
        test_module="tests.test_standards14.StandardsV14Tests.test_standards_catalog_coverage",
        limitations_and_gaps="Explicitly documents that empirical test success does not constitute formal FIPS 140-3 Cryptographic Module Validation Program (CMVP) validation.",
    ),
]


TRACEABILITY_CATALOG = STANDARDS_CATALOG


def get_standards_catalog() -> list[dict[str, Any]]:
    """Return catalog of all standards requirements and implementation links."""
    return [asdict(r) for r in STANDARDS_CATALOG]


get_traceability_matrix = get_standards_catalog

