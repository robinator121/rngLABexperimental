"""RNGLab v0.13+ Standards Integration (NIST SP 800-90A/B/C, SP 800-22)."""
from .drbg_90a import HMAC_DRBG, NIST_WORKED_EXAMPLE, verify_kat_hmac_drbg
from .entropy_90b import (
    MinEntropyReport,
    assess_min_entropy,
    estimate_collision,
    estimate_compression,
    estimate_markov,
    estimate_mcv,
)
from .health_90b import AdaptiveProportionTest, HealthTestResult, RepetitionCountTest
from .sp800_22 import (
    SP800_22_Result,
    block_frequency_test,
    frequency_monobit_test,
    run_sp800_22_battery,
    runs_test,
)
from .traceability import TRACEABILITY_CATALOG, StandardsRequirement, get_traceability_matrix

__all__ = [
    "RepetitionCountTest",
    "AdaptiveProportionTest",
    "HealthTestResult",
    "assess_min_entropy",
    "estimate_mcv",
    "estimate_collision",
    "estimate_markov",
    "estimate_compression",
    "MinEntropyReport",
    "HMAC_DRBG",
    "verify_kat_hmac_drbg",
    "NIST_WORKED_EXAMPLE",
    "frequency_monobit_test",
    "block_frequency_test",
    "runs_test",
    "run_sp800_22_battery",
    "SP800_22_Result",
    "TRACEABILITY_CATALOG",
    "StandardsRequirement",
    "get_traceability_matrix",
]
