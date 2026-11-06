"""RNGLab v0.18 namespace: monitors and their calibration.

Re-exports, unchanged:

- FrozenRankCalibration (rnglab14.monitoring.calibration): conservative
  marginal rank p-values for whole-trajectory scores under exchangeability.
- rank_threshold (rnglab17.endpoints): the healthy rank threshold used by the
  v0.17 whole-path calibration.
- RepetitionCountTest, AdaptiveProportionTest (rnglab14.standards.health_90b):
  SP 800-90B-style continuous health tests with binomial cutoffs.

Not re-exported: rnglab14.monitoring.DualIntegrityMonitor, whose v0.14
"0.0% hazardous gap" result is superseded (docs/v018/AUDIT.md section 2.2).
The Phase II marginal, joint-channel and application-shape monitors exist only
in research/phase2-20261030/experiment/fault_study.py; the Phase 3 monitors
will be added here as new v0.18 code.
"""
from rnglab14.monitoring.calibration import FrozenRankCalibration
from rnglab14.standards.health_90b import AdaptiveProportionTest, HealthTestResult, RepetitionCountTest
from rnglab17.endpoints import rank_threshold

__all__ = [
    "FrozenRankCalibration",
    "rank_threshold",
    "RepetitionCountTest",
    "AdaptiveProportionTest",
    "HealthTestResult",
]
