"""RNGLab v0.18 namespace: statistics.

Re-exports, unchanged:

- summarize_series (rnglab12.uncertainty): non-overlapping batch means at
  chosen batch sizes plus the Geyer initial-positive/monotone estimator.
- batch_intervals, rate, causal_loss, window_metrics (rnglab17.endpoints):
  vectorized batch-means intervals across replicas, Clopper-Pearson rates,
  paired excess squared loss, rolling-window endpoints.
- summarize_independent_reference (rnglab14.models.pimc_harmonic): Student-t
  interval for independently sampled estimators.
- Bonferroni, Holm and Benjamini-Hochberg (rnglab14.statistics.multitest).

The v0.18 Phase 1 estimators (overlapping batch means, spectral variance,
Geyer, ESS and the batch-length diagnostic) will be added to this module.
"""
from rnglab12.uncertainty import summarize_series
from rnglab14.models.pimc_harmonic import summarize_independent_reference
from rnglab14.statistics.multitest import (
    MultiTestingResult,
    benjamini_hochberg_fdr,
    bonferroni_correction,
    holm_bonferroni_correction,
)
from rnglab17.endpoints import batch_intervals, causal_loss, rate, window_metrics

__all__ = [
    "summarize_series",
    "summarize_independent_reference",
    "batch_intervals",
    "rate",
    "causal_loss",
    "window_metrics",
    "bonferroni_correction",
    "holm_bonferroni_correction",
    "benjamini_hochberg_fdr",
    "MultiTestingResult",
]
