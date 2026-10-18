"""RNGLab v0.14: Statistics Package."""

from .multitest import (
    MultiTestingResult,
    benjamini_hochberg_fdr,
    bonferroni_correction,
    holm_bonferroni_correction,
)

__all__ = [
    "MultiTestingResult",
    "bonferroni_correction",
    "holm_bonferroni_correction",
    "benjamini_hochberg_fdr",
]
