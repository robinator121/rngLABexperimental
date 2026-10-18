"""RNGLab v0.14: Multi-Testing Significance and Error Rate Control.

Remedies multi-testing inflation arising from repeated campaign cycles:
1. Bonferroni Single-Step FWER Correction
2. Holm-Bonferroni Step-Down FWER Procedure
3. Benjamini-Hochberg False Discovery Rate (FDR) Procedure
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple


@dataclass(frozen=True)
class MultiTestingResult:
    method: str
    nominal_alpha: float
    total_tests: int
    rejected_count: int
    adjusted_p_values: List[float]
    reject_decisions: List[bool]


def bonferroni_correction(
    p_values: Sequence[float],
    alpha: float = 0.05,
) -> MultiTestingResult:
    """Bonferroni FWER control: rejects if p_i <= alpha / m."""
    m = len(p_values)
    if m == 0:
        return MultiTestingResult("bonferroni", alpha, 0, 0, [], [])

    adj_p = [min(1.0, p * m) for p in p_values]
    decisions = [adj <= alpha for adj in adj_p]
    return MultiTestingResult(
        method="bonferroni",
        nominal_alpha=alpha,
        total_tests=m,
        rejected_count=sum(decisions),
        adjusted_p_values=adj_p,
        reject_decisions=decisions,
    )


def holm_bonferroni_correction(
    p_values: Sequence[float],
    alpha: float = 0.05,
) -> MultiTestingResult:
    """Holm-Bonferroni step-down procedure: uniformly more powerful than Bonferroni."""
    m = len(p_values)
    if m == 0:
        return MultiTestingResult("holm", alpha, 0, 0, [], [])

    # Sort p-values with original indices
    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    decisions = [False] * m
    adj_p = [0.0] * m

    running_max = 0.0
    for rank, (orig_idx, p_val) in enumerate(indexed):
        k = rank + 1
        adj = min(1.0, p_val * (m - k + 1))
        running_max = max(running_max, adj)
        adj_p[orig_idx] = running_max

    for i in range(m):
        decisions[i] = (adj_p[i] <= alpha)

    return MultiTestingResult(
        method="holm_bonferroni",
        nominal_alpha=alpha,
        total_tests=m,
        rejected_count=sum(decisions),
        adjusted_p_values=adj_p,
        reject_decisions=decisions,
    )


def benjamini_hochberg_fdr(
    p_values: Sequence[float],
    q_fdr: float = 0.05,
) -> MultiTestingResult:
    """Benjamini-Hochberg procedure controlling the False Discovery Rate at level q."""
    m = len(p_values)
    if m == 0:
        return MultiTestingResult("benjamini_hochberg", q_fdr, 0, 0, [], [])

    indexed = sorted(enumerate(p_values), key=lambda x: x[1])
    decisions = [False] * m
    adj_p = [0.0] * m

    # Calculate step-up critical threshold
    # P_(k) <= (k / m) * q
    max_k_to_reject = -1
    for rank, (orig_idx, p_val) in enumerate(indexed):
        k = rank + 1
        threshold = (k / m) * q_fdr
        if p_val <= threshold:
            max_k_to_reject = rank

    if max_k_to_reject >= 0:
        for rank in range(max_k_to_reject + 1):
            orig_idx = indexed[rank][0]
            decisions[orig_idx] = True

    # Compute adjusted p-values (q-values)
    running_min = 1.0
    for rank in range(m - 1, -1, -1):
        k = rank + 1
        orig_idx, p_val = indexed[rank]
        q_val = min(1.0, (m / k) * p_val)
        running_min = min(running_min, q_val)
        adj_p[orig_idx] = running_min

    return MultiTestingResult(
        method="benjamini_hochberg",
        nominal_alpha=q_fdr,
        total_tests=m,
        rejected_count=sum(decisions),
        adjusted_p_values=adj_p,
        reject_decisions=decisions,
    )
