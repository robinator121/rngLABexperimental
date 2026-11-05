"""RNGLab v0.18: multi-chain convergence diagnostic.

Rank-normalized split-R-hat (Vehtari, Gelman, Simpson, Carpenter & Buerkner,
Bayesian Analysis 16, 667-718, 2021), the "bulk" version without the folded
variant. Each chain is split into halves; all draws of one ensemble are ranked
together (average ranks for ties), mapped to normal scores
Phi^-1((r - 3/8)/(S - 1/4)), and the classic potential scale reduction is
computed on the scores. A large value says the chains disagree; a value near 1
does not prove convergence.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm, rankdata

NOT_CONVERGED_THRESHOLD = 1.01


def split_rhat(chains) -> np.ndarray:
    """chains: (n_chains, n_draws) or (ensembles, n_chains, n_draws). Returns one R-hat per ensemble.

    An odd number of draws drops the middle draw of each chain. NaN means the
    within-chain variance of the scores is zero (all draws tied).
    """
    a = np.asarray(chains, dtype=float)
    if a.ndim == 2:
        a = a[None]
    if a.ndim != 3 or a.shape[1] < 1 or a.shape[2] < 4 or not np.isfinite(a).all():
        raise ValueError("finite chains with shape (chains, draws) or (ensembles, chains, draws), draws >= 4")
    ensembles, n_chains, n_draws = a.shape
    half = n_draws // 2
    halves = np.concatenate([a[:, :, :half], a[:, :, n_draws - half:]], axis=1)
    flat = halves.reshape(ensembles, -1)
    total = flat.shape[1]
    scores = norm.ppf((rankdata(flat, method="average", axis=1) - 0.375) / (total - 0.25))
    scores = scores.reshape(ensembles, 2 * n_chains, half)
    within = scores.var(axis=2, ddof=1).mean(axis=1)
    between = half * scores.mean(axis=2).var(axis=1, ddof=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        pooled = (half - 1) / half * within + between / half
        return np.where(within > 0, np.sqrt(pooled / within), np.nan)
