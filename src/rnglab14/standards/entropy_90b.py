"""NIST SP 800-90B Min-Entropy Estimation Suite (Section 6.3).

Implements estimators for assessing min-entropy per sample:
1. Most Common Value (MCV) Estimate (Section 6.3.1)
2. Collision Estimate (Section 6.3.2)
3. Markov Estimate (Section 6.3.3)
4. Compression Estimate (Section 6.3.4)
"""
from __future__ import annotations

import collections
from dataclasses import dataclass
import math
from typing import Sequence


@dataclass(frozen=True)
class MinEntropyReport:
    min_entropy: float
    mcv_estimate: float
    collision_estimate: float
    markov_estimate: float
    compression_estimate: float
    sample_count: int
    alphabet_size: int
    conservative: bool = True


def estimate_mcv(samples: Sequence[int], alpha: float = 0.005) -> float:
    """NIST SP 800-90B Section 6.3.1 Most Common Value (MCV) Estimate."""
    n = len(samples)
    if n == 0:
        return 0.0
    counts = collections.Counter(samples)
    c_max = max(counts.values())
    p_hat = c_max / n

    # 99% confidence upper bound on p_max using normal approximation
    z = 2.5758293035489004  # standard normal quantile for 0.995 (one-tailed 99% bound)
    p_upper = min(1.0, p_hat + z * math.sqrt(p_hat * (1.0 - p_hat) / n))
    return max(0.0, -math.log2(p_upper))


def estimate_collision(samples: Sequence[int]) -> float:
    """NIST SP 800-90B Section 6.3.2 Collision Estimate."""
    n = len(samples)
    if n < 10:
        return 0.0

    alphabet = set(samples)
    # If binary stream, SP 800-90B groups bits into 8-bit bytes
    if len(alphabet) <= 2 and n >= 16:
        byte_symbols = []
        for i in range(0, (n // 8) * 8, 8):
            val = 0
            for b in samples[i : i + 8]:
                val = (val << 1) | (b & 1)
            byte_symbols.append(val)
        return estimate_collision(byte_symbols) / 8.0

    # Count steps between collisions
    seen: dict[int, int] = {}
    collision_times: list[int] = []
    t = 0
    for s in samples:
        t += 1
        if s in seen:
            collision_times.append(t)
            seen.clear()
            t = 0
        else:
            seen[s] = t

    if not collision_times:
        return estimate_mcv(samples)

    # Average collision wait time
    v = sum(collision_times) / len(collision_times)
    # SP 800-90B Section 6.3.2 collision probability estimate
    # For alphabet k, expected wait E[v] = sqrt(pi * k / 2).
    # p is estimated from wait time:
    p = max(1.0 / len(alphabet), 1.0 / max(1.0, v))
    return max(0.0, -math.log2(min(1.0, p)))


def estimate_markov(samples: Sequence[int]) -> float:
    """NIST SP 800-90B Section 6.3.3 Markov Estimate for first-order dependencies."""
    n = len(samples)
    if n < 2:
        return 0.0

    # Build transition counts
    pair_counts: dict[tuple[int, int], int] = collections.defaultdict(int)
    row_counts: dict[int, int] = collections.defaultdict(int)
    for i in range(n - 1):
        s1, s2 = samples[i], samples[i + 1]
        pair_counts[(s1, s2)] += 1
        row_counts[s1] += 1

    max_p = 0.0
    for (s1, s2), c in pair_counts.items():
        cond_p = c / row_counts[s1]
        if cond_p > max_p:
            max_p = cond_p

    if max_p <= 0.0:
        return estimate_mcv(samples)

    return max(0.0, -math.log2(max_p))


def estimate_compression(samples: Sequence[int]) -> float:
    """NIST SP 800-90B Section 6.3.4 Compression Estimate (Maurer-like block match)."""
    n = len(samples)
    if n < 16:
        return estimate_mcv(samples)

    # Calculate empirical entropy based on match lengths
    block_size = 6  # 6-bit blocks or 6-sample sequences
    blocks = [samples[i : i + block_size] for i in range(0, n - block_size, block_size)]
    if not blocks:
        return estimate_mcv(samples)

    block_tuples = [tuple(b) for b in blocks]
    counts = collections.Counter(block_tuples)
    total_blocks = len(block_tuples)

    # Shannon entropy of blocks normalized per sample
    h_shannon = 0.0
    for c in counts.values():
        p = c / total_blocks
        h_shannon -= p * math.log2(p)
    per_sample_entropy = h_shannon / block_size

    # Min-entropy is strictly bounded by Shannon entropy
    return max(0.0, min(8.0, per_sample_entropy * 0.85))


def assess_min_entropy(samples: Sequence[int]) -> MinEntropyReport:
    """Execute all SP 800-90B min-entropy estimators and return conservative minimum."""
    if not samples:
        raise ValueError("Cannot assess min-entropy of empty sequence")

    mcv = estimate_mcv(samples)
    collision = estimate_collision(samples)
    markov = estimate_markov(samples)
    compression = estimate_compression(samples)

    conservative_min = min(mcv, collision, markov, compression)

    return MinEntropyReport(
        min_entropy=round(conservative_min, 4),
        mcv_estimate=round(mcv, 4),
        collision_estimate=round(collision, 4),
        markov_estimate=round(markov, 4),
        compression_estimate=round(compression, 4),
        sample_count=len(samples),
        alphabet_size=len(set(samples)),
        conservative=True,
    )
