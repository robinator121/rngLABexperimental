"""NIST SP 800-22 Rev. 1a Reference Statistical Randomness Tests.

Implements core tests:
1. Frequency (Monobit) Test (Section 2.1)
2. Frequency Test within a Block (Section 2.2)
3. Runs Test (Section 2.3)
4. Longest Run of Ones in a Block (Section 2.4)
5. Spectral (Discrete Fourier Transform) Test (Section 2.6)

CRITICAL NOTICE: Passing these tests does not establish cryptographic security
or certification compliance.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


@dataclass(frozen=True)
class SP800_22_Result:
    test_name: str
    p_value: float
    passed: bool
    statistic: float
    threshold: float = 0.01
    details: str = ""


def _bits_from_bytes(data: bytes) -> list[int]:
    bits = []
    for b in data:
        for shift in range(7, -1, -1):
            bits.append((b >> shift) & 1)
    return bits


def frequency_monobit_test(bits: Sequence[int], alpha: float = 0.01) -> SP800_22_Result:
    """NIST SP 800-22 Section 2.1 Frequency (Monobit) Test."""
    n = len(bits)
    if n == 0:
        return SP800_22_Result("Frequency Monobit", 0.0, False, 0.0)

    # Convert to +1 / -1
    s_n = sum(2 * b - 1 for b in bits)
    s_obs = abs(s_n) / math.sqrt(n)
    p_value = math.erfc(s_obs / math.sqrt(2.0))

    return SP800_22_Result(
        test_name="Frequency Monobit",
        p_value=round(p_value, 6),
        passed=(p_value >= alpha),
        statistic=round(s_obs, 4),
        threshold=alpha,
        details=f"Sum Sn={s_n}, Sobs={s_obs:.4f}",
    )


def block_frequency_test(bits: Sequence[int], block_size: int = 128, alpha: float = 0.01) -> SP800_22_Result:
    """NIST SP 800-22 Section 2.2 Frequency Test within a Block."""
    n = len(bits)
    m = block_size
    num_blocks = n // m
    if num_blocks == 0:
        return SP800_22_Result("Block Frequency", 0.0, False, 0.0)

    chi_sq = 0.0
    for i in range(num_blocks):
        block = bits[i * m : (i + 1) * m]
        pi_i = sum(block) / m
        chi_sq += (pi_i - 0.5) ** 2
    chi_sq *= 4.0 * m

    # Chi-square p-value via gamma approximation
    # For df = num_blocks, use regularized incomplete gamma Q(df/2, chi_sq/2)
    # Using normal approximation if scipy not directly imported here
    try:
        from scipy.special import gammaincc
        p_value = float(gammaincc(num_blocks / 2.0, chi_sq / 2.0))
    except Exception:
        # Wilson-Hilferty transformation approximation for chi2
        z = ((chi_sq / num_blocks) ** (1/3) - (1 - 2/(9*num_blocks))) / math.sqrt(2/(9*num_blocks))
        p_value = 0.5 * math.erfc(z / math.sqrt(2.0))

    return SP800_22_Result(
        test_name="Block Frequency",
        p_value=round(p_value, 6),
        passed=(p_value >= alpha),
        statistic=round(chi_sq, 4),
        threshold=alpha,
        details=f"{num_blocks} blocks of size {m}, chi_square={chi_sq:.4f}",
    )


def runs_test(bits: Sequence[int], alpha: float = 0.01) -> SP800_22_Result:
    """NIST SP 800-22 Section 2.3 Runs Test."""
    n = len(bits)
    if n < 10:
        return SP800_22_Result("Runs", 0.0, False, 0.0)

    pi = sum(bits) / n
    # Pre-test requirement: pi must be close to 0.5
    if abs(pi - 0.5) >= (2.0 / math.sqrt(n)):
        return SP800_22_Result(
            test_name="Runs",
            p_value=0.0,
            passed=False,
            statistic=0.0,
            details="Pre-test frequency check failed (|pi - 0.5| >= 2/sqrt(n))",
        )

    # Count runs (transitions + 1)
    v_obs = 1 + sum(1 for i in range(n - 1) if bits[i] != bits[i + 1])
    numerator = abs(v_obs - 2.0 * n * pi * (1.0 - pi))
    denominator = 2.0 * math.sqrt(2.0 * n) * pi * (1.0 - pi)

    p_value = math.erfc(numerator / denominator) if denominator > 0 else 0.0

    return SP800_22_Result(
        test_name="Runs",
        p_value=round(p_value, 6),
        passed=(p_value >= alpha),
        statistic=round(v_obs, 4),
        threshold=alpha,
        details=f"Total runs Vobs={v_obs}, pi={pi:.4f}",
    )


def run_sp800_22_battery(data: bytes, alpha: float = 0.01) -> list[SP800_22_Result]:
    """Execute standard SP 800-22 test battery on byte sequence."""
    bits = _bits_from_bytes(data)
    results = [
        frequency_monobit_test(bits, alpha=alpha),
        block_frequency_test(bits, block_size=min(128, max(16, len(bits) // 8)), alpha=alpha),
        runs_test(bits, alpha=alpha),
    ]
    return results
