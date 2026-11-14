"""Independent recount library for the v0.19 Stage D pilot (H2-R).

Independence rule. This file imports only the standard library, NumPy and SciPy.
It imports nothing from research/v019-d-pilot, research/v018-h1-detection,
rnglab12.uncertainty or rnglab11.statistics. Every quantity is recomputed from
its definition in research/v019-d-pilot/protocol.json ("P" below, with the JSON
key), through a different formula or code path wherever one exists:

    quantity                     primary pipeline              this file
    base stream                  pure-Python PCG64Source       numpy.random.PCG64 with explicit state
    sweep of a word              searchsorted on the marks     np.repeat over words per sweep
    fault transforms             per-word scalar adapter       vectorised (SHIFT, DUP), own loop (LAG1)
    reference energy E_M         closed form (tanh/asinh)      normal-mode sum
    SP 800-22 p-values           erfc / gammaincc              normal and chi-square survival functions
    serial / ApEn counts         order-9 counts, marginalised  one direct count per order
    longest run of ones          padded differences            running reset index
    duplicate p-value            Poisson survival function     regularised incomplete gamma
    batch-means interval         NumPy float means             exactly rounded sums (math.fsum)
    Clopper-Pearson              beta quantiles                root of the binomial tail equations
    order-statistic threshold    full sort                     partial selection
    least-squares slope          float sums                    exact rational arithmetic
    logistic optimum             L-BFGS-B                      damped Newton on the same objective
    validity probability         scipy.stats.nct               quadrature over the chi-square mixture

Where bit-identical agreement is required and only one implementation exists
(statistics.NormalDist for LAG1, the NumPy generator calls that define the
profile draws and the bootstrap indices), the same library call is used and
the surrounding driver code is independent. Those shared calls are listed in
docs/v019/IMPLEMENTATION_AUDIT.md.
"""
from __future__ import annotations

import calendar
import datetime
import hashlib
import json
import math
import re
import zoneinfo
from fractions import Fraction
from pathlib import Path
from statistics import NormalDist
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from numpy.typing import NDArray
from scipy import integrate, optimize, special, stats

HERE = Path(__file__).resolve().parent
PROTOCOL_PATH = HERE.parent / "v019-d-pilot" / "protocol.json"

H = 40960
LENGTHS = (10240, 20480, 40960)
BATTERY_WINDOW, APPLICATION_WINDOW, STRIDE = 512, 1024, 256
BATTERY_ENDS = np.arange(BATTERY_WINDOW, H + 1, STRIDE, dtype=np.int64)
APPLICATION_ENDS = np.arange(APPLICATION_WINDOW, H + 1, STRIDE, dtype=np.int64)
P_KEYS = ("frequency", "block_frequency", "runs", "serial_1", "serial_2", "approximate_entropy",
          "cumulative_sums", "longest_run", "lag1_correlation", "adjacent_duplicates")
BIT_CAP_BITS = 1_000_000
P_FLOOR = 1e-300
TOP = 1 - 2.0**-53
FAMILIES = ("SHIFT", "LAG1", "DUP")
TARGETS = {"SHIFT": "frequency_signed_z", "LAG1": "lag1_feature", "DUP": "duplicate_rate"}
MODELS = ("F0", "F0+", "F-B-lite", "F-AB-lite")
FloatArray = NDArray[np.float64]
_MASK128 = (1 << 128) - 1
# NIST SP 800-22 Rev 1a, sections 2.4 and 3.4: block size -> (K, lowest class, class probabilities).
_LONGEST_RUN = {8: (3, 1, (0.2148, 0.3672, 0.2305, 0.1875)),
                128: (5, 4, (0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124)),
                10000: (6, 10, (0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727))}


class RecountError(ValueError):
    """An input violates a structural requirement of the frozen protocol."""


def load_protocol(path: Path = PROTOCOL_PATH) -> dict[str, Any]:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def protocol_caps(protocol: Mapping[str, Any] | None = None) -> dict[str, float]:
    """P profile.caps: read from the frozen document, never from pipeline code."""
    protocol = load_protocol() if protocol is None else protocol
    return {family: float(protocol["profile"]["caps"][family]["s_max"]) for family in FAMILIES}


# ---------------------------------------------------------------------------
# Seeds, profile draws, severity (P data.seeds, profile)
# ---------------------------------------------------------------------------
def namespace_seed(namespace: str, condition: str, role: str, replica: int) -> int:
    label = "/".join([namespace, condition, role, f"{replica:05d}"])
    return int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")


def analysis_seed(namespace: str, condition: str, purpose: str) -> int:
    """P data.seeds.analysis: the role is 'analysis', the unit is the purpose string."""
    label = "/".join([namespace, condition, "analysis", purpose])
    return int.from_bytes(hashlib.sha256(label.encode()).digest()[:8], "big")


def run_seeds(namespace: str, condition: str, replica: int) -> dict[str, int]:
    return {role: namespace_seed(namespace, condition, role, replica)
            for role in ("local", "centroid", "start", "fault", "profile")}


# The author's decision of 2026-11-13 on interrupted runs (not part of the frozen protocol, which is
# silent): the same identity is run again, at most three times; after that the run is unresolved.
INTERRUPTED_RUN_POLICY = "regenerate_same_id"
MAX_INFRASTRUCTURE_RETRIES = 3
INFRASTRUCTURE_REASONS = ("process_ended_without_record", "operating_system_error", "memory_error")
INFRASTRUCTURE_PHASE = "infrastructure"
UNRESOLVED_INFRASTRUCTURE = "UNRESOLVED_INFRASTRUCTURE"
RNG_INPUT_WORDS = 128
# Stop criteria that use the runs of a condition; an unresolved run keeps them from PASS.
CRITERIA_USING = {"calibration": (2, 3, 4), "healthy": (1, 2, 3, 4), "fault": (2, 3, 4)}


def rng_input_sha256(seeds: Mapping[str, int]) -> str:
    """Hash of what a run's five generators emit first, from generators built here.

    Per role the SHA-256 of the first 128 outputs as big-endian 64-bit words: the repository's
    PCG64 source (reproduced by `pcg64_words`) for local, centroid and start, and the bit
    generator of numpy.random.default_rng(seed) for fault and profile. The result is the
    SHA-256 of the compact, key-sorted JSON of the five hashes. Integer arithmetic only.
    """
    roles = {}
    for role in ("local", "centroid", "start", "fault", "profile"):
        if role in ("local", "centroid", "start"):
            words = pcg64_words(seeds[role], RNG_INPUT_WORDS)
        else:
            words = np.random.default_rng(seeds[role]).bit_generator.random_raw(RNG_INPUT_WORDS)
        roles[role] = hashlib.sha256(np.asarray(words, dtype=np.uint64).astype(">u8").tobytes()).hexdigest()
    return hashlib.sha256(json.dumps(roles, sort_keys=True, separators=(",", ":")).encode("ascii")).hexdigest()


def criteria_blocked_by(conditions: Iterable[str]) -> list[int]:
    """Stop criteria that cannot be PASS when runs of these conditions are unresolved."""
    return sorted({key for condition in conditions
                   for key in CRITERIA_USING["fault" if condition.startswith("fault/") else condition]})


def draw_profile(profile_seed: int, s_max: float) -> dict[str, Any]:
    """P profile.draw_order. These generator calls are the definition of the draw."""
    rng = np.random.default_rng(profile_seed)
    t0 = int(rng.integers(2048, 8193))
    t_cap = np.float64(np.exp(rng.uniform(np.log(2048.0), np.log(163840.0))))
    return {"t0": t0, "T_cap": t_cap, "r": np.float64(np.float64(s_max) / t_cap), "s_max": np.float64(s_max)}


def severity_path(t0: int | None, t_cap: float | None, s_max: float | None, sweeps: int = H) -> FloatArray:
    """P profile.severity, vectorised: 0 before t0, min(s_max, r*(t-t0)) from t0, r = s_max/T_cap."""
    out = np.zeros(sweeps, dtype=np.float64)
    if t0 is None:
        return out
    r = np.float64(s_max) / np.float64(t_cap)
    t = np.arange(sweeps, dtype=np.int64)
    active = t >= t0
    out[active] = np.minimum(np.float64(s_max), r * (t[active] - t0).astype(np.float64))
    return out


def s_bar_closed_form(length: int, t0: int, t_cap: float, s_max: float) -> float:
    """P oracle.S_bar, algebraic form: mean severity over sweeps [0, length)."""
    m = max(length - 1 - t0, 0)
    q = min(m, math.floor(t_cap))
    return float((s_max / t_cap * q * (q + 1) / 2 + s_max * (m - q)) / length)


# ---------------------------------------------------------------------------
# Byte ledger and streams (P sweep_clock, fault_transforms)
# ---------------------------------------------------------------------------
def check_ledger(marks: NDArray[np.int64], sweeps: int | None = None) -> list[str]:
    """I-1 conditions on the recorded marks; returns the list of violations (empty when clean)."""
    problems: list[str] = []
    marks = np.asarray(marks)
    if marks.dtype != np.dtype("int64"):
        problems.append(f"marks dtype {marks.dtype}, expected int64")
    if marks.ndim != 1:
        return problems + [f"marks has {marks.ndim} dimensions"]
    if sweeps is not None and len(marks) != sweeps:
        problems.append(f"{len(marks)} marks, expected {sweeps}")
    step = np.diff(marks.astype(np.int64), prepend=0)
    bad = np.flatnonzero((step < 128) | (step > 256) | (step % 8 != 0))
    if len(bad):
        problems.append(f"{len(bad)} sweeps outside 128+8k, k in 0..16; first at sweep {int(bad[0])} ({int(step[bad[0]])} bytes)")
    return problems


def words_per_sweep(marks: NDArray[np.int64]) -> NDArray[np.int64]:
    step = np.diff(np.asarray(marks, dtype=np.int64), prepend=0)
    if np.any(step <= 0) or np.any(step % 8):
        raise RecountError("marks must be increasing and word aligned")
    return step // 8


def sweep_of_word(marks: NDArray[np.int64]) -> NDArray[np.int64]:
    """Sweep index of every consumed word: word w lies in sweep t iff B(t-1) <= 8w < B(t)."""
    return np.repeat(np.arange(len(marks), dtype=np.int64), words_per_sweep(marks))


def pcg64_words(seed: int, count: int) -> NDArray[np.uint64]:
    """Base stream of the repository's PCG64 source, from NumPy's own PCG64.

    The source uses state = seed, increment = 3 (stream id 1), XSL-RR 128/64 with
    the default multiplier, and discards one output at construction.
    """
    generator = np.random.PCG64(0)
    state = generator.state
    state["state"] = {"state": int(seed) & _MASK128, "inc": 3}
    state["has_uint32"], state["uinteger"] = 0, 0
    generator.state = state
    generator.random_raw(1)
    return generator.random_raw(count)


def words_to_bytes(words: NDArray[np.uint64]) -> bytes:
    return np.asarray(words, dtype=np.uint64).astype(">u8").tobytes()


def bytes_to_words(data: bytes) -> NDArray[np.uint64]:
    if len(data) % 8:
        raise RecountError("stream is not word aligned")
    return np.frombuffer(data, dtype=">u8").astype(np.uint64)


def apply_fault(family: str, base: NDArray[np.uint64], word_sweep: NDArray[np.int64], severity: FloatArray,
                t0: int | None, fault_seed: int) -> NDArray[np.uint64]:
    """Emitted words from base words, the sweep of each word and the per-sweep severity."""
    base = np.asarray(base, dtype=np.uint64)
    if family == "healthy":
        return base.copy()
    if family not in FAMILIES:
        raise RecountError(f"unknown family {family!r}")
    if len(base) != len(word_sweep):
        raise RecountError("one sweep index per word is required")
    post = word_sweep >= t0
    s_word = np.asarray(severity, dtype=np.float64)[word_sweep]
    low, mantissa = base & np.uint64(0x7FF), base >> np.uint64(11)
    u = mantissa.astype(np.float64) * 2.0**-53
    out = base.copy()
    if family == "SHIFT":
        shifted = np.minimum(u + s_word, TOP)
        out[post] = ((shifted * 2.0**53).astype(np.uint64) << np.uint64(11) | low)[post]
        return out
    index = np.flatnonzero(post)
    if family == "DUP":
        duplicate = np.zeros(len(base), dtype=bool)
        duplicate[index] = np.random.default_rng(fault_seed).random(len(index)) < s_word[index]
        if len(base):
            duplicate[0] = False  # no previous word exists
        source = np.maximum.accumulate(np.where(duplicate, -1, np.arange(len(base))))
        return base[source]
    std, z, mapped = NormalDist(), None, np.empty(len(index), dtype=np.float64)
    for i, (p, rho) in enumerate(zip((u[index] + 2.0**-54).tolist(), s_word[index].tolist())):
        e = std.inv_cdf(p)  # a domain error is a run failure under P fault_transforms.LAG1
        z = e if z is None else rho * z + math.sqrt(1 - rho * rho) * e
        mapped[i] = min(std.cdf(z), TOP)
    out[index] = (mapped * 2.0**53).astype(np.uint64) << np.uint64(11) | low[index]
    return out


def window_words(marks: NDArray[np.int64], start: int, end: int) -> tuple[int, int]:
    """Word range of sweeps [start, end)."""
    return (0 if start == 0 else int(marks[start - 1]) // 8), int(marks[end - 1]) // 8


# ---------------------------------------------------------------------------
# Conventional battery (P operational_alarm, signal_schema)
# ---------------------------------------------------------------------------
def cyclic_pattern_values(bits: NDArray[np.uint8], order: int) -> NDArray[np.int32]:
    """Integer value of the order-bit pattern starting at each position of the cyclic sequence."""
    n = len(bits)
    extended = np.concatenate([bits, bits[:order - 1]]).astype(np.int32)
    value = np.zeros(n, dtype=np.int32)
    for k in range(order):
        value += extended[k:k + n] << (order - 1 - k)
    return value


def cyclic_pattern_counts(bits: NDArray[np.uint8], order: int,
                          values: tuple[int, NDArray[np.int32]] | None = None) -> NDArray[np.int64]:
    """Counts of the n overlapping order-bit patterns of the cyclically extended sequence.

    With `values` = (k, order-k pattern values), k >= order, the order-bit pattern at a
    position is the leading bits of the longer pattern at that position (a right shift
    of the value), so one pass serves several orders.
    """
    if values is None:
        return np.bincount(cyclic_pattern_values(bits, order), minlength=2**order)
    longest, value = values
    return np.bincount(value >> (longest - order), minlength=2**order)


def _psi2(counts: NDArray[np.int64], order: int, n: int) -> float:
    c = counts.astype(np.float64)
    return float(2**order / n * np.sum(c * c) - n)


def _phi(counts: NDArray[np.int64], n: int) -> float:
    c = counts / n
    c = c[c > 0]
    return float(np.sum(c * np.log(c)))


def longest_runs_per_block(bits: NDArray[np.uint8], block: int) -> NDArray[np.int64]:
    """Longest run of ones in each complete block, from the index of the last reset."""
    count = len(bits) // block
    used = np.asarray(bits[:count * block])
    position = np.arange(count * block, dtype=np.int64)
    last_zero = np.maximum.accumulate(np.where(used == 0, position, -1))
    block_start = (position // block) * block
    run_length = position - np.maximum(last_zero, block_start - 1)
    return run_length.reshape(count, block).max(axis=1) if count else np.zeros(0, dtype=np.int64)


def frequency_test(bits: NDArray[np.uint8]) -> tuple[float, float]:
    """SP 800-22 section 2.1: (signed standardised imbalance, p)."""
    n, ones = len(bits), int(np.sum(bits, dtype=np.int64))
    z = (2 * ones - n) / math.sqrt(n)
    return z, float(2 * stats.norm.sf(abs(z)))


def block_frequency_test(bits: NDArray[np.uint8], block: int) -> tuple[float, float]:
    """SP 800-22 section 2.2: (chi-square, p) over the complete blocks."""
    blocks = len(bits) // block
    proportion = np.asarray(bits[:blocks * block]).reshape(blocks, block).sum(axis=1, dtype=np.int64) / block
    chi2 = float(4 * block * np.sum((proportion - 0.5)**2))
    return chi2, float(stats.chi2.sf(chi2, blocks))


def runs_test(bits: NDArray[np.uint8]) -> tuple[int, float, float | None, float | None]:
    """SP 800-22 section 2.3 without the prerequisite: (V, residual, signed z, p)."""
    n = len(bits)
    pi = int(np.sum(bits, dtype=np.int64)) / n
    transitions = 1 + int(np.count_nonzero(bits[1:] != bits[:-1]))
    residual = transitions - 2 * n * pi * (1 - pi)
    if not 0 < pi < 1:
        return transitions, residual, None, None
    z = residual / (2 * math.sqrt(n) * pi * (1 - pi))
    return transitions, residual, z, float(2 * stats.norm.sf(abs(z)))


def serial_test(bits: NDArray[np.uint8], m: int) -> dict[str, Any]:
    """SP 800-22 section 2.11: psi-squared at orders m-2, m-1, m and both p-values."""
    n = len(bits)
    values = (m, cyclic_pattern_values(bits, m))
    psi = [_psi2(cyclic_pattern_counts(bits, k, values), k, n) if k > 0 else 0.0 for k in (m - 2, m - 1, m)]
    d1, d2 = psi[2] - psi[1], psi[2] - 2 * psi[1] + psi[0]
    return {"psi2": psi, "delta1": d1, "delta2": d2,
            "p1": float(stats.chi2.sf(d1, 2**(m - 1))), "p2": float(stats.chi2.sf(d2, 2**(m - 2)))}


def approximate_entropy_test(bits: NDArray[np.uint8], m: int) -> dict[str, Any]:
    """SP 800-22 section 2.12: phi at orders m and m+1, chi-square and p."""
    n = len(bits)
    values = (m + 1, cyclic_pattern_values(bits, m + 1))
    phi = [_phi(cyclic_pattern_counts(bits, k, values), n) for k in (m, m + 1)]
    chi2 = 2 * n * (math.log(2) - (phi[0] - phi[1]))
    return {"phi": phi, "chi2": chi2, "p": float(stats.chi2.sf(chi2, 2**m))}


def cumulative_sums_forward_test(bits: NDArray[np.uint8]) -> tuple[int, float | None]:
    """SP 800-22 section 2.13, forward mode: (largest absolute partial sum, p).

    The summation limits use integer division that truncates toward zero, as the
    reference implementation does; math.trunc of the exact small quotient gives it.
    """
    n = len(bits)
    z = int(np.max(np.abs(np.cumsum(2 * np.asarray(bits, dtype=np.int64) - 1)))) if n else 0
    if z == 0:
        return z, None
    root = math.sqrt(n)
    upper = math.trunc((math.trunc(n / z) - 1) / 4)
    k1 = np.arange(math.trunc((math.trunc(-n / z) + 1) / 4), upper + 1)
    k2 = np.arange(math.trunc((math.trunc(-n / z) - 3) / 4), upper + 1)
    s1 = np.sum(special.ndtr((4 * k1 + 1) * z / root) - special.ndtr((4 * k1 - 1) * z / root))
    s2 = np.sum(special.ndtr((4 * k2 + 3) * z / root) - special.ndtr((4 * k2 + 1) * z / root))
    return z, float(1 - s1 + s2)


def longest_run_test(bits: NDArray[np.uint8]) -> dict[str, Any] | None:
    """SP 800-22 section 2.4 with the block size chosen from n as in section 2.4.2."""
    n = len(bits)
    block = 10000 if n >= 750000 else 128 if n >= 6272 else 8
    classes, lowest, probabilities = _LONGEST_RUN[block]
    blocks = n // block
    if blocks == 0:
        return None
    category = np.clip(longest_runs_per_block(bits, block), lowest, lowest + classes) - lowest
    nu = np.bincount(category, minlength=classes + 1)
    expected = blocks * np.asarray(probabilities)
    chi2 = float(np.sum((nu - expected)**2 / expected))
    return {"block": block, "blocks": blocks, "nu": nu.tolist(), "chi2": chi2, "p": float(stats.chi2.sf(chi2, classes))}


def bit_statistics(bits: NDArray[np.uint8]) -> dict[str, Any]:
    """The eight bit-level tests with the pilot's parameters and applicability rules.

    Block frequency uses M = max(128, n//99 + 1); serial and approximate entropy use
    m = 8. A p-value is None when its test is not applicable (n < 100, the block
    conditions, or the runs prerequisite |pi - 1/2| < 2/sqrt(n)).
    """
    b = np.asarray(bits, dtype=np.uint8)
    n = len(b)
    out: dict[str, Any] = {"n": n, "ones": int(b.sum(dtype=np.int64))}
    for key in ("frequency", "block_frequency", "runs", "serial_1", "serial_2", "approximate_entropy",
                "cumulative_sums", "longest_run", "frequency_signed_z", "block_frequency_chi2",
                "runs_residual", "runs_signed_z"):
        out[key] = None
    if n == 0:
        return out
    enough = n >= 100
    z, p = frequency_test(b)
    out["frequency_signed_z"], out["frequency"] = z, (p if enough else None)
    block = max(128, n // 99 + 1)
    if enough and block <= n and block > 0.01 * n and n // block < 100:
        out["block_frequency_chi2"], out["block_frequency"] = block_frequency_test(b, block)
    transitions, residual, runs_z, runs_p = runs_test(b)
    pi = out["ones"] / n
    out.update(V=transitions, pi=pi, runs_residual=residual, runs_signed_z=runs_z,
               prerequisite_failed=bool(abs(pi - 0.5) >= 2 / math.sqrt(n)))
    if enough and not out["prerequisite_failed"]:
        out["runs"] = runs_p
    if n >= 9:
        values = (9, cyclic_pattern_values(b, 9))
        counts = {k: cyclic_pattern_counts(b, k, values) for k in (6, 7, 8, 9)}
        psi = [_psi2(counts[k], k, n) for k in (6, 7, 8)]
        d1, d2 = psi[2] - psi[1], psi[2] - 2 * psi[1] + psi[0]
        phi8, phi9 = _phi(counts[8], n), _phi(counts[9], n)
        apen_chi2 = 2 * n * (math.log(2) - (phi8 - phi9))
        out.update(serial_psi2=psi, serial_delta1=d1, serial_delta2=d2,
                   serial_1=float(stats.chi2.sf(d1, 2**7)), serial_2=float(stats.chi2.sf(d2, 2**6)),
                   approximate_entropy_phi8=phi8, approximate_entropy_phi9=phi9,
                   approximate_entropy_chi2=apen_chi2, approximate_entropy=float(stats.chi2.sf(apen_chi2, 2**8)))
    out["cumulative_sums_max_abs"], out["cumulative_sums"] = cumulative_sums_forward_test(b)
    longest = longest_run_test(b)
    if longest is not None:
        out.update(longest_run_categories=longest["nu"], longest_run_chi2=longest["chi2"],
                   longest_run_block_size=longest["block"], longest_run=longest["p"])
    return out


def word_statistics(words: NDArray[np.uint64]) -> dict[str, Any]:
    """Lag-1 correlation of the 53-bit uniforms and adjacent identical words, full window."""
    w = np.asarray(words, dtype=np.uint64)
    count, pairs = len(w), len(w) - 1
    out: dict[str, Any] = {"word_count": count, "pair_count": pairs, "lag1_correlation": None,
                           "adjacent_duplicates": None, "lag1_r": None, "lag1_feature": None,
                           "duplicate_count": 0, "duplicate_rate": None}
    if pairs < 1:
        return out
    u = (w >> np.uint64(11)).astype(np.float64) * 2.0**-53
    a, c = u[:-1] - u[:-1].mean(), u[1:] - u[1:].mean()
    denominator = math.sqrt(float(np.dot(a, a)) * float(np.dot(c, c)))
    # P signal_schema lag1: unavailable when either paired vector has zero variance. The test is on
    # the values themselves: a rounding residue in the mean of a constant vector must not pass as variance.
    constant = np.ptp(u[:-1]) == 0 or np.ptp(u[1:]) == 0
    if pairs >= 2 and not constant and denominator > 0:
        r = float(np.dot(a, c)) / denominator
        out.update(lag1_r=r, lag1_test_z=math.sqrt(count) * r, lag1_feature=math.sqrt(pairs) * r,
                   lag1_correlation=float(2 * stats.norm.sf(abs(math.sqrt(count) * r))))
    duplicates = int(np.count_nonzero(np.diff(w) == 0))
    lam = pairs * 2.0**-64
    out.update(duplicate_count=duplicates, duplicate_rate=duplicates / pairs, duplicate_lambda=lam,
               adjacent_duplicates=1.0 if duplicates == 0 else max(float(special.gammainc(duplicates, lam)), P_FLOOR))
    return out


def battery_checkpoint(window: NDArray[np.uint64]) -> dict[str, Any]:
    """P operational_alarm.bit_cap: bit tests on the first 1,000,000 bits, word tests on all words."""
    data = words_to_bytes(window)
    bits = np.unpackbits(np.frombuffer(data, dtype=np.uint8), bitorder="big")[:BIT_CAP_BITS]
    out = {**bit_statistics(bits), **word_statistics(window)}
    out.update(original_bits=8 * len(data), bit_used_bits=len(bits), was_cut=8 * len(data) > BIT_CAP_BITS)
    out["p_values"] = [out[key] for key in P_KEYS]
    out["S"], out["available_p_count"] = alarm_statistic(out["p_values"])
    out["S_prime"], out["available_p_count_prime"] = alarm_statistic(out["p_values"][:8])
    return out


def alarm_statistic(p_values: Iterable[float | None]) -> tuple[float, int]:
    """P operational_alarm.statistic: -log10(max(min available p, floor)); 0 when none is available."""
    available = [float(p) for p in p_values if p is not None and math.isfinite(p) and 0.0 <= p <= 1.0]
    return (-math.log10(max(min(available), P_FLOOR)) if available else 0.0), len(available)


# ---------------------------------------------------------------------------
# Alarm calibration and issuance (P operational_alarm)
# ---------------------------------------------------------------------------
def order_statistic_threshold(maxima: Sequence[float], order: int = 950, required: int = 999) -> float:
    values = np.asarray(maxima, dtype=np.float64)
    if values.shape != (required,) or not np.isfinite(values).all():
        raise RecountError(f"exactly {required} finite calibration maxima are required")
    return float(np.partition(values, order - 1)[order - 1])


def first_alarm(statistic: Sequence[float], threshold: float, ends: NDArray[np.int64] = BATTERY_ENDS) -> int | None:
    """First strict exceedance in chronological order; uses no onset information."""
    s = np.asarray(statistic, dtype=np.float64)
    if s.shape != ends.shape:
        raise RecountError("one statistic per checkpoint is required")
    if not np.isfinite(s).all():
        raise RecountError("a non-finite checkpoint statistic leaves the first alarm unknown")
    hits = np.flatnonzero(s > threshold)
    return int(ends[hits[0]]) if len(hits) else None


def onset_class(alarm_sweep: int | None, t0: int | None) -> str | None:
    if alarm_sweep is None:
        return "no_alarm"
    if t0 is None:
        return None
    return "pre_onset" if alarm_sweep <= t0 else "post_onset"


def eligible(alarm_sweep: int | None, length: int) -> bool:
    return alarm_sweep is not None and alarm_sweep < length


def clopper_pearson(k: int, n: int, level: float = 0.95) -> list[float] | None:
    """Two-sided Clopper-Pearson interval from the binomial tail equations."""
    if n == 0:
        return None
    if not 0 <= k <= n:
        raise RecountError("0 <= k <= n required")
    tail = (1 - level) / 2
    solve = lambda f: float(optimize.brentq(f, 0.0, 1.0, xtol=1e-15, rtol=8.9e-16, maxiter=500))
    lower = 0.0 if k == 0 else solve(lambda p: stats.binom.sf(k - 1, n, p) - tail)
    upper = 1.0 if k == n else solve(lambda p: stats.binom.cdf(k, n, p) - tail)
    return [lower, upper]


# ---------------------------------------------------------------------------
# Endpoint (P endpoint)
# ---------------------------------------------------------------------------
def reference_energy(beads: int = 16, beta: float = 2.0) -> float:
    """Finite-bead primitive-action harmonic energy as a normal-mode sum.

    ln Z_M = -(1/2) sum_k ln(4 sin^2(pi k / M) + (beta/M)^2), so
    E_M = -d ln Z_M / d beta = sum_k (beta/M^2) / (4 sin^2(pi k / M) + (beta/M)^2).
    """
    return math.fsum((beta / beads**2) / (4 * math.sin(math.pi * k / beads)**2 + (beta / beads)**2)
                     for k in range(beads))


T975_15 = float(special.stdtrit(15, 0.975))


def endpoint(cv: FloatArray, n: int | None, reference: float) -> dict[str, Any]:
    """Sixteen-batch batch-means interval of cv[0:n], with exactly rounded sums.

    reason is None for a valid endpoint, otherwise one of NONFINITE_INPUT,
    INCOMPLETE_TRACE, EXCEPTION (fewer than 32 samples), DEGENERATE_BATCH_VARIANCE,
    UNCOVERED. margin = half_width - |estimate - reference|: positive inside the
    interval, and the distance to the decision boundary either way.
    """
    out: dict[str, Any] = {"n": n, "valid": False, "reason": None, "estimate": None, "signed_error": None,
                           "mcse": None, "half_width": None, "interval_95": None, "covered": None,
                           "batch_means": None, "margin": None}
    if n is None:
        return {**out, "valid": None, "reason": "NO_ALARM"}
    if n <= 0 or n % 16:
        raise RecountError("endpoint length must be a positive multiple of 16")
    values = np.asarray(cv, dtype=np.float64)
    prefix = values[:n] if values.ndim == 1 else values.reshape(-1)[:n]
    if not np.isfinite(prefix).all():
        return {**out, "reason": "NONFINITE_INPUT"}
    if values.ndim != 1 or len(prefix) != n:
        return {**out, "reason": "INCOMPLETE_TRACE"}
    size = n // 16
    data = prefix.tolist()
    means = [math.fsum(data[j * size:(j + 1) * size]) / size for j in range(16)]
    estimate = math.fsum(means) / 16
    out.update(batch_means=means, estimate=estimate, signed_error=estimate - reference)
    if n < 32:
        return {**out, "reason": "EXCEPTION"}
    variance = math.fsum((m - estimate)**2 for m in means) / 15
    if not math.isfinite(variance) or variance <= 0:
        return {**out, "reason": "DEGENERATE_BATCH_VARIANCE"}
    mcse = math.sqrt(variance / 16)
    half = T975_15 * mcse
    bounds = [estimate - half, estimate + half]
    covered = bounds[0] <= reference <= bounds[1]
    out.update(mcse=mcse, half_width=half, interval_95=bounds, covered=covered,
               margin=half - abs(estimate - reference), valid=covered, reason=None if covered else "UNCOVERED")
    return out


def issuance_length(alarm_sweep: int | None) -> int | None:
    return None if alarm_sweep is None else 16 * (alarm_sweep // 16)


# ---------------------------------------------------------------------------
# Forecaster features and models (P forecasters)
# ---------------------------------------------------------------------------
def trailing_slope(ends: Sequence[int], values: Sequence[float | None], alarm_index: int) -> tuple[float | None, int]:
    """P forecasters.slope: least squares over the last min(4, j+1) checkpoints, exact rationals."""
    first = max(0, alarm_index - 3)
    x = [int(e) for e in ends[first:alarm_index + 1]]
    v = list(values[first:alarm_index + 1])
    count = len(x)
    if any(value is None or not math.isfinite(value) for value in v):
        return None, count
    if count == 1:
        return 0.0, 1
    xf, vf = [Fraction(e) for e in x], [Fraction(float(value)) for value in v]
    xm, vm = sum(xf) / count, sum(vf) / count
    return float(sum((a - xm) * (b - vm) for a, b in zip(xf, vf)) / sum((a - xm)**2 for a in xf)), count


def application_index(alarm_sweep: int, ends: NDArray[np.int64] = APPLICATION_ENDS) -> int | None:
    """Last application checkpoint with end <= alarm sweep, or None (alarms at 512 and 768)."""
    position = int(np.searchsorted(ends, alarm_sweep, side="right")) - 1
    return position if position >= 0 else None


def standardise(x_fit: FloatArray) -> tuple[FloatArray, FloatArray]:
    """P forecasters.standardisation: fitting-row mean and population SD; scale 1 for constants."""
    x = np.asarray(x_fit, dtype=np.float64)
    mu = x.mean(axis=0) if x.size else np.zeros(x.shape[1])
    sd = x.std(axis=0, ddof=0) if x.size else np.ones(x.shape[1])
    return mu, np.where(sd == 0, 1.0, sd)


def logistic_objective(theta: FloatArray, z: FloatArray, y: FloatArray) -> tuple[float, FloatArray]:
    """P forecasters.logistic.objective and gradient: summed log-loss + ||w||^2/2, free intercept."""
    w, b = theta[:-1], theta[-1]
    eta = z @ w + b
    value = float(np.sum(np.logaddexp(0.0, eta) - y * eta) + 0.5 * np.dot(w, w))
    residual = special.expit(eta) - y
    return value, np.concatenate([z.T @ residual + w, [np.sum(residual)]])


def fit_logistic_newton(z: FloatArray, y: FloatArray, tolerance: float = 1e-10, iterations: int = 200) -> dict[str, Any]:
    """Minimiser of the frozen objective by damped Newton steps (reference optimum, not the pilot fit)."""
    z, y = np.asarray(z, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if z.ndim != 2 or len(z) != len(y) or len(y) == 0 or y.min() == y.max():
        raise RecountError("both classes and a 2-D design matrix are required")
    d = z.shape[1]
    design = np.column_stack([z, np.ones(len(y))])
    penalty = np.diag(np.r_[np.ones(d), 0.0])
    theta = np.zeros(d + 1)
    value, gradient = logistic_objective(theta, z, y)
    for iteration in range(iterations):
        if np.max(np.abs(gradient)) <= tolerance:
            return {"theta": theta, "objective": value, "gradient_inf": float(np.max(np.abs(gradient))),
                    "iterations": iteration, "converged": True}
        p = special.expit(design @ theta)
        hessian = design.T @ (design * (p * (1 - p))[:, None]) + penalty
        step = np.linalg.solve(hessian + 1e-12 * np.eye(d + 1), gradient)
        scale = 1.0
        while True:
            candidate = theta - scale * step
            new_value, new_gradient = logistic_objective(candidate, z, y)
            if new_value <= value + 1e-12 * abs(value) or scale < 1e-12:
                break
            scale /= 2
        theta, value, gradient = candidate, new_value, new_gradient
    return {"theta": theta, "objective": value, "gradient_inf": float(np.max(np.abs(gradient))),
            "iterations": iterations, "converged": bool(np.max(np.abs(gradient)) <= tolerance)}


def logistic_predict(theta: Sequence[float], mu: Sequence[float], scale: Sequence[float], x: FloatArray) -> FloatArray:
    theta = np.asarray(theta, dtype=np.float64)
    z = (np.asarray(x, dtype=np.float64) - np.asarray(mu)) / np.asarray(scale)
    return special.expit(z @ theta[:-1] + theta[-1])


# ---------------------------------------------------------------------------
# Scores (P scores)
# ---------------------------------------------------------------------------
def brier(q: Sequence[float], y: Sequence[int]) -> float | None:
    q, y = np.asarray(q, dtype=np.float64), np.asarray(y, dtype=np.float64)
    if len(y) == 0 or not np.isfinite(q).all():
        return None
    return float(np.mean((q - y)**2))


def skill(score: float | None, baseline: float | None) -> float | None:
    if score is None or baseline is None or baseline == 0:
        return None
    return 1 - score / baseline


def calibration_bins(q: Sequence[float], y: Sequence[int]) -> list[dict[str, Any]]:
    """Five bins [0,.2), [.2,.4), [.4,.6), [.6,.8), [.8,1]; q = 1 falls in the last bin."""
    q, y = np.asarray(q, dtype=np.float64), np.asarray(y, dtype=np.int64)
    if len(q) and (not np.isfinite(q).all() or q.min() < 0 or q.max() > 1):
        raise RecountError("calibration needs finite probabilities in [0, 1]")
    # Bin index = number of inner edges at or below q, compared against the stated edge values.
    index = np.array([sum(1 for edge in (0.2, 0.4, 0.6, 0.8) if value >= edge) for value in q], dtype=np.int64)
    rows = []
    for b in range(5):
        inside = index == b
        count, hits = int(inside.sum()), int(y[inside].sum())
        rows.append({"bin": b, "count": count, "mean_predicted": float(q[inside].mean()) if count else None,
                     "observed_fraction": hits / count if count else None, "cp95": clopper_pearson(hits, count)})
    return rows


def bootstrap_indices(seed: int, n: int, resamples: int = 10000) -> NDArray[np.int64]:
    """P scores.bootstrap.exact_draws. One block draw equals `resamples` successive draws of size n
    for this NumPy generator; tests/v019/test_audit_recount.py asserts that equality."""
    return np.random.default_rng(seed).integers(0, n, size=(resamples, n))


def paired_bootstrap(y: Sequence[int], predictions: Mapping[str, Sequence[float]], seed: int,
                     baseline: str = "F0+", resamples: int = 10000,
                     criterion_level: float = 1 - 0.05 / 3, other_level: float = 0.95) -> dict[str, Any]:
    """Percentile intervals over trajectories for every model's Brier score and skill.

    Skill intervals are null as soon as one resample has a zero baseline score
    (P scores.bootstrap.zero_denominator); the count of such resamples is reported.
    """
    y = np.asarray(y, dtype=np.float64)
    n = len(y)
    if n == 0:
        return {"n": 0, "models": {}}
    index = bootstrap_indices(seed, n, resamples)
    loss = {m: (np.asarray(q, dtype=np.float64) - y)**2 for m, q in predictions.items()}
    score = {m: values[index].mean(axis=1) for m, values in loss.items()}
    undefined = int(np.count_nonzero(score[baseline] == 0)) if baseline in score else resamples

    def interval(values: FloatArray, level: float) -> list[float]:
        return np.quantile(values, [(1 - level) / 2, 1 - (1 - level) / 2], method="linear").tolist()

    out: dict[str, Any] = {"n": n, "seed": seed, "resamples": resamples, "undefined_skill_resamples": undefined,
                           "models": {}}
    for m in predictions:
        row = {"brier": float(loss[m].mean()), "brier_ci95": interval(score[m], other_level)}
        row["skill"] = skill(row["brier"], float(loss[baseline].mean()) if baseline in loss else None)
        if undefined == 0:
            values = 1 - score[m] / score[baseline]
            row["skill_ci95"], row["skill_ci_criterion"] = interval(values, other_level), interval(values, criterion_level)
        else:
            row["skill_ci95"], row["skill_ci_criterion"] = None, None
        out["models"][m] = row
    return out


# ---------------------------------------------------------------------------
# Oracle reference (P oracle)
# ---------------------------------------------------------------------------
def slope_through_origin(x: Sequence[float], e: Sequence[float]) -> float | None:
    x, e = np.asarray(x, dtype=np.float64), np.asarray(e, dtype=np.float64)
    denominator = math.fsum((x * x).tolist())
    return math.fsum((x * e).tolist()) / denominator if denominator > 0 else None


def sigma_eff(length: int, mcse: Sequence[float]) -> float:
    return math.sqrt(length) * math.fsum(mcse) / len(mcse)


def validity_probability_frozen(delta: float) -> float:
    """P oracle.probability exactly as frozen. May be NaN: the caller must then report it unavailable."""
    c = stats.t.ppf(0.975, df=15)
    return float(stats.nct.cdf(c, df=15, nc=delta) - stats.nct.cdf(-c, df=15, nc=delta))


_GL_NODES, _GL_WEIGHTS = np.polynomial.legendre.leggauss(200)


def validity_probability_stable(delta: float) -> float:
    """The same probability, P(|T| <= c) for T ~ noncentral t(15, delta), by quadrature.

    T = (Z + delta) / W with Z standard normal and W = sqrt(V/15), V chi-square(15), hence
    P = E_W[Phi(c W - |delta|) - Phi(-c W - |delta|)] (the probability is even in delta).
    Gauss-Legendre panels of width at most 0.75 cover [0, |delta|/c + 6]; the density of W
    beyond that is below exp(-270) of its value at |delta|/c. Against a 70-digit reference
    (research/v019-audit/f01_oracle.py) the absolute error is below 3e-15 everywhere and
    the relative error below 1e-13 down to 1e-306.
    """
    c, d = T975_15, abs(float(delta))
    top = d / c + 6.0
    edges = np.linspace(0.0, top, max(8, int(math.ceil(top / 0.75))) + 1)
    log_constant = math.log(2.0) + 7.5 * math.log(7.5) - math.lgamma(7.5)
    total = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        w = (a + b) / 2 + (b - a) / 2 * _GL_NODES
        density = np.exp(log_constant + 14 * np.log(w) - 7.5 * w * w)
        inside = special.ndtr(c * w - d) - special.ndtr(-c * w - d)
        total += float(np.sum(_GL_WEIGHTS * inside * density)) * (b - a) / 2
    return float(min(max(total, 0.0), 1.0))


# ---------------------------------------------------------------------------
# Stop criteria (P stop_criteria)
# ---------------------------------------------------------------------------
def criterion_1(covered_by_length: Mapping[int, int], runs: int = 2000) -> dict[str, Any]:
    rows, labels = {}, []
    for length, k in covered_by_length.items():
        lower, upper = clopper_pearson(k, runs, 0.95)
        label = ("PASS" if lower >= 0.92 and upper <= 0.98 else
                 "OUTSIDE_BAND" if upper < 0.92 or lower > 0.98 else "STRADDLES_BAND")
        rows[length] = {"covered": k, "runs": runs, "cp95": [lower, upper], "label": label}
        labels.append(label)
    overall = ("PASS" if all(v == "PASS" for v in labels) else
               "FAIL" if "OUTSIDE_BAND" in labels else "INSUFFICIENT_INFORMATION")
    return {"per_length": rows, "result": overall, "stops_D": "OUTSIDE_BAND" in labels}


def _family_vote(passes: Mapping[str, bool | None], needed: int) -> str:
    """PASS with `needed` passing families; INSUFFICIENT if unresolved ones could still reach it."""
    passed = sum(1 for v in passes.values() if v is True)
    unresolved = sum(1 for v in passes.values() if v is None)
    return "PASS" if passed >= needed else "INSUFFICIENT_INFORMATION" if passed + unresolved >= needed else "FAIL"


def criteria_2_3(per_family: Mapping[str, Mapping[str, Any]], length: int = 20480) -> dict[str, Any]:
    """per_family[f] = {alarm: [...], valid: [...], already_invalid: [...]} over all planned fault runs.

    alarm entries are sweeps or None (no alarm); an entry "UNKNOWN" marks unknown issuance.
    """
    c2, c3, rows = {}, {}, {}
    for family, data in per_family.items():
        alarm = list(data["alarm"])
        unknown = sum(1 for a in alarm if a == "UNKNOWN")
        chosen = [i for i, a in enumerate(alarm) if a not in (None, "UNKNOWN") and a < length]
        n = len(chosen)
        p = sum(bool(data["valid"][i]) for i in chosen) / n if n else None
        a = sum(bool(data["already_invalid"][i]) for i in chosen) / n if n else None
        median = float(np.median([length - alarm[i] for i in chosen])) if n else None
        sufficient = n >= 300 and unknown == 0
        c2[family] = (0.2 <= p <= 0.8 and a < 0.5) if sufficient else None
        c3[family] = (median >= 5120) if sufficient else None
        rows[family] = {"eligible": n, "unknown_issuance": unknown, "fraction_valid": p,
                        "fraction_already_invalid": a, "median_remaining": median,
                        "criterion_2": c2[family], "criterion_3": c3[family]}
    return {"per_family": rows, "criterion_2": _family_vote(c2, 2), "criterion_3": _family_vote(c3, 2)}


def criterion_4(per_family: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """per_family[f] = paired_bootstrap(...) output for eligible evaluation runs at L = 20480."""
    votes, rows = {}, {}
    for family, result in per_family.items():
        model = result.get("models", {}).get("F-B-lite")
        base = result.get("models", {}).get("F0+")
        n = result.get("n", 0)
        if n < 240 or model is None or base is None or not base["brier"] or model["skill_ci_criterion"] is None:
            votes[family] = None
        else:
            votes[family] = bool(model["skill"] >= 0.10 and model["skill_ci_criterion"][0] > 0)
        rows[family] = {"eligible_evaluation": n, "skill": None if model is None else model["skill"],
                        "interval": None if model is None else model["skill_ci_criterion"], "pass": votes[family]}
    return {"per_family": rows, "result": _family_vote(votes, 1)}


def gate_v2d(alarms: int, runs: int = 2000) -> dict[str, Any]:
    interval = clopper_pearson(alarms, runs, 0.95)
    return {"alarms": alarms, "runs": runs, "cp95": interval, "pass": bool(interval[1] <= 0.10)}


def _calendar_date(value: Any) -> datetime.date | None:
    """A date written as YYYY-MM-DD, or None when the text is not one."""
    if not isinstance(value, str) or not re.fullmatch(r"\d{4}-\d{2}-\d{2}", value, flags=re.ASCII):
        return None
    try:
        return datetime.date.fromisoformat(value)
    except ValueError:
        return None


def months_earlier(day: datetime.date, months: int) -> datetime.date:
    """Calendar months back, clamped to the last day of the target month (P stop_criteria.5.calculation)."""
    year, month0 = divmod(day.year * 12 + day.month - 1 - months, 12)
    return datetime.date(year, month0 + 1, min(day.day, calendar.monthrange(year, month0 + 1)[1]))


def manifest_date(written_at: Any, zone: str) -> datetime.date | None:
    """Local calendar date, in the protocol's time zone, of a manifest's written_at stamp."""
    if not isinstance(written_at, str):
        return None
    try:
        stamp = datetime.datetime.fromisoformat(written_at)
    except ValueError:
        return None
    return stamp.astimezone(zoneinfo.ZoneInfo(zone)).date() if stamp.tzinfo else stamp.date()


def criterion_5(tdk_deadline: Any, engineering_complete_on: Any, results_available_on: Any,
                manifest_written_at: Any = None, protocol: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """P stop_criteria.5 and the last sentence of implementation_checks.recount.

    The label follows the two deadlines only: engineering on or before the approval date plus six
    weeks, results on or before three calendar months before the declared deadline. A missing date
    (None or an empty text) gives INSUFFICIENT_INFORMATION. A date that is not a calendar date gives no label. The order
    of the freeze steps (engineering, then the pre-generation manifest, then generation) is checked
    against the manifest's date and reported in date_problems; it does not change the label.
    """
    protocol = load_protocol() if protocol is None else protocol
    rule = protocol["stop_criteria"]["5"]
    engineering_deadline = datetime.date.fromisoformat(protocol["approval"]["date"]) + datetime.timedelta(weeks=6)
    given = {"tdk_deadline": tdk_deadline, "engineering_complete_on": engineering_complete_on,
             "results_available_on": results_available_on}
    dates = {name: _calendar_date(value) for name, value in given.items()}
    problems = [f"{name} is not a calendar date (YYYY-MM-DD): {value!r}"
                for name, value in given.items() if value not in (None, "") and dates[name] is None]
    written = manifest_date(manifest_written_at, rule["timezone"])
    if manifest_written_at is not None and written is None:
        problems.append(f"the manifest's written_at is not a time stamp: {manifest_written_at!r}")
    out: dict[str, Any] = {**given, "engineering_deadline": engineering_deadline.isoformat(), "results_deadline": None,
                           "manifest_date": None if written is None else written.isoformat(), "timezone": rule["timezone"]}
    if dates["tdk_deadline"] is not None:
        out["results_deadline"] = months_earlier(dates["tdk_deadline"], 3).isoformat()
    if written is not None:
        if dates["engineering_complete_on"] is not None and dates["engineering_complete_on"] > written:
            problems.append("engineering_complete_on is later than the pre-generation manifest (freeze steps 4 and 6)")
        if dates["results_available_on"] is not None and dates["results_available_on"] < written:
            problems.append("results_available_on is earlier than the pre-generation manifest (freeze steps 6 and 7)")
    if any(value in (None, "") for value in given.values()):
        result: str | None = "INSUFFICIENT_INFORMATION"
    elif any(day is None for day in dates.values()):
        result = None
    else:
        on_time = (dates["engineering_complete_on"] <= engineering_deadline
                   and dates["results_available_on"] <= months_earlier(dates["tdk_deadline"], 3))
        result = "PASS" if on_time else "FAIL"
    return {**out, "result": result, "date_problems": problems}


def compare(primary: Any, independent: Any, rtol: float = 1e-12, atol: float = 1e-12) -> bool:
    """P implementation_checks.recount: exact for flags, counts, labels and null patterns; isclose for floats."""
    if isinstance(primary, (bool, np.bool_)) or isinstance(independent, (bool, np.bool_)):
        return isinstance(primary, (bool, np.bool_)) and isinstance(independent, (bool, np.bool_)) and bool(primary) == bool(independent)
    if primary is None or independent is None:
        return primary is None and independent is None
    if isinstance(primary, (str, bytes)) or isinstance(independent, (str, bytes)):
        return primary == independent
    if isinstance(primary, Mapping) and isinstance(independent, Mapping):
        return primary.keys() == independent.keys() and all(compare(primary[k], independent[k], rtol, atol) for k in primary)
    if isinstance(primary, (list, tuple, np.ndarray)) or isinstance(independent, (list, tuple, np.ndarray)):
        a, b = list(primary), list(independent)
        return len(a) == len(b) and all(compare(x, y, rtol, atol) for x, y in zip(a, b))
    if isinstance(primary, (int, np.integer)) and isinstance(independent, (int, np.integer)):
        return int(primary) == int(independent)
    return bool(np.isclose(float(primary), float(independent), rtol=rtol, atol=atol, equal_nan=False))
