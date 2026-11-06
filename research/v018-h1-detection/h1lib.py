"""Stage C (H1) library: conventional tests, fault adapters, sweep clock, runs and arm statistics.

RNGLab v0.18 (written by Claude Code, 2026-11-06) for the frozen protocol
rnglab-v018-h1-detection-v1 (research/v018-h1-detection/protocol.json). Section
numbers below refer to PROTOCOL.md and to NIST SP 800-22 Rev 1a.
"""
from __future__ import annotations

import hashlib
import math
import struct
import sys
from pathlib import Path
from statistics import NormalDist

import numpy as np
from scipy.special import gammaincc
from scipy.stats import norm, poisson

REPO = Path(__file__).resolve().parents[2]
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))

from rnglab.seeds import namespace_seed  # noqa: E402
from rnglab11.statistics import HealthTests, evaluate_bits  # noqa: E402
from rnglab12.uncertainty import summarize_series  # noqa: E402
from rnglab14.models.pimc_harmonic import (  # noqa: E402
    PathIntegralHarmonicOscillator, exact_finite_bead_energy, sample_exact_finite_bead_paths,
)
from rnglab14.sources import PCG64Source  # noqa: E402
from rnglab14.sources.base import RandomSource, SourceManifest  # noqa: E402

M, BETA, DELTA, SWEEPS = 16, 2.0, 0.4, 10240
ONSET_BYTES = 2**19
HARM_START, HARM_END, HARM_BATCH = 2560, 10240, 480
CONV_WINDOW, APP_WINDOW, STRIDE = 512, 1024, 256
PROPOSAL_BYTES_PER_SWEEP = 8 * M
FAMILIES = ("healthy", "SHIFT", "LAG1", "DUP", "GFSR")
GFSR_PAIRS = ((31, 13), (89, 38), (127, 63), (250, 103), (521, 168), (1279, 418))
TOP = 1 - 2.0**-53
STD = NormalDist()
P_FLOOR = 1e-300

# --------------------------------------------------------------------------
# SP 800-22 Rev 1a tests added for Stage C (sections 2.4, 2.11, 2.12, 2.13).
# Pi values for the longest-run test are the published section 3.4 tables.
# --------------------------------------------------------------------------
LONGEST_RUN_TABLES = {
    8: (3, (1, 2, 3, 4), (0.2148, 0.3672, 0.2305, 0.1875)),
    128: (5, (4, 5, 6, 7, 8, 9), (0.1174, 0.2430, 0.2493, 0.1752, 0.1027, 0.1124)),
    10000: (6, (10, 11, 12, 13, 14, 15, 16), (0.0882, 0.2092, 0.2483, 0.1933, 0.1208, 0.0675, 0.0727)),
}


def pattern_counts(bits, m_max: int) -> dict[int, np.ndarray]:
    """Counts of the n overlapping m-bit patterns of the cyclically augmented sequence, m = 0..m_max.

    Computed once at m_max; with cyclic augmentation the (m-1)-bit pattern at
    position i is the prefix of the m-bit pattern at i, so lower orders follow
    exactly by summing over the last bit.
    """
    if not 1 <= m_max <= 16:
        raise ValueError("1 <= m_max <= 16 required")
    b = np.asarray(bits, dtype=np.uint16)
    n = len(b)
    extended = np.concatenate([b, b[:m_max - 1]])
    values = np.zeros(n, dtype=np.uint16)
    for k in range(m_max):
        values = (values << 1) | extended[k:k + n]
    counts = {m_max: np.bincount(values, minlength=2**m_max)}
    for m in range(m_max, 0, -1):
        counts[m - 1] = counts[m][0::2] + counts[m][1::2]
    return counts


def serial_test(bits, m: int = 8, counts: dict | None = None) -> tuple[float, float]:
    """Section 2.11: (P-value1, P-value2). psi^2 of order <= 0 is zero."""
    n = len(bits)
    counts = counts or pattern_counts(bits, m)

    def psi2(k):
        return float(2**k / n * np.sum(counts[k].astype(float)**2) - n) if k > 0 else 0.0

    p_m, p_m1, p_m2 = psi2(m), psi2(m - 1), psi2(m - 2)
    return (float(gammaincc(2**(m - 2), (p_m - p_m1) / 2)),
            float(gammaincc(2**(m - 3), (p_m - 2 * p_m1 + p_m2) / 2)))


def approximate_entropy_test(bits, m: int = 8, counts: dict | None = None) -> float:
    """Section 2.12."""
    n = len(bits)
    counts = counts or pattern_counts(bits, m + 1)

    def phi(k):
        c = counts[k] / n
        c = c[c > 0]
        return float(np.sum(c * np.log(c)))

    apen = phi(m) - phi(m + 1)
    chi2 = 2 * n * (math.log(2) - apen)
    return float(gammaincc(2**(m - 1), chi2 / 2))


def _tdiv(a: int, b: int) -> int:
    """C integer division (truncation toward zero), as in the NIST reference code."""
    q = abs(a) // abs(b)
    return q if (a >= 0) == (b > 0) else -q


def cumulative_sums_forward(bits) -> float:
    """Section 2.13, mode 0. Summation limits use the reference code's integer
    arithmetic; with floor instead, the 10-bit worked example gives 0.4115847
    rather than the published 0.4116588."""
    x = 2 * np.asarray(bits, dtype=np.int64) - 1
    n = len(x)
    z = int(np.max(np.abs(np.cumsum(x))))
    sqn = math.sqrt(n)
    k1 = np.arange(_tdiv(_tdiv(-n, z) + 1, 4), _tdiv(_tdiv(n, z) - 1, 4) + 1)
    k2 = np.arange(_tdiv(_tdiv(-n, z) - 3, 4), _tdiv(_tdiv(n, z) - 1, 4) + 1)
    s1 = np.sum(norm.cdf((4 * k1 + 1) * z / sqn) - norm.cdf((4 * k1 - 1) * z / sqn))
    s2 = np.sum(norm.cdf((4 * k2 + 3) * z / sqn) - norm.cdf((4 * k2 + 1) * z / sqn))
    return float(1 - s1 + s2)


def longest_run_block_size(n: int) -> int:
    return 10000 if n >= 750000 else 128 if n >= 6272 else 8


def longest_run_test(bits) -> tuple[float, list[int]]:
    """Section 2.4: (P-value, category counts nu)."""
    b = np.asarray(bits, dtype=np.int8)
    block = longest_run_block_size(len(b))
    k, edges, pi = LONGEST_RUN_TABLES[block]
    count = len(b) // block
    rows = b[:count * block].reshape(count, block)
    padded = np.concatenate([np.zeros((count, 1), np.int8), rows, np.zeros((count, 1), np.int8)], axis=1).ravel()
    change = np.diff(padded.astype(np.int16))
    starts, ends = np.flatnonzero(change == 1), np.flatnonzero(change == -1)
    longest = np.zeros(count, dtype=np.int64)
    np.maximum.at(longest, starts // (block + 2), ends - starts)
    categories = np.clip(longest, edges[0], edges[-1])
    nu = [int(np.sum(categories == e)) for e in edges]
    chi2 = sum((nu[i] - count * pi[i])**2 / (count * pi[i]) for i in range(k + 1))
    return float(gammaincc(k / 2, chi2 / 2)), nu


def word_tests(data: bytes) -> dict[str, float]:
    """Lag-1 serial correlation of the 53-bit uniforms and adjacent identical 64-bit words."""
    words = np.frombuffer(data, dtype=">u8")
    u = (words >> np.uint64(11)).astype(float) * 2.0**-53
    n = len(u)
    r = float(np.corrcoef(u[:-1], u[1:])[0, 1])
    duplicates = int(np.sum(words[1:] == words[:-1]))
    lam = (n - 1) * 2.0**-64
    return {"lag1_correlation": float(math.erfc(abs(math.sqrt(n) * r) / math.sqrt(2))),
            "adjacent_duplicates": 1.0 if duplicates == 0 else max(float(poisson.sf(duplicates - 1, lam)), P_FLOOR),
            "lag1_r": r, "duplicate_count": duplicates}


def conventional_window(data: bytes) -> dict[str, float | None]:
    """All p-values of arm C for one window (PROTOCOL.md section 2.5)."""
    bits = np.unpackbits(np.frombuffer(data, dtype=np.uint8))
    n = len(bits)
    base = evaluate_bits(bits, block_size=max(128, n // 99 + 1))["tests"]
    counts = pattern_counts(bits, 9)
    serial = serial_test(bits, 8, counts)
    words = word_tests(data)
    return {"frequency": base["frequency"]["p_value"], "block_frequency": base["block_frequency"]["p_value"],
            "runs": base["runs"]["p_value"], "serial_1": serial[0], "serial_2": serial[1],
            "approximate_entropy": approximate_entropy_test(bits, 8, counts), "cumulative_sums": cumulative_sums_forward(bits),
            "longest_run": longest_run_test(bits)[0], "lag1_correlation": words["lag1_correlation"],
            "adjacent_duplicates": words["adjacent_duplicates"]}


CONV_P_KEYS = ("frequency", "block_frequency", "runs", "serial_1", "serial_2", "approximate_entropy",
               "cumulative_sums", "longest_run", "lag1_correlation", "adjacent_duplicates")


def health_test_events(symbols: np.ndarray, tests: HealthTests) -> tuple[np.ndarray, np.ndarray]:
    """Vectorized RCT and APT event sample numbers (1-based), matching HealthTests.update."""
    x = np.asarray(symbols, dtype=np.int64)
    change = np.flatnonzero(np.diff(x) != 0) + 1
    run_starts = np.concatenate([[0], change])
    run_lengths = np.diff(np.concatenate([run_starts, [len(x)]]))
    long = run_lengths >= tests.rct_cutoff
    rct = run_starts[long] + tests.rct_cutoff
    w = tests.window
    apt = []
    for start in range(0, len(x), w):
        window = x[start:start + w]
        hits = np.cumsum(window == window[0])
        reached = np.flatnonzero(hits >= tests.apt_cutoff)
        if len(reached):
            apt.append(start + int(reached[0]) + 1)
    return rct, np.asarray(apt, dtype=np.int64)


# --------------------------------------------------------------------------
# GFSR and GF(2) irreducibility (PROTOCOL.md section 2.3)
# --------------------------------------------------------------------------
def _gf2_square(a: int) -> int:
    return int(bin(a)[2:].replace("", "0")[1:] or "0", 2) if a else 0


def _gf2_mod_trinomial(a: int, p: int, q: int) -> int:
    mask = (1 << p) - 1
    while a >> p:
        high = a >> p
        a = (a & mask) ^ (high << q) ^ high
    return a


def _gf2_mod(a: int, f: int) -> int:
    df = f.bit_length() - 1
    while a and a.bit_length() - 1 >= df:
        a ^= f << (a.bit_length() - 1 - df)
    return a


def _gf2_gcd(a: int, b: int) -> int:
    while b:
        a, b = b, _gf2_mod(a, b)
    return a


def _prime_factors(n: int) -> list[int]:
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    return out + ([n] if n > 1 else [])


def trinomial_irreducible(p: int, q: int) -> bool:
    """Rabin's test for x^p + x^q + 1 over GF(2)."""
    f = (1 << p) | (1 << q) | 1
    x = 2

    def x_power_2k(k):
        a = x
        for _ in range(k):
            a = _gf2_mod_trinomial(_gf2_square(a), p, q)
        return a

    if x_power_2k(p) != x:
        return False
    return all(_gf2_gcd(f, x_power_2k(p // r) ^ x) == 1 for r in _prime_factors(p))


class GFSR:
    """x_n = x_{n-p} XOR x_{n-q} on 32-bit words; state from numpy default_rng(seed); 10p warm-up outputs."""

    def __init__(self, p: int, q: int, seed: int):
        self.p, self.q = p, q
        self.state = [int(v) for v in np.random.default_rng(seed).integers(0, 2**32, size=p, dtype=np.uint64)]
        self.i = 0
        for _ in range(10 * p):
            self.next()

    def next(self) -> int:
        p, q, i = self.p, self.q, self.i
        value = self.state[i] ^ self.state[(i + p - q) % p]
        self.state[i] = value
        self.i = (i + 1) % p
        return value

    def word(self) -> bytes:
        return struct.pack(">Q", (self.next() << 32) | self.next())


# --------------------------------------------------------------------------
# Fault adapter and sweep clock (PROTOCOL.md sections 2.2-2.3)
# --------------------------------------------------------------------------
class FaultyWords(RandomSource):
    """Local source wrapper: records every emitted byte; transforms 8-byte words from the onset offset."""

    def __init__(self, base: RandomSource, family: str, severity, fault_seed: int, gfsr_seed: int,
                 onset_bytes: int = ONSET_BYTES):
        if family not in FAMILIES:
            raise ValueError(f"unknown family {family!r}")
        super().__init__(SourceManifest(source_id=f"h1_{family}_{severity}", family="synthetic_fault" if family != "GFSR" else "weak_generator",
                                        is_deterministic=True, is_physical=False,
                                        metadata={"base": base.manifest.source_id, "family": family, "severity": str(severity),
                                                  "onset_bytes": onset_bytes}))
        self.base, self.family, self.severity, self.onset = base, family, severity, onset_bytes
        self.emitted = bytearray()
        self.fault_rng = np.random.default_rng(fault_seed)
        self.gfsr_seed = gfsr_seed
        self.gfsr = None
        self.prev: bytes | None = None
        self.z: float | None = None

    def _transform(self, raw: bytes) -> bytes:
        if self.family == "DUP":
            duplicate = self.fault_rng.random() < self.severity
            return self.prev if duplicate and self.prev is not None else raw
        word = struct.unpack(">Q", raw)[0]
        u = (word >> 11) * 2.0**-53
        if self.family == "SHIFT":
            u = min(u + self.severity, TOP)
        else:  # LAG1
            e = STD.inv_cdf(u + 2.0**-54)
            rho = self.severity
            self.z = e if self.z is None else rho * self.z + math.sqrt(1 - rho * rho) * e
            u = min(STD.cdf(self.z), TOP)
        return struct.pack(">Q", (int(u * 2.0**53) << 11) | (word & 0x7FF))

    def _read_bytes(self, n: int) -> bytes:
        post = self.family != "healthy" and len(self.emitted) >= self.onset
        if post and n != 8:
            raise ValueError("post-onset reads must be single 8-byte words")
        if post and self.family == "GFSR":
            if self.gfsr is None:
                self.gfsr = GFSR(*self.severity, self.gfsr_seed)
            raw = self.gfsr.word()
        else:
            raw = self.base.read_bytes(n)
            if post:
                raw = self._transform(raw)
        if n == 8:
            self.prev = raw
        self.emitted += raw
        return raw


class SweepClock(RandomSource):
    """Healthy centroid source wrapper; records the local consumed bytes at each sweep boundary."""

    def __init__(self, base: RandomSource, local: RandomSource):
        super().__init__(SourceManifest(source_id=f"clock_{base.manifest.source_id}", family="reference_prng",
                                        is_deterministic=True, is_physical=False))
        self.base, self.local, self.marks = base, local, []

    def _read_bytes(self, n: int) -> bytes:
        self.marks.append(self.local.consumed_bytes)
        return self.base.read_bytes(n)


# --------------------------------------------------------------------------
# One run (PROTOCOL.md sections 2-3)
# --------------------------------------------------------------------------
def seeds_for(ns: str, condition: str, i: int) -> dict[str, int]:
    unit = f"{i:05d}"
    return {role: namespace_seed(ns, condition, role, unit) for role in ("local", "centroid", "start", "fault", "gfsr")}


def window_bytes(marks: np.ndarray, start: int, end: int) -> tuple[int, int]:
    """Byte range of sweeps [start, end): [B(start-1), B(end-1)) with B(-1) = 0."""
    return (0 if start == 0 else int(marks[start - 1])), int(marks[end - 1])


def run_one(ns: str, condition: str, family: str, severity, i: int, sweeps: int = SWEEPS,
            onset_bytes: int = ONSET_BYTES) -> tuple[dict, dict]:
    seeds = seeds_for(ns, condition, i)
    local = FaultyWords(PCG64Source(seeds["local"]), family, severity, seeds["fault"], seeds["gfsr"], onset_bytes)
    clock = SweepClock(PCG64Source(seeds["centroid"]), local)
    start = sample_exact_finite_bead_paths(PCG64Source(seeds["start"]), M, BETA, 2)[0]
    model = PathIntegralHarmonicOscillator(beads=M, temperature=1 / BETA, displacement_delta=DELTA)
    result = model.run_simulation(local, sweeps=sweeps, burn_in=0, initial_path=start, centroid_source=clock,
                                  initial_path_provenance={"namespace": ns, "condition": condition, "role": "start",
                                                           "replica": i, "seed": seeds["start"]},
                                  uncertainty_batch_sizes=(512,), record_trace=True)
    trace = result.trace
    marks = np.asarray(clock.marks, dtype=np.int64)
    stream = bytes(local.emitted)
    onset_sweep = int(np.searchsorted(marks, onset_bytes, side="right")) if marks[-1] > onset_bytes else None
    record = {"namespace": ns, "condition": condition, "family": family, "severity": severity, "replica": i,
              "seeds": seeds, "sweeps": sweeps, "local_bytes": local.consumed_bytes,
              "local_sha256": hashlib.sha256(stream).hexdigest(),
              "prefix_sha256": hashlib.sha256(stream[:onset_bytes]).hexdigest(),
              "onset_sweep": onset_sweep, "acceptance_rate_4dp": result.acceptance_rate,
              "start_sha256": result.initialization_protocol["initial_path_sha256"]}
    if sweeps >= HARM_END:
        bm = summarize_series(trace["CV"][HARM_START:HARM_END], reference=exact_finite_bead_energy(M, BETA),
                              batch_sizes=(HARM_BATCH,))["batch_means"][str(HARM_BATCH)]
        record["harm"] = {"mean": bm["prefix_mean"], "mcse": bm["mcse"], "interval": bm["interval_95"],
                          "covered": bm["covered"], "status": bm["status"],
                          "wrong": not (bm["status"] == "ESTIMATED" and bm["covered"])}
    conv = []
    for end in range(CONV_WINDOW, sweeps + 1, STRIDE):
        lo, hi = window_bytes(marks, end - CONV_WINDOW, end)
        conv.append({"end": end, **conventional_window(stream[lo:hi])})
    d = trace["TD"] - trace["CV"]
    per_sweep_bytes = np.diff(np.concatenate([[0], marks]))
    uphill = (per_sweep_bytes - PROPOSAL_BYTES_PER_SWEEP) / (8 * M)
    app = [{"end": end, "mean_D": float(d[end - APP_WINDOW:end].mean()),
            "uphill_fraction": float(uphill[end - APP_WINDOW:end].mean()),
            "mean_rg2": float(trace["rg2"][end - APP_WINDOW:end].mean())}
           for end in range(APP_WINDOW, sweeps + 1, STRIDE)]
    tests = HealthTests(min_entropy=7, alpha=2**-20, alphabet_size=256)
    symbols = np.frombuffer(stream, dtype=np.uint8)
    rct, apt = health_test_events(symbols, tests)
    record["health_tests"] = {"rct_cutoff": tests.rct_cutoff, "apt_cutoff": tests.apt_cutoff, "window": tests.window,
                              "rct_pre": int(np.sum(rct <= onset_bytes)), "rct_post": int(np.sum(rct > onset_bytes)),
                              "apt_pre": int(np.sum(apt <= onset_bytes)), "apt_post": int(np.sum(apt > onset_bytes))}
    record["conventional"] = conv
    record["application"] = app
    record["per_sweep_bytes_ok"] = bool(np.all(per_sweep_bytes >= PROPOSAL_BYTES_PER_SWEEP)
                                        and np.all((per_sweep_bytes - PROPOSAL_BYTES_PER_SWEEP) % 8 == 0)
                                        and len(marks) == sweeps and np.all(np.diff(marks) > 0))
    arrays = {"CV": trace["CV"], "TD": trace["TD"], "rg2": trace["rg2"], "marks": marks}
    return record, arrays


# --------------------------------------------------------------------------
# Arm statistics (PROTOCOL.md sections 2.5-2.7)
# --------------------------------------------------------------------------
def conventional_statistic(checkpoint: dict) -> float:
    ps = [checkpoint[k] for k in CONV_P_KEYS if checkpoint.get(k) is not None]
    return -math.log10(max(min(ps), P_FLOOR))


def application_constants(calibration_records: list[dict]) -> dict[str, float]:
    d = np.array([c["mean_D"] for r in calibration_records for c in r["application"]])
    f = np.array([c["uphill_fraction"] for r in calibration_records for c in r["application"]])
    return {"sigma_D": float(d.std(ddof=1)), "mu_f": float(f.mean()), "sigma_f": float(f.std(ddof=1))}


def application_statistic(checkpoint: dict, constants: dict) -> float:
    return max(abs(checkpoint["mean_D"]) / constants["sigma_D"],
               abs(checkpoint["uphill_fraction"] - constants["mu_f"]) / constants["sigma_f"])


def arm_series(record: dict, arm: str, constants: dict | None = None) -> list[tuple[int, float]]:
    if arm == "C":
        return [(c["end"], conventional_statistic(c)) for c in record["conventional"]]
    return [(c["end"], application_statistic(c, constants)) for c in record["application"]]


def threshold(calibration_maxima, order_statistic: int) -> float:
    return float(np.sort(np.asarray(calibration_maxima))[order_statistic - 1])


def detection(series: list[tuple[int, float]], tau: float, onset_sweep: int | None,
              marks: np.ndarray | None = None, onset_bytes: int = ONSET_BYTES) -> dict:
    whole_run = any(s > tau for _, s in series)
    post = [(end, s) for end, s in series if onset_sweep is not None and end > onset_sweep and s > tau]
    out = {"whole_run_alarm": whole_run, "post_onset_alarm": bool(post), "alarm_end": post[0][0] if post else None}
    if post:
        out["latency_sweeps"] = post[0][0] - onset_sweep
        if marks is not None:
            out["latency_bytes"] = int(marks[post[0][0] - 1]) - onset_bytes
    return out
