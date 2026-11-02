"""Frozen whole-trajectory source diagnostics and explicit entropy policies."""
from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
from scipy.special import erfc, gammaincc

from rnglab14.statistics.multitest import holm_bonferroni_correction
from rnglab15.assessment import Health
from rnglab17.endpoints import rank_threshold


def diagnostics(data: bytes, alpha=.01):
    """Three SP800-22-inspired diagnostics; not the full STS battery.

    Unrounded p-values; Holm correction WITHIN this stream. Runs uses its
    frequency prerequisite. No min-entropy lower bound is estimated here.
    """
    bits = np.unpackbits(np.frombuffer(data, dtype=np.uint8)).astype(float)
    n = len(bits)
    if n < 1024 or not 0 < alpha < 1:
        raise ValueError("at least 1024 bits and valid alpha required")
    p = float(bits.mean())
    monobit = float(erfc(abs(2*bits.sum()-n)/math.sqrt(2*n)))
    blocks = bits[:n//128*128].reshape(-1, 128).mean(1)
    block = float(gammaincc(len(blocks)/2, 4*128*np.sum((blocks-.5)**2)/2))
    runs = 1 + np.count_nonzero(bits[1:] != bits[:-1])
    run_p = 0.0 if abs(p-.5) >= 2/math.sqrt(n) else float(
        erfc(abs(runs-2*n*p*(1-p))/(2*math.sqrt(2*n)*p*(1-p))))
    p_values = [monobit, block, run_p]
    adjusted = holm_bonferroni_correction(p_values, alpha)
    return {"p_values": dict(zip(("monobit", "block_frequency", "runs"), p_values)),
            "holm_adjusted": adjusted.adjusted_p_values,
            "rejected": any(adjusted.reject_decisions), "alpha_per_stream": alpha,
            "bits": n, "security_entropy_credit_bits": 0}


def source_scores(blocks: list[bytes]):
    """Causal score per fixed-size source block; no future input in a score."""
    if not blocks or len({len(b) for b in blocks}) != 1 or len(blocks[0]) < 2:
        raise ValueError("nonempty equal-sized source blocks required")
    bits = np.unpackbits(np.frombuffer(b"".join(blocks), dtype=np.uint8)).reshape(len(blocks), -1)
    n = bits.shape[1]
    bias = 2*math.sqrt(n)*abs(bits.mean(1)-.5)
    equal = bits[:, 1:] == bits[:, :-1]
    serial = 2*math.sqrt(n-1)*abs(equal.mean(1)-.5)
    repeated = np.array([False]+[a == b for a, b in zip(blocks[:-1], blocks[1:])])
    return np.maximum.reduce((bias, serial, repeated*math.sqrt(n)))


def calibrate(paths: list[list[bytes]], alpha=.05):
    """Joint max over scores AND all times; strict rank cutoff on separate paths."""
    scores = [float(source_scores(p).max()) for p in paths]
    return {"threshold": rank_threshold(scores, alpha), "alpha_whole_trajectory": alpha,
            "replicas": len(paths), "maximum_scores": scores,
            "rule": "strict > ceil((n+1)*(1-alpha))-th calibration order statistic",
            "assumption": "exchangeable healthy synthetic trajectories; no hardware claim"}


def health_screen(data: bytes):
    h = Health(h=1., alpha=2**-20, window=512)
    for bit in np.unpackbits(np.frombuffer(data, dtype=np.uint8)):
        h.step(int(bit))
    return h.report()


@dataclass(frozen=True)
class SeedEvidence:
    conditional_entropy_bits: float | None
    basis: str

    def __post_init__(self):
        value = self.conditional_entropy_bits
        if value is not None and (type(value) not in (float, int) or not math.isfinite(value) or not 0 <= value <= 256):
            raise ValueError("invalid conditional seed uncertainty")
        if self.basis not in ("synthetic_model", "external_assessment", "unvalidated"):
            raise ValueError("unknown evidence basis")


class EntropyGate:
    """Fail-closed model policy. Caller evidence is an assumption, not validation.

    A health alarm latches quarantine. Recovery requires a fresh source check
    AND adequate independently justified seed uncertainty, followed by reseeding
    in the caller. Source restoration alone cannot recover historical secrets.
    """
    def __init__(self, minimum_bits=128):
        if type(minimum_bits) is not int or not 1 <= minimum_bits <= 256:
            raise ValueError("invalid startup budget")
        self.minimum_bits = minimum_bits
        self.quarantined = False

    def check(self, evidence: SeedEvidence, *, healthy=True, recovery=False):
        sufficient = (evidence.basis != "unvalidated" and evidence.conditional_entropy_bits is not None
                      and evidence.conditional_entropy_bits >= self.minimum_bits)
        if not healthy or not sufficient:
            self.quarantined = True
            raise RuntimeError("entropy policy blocked operation")
        if self.quarantined and not recovery:
            raise RuntimeError("quarantine requires explicit assessed recovery")
        self.quarantined = False
