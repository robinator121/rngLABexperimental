"""A documented SP 800-22 subset, streaming health checks, and diagnostics.

No result is an entropy assessment, certification, or cryptographic guarantee.
Health-test entropy bounds must be supplied by the caller; they are not inferred
from the same output being checked.
"""
from __future__ import annotations

import math
from typing import Any, Iterable

import numpy as np
from scipy.special import gammaincc
from scipy.stats import binom


MAX_OFFLINE_SAMPLES = 1_000_000


def _alpha(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float, np.number)):
        raise ValueError("alpha must be a finite number in (0,1)")
    value = float(value)
    if not math.isfinite(value) or not 0 < value < 1:
        raise ValueError("alpha must be a finite number in (0,1)")
    return value


def _positive_int(value: int, name: str, minimum: int = 1) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def bit_array(bits: Any, maximum: int = MAX_OFFLINE_SAMPLES) -> np.ndarray:
    """Accept a one-dimensional 0/1 sequence, never implicit packed-byte decoding."""
    if isinstance(bits, str):
        if set(bits) - {"0", "1"}:
            raise ValueError("bit string contains characters other than 0 and 1")
        data = np.frombuffer(bits.encode("ascii"), dtype=np.uint8) - 48
    else:
        data = np.asarray(bits)
    if data.ndim != 1 or data.size > maximum:
        raise ValueError(f"one-dimensional bits required, at most {maximum} per call")
    if data.dtype.kind not in "buif" or not np.isfinite(data).all() or np.any((data != 0) & (data != 1)):
        raise ValueError("bits must contain only finite numeric 0/1 values")
    return data.astype(np.uint8, copy=False)


def _result(pvalue: float | None, alpha: float, reason: str | None = None, **data: Any) -> dict[str, Any]:
    if pvalue is None:
        return {"status": "NOT_APPLICABLE", "p_value": None, "alpha": alpha, "reason": reason, **data}
    p = float(pvalue)
    return {"status": "FAIL" if p < alpha else "PASS", "p_value": p, "alpha": alpha, **data}


def evaluate_bits(bits: Any, alpha: float = .01, block_size: int = 128) -> dict[str, Any]:
    """Frequency, block frequency, and runs using SP 800-22 Rev 1a formulas.

    Input recommendations are enforced, including runs' prerequisite and block
    frequency M>=20, M>0.01*n, N<100. Alpha is pointwise per test/window.
    ``PASS`` means failure to reject the test null, not proven randomness.
    """
    a = _alpha(alpha)
    m = _positive_int(block_size, "block_size", 20)
    b = bit_array(bits)
    n = len(b)
    tests: dict[str, dict[str, Any]] = {}
    if n < 100:
        tests = {name: _result(None, a, "INSUFFICIENT_DATA: n >= 100 recommended", n=n)
                 for name in ("frequency", "block_frequency", "runs")}
    else:
        ones = int(b.sum(dtype=np.int64))
        pi = ones / n
        statistic = abs(2 * ones - n) / math.sqrt(n)
        tests["frequency"] = _result(math.erfc(statistic / math.sqrt(2)), a,
                                     n=n, ones=ones, statistic=statistic)
        count = n // m
        if m > n or m <= .01 * n or count >= 100:
            tests["block_frequency"] = _result(None, a,
                "INVALID_PARAMETERS: require M<=n, M>0.01*n and floor(n/M)<100",
                n=n, block_size=m, blocks=count)
        else:
            proportions = b[:count*m].reshape(count, m).mean(axis=1)
            chi2 = float(4*m*np.sum((proportions-.5)**2))
            tests["block_frequency"] = _result(float(gammaincc(count/2, chi2/2)), a,
                n=n, block_size=m, blocks=count, unused_bits=n-count*m, statistic=chi2)
        runs = 1 + int(np.count_nonzero(b[1:] != b[:-1]))
        tau = 2 / math.sqrt(n)
        if abs(pi - .5) >= tau:
            tests["runs"] = _result(None, a,
                "PREREQUISITE_FAILED: frequency deviation >= 2/sqrt(n); inspect frequency result",
                n=n, runs=runs, proportion_ones=pi)
        else:
            numerator = abs(runs - 2*n*pi*(1-pi))
            denominator = 2*math.sqrt(2*n)*pi*(1-pi)
            tests["runs"] = _result(math.erfc(numerator/denominator), a,
                                     n=n, runs=runs, proportion_ones=pi)
    return {"standard": "NIST SP 800-22 Rev 1a: three-test subset",
            "full_suite": False, "certified": False, "n_bits": n, "tests": tests,
            "multiplicity": "Pointwise alpha only; repeated windows/tests need prospective joint calibration",
            "interpretation": "A pass does not establish physical entropy or cryptographic security"}


class HealthTests:
    """SP 800-90B-informed RCT/APT on integer source symbols.

    Uses the NIST cutoff expression (section 4.4, APT footnote 10) with
    Binomial(W, 2**-H). Chunk boundaries have no statistical meaning. APT has
    nonoverlapping W-symbol windows anchored by the first symbol. Errors latch
    within a run/window to keep reporting bounded; checks continue afterward.
    """

    def __init__(self, min_entropy: float, alpha: float = 2**-20,
                 window: int | None = None, alphabet_size: int = 2):
        self.alphabet_size = _positive_int(alphabet_size, "alphabet_size", 2)
        self.alpha = _alpha(alpha)
        if (isinstance(min_entropy, bool) or not isinstance(min_entropy, (int, float, np.number))
                or not math.isfinite(float(min_entropy)) or not 0 < min_entropy <= math.log2(self.alphabet_size)):
            raise ValueError("an asserted positive min-entropy bound <= log2(alphabet_size) is required")
        self.min_entropy = float(min_entropy)
        self.window = _positive_int(window if window is not None else (1024 if alphabet_size == 2 else 512), "window", 2)
        self.rct_cutoff = 1 + math.ceil(-math.log2(self.alpha)/self.min_entropy)
        p = 2**(-self.min_entropy)
        # Search survival probabilities to avoid catastrophic cancellation in 1-alpha.
        lo, hi = 1, self.window + 1
        while lo < hi:
            mid = (lo+hi)//2
            if float(binom.sf(mid-1, self.window, p)) <= self.alpha:
                hi = mid
            else:
                lo = mid+1
        self.apt_cutoff = lo
        self.samples = 0
        self._last: int | None = None
        self._run = 0
        self._run_alerted = False
        self._anchor: int | None = None
        self._position = self._matches = 0
        self._apt_alerted = False

    def update(self, symbols: Iterable[int]) -> list[dict[str, Any]]:
        # A bounded chunk avoids an accidental unbounded generator/materialization.
        x = np.asarray(symbols)
        if x.ndim != 1 or len(x) > MAX_OFFLINE_SAMPLES or x.dtype.kind not in "bui":
            raise ValueError("bounded one-dimensional integer source symbols required")
        if np.any(x < 0) or np.any(x >= self.alphabet_size):
            raise ValueError("source symbol outside declared alphabet")
        events = []
        for raw in x:
            value = int(raw)
            self.samples += 1
            if value == self._last:
                self._run += 1
            else:
                self._last, self._run, self._run_alerted = value, 1, False
            if self._run >= self.rct_cutoff and not self._run_alerted:
                events.append({"test": "repetition_count", "status": "FAIL", "sample": self.samples,
                               "count": self._run, "cutoff": self.rct_cutoff})
                self._run_alerted = True
            if self._position == 0:
                self._anchor, self._matches, self._apt_alerted = value, 0, False
            self._position += 1
            self._matches += int(value == self._anchor)
            if self._matches >= self.apt_cutoff and not self._apt_alerted:
                events.append({"test": "adaptive_proportion", "status": "FAIL", "sample": self.samples,
                               "window_start": self.samples-self._position+1,
                               "count": self._matches, "cutoff": self.apt_cutoff})
                self._apt_alerted = True
            if self._position == self.window:
                self._position = 0
        return events

    def parameters(self) -> dict[str, Any]:
        return {"standard": "NIST SP 800-90B-informed health tests", "certified": False,
                "asserted_min_entropy_per_symbol": self.min_entropy, "alpha_per_opportunity": self.alpha,
                "alphabet_size": self.alphabet_size, "window": self.window,
                "standard_window": self.window == (1024 if self.alphabet_size == 2 else 512),
                "repetition_count_cutoff": self.rct_cutoff, "adaptive_proportion_cutoff": self.apt_cutoff,
                "entropy_assessment": "NOT_IMPLEMENTED: externally justified bound required",
                "assumptions": "Threshold probabilities require the source model behind the asserted bound; marginal entropy alone does not bound arbitrary dependent processes",
                "false_alarm_scope": "Per opportunity/window; neither lifetime nor combined-test alpha"}

    def state_dict(self) -> dict[str, Any]:
        return {"version": 1, "parameters": {"min_entropy": self.min_entropy, "alpha": self.alpha,
                "window": self.window, "alphabet_size": self.alphabet_size},
                "samples": self.samples, "last": self._last, "run": self._run,
                "run_alerted": self._run_alerted, "anchor": self._anchor,
                "position": self._position, "matches": self._matches, "apt_alerted": self._apt_alerted}

    @classmethod
    def from_state_dict(cls, state: dict[str, Any]) -> "HealthTests":
        if state.get("version") != 1:
            raise ValueError("unknown health-test state version")
        result = cls(**state["parameters"])
        for field in ("samples", "run", "position", "matches"):
            value = state[field]
            if isinstance(value, bool) or not isinstance(value, int) or value < 0:
                raise ValueError("invalid health-test state counters")
        if state["position"] >= result.window or state["matches"] > result.window:
            raise ValueError("invalid APT checkpoint")
        for name in ("last", "anchor"):
            value = state[name]
            if value is not None and (not isinstance(value, int) or not 0 <= value < result.alphabet_size):
                raise ValueError("invalid checkpoint symbol")
        result.samples, result._run = state["samples"], state["run"]
        result._last, result._run_alerted = state["last"], bool(state["run_alerted"])
        result._anchor, result._position = state["anchor"], state["position"]
        result._matches, result._apt_alerted = state["matches"], bool(state["apt_alerted"])
        return result


def uniform_array(samples: Any, maximum: int = MAX_OFFLINE_SAMPLES) -> np.ndarray:
    u = np.asarray(samples)
    if u.ndim == 1:
        u = u[:, None]
    if u.ndim != 2 or not u.shape[0] or not u.shape[1] or u.size > maximum:
        raise ValueError("a nonempty bounded samples-by-channels matrix is required")
    if u.dtype.kind not in "uif" or not np.isfinite(u).all() or np.any(u < 0) or np.any(u >= 1):
        raise ValueError("uniform values must be finite numeric values in [0,1)")
    return np.asarray(u, dtype=float)


def _pearson(a: np.ndarray, b: np.ndarray) -> float | None:
    a, b = a-a.mean(), b-b.mean()
    scale = math.sqrt(float(np.dot(a, a)*np.dot(b, b)))
    return None if not scale else float(np.clip(np.dot(a, b)/scale, -1, 1))


def diagnostic_window(uniforms: Any, lags: tuple[int, ...] = (1, 2, 8, 32, 128, 256),
                      word_bits: int = 8) -> dict[str, Any]:
    """Descriptive offline/sliding-window analyses, without invented p-values.

    All sequence diagnostics use the declared numerical uniform representation.
    Repeated words use nonoverlapping thresholded (u>=.5) bit words, not source
    bytes. Constant-column correlations are unavailable, never evidence of zero
    dependence. Sliding callers must independently calibrate decisions.
    """
    u = uniform_array(uniforms)
    word_bits = _positive_int(word_bits, "word_bits", 2)
    if word_bits > 20 or any(isinstance(lag, bool) or not isinstance(lag, int) or lag < 1 for lag in lags):
        raise ValueError("word_bits<=20 and positive integer lags required")
    n, d = u.shape
    side = u >= .5
    result: dict[str, Any] = {"status": "DESCRIPTIVE_ONLY", "rows": n, "channels": d,
         "units": "numerical_uniforms", "p_values": None, "channels_detail": [], "cross_channel": []}
    for j in range(d):
        x, b = u[:, j], side[:, j]
        runs = 1 + int(np.count_nonzero(b[1:] != b[:-1]))
        lengths = np.diff(np.r_[0, np.flatnonzero(b[1:] != b[:-1])+1, n])
        split = n//2
        chunks = [chunk for chunk in np.array_split(x, min(8, n)) if len(chunk)]
        correlations = {str(lag): _pearson(x[:-lag], x[lag:]) if n-lag >= 3 else None for lag in lags}
        centered = x-x.mean()
        power = np.abs(np.fft.rfft(centered))**2
        if len(power) > 1 and power[1:].sum() > 0:
            peak = 1 + int(np.argmax(power[1:]))
            periodicity = {"peak_cycles_per_sample": peak/n,
                           "peak_power_fraction": float(power[peak]/power[1:].sum())}
        else:
            periodicity = {"peak_cycles_per_sample": None, "peak_power_fraction": None}
        count = n//word_bits
        if count:
            words = b[:count*word_bits].reshape(count, word_bits).astype(np.uint64) @ (2**np.arange(word_bits, dtype=np.uint64))
            _, multiplicities = np.unique(words, return_counts=True)
            repeated = {"word_bits": word_bits, "words": count,
                        "duplicate_word_occurrences": count-len(multiplicities),
                        "maximum_multiplicity": int(multiplicities.max()),
                        "unused_threshold_bits": n-count*word_bits}
        else:
            repeated = {"status": "NOT_APPLICABLE", "reason": "insufficient bits for a word"}
        result["channels_detail"].append({"channel": j, "mean": float(x.mean()),
            "bias_from_half": float(x.mean()-.5), "threshold_bit_fraction": float(b.mean()),
            "runs": runs, "longest_threshold_run": int(lengths.max()), "autocorrelation": correlations,
            "serial_same_side_fraction": float(np.mean(b[1:] == b[:-1])) if n>1 else None,
            "half_mean_shift": float(x[split:].mean()-x[:split].mean()) if split else None,
            "block_mean_range": float(np.ptp([v.mean() for v in chunks])),
            "extreme_value_fraction": float(np.mean((x<.01)|(x>=.99))),
            "constant": bool(np.ptp(x)==0), "periodicity": periodicity, "repeated_words": repeated})
    for j in range(d):
        for k in range(j+1, d):
            result["cross_channel"].append({"channels": [j, k],
                "pearson": _pearson(u[:, j], u[:, k]) if n>=3 else None,
                "threshold_agreement": float(np.mean(side[:, j] == side[:, k]))})
    return result
