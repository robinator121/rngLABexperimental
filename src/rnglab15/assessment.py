"""Streaming descriptive screens; these are NOT SP 800-90B estimators."""
from __future__ import annotations
from collections import Counter
import math
from scipy.stats import binom


class Moments:
    def __init__(self):
        self.n = 0
        self.mean = self.m2 = 0.0
        self.minimum = self.maximum = None

    def add(self, x):
        self.n += 1
        delta = x - self.mean
        self.mean += delta / self.n
        self.m2 += delta * (x - self.mean)
        self.minimum = x if self.minimum is None else min(self.minimum, x)
        self.maximum = x if self.maximum is None else max(self.maximum, x)

    def report(self):
        return {"n": self.n, "mean": self.mean if self.n else None,
                "sd": math.sqrt(max(0, self.m2) / (self.n - 1)) if self.n > 1 else None,
                "min": self.minimum, "max": self.maximum}


class Timing:
    def __init__(self):
        self.previous = None
        self.count = self.zero = self.backward = self.lattice = 0
        self.gaps = Moments()
        self.previous_gap = None
        self.xy = self.xx = self.yy = self.sx = self.sy = 0.0
        self.pairs = 0

    def add(self, tick):
        self.count += 1
        gap = None if self.previous is None else tick - self.previous
        self.previous = tick
        if gap is None:
            return None
        self.zero += gap == 0
        self.backward += gap < 0
        if gap <= 0:
            self.previous_gap = None
            return gap
        self.lattice = math.gcd(self.lattice, gap)
        self.gaps.add(gap)
        if self.previous_gap is not None:
            x, y = self.previous_gap, gap
            self.pairs += 1
            self.sx += x; self.sy += y
            self.xx += x*x; self.yy += y*y; self.xy += x*y
        self.previous_gap = gap
        return gap

    def report(self, scale=None):
        n = self.pairs
        den = (n*self.xx-self.sx**2)*(n*self.yy-self.sy**2)
        corr = (n*self.xy-self.sx*self.sy)/math.sqrt(den) if den > 0 else None
        g = self.gaps.report()
        return {"records": self.count, "zero_gaps": self.zero, "backward_gaps": self.backward,
                "positive_gap_gcd_raw_units": self.lattice, "positive_gaps": g,
                "gap_cv": g["sd"]/g["mean"] if g["n"] > 1 and g["mean"] else None,
                "lag1_positive_gap_correlation": corr,
                "positive_gap_rate_hz": 1/(float(scale)*g["mean"]) if scale and g["mean"] else None,
                "dead_time": "unknown; minimum recorded gap is not a dead-time estimate",
                "clock_scale": str(scale) if scale else None}


class Health:
    """State survives chunks. Threshold assumptions are explicit, alarms latch.

    APT counts the first symbol and Binomial(W-1,p) matches. The cutoff is the
    smallest total count C with Pr[1+Binomial(W-1,p)>=C] <= alpha.
    Binary input, fixed non-overlapping windows, no production validation claim.
    """
    def __init__(self, h=0.5, alpha=2**-20, window=512):
        if not (math.isfinite(h) and 0 < h <= 1 and 0 < alpha < 1 and window >= 2):
            raise ValueError("invalid binary health-test assumptions")
        self.h, self.alpha, self.window = h, alpha, window
        self.rct_cutoff = 1 + math.ceil(-math.log2(alpha)/h)
        p = 2**-h
        self.apt_cutoff = next((c for c in range(2, window+2)
                                if binom.sf(c-2, window-1, p) <= alpha), window+1)
        self.n = self.run = self.pos = self.matches = self.failures = 0
        self.last = self.target = None
        self.first_alarm = None

    def step(self, bit):
        if bit not in (0, 1):
            raise ValueError("binary health test requires 0/1 symbols")
        self.n += 1
        self.run = self.run + 1 if bit == self.last else 1
        self.last = bit
        if self.pos == 0:
            self.target, self.matches = bit, 0
        self.matches += bit == self.target
        self.pos = (self.pos + 1) % self.window
        fail = self.run >= self.rct_cutoff or self.matches >= self.apt_cutoff
        self.failures += fail
        if fail and self.first_alarm is None:
            self.first_alarm = self.n - 1
        return self.first_alarm is None

    def report(self):
        return {"assumed_h_bits_per_binary_symbol": self.h, "alpha_per_test_opportunity": self.alpha,
                "window": self.window, "rct_cutoff": self.rct_cutoff, "apt_cutoff": self.apt_cutoff,
                "symbols": self.n, "alarm_positions": self.failures, "first_alarm_zero_based": self.first_alarm,
                "incomplete_apt_window_symbols": self.pos, "latched_failure": self.first_alarm is not None,
                "qualification": "exploratory thresholds; source H unvalidated; no whole-run alpha guarantee"}


def screen(bits):
    """IID-only simultaneous Hoeffding bound, plus frozen held-out predictors.

    Predictor scores cannot supply a conditional min-entropy lower bound for
    arbitrary dependent data. Security entropy credit remains ZERO.
    """
    n = len(bits)
    if n == 0 or any(x not in (0, 1) for x in bits):
        raise ValueError("nonempty binary symbols required")
    count = Counter(bits)
    delta = 0.01
    p_upper = min(1.0, max(count.values())/n + math.sqrt(math.log(2/delta)/(2*n)))
    cut = n//2
    train, test = bits[:cut], bits[cut:]
    predictors = {}
    majority = int(sum(train) > len(train)/2)
    for order in (0, 1, 2, 8):
        table = {}
        for i in range(order, len(train)):
            context = tuple(train[i-order:i])
            if context not in table:
                table[context] = [0, 0]
            table[context][train[i]] += 1
        correct = 0
        history = list(train[-order:]) if order else []
        for b in test:
            counts = table.get(tuple(history[-order:]) if order else (), [1-majority, majority])
            prediction = int(counts[1] > counts[0])
            correct += prediction == b
            if order:
                history.append(b)
                del history[:-order]
        predictors[str(order)] = {"correct": correct, "tested": len(test),
                                  "accuracy": correct/len(test) if test else None}
    return {"symbols": n, "ones_fraction": count[1]/n,
            "iid_only_99pct_marginal_min_entropy_lower_bits": -math.log2(p_upper),
            "iid_assumption": "not established; bound does not apply to dependent samples",
            "heldout_predictors_context_lengths": predictors,
            "non_iid_min_entropy_lower_bits": 0.0, "validated_entropy_credit_bits": 0,
            "nist_90b_status": "not run; use official IID/non-IID and restart tools on a defined source"}
