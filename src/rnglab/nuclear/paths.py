"""RNGLab v0.20 nuclear module, milestone B2: decay-path probabilities over one interval.

A nucleus that starts an interval of length delta in nuclide s0 follows one
decay path within it: a (possibly empty) sequence of channels s0 -> s1 -> ... -> sm,
after which it does not decay again before the interval ends. In an acyclic
scheme the path prefixes from s0 form a finite set of mutually exclusive
categories whose probabilities sum to 1:

    P(path) = w * B(lambda_0, ..., lambda_m; delta),
    w = product of the channel probabilities along the path (``exact_channel_probabilities``),
    B = (prod_{q<m} lambda_q) * sum_q exp(-lambda_q delta) / prod_{r != q} (lambda_r - lambda_q),

lambda_m being the decay constant of the last nuclide (0 if stable). This is
the closed-form Bateman solution of Layer A, used per path.

Precision
---------
The alternating sum cancels when rates are close, and probabilities of long
paths over short intervals are tiny. Both are handled by evaluating B in
decimal arithmetic (Python's ``decimal``) with adaptive precision: the number
of digits lost to cancellation, log10(sum |terms| / |sum|), is measured and
the evaluation is repeated with enough digits to keep ``EXTRA_DIGITS``
significant digits. Results therefore carry relative accuracy, also for tiny
probabilities, as long as they are representable after conversion to float.
Rates are ln 2 / T_half computed in decimal from the exact binary value of T_half.

Limits: equal decay constants along a path are refused (the closed form
divides by their difference; the confluent case is not implemented), as is any
evaluation needing more than ``MAX_DIGITS`` digits.

Independent checks (no shared code with the decimal route):
- ``tree_probabilities``: the matrix exponential of the generator on the tree of
  path prefixes (float64, absolute accuracy);
- ``mean_value_bounds``: by the mean value theorem for divided differences,
  P(path) = w * prod_{q<m}(lambda_q delta) * exp(-xi delta) / m! for some xi between
  the smallest and largest rate on the path, which brackets every probability.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from decimal import Context, Decimal, DivisionByZero, InvalidOperation, localcontext
from typing import Sequence

import numpy as np
from scipy.linalg import expm

from .decay import DecayScheme

MIN_DIGITS = 50
EXTRA_DIGITS = 30
MAX_DIGITS = 1000


@dataclass(frozen=True)
class PathCategory:
    """One path prefix from ``nuclides[0]``: the nuclides visited and the channel index of each step."""

    nuclides: tuple[int, ...]
    channels: tuple[int, ...]

    @property
    def end(self) -> int:
        return self.nuclides[-1]

    @property
    def steps(self) -> tuple[tuple[int, int], ...]:
        """(parent nuclide, channel index) of every decay on the path."""
        return tuple(zip(self.nuclides[:-1], self.channels))


def path_categories(scheme: DecayScheme, start: int) -> list[PathCategory]:
    """All path prefixes from ``start``, depth first, the empty path first."""
    out: list[PathCategory] = []

    def walk(nuclides: tuple[int, ...], channels: tuple[int, ...]):
        out.append(PathCategory(nuclides, channels))
        for k, channel in enumerate(scheme.nuclides[nuclides[-1]].channels):
            walk(nuclides + (scheme.index[channel.daughter],), channels + (k,))

    walk((int(start),), ())
    return out


def _context(digits: int) -> Context:
    return Context(prec=digits, Emin=-10 ** 9, Emax=10 ** 9, traps=[DivisionByZero, InvalidOperation])


def exact_channel_probabilities(nuclide) -> list[Decimal]:
    """The sampler's branching law as exact decimals: differences of its float thresholds, the last edge 1.

    The thresholds are the running float sums of the declared ratios, as in
    DecayScheme and Nuclide.channel_probabilities; their exact differences sum
    to 1 exactly, whereas the rounded floats of channel_probabilities need not.
    """
    running, edges = 0.0, [Decimal(0)]
    for channel in nuclide.channels[:-1]:
        running += channel.branching
        edges.append(Decimal(running))
    edges.append(Decimal(1))
    with localcontext(_context(1100)):  # enough digits for the exact difference of any two doubles
        return [high - low for low, high in zip(edges, edges[1:])]


def _rates_from_half_lives(half_lives: Sequence[float], digits: int) -> list[Decimal]:
    """ln 2 / T_half in decimal at ``digits`` digits (0 for a stable nuclide, T_half = inf)."""
    with localcontext(_context(digits)):
        ln2 = Decimal(2).ln()
        return [Decimal(0) if math.isinf(h) else ln2 / Decimal(h) for h in half_lives]


def decimal_rates(scheme: DecayScheme, digits: int) -> list[Decimal]:
    return _rates_from_half_lives([n.half_life for n in scheme.nuclides], digits)


def _bateman_decimal(rates: Sequence[Decimal], delta: Decimal, digits: int) -> tuple[Decimal, Decimal]:
    """(B, prefactor * sum |terms|) at the given precision."""
    with localcontext(_context(digits)):
        prefactor = Decimal(1)
        for rate in rates[:-1]:
            prefactor *= rate
        total, magnitude = Decimal(0), Decimal(0)
        for q, rate in enumerate(rates):
            denominator = Decimal(1)
            for r, other in enumerate(rates):
                if r != q:
                    denominator *= other - rate
            term = (-rate * delta).exp() / denominator
            total += term
            magnitude += abs(term)
        return prefactor * total, prefactor * magnitude


def bateman_decimal(half_lives: Sequence[float], delta: Decimal, *, min_digits: int = MIN_DIGITS,
                    extra_digits: int = EXTRA_DIGITS) -> tuple[Decimal, int]:
    """Closed-form Bateman probability with adaptive precision: (value, digits used).

    The rates are recomputed from the exact half-lives at every precision, so
    their differences gain digits together with the sum.
    """
    digits = min_digits
    while True:
        rates = _rates_from_half_lives(half_lives, digits + 20)
        try:
            value, magnitude = _bateman_decimal(rates, delta, digits)
        except (DivisionByZero, InvalidOperation) as error:
            raise ValueError("equal decay constants along a path are not supported by the closed form") from error
        if magnitude == 0:
            return Decimal(0), digits
        if value > 0:
            with localcontext(_context(digits)):
                lost = (magnitude / value).log10()
            if lost + extra_digits <= digits:
                return value, digits
            needed = int(lost) + extra_digits + 10
        else:
            needed = 2 * digits
        if needed > MAX_DIGITS:
            if digits >= MAX_DIGITS:
                raise ValueError(f"closed form needs more than {MAX_DIGITS} digits here; rates too close")
            needed = MAX_DIGITS  # the growth estimate can overshoot: try the maximum before refusing
        digits = min(MAX_DIGITS, max(needed, digits + 10))


@dataclass(frozen=True)
class IntervalLaw:
    """Path categories and their probabilities for one start nuclide and one interval length."""

    start: int
    delta: float
    categories: tuple[PathCategory, ...]
    probabilities: tuple[Decimal, ...]
    digits_used: int

    def as_float(self) -> np.ndarray:
        return np.array([float(p) for p in self.probabilities])

    def conditional(self) -> tuple[np.ndarray, np.ndarray]:
        """For k < K - 1: p_k = q_k / T_k and its complement T_{k+1} / T_k, T_k = sum_{j >= k} q_j.

        Computed in decimal and rounded once, so both keep relative accuracy.
        A category whose tail mass is zero gets p = 0 and complement 1.
        """
        with localcontext(_context(max(self.digits_used, MIN_DIGITS))):
            tails = [Decimal(0)] * (len(self.probabilities) + 1)
            for k in range(len(self.probabilities) - 1, -1, -1):
                tails[k] = tails[k + 1] + self.probabilities[k]
            p, r = [], []
            for k in range(len(self.probabilities) - 1):
                if tails[k] == 0:
                    p.append(0.0)
                    r.append(1.0)
                else:
                    p.append(float(self.probabilities[k] / tails[k]))
                    r.append(float(tails[k + 1] / tails[k]))
        return np.array(p), np.array(r)


def interval_law(scheme: DecayScheme, start: int, delta: float, *, min_digits: int = MIN_DIGITS,
                 extra_digits: int = EXTRA_DIGITS) -> IntervalLaw:
    """The exact law of the path followed within ``delta`` by a nucleus starting in ``start``.

    ``min_digits`` and ``extra_digits`` only change the working precision (for
    checks of the precision itself); the engine uses the defaults.
    """
    if not (isinstance(delta, (int, float)) and math.isfinite(delta) and delta > 0):
        raise ValueError("delta must be a finite positive number")
    categories = path_categories(scheme, start)
    digits_used = min_digits
    probabilities = []
    work_digits = min_digits + extra_digits
    half_lives = [n.half_life for n in scheme.nuclides]
    d = Decimal(delta)
    for category in categories:
        nuclide_weights = [exact_channel_probabilities(scheme.nuclides[parent])[channel]
                           for parent, channel in category.steps]
        value, digits = bateman_decimal([half_lives[i] for i in category.nuclides], d, min_digits=min_digits,
                                        extra_digits=extra_digits)
        digits_used = max(digits_used, digits)
        with localcontext(_context(max(digits, work_digits))):
            weight = Decimal(1)
            for w in nuclide_weights:
                weight *= w
            probabilities.append(weight * value)
    return IntervalLaw(int(start), float(delta), tuple(categories), tuple(probabilities), digits_used)


def tree_probabilities(scheme: DecayScheme, start: int, delta: float) -> np.ndarray:
    """Independent float64 route: matrix exponential on the tree of path prefixes."""
    categories = path_categories(scheme, start)
    index = {(c.nuclides, c.channels): i for i, c in enumerate(categories)}
    generator = np.zeros((len(categories), len(categories)))
    for i, category in enumerate(categories):
        nuclide = scheme.nuclides[category.end]
        if nuclide.stable:
            continue
        rate = math.log(2.0) / nuclide.half_life
        generator[i, i] = -rate
        for k, (channel, probability) in enumerate(zip(nuclide.channels, nuclide.channel_probabilities())):
            child = index[(category.nuclides + (scheme.index[channel.daughter],), category.channels + (k,))]
            generator[child, i] += rate * probability
    return expm(generator * delta)[:, 0]


def bounds_hold(scheme: DecayScheme, law: IntervalLaw, rtol: str = "1e-25") -> bool:
    """Whether every probability of ``law`` lies in its mean-value bracket, compared at high precision.

    Decimal comparisons must not run in the default 28-digit context, where a
    tolerance such as 1 - 1e-30 rounds to 1.
    """
    bounds = mean_value_bounds(scheme, law)
    with localcontext(_context(MAX_DIGITS)):
        tol = Decimal(rtol)
        return all(low * (1 - tol) <= value <= high * (1 + tol)
                   for (low, high), value in zip(bounds, law.probabilities))


def mean_value_bounds(scheme: DecayScheme, law: IntervalLaw) -> list[tuple[Decimal, Decimal]]:
    """Rigorous bracket of every path probability (mean value theorem for divided differences)."""
    rates = decimal_rates(scheme, MIN_DIGITS)
    bounds = []
    with localcontext(_context(MIN_DIGITS)):
        d = Decimal(law.delta)
        for category in law.categories:
            path_rates = [rates[i] for i in category.nuclides]
            weight = Decimal(1)
            for parent, channel in category.steps:
                weight *= exact_channel_probabilities(scheme.nuclides[parent])[channel]
            scale = weight
            for rate in path_rates[:-1]:
                scale *= rate * d
            scale /= math.factorial(len(path_rates) - 1)
            bounds.append((scale * (-max(path_rates) * d).exp(), scale * (-min(path_rates) * d).exp()))
    return bounds
