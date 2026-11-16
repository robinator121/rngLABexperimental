"""RNGLab v0.20 nuclear module, milestone B2: exact binomial sampling by inversion.

One uniform per draw. A word gives u = k * 2**-53 in [0, 1), which stands for
the cell [u, u + 2**-53) of a continuous uniform; the draw inverts at the cell's
midpoint v = u + 2**-54 (RNGLab's convention since Stage C, Phi^-1(u + 2**-54)):

    X = min{x in 0..n : F(x) > v},   F the Binomial(n, p) distribution function,

so X has the Binomial(n, p) law up to the accuracy with which F is evaluated
and the 2**-53 resolution of u. The midpoint matters at the ends: inverting at
u = 0 itself would return x = 0 with probability 2**-53 even when that outcome
is astronomically unlikely (at n = 2**53, about 7e7 standard deviations below
the mean, enough to double the variance); the midpoint keeps the extreme cells
at the 2**-54 and 1 - 2**-54 quantiles. v itself is not representable above
1/2, so the comparisons use exact forms: F(x) - u > 2**-54 (the subtraction is
exact for F within a factor 2 of u) and S(x) < (1 - u) - 2**-54 (exact).

Discretisation bound: with exact F, the sampled distribution function equals
(number of midpoints below F(x)) * 2**-53, within 2**-54 of F(x) at every x,
and the support is [Q(2**-54), Q(1 - 2**-54)], about +-8.3 standard deviations
for large n. Summation by parts then bounds the error of the variance by
2**-54 * sum_x |2 (x - mean) + 1| over that range, a relative error of order
1e-14 (research/v020-nuclear-b2-numerics, midpoint_convention).

Evaluation of F, in two regimes (research/v020-nuclear-b2-numerics measures both)
-------------------------------------------------------------------------------
The predicate F(x) > v is evaluated either on the lower tail F or as
S(x) < 1 - v on the upper tail S = 1 - F; ``uses_upper_tail`` says which.
Extreme u needs the upper tail: placing u = 1 - 2**-53 correctly needs S near
1e-16 to relative, not absolute, accuracy.
- Small mean, n p <= TABLE_MEAN_MAX: F decides u < 1/2 and S decides u >= 1/2.
  Both tails are built from the exact
  recurrence w_0 = 1, w_{x+1} = w_x (n - x) / (x + 1) * p / (1 - p), normalised by
  its own sum (no exp or log of n p is needed), up to n p + 40 sd + 60 (the mass
  beyond is below 1e-300): F by summing from the left, S from the right.
- Otherwise: F(k) = betaincc(k + 1, n - k, p) and S(k) = betainc(k + 1, n - k, p)
  (scipy.special, Boost), which take p itself, never 1 - p.
Why two regimes: Boost's incomplete beta loses accuracy for small means at
1e6 <= n <= 1e9 (absolute error up to 7.4e-12 for F and 2.3e-8 for S at
n = 1e9, mean 10-30), while the recurrence table is accurate to about 1e-15
there. Above mean 256, Boost's F is accurate to about 1e-16 absolute up to
n = 1e9 and about 1e-13 at n = 1e15; its S is weaker near the centre (up to
4.8e-13 at n = 1e9, 7.8e-12 at n = 1e12, 2.8e-10 at n = 2**53) but accurate to
about 4e-13 relative in the far upper tail up to n = 1e9. So in the Boost
regime F decides every u except the far upper tail, 1 - u < UPPER_TAIL_SWITCH
= 2**-10, where S is consulted up to n = UPPER_TAIL_MAX_N = 1e9. Above 1e9, F
decides every u, and uniforms within about 1e-13 of 1 are resolved only to F's
absolute accuracy (a total-variation error below about 1e-13).
scipy.special.bdtr is not used at all: it is wrong by 0.34 at n = 1e9,
p = 1/2 and returns NaN from n = 1e12.

Search (Boost regime)
---------------------
A quantile guess (normal approximation with skewness and continuity
corrections) gives a bracket that is verified by evaluating the predicate at
both ends; a side that fails falls back to the global bound (-1 or n). Integer
bisection then finds X. The guess changes the amount of work, never the answer:
any valid bracket leads to the same X, because the predicate is monotone in k.

Domain
------
0 <= n <= PART_MAX = 2**53 per draw (every integer up to it is exact in double
precision) and 0 <= p <= 1/2. Callers sample min(p, 1 - p) and flip
(x -> n - x) when p > 1/2, with the complement computed accurately by them.
Larger n is handled by ``split_counts``: Binomial(n, p) is the sum of
independent binomials over parts of at most PART_MAX trials, one uniform each.
"""
from __future__ import annotations

import math

import numpy as np
from scipy import special

PART_MAX = 2 ** 53
TABLE_MEAN_MAX = 256.0
UPPER_TAIL_MAX_N = 10 ** 9  # above it Boost's upper tail S loses accuracy; F is used for every u
UPPER_TAIL_SWITCH = 2.0 ** -10  # Boost regime: S is consulted only when 1 - u < this (far upper tail)
_INITIAL_HALF_WIDTH = 4


def table_tails(n: int, p: float) -> tuple[np.ndarray, np.ndarray]:
    """(F(0..L), S(0..L)) of Binomial(n, p) from the normalised exact recurrence (small means only).

    F is summed from the left and S = P(X > k) from the right, so each tail
    keeps relative accuracy; beyond L both are 1 and 0 to within 1e-300.
    """
    mean = n * p
    if not mean <= TABLE_MEAN_MAX:
        raise ValueError(f"the recurrence table serves means up to {TABLE_MEAN_MAX} only (weights would overflow)")
    length = int(min(n, math.ceil(mean + 40.0 * math.sqrt(mean) + 60.0)))
    x = np.arange(length, dtype=np.float64)
    weights = np.empty(length + 1)
    weights[0] = 1.0
    weights[1:] = np.cumprod((n - x) / (x + 1.0) * (p / (1.0 - p)))
    total = weights.sum()
    lower = np.cumsum(weights) / total
    upper = np.zeros(length + 1)
    upper[:-1] = np.cumsum(weights[::-1])[::-1][1:] / total
    return lower, upper


def table_cdf(n: int, p: float) -> np.ndarray:
    """F(0..L) of the table regime (see ``table_tails``)."""
    return table_tails(n, p)[0]


def sampler_tails(k: int, n: int, p: float) -> tuple[float, float]:
    """(F(k), S(k)) exactly as the sampler evaluates them (for verification)."""
    if k < 0:
        return 0.0, 1.0
    if k >= n:
        return 1.0, 0.0
    if n * p <= TABLE_MEAN_MAX:
        lower, upper = table_tails(n, p)
        return (float(lower[k]), float(upper[k])) if k < lower.size else (1.0, 0.0)
    return float(special.betaincc(k + 1.0, n - k, p)), float(special.betainc(k + 1.0, n - k, p))


def uses_upper_tail(n: int, p: float, u: float) -> bool:
    """Whether the sampler decides u through the upper tail S (else through F)."""
    if n * p <= TABLE_MEAN_MAX:
        return u >= 0.5
    return n <= UPPER_TAIL_MAX_N and 1.0 - u < UPPER_TAIL_SWITCH


HALF_CELL = 2.0 ** -54


def _table_inverse(n: int, p: float, u: float) -> int:
    lower, upper = table_tails(n, p)
    if u < 0.5:
        return int(np.searchsorted(lower, u + HALF_CELL, side="right"))  # u + 2**-54 < 1/2 is exact
    below = upper < (1.0 - u) - HALF_CELL  # exact: 1 - u <= 1/2
    return int(np.argmax(below)) if below.any() else n


def _predicate(k: np.ndarray, n: np.ndarray, p: np.ndarray, u: np.ndarray) -> np.ndarray:
    """F(k) > u + 2**-54 for integer-valued k in [0, n - 1] (float arrays), Boost regime, tail chosen by u."""
    out = np.empty(k.shape, dtype=bool)
    lower = (1.0 - u >= UPPER_TAIL_SWITCH) | (n > UPPER_TAIL_MAX_N)
    if lower.any():
        kl = k[lower]
        f = special.betaincc(kl + 1.0, n[lower] - kl, p[lower])
        out[lower] = f - u[lower] > HALF_CELL  # exact when f is within a factor 2 of u; decisive otherwise
    upper = ~lower
    if upper.any():
        ku = k[upper]
        out[upper] = special.betainc(ku + 1.0, n[upper] - ku, p[upper]) < (1.0 - u[upper]) - HALF_CELL
    return out


def _predicate_full(k: np.ndarray, n: np.ndarray, p: np.ndarray, u: np.ndarray) -> np.ndarray:
    """The predicate on int64 k and n, with F(k) = 0 below 0 and F(k) = 1 from n on.

    k and n are converted to float only inside [0, n - 1], where they are exact (n <= 2**53).
    """
    out = k >= n
    inside = (k >= 0) & (k < n)
    if inside.any():
        out[inside] = _predicate(k[inside].astype(np.float64), n[inside].astype(np.float64), p[inside], u[inside])
    return out


def binomial_inverse(n, p, u, *, half_width: int = _INITIAL_HALF_WIDTH) -> np.ndarray:
    """Exact inversion of Binomial(n, p) at uniforms u (vectorised; int64 result).

    ``half_width`` only sets the initial bracket around the quantile guess; it
    never changes the result (tests rely on that).
    """
    n = np.asarray(n, dtype=np.int64)
    p = np.asarray(p, dtype=np.float64)
    u = np.asarray(u, dtype=np.float64)
    n, p, u = np.broadcast_arrays(n, p, u)
    shape = n.shape
    n, p, u = n.ravel(), p.ravel(), u.ravel()
    if np.any(n < 0) or np.any(n > PART_MAX):
        raise ValueError(f"n must lie in [0, {PART_MAX}] per draw; split larger counts")
    if not np.all(np.isfinite(p)) or np.any(p < 0.0) or np.any(p > 0.5):
        raise ValueError("p must lie in [0, 1/2]; sample the complement and flip")
    if not np.all(np.isfinite(u)) or np.any(u < 0.0) or np.any(u >= 1.0):
        raise ValueError("u must lie in [0, 1)")
    out = np.zeros(n.shape, dtype=np.int64)
    active = (n > 0) & (p > 0.0)
    small = active & (n.astype(np.float64) * p <= TABLE_MEAN_MAX)
    for index in np.flatnonzero(small):
        out[index] = _table_inverse(int(n[index]), float(p[index]), float(u[index]))
    active &= ~small
    if not active.any():
        return out.reshape(shape)
    na, pa, ua = n[active], p[active], u[active]
    nf = na.astype(np.float64)
    mean = nf * pa
    sd = np.sqrt(nf * pa * (1.0 - pa))
    with np.errstate(divide="ignore", invalid="ignore"):
        z = special.ndtri(ua + HALF_CELL)
    z = np.nan_to_num(z, nan=0.0, posinf=40.0, neginf=-40.0)
    guess_f = np.floor(mean + sd * z + (1.0 - 2.0 * pa) * (z * z - 1.0) / 6.0 - 0.5)
    guess = np.clip(np.nan_to_num(guess_f, nan=0.0), 0.0, nf).astype(np.int64)
    width = max(int(half_width), 1)
    low_try = np.maximum(guess - width, -1)
    high_try = np.minimum(guess + width, na)
    low_ok = ~_predicate_full(low_try, na, pa, ua)      # F(low_try) <= u: a valid lower end
    high_ok = _predicate_full(high_try, na, pa, ua)     # F(high_try) > u: a valid upper end
    lo = np.where(low_ok, low_try, -1).astype(np.int64)
    hi = np.where(high_ok, high_try, na).astype(np.int64)
    # Invariant: predicate(lo) is False and predicate(hi) is True (integers, int64 arithmetic).
    while True:
        open_ = hi - lo > 1
        if not open_.any():
            break
        lo_open, hi_open = lo[open_], hi[open_]
        mid = lo_open + (hi_open - lo_open) // 2
        result = _predicate_full(mid, na[open_], pa[open_], ua[open_])
        hi[open_] = np.where(result, mid, hi_open)
        lo[open_] = np.where(result, lo_open, mid)
    out[active] = hi
    return out.reshape(shape)


def split_counts(n: int, parts: int, part_max: int = PART_MAX) -> list[int]:
    """Deterministic split of n trials into ``parts`` parts of at most part_max each."""
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or n < 0:
        raise ValueError("n must be a nonnegative integer")
    if n > parts * part_max:
        raise ValueError("n exceeds parts * part_max")
    n = int(n)
    return [min(part_max, max(0, n - s * part_max)) for s in range(parts)]


def binomial_draw(n: int, p_small: float, flip: bool, uniforms: np.ndarray, *, part_max: int = PART_MAX) -> int:
    """Binomial(n, p) from len(uniforms) parts; p = p_small, or 1 - p_small when flip (scalar reference)."""
    sizes = split_counts(n, len(uniforms), part_max)
    draws = binomial_inverse(np.array(sizes, dtype=np.int64), np.full(len(sizes), p_small), np.asarray(uniforms))
    total = sum(int(d) for d in draws)
    return n - total if flip else total
