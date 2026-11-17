"""B2 binomial sampler: exact inversion against rational arithmetic, domain, splitting."""
import math
from fractions import Fraction

import numpy as np
import pytest

from rnglab.nuclear.binomial import PART_MAX, binomial_draw, binomial_inverse, split_counts


HALF_CELL = Fraction(1, 2 ** 54)


def _exact_quantiles(n, p, us):
    """min{x : F(x) > u + 2**-54} (the midpoint of u's cell) with F computed exactly in rational arithmetic."""
    P = Fraction(p)
    cdf, total = [], Fraction(0)
    for x in range(n + 1):
        total += math.comb(n, x) * P ** x * (1 - P) ** (n - x)
        cdf.append(total)
    return np.array([next(x for x in range(n + 1) if cdf[x] > Fraction(u) + HALF_CELL) for u in us])


@pytest.mark.parametrize("n", [1, 2, 5, 17, 64])
@pytest.mark.parametrize("p", [1e-6, 0.01, 0.2, 0.5])
def test_inversion_matches_exact_rational_inversion(n, p):
    grid = np.linspace(0.0, 1.0 - 2.0 ** -53, 301)
    edges = np.array([0.0, 2.0 ** -53, 0.5, 0.5 - 2.0 ** -53, 1.0 - 2.0 ** -53])
    us = np.concatenate([grid, edges])
    got = binomial_inverse(np.full(us.size, n), np.full(us.size, p), us)
    assert np.array_equal(got, _exact_quantiles(n, p, us))


def test_inversion_just_either_side_of_every_cdf_step():
    """u within 1e-12 of each exact CDF step lands on the correct side.

    Exactly at a step the answer depends on the last bits of F, which no
    double-precision evaluation resolves; that limit is documented, not tested.
    """
    n, p = 12, 0.25
    P = Fraction(p)
    cdf, total = [], Fraction(0)
    for x in range(n):
        total += math.comb(n, x) * P ** x * (1 - P) ** (n - x)
        cdf.append(float(total))
    us = np.array([c + d for c in cdf for d in (-1e-12, 1e-12) if 0.0 <= c + d < 1.0])
    assert np.array_equal(binomial_inverse(np.full(us.size, n), np.full(us.size, p), us),
                          _exact_quantiles(n, p, us))


def test_result_does_not_depend_on_the_bracket():
    rng = np.random.default_rng(7)
    n = rng.integers(0, 10 ** 9, 200)
    p = rng.uniform(0.0, 0.5, 200)
    u = rng.random(200)
    reference = binomial_inverse(n, p, u)
    for width in (1, 3, 50, 10 ** 6):
        assert np.array_equal(binomial_inverse(n, p, u, half_width=width), reference)


def test_degenerate_cases_and_domain():
    assert binomial_inverse([0, 10, 10], [0.3, 0.0, 0.5], [0.9, 0.9, 0.0]).tolist() == [0, 0, 0]
    assert binomial_inverse([10], [0.5], [1.0 - 2.0 ** -53]).tolist() == [10]
    for bad in (dict(n=[-1], p=[0.1], u=[0.1]), dict(n=[PART_MAX + 1], p=[0.1], u=[0.1]),
                dict(n=[5], p=[0.6], u=[0.1]), dict(n=[5], p=[float("nan")], u=[0.1]),
                dict(n=[5], p=[0.1], u=[1.0]), dict(n=[5], p=[0.1], u=[-0.1])):
        with pytest.raises(ValueError):
            binomial_inverse(**bad)


def test_largest_supported_n_runs_and_stays_in_range():
    u = np.array([0.0, 1e-12, 0.3, 0.5, 0.9, 1.0 - 2.0 ** -53])
    x = binomial_inverse(np.full(u.size, PART_MAX), np.full(u.size, 0.5), u)
    assert np.all((x >= 0) & (x <= PART_MAX)) and np.all(np.diff(x) >= 0)
    z = (x - PART_MAX / 2) / math.sqrt(PART_MAX / 4)
    assert np.all(np.abs(z[2:5]) < 3)
    # The extreme cells sit at the 2**-54 and 1 - 2**-54 quantiles, about 8.3 sd out, not at 0 or n.
    assert -9.0 < z[0] < -7.5 and 7.0 < z[-1] < 9.5


def test_split_counts():
    assert split_counts(10, 3, 4) == [4, 4, 2]
    assert split_counts(0, 2, 4) == [0, 0]
    with pytest.raises(ValueError):
        split_counts(13, 3, 4)


def test_split_and_flipped_draws_follow_the_exact_law():
    """Many parts and flipped probabilities: the pooled draws match Binomial(n, p) (fixed-seed regression)."""
    from scipy import stats
    rng = np.random.default_rng(11)
    n, p, parts, part_max = 23, 0.7, 6, 4
    draws = np.array([binomial_draw(n, 1.0 - p, True, rng.random(parts), part_max=part_max) for _ in range(4000)])
    observed = np.bincount(draws, minlength=n + 1)
    expected = stats.binom.pmf(np.arange(n + 1), n, p) * draws.size
    keep = expected >= 5
    chi2 = (((observed[keep] - expected[keep]) ** 2) / expected[keep]).sum() + \
        (observed[~keep].sum() - expected[~keep].sum()) ** 2 / expected[~keep].sum()
    assert stats.chi2.sf(chi2, keep.sum()) > 1e-4
    assert draws.min() >= 0 and draws.max() <= n


def test_table_regime_cdf_is_exact_for_small_n():
    from rnglab.nuclear.binomial import table_cdf
    for n, p in ((1, 0.3), (9, 0.5), (40, 1e-3), (200, 0.2)):
        P = Fraction(p)
        total, exact = Fraction(0), []
        for x in range(n + 1):
            total += math.comb(n, x) * P ** x * (1 - P) ** (n - x)
            exact.append(float(total))
        table = table_cdf(n, p)
        assert np.max(np.abs(table - np.array(exact[:table.size]))) < 1e-14


@pytest.mark.parametrize("p", [0.42, 0.426, 0.43, 0.44])
def test_both_sides_of_the_regime_boundary_match_exact_inversion(p):
    """n = 600: means 252 and 255.6 use the table, 258 and 264 use the incomplete beta."""
    from rnglab.nuclear.binomial import TABLE_MEAN_MAX
    n = 600
    us = np.linspace(0.0, 1.0 - 2.0 ** -53, 97)
    got = binomial_inverse(np.full(us.size, n), np.full(us.size, p), us)
    assert np.array_equal(got, _exact_quantiles(n, p, us))
    assert (n * p <= TABLE_MEAN_MAX) == (p <= 0.426)


def test_scalar_and_shaped_inputs_keep_their_shape():
    assert binomial_inverse(10, 0.3, 0.5).shape == ()
    assert int(binomial_inverse(10, 0.3, 0.5)) == int(binomial_inverse([10], [0.3], [0.5])[0])
    assert int(binomial_inverse(10 ** 12, 0.3, 0.5)) == int(binomial_inverse([10 ** 12], [0.3], [0.5])[0])
    grid = binomial_inverse(np.full((2, 3), 50), 0.2, np.array([[0.1, 0.5, 0.9], [0.0, 0.3, 0.99]]))
    assert grid.shape == (2, 3)
    assert int(binomial_inverse(0, 0.2, 0.7)) == 0 and int(binomial_inverse(7, 0.2, 0.0)) == 0


def test_a_zero_uniform_is_not_an_outlier():
    """u = 0 has probability 2**-53. Inverting at u itself would return x = 0, which at n = 2**53, p = 1/2
    lies 9.5e7 sd below the mean and would by itself double the variance (an independent review's finding;
    at the present part limit 2**50 it would still add 12.5%). The midpoint convention returns the 2**-54
    quantile instead."""
    for n, p in ((PART_MAX, 0.5), (10 ** 12, 0.25), (10 ** 6, 0.5), (5000, 0.4)):
        sd = math.sqrt(n * p * (1 - p))
        low = (int(binomial_inverse(n, p, 0.0)) - n * p) / sd
        high = (int(binomial_inverse(n, p, 1.0 - 2.0 ** -53)) - n * p) / sd
        assert -9.5 < low < -7.0 and 7.0 < high < 9.5, (n, p, low, high)


def test_midpoint_convention_on_a_tiny_case():
    """n = 1: X = 1 iff u + 2**-54 >= 1 - p, checked on both sides of the step for several p."""
    for p in (0.5, 0.25, 2.0 ** -53, 3 * 2.0 ** -53):
        step = 1.0 - p                      # exact for these p
        below, at = step - 2.0 ** -53, step
        assert int(binomial_inverse(1, p, below)) == (1 if Fraction(below) + HALF_CELL >= Fraction(step) else 0)
        assert int(binomial_inverse(1, p, at)) == 1
        assert int(binomial_inverse(1, p, 0.0)) == 0
