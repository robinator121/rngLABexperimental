"""B2 interval laws: completeness, two independent routes, bounds, precision limits."""
import math
from decimal import Context, Decimal, localcontext

import numpy as np
import pytest

from rnglab.nuclear import Channel, DecayScheme, Mode, Nuclide, state_probabilities
from rnglab.nuclear.paths import (
    bounds_hold, interval_law, path_categories, tree_probabilities,
)


def _chain(half_lives):
    """A linear beta- chain ending in a stable nuclide."""
    nuclides = [Nuclide(f"C{i}", 150 + i, 400, h, (Channel(f"C{i + 1}", Mode.BETA_MINUS, 1.0),))
                for i, h in enumerate(half_lives)]
    nuclides.append(Nuclide(f"C{len(half_lives)}", 150 + len(half_lives), 400, math.inf))
    return DecayScheme("chain", nuclides)


def _sum(values):
    with localcontext(Context(prec=300)):
        return sum(values, Decimal(0))


def test_categories_are_all_path_prefixes(branching, competing):
    assert len(path_categories(branching, 0)) == 9
    assert len(path_categories(competing, 0)) == 1 + 2 * (1 + 2 + 1 + 1)
    for scheme in (branching, competing):
        for i, nuclide in enumerate(scheme.nuclides):
            cats = path_categories(scheme, i)
            assert cats[0].nuclides == (i,) and cats[0].channels == ()
            assert len(set((c.nuclides, c.channels) for c in cats)) == len(cats)
            if nuclide.stable:
                assert len(cats) == 1


@pytest.mark.parametrize("delta", [1e-9, 0.3, 2.0, 50.0, 1e4])
def test_two_routes_and_layer_a_agree(branching, competing, delta):
    for scheme in (branching, competing):
        for i in range(len(scheme)):
            law = interval_law(scheme, i, delta)
            assert abs(_sum(law.probabilities) - 1) < Decimal("1e-40")
            assert np.max(np.abs(law.as_float() - tree_probabilities(scheme, i, delta))) < 1e-13
            ends = np.zeros(len(scheme))
            for category, q in zip(law.categories, law.as_float()):
                ends[category.end] += q
            assert np.max(np.abs(ends - state_probabilities(scheme, scheme.nuclides[i].name, delta)[0])) < 1e-13
            assert bounds_hold(scheme, law)


def test_tiny_intervals_keep_relative_accuracy(branching):
    """For delta -> 0 a path of m decays has probability w * prod(lambda delta) / m! to leading order."""
    delta = 1e-12
    law = interval_law(branching, 0, delta)
    rates = [math.log(2) / n.half_life if not n.stable else 0.0 for n in branching.nuclides]
    for category, q in zip(law.categories, law.probabilities):
        m = len(category.steps)
        weight = math.prod(branching.nuclides[a].channel_probabilities()[k] for a, k in category.steps)
        leading = weight * math.prod(rates[a] * delta for a, _ in category.steps) / math.factorial(m)
        assert float(q) == pytest.approx(leading, rel=1e-9)
    assert float(law.probabilities[3]) < 1e-30


@pytest.mark.parametrize("spacing", [1e-3, 1e-6, 1e-9, 1e-12])
def test_close_rates_use_more_digits_and_stay_accurate(spacing):
    scheme = _chain([1.0 * (1 + i * spacing) for i in range(4)])
    law = interval_law(scheme, 0, 1.0)
    assert np.max(np.abs(law.as_float() - tree_probabilities(scheme, 0, 1.0))) < 1e-14
    assert bounds_hold(scheme, law)
    if spacing <= 1e-9:
        assert law.digits_used > 50


def test_equal_rates_are_refused():
    scheme = _chain([1.0, 1.0])
    with pytest.raises(ValueError, match="equal decay constants"):
        interval_law(scheme, 0, 1.0)


@pytest.mark.parametrize("half_lives", [[1e-9, 1e9], [1e9, 1e-9, 1.0], [1e-6, 1e3, 1e-3]])
def test_extreme_rate_ratios(half_lives):
    scheme = _chain(half_lives)
    for delta in (1e-10, 1.0, 1e10):
        law = interval_law(scheme, 0, delta)
        assert abs(_sum(law.probabilities) - 1) < Decimal("1e-30")
        assert np.max(np.abs(law.as_float() - tree_probabilities(scheme, 0, delta))) < 1e-12
        assert bounds_hold(scheme, law)


def test_long_intervals_put_all_mass_on_stable_ends(branching):
    law = interval_law(branching, 0, 1e6)
    stable = {i for i, n in enumerate(branching.nuclides) if n.stable}
    mass = _sum([q for c, q in zip(law.categories, law.probabilities) if c.end in stable])
    assert abs(mass - 1) < Decimal("1e-40")


def test_conditional_probabilities_and_complements(branching):
    law = interval_law(branching, 0, 2.0)
    p, r = law.conditional()
    assert p.shape == r.shape == (len(law.categories) - 1,)
    assert np.all((p >= 0) & (p <= 1) & (r >= 0) & (r <= 1))
    assert np.max(np.abs(p + r - 1)) < 1e-15
    q = law.as_float()
    tails = np.cumsum(q[::-1])[::-1]
    assert np.allclose(p, q[:-1] / tails[:-1], rtol=1e-12, atol=0)


def test_bad_interval_lengths(branching):
    for delta in (0.0, -1.0, float("inf"), float("nan")):
        with pytest.raises(ValueError):
            interval_law(branching, 0, delta)
