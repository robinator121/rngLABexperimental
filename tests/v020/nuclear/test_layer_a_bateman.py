"""Layer A analytic references: two independent routes, limits and hand-derived cases."""
import math

import numpy as np
import pytest

from rnglab.nuclear import (
    bateman_chain, decay_probabilities, expected_populations, generator_matrix, path_state_probabilities,
    state_probabilities, survival_probability, transient_equilibrium_ratio, truncated_exponential_cdf,
)

TIMES = [0.0, 0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 40.0]


def test_generator_conserves_probability(branching, competing):
    for scheme in (branching, competing):
        matrix = generator_matrix(scheme)
        assert np.max(np.abs(matrix.sum(axis=0))) <= 1e-15
        assert np.all(matrix - np.diag(np.diag(matrix)) >= 0)


@pytest.mark.parametrize("scheme_name, start", [("branching", "Xa-400"), ("branching", "Xc-400m"),
                                                ("competing", "Pa-300"), ("competing", "Pb-300")])
def test_matrix_exponential_agrees_with_closed_form_bateman(request, scheme_name, start):
    scheme = request.getfixturevalue(scheme_name)
    by_expm = state_probabilities(scheme, start, TIMES)
    by_paths = path_state_probabilities(scheme, start, TIMES)
    assert np.max(np.abs(by_expm - by_paths)) < 1e-12
    assert np.allclose(by_expm.sum(axis=1), 1.0, atol=1e-13, rtol=0)
    assert np.all(by_expm > -1e-15)


def test_two_member_chain_against_the_textbook_formula():
    l1, l2 = math.log(2) / 3.0, math.log(2) / 0.7
    t = np.array([0.0, 0.4, 1.0, 5.0])
    expected = l1 / (l2 - l1) * (np.exp(-l1 * t) - np.exp(-l2 * t))
    assert np.allclose(bateman_chain([l1, l2], t), expected, rtol=1e-13, atol=0)
    assert np.allclose(bateman_chain([l1], t), np.exp(-l1 * t), rtol=1e-15, atol=0)
    assert bateman_chain([l1, l2, 0.0], [1e6])[0] == pytest.approx(1.0, abs=1e-12)


def test_degenerate_or_invalid_rates_are_refused():
    with pytest.raises(ValueError, match="distinct"):
        bateman_chain([0.5, 0.5 * (1 + 1e-9)], 1.0)
    with pytest.raises(ValueError):
        bateman_chain([0.0, 0.5], 1.0)
    with pytest.raises(ValueError):
        bateman_chain([], 1.0)


def test_survival_and_single_nuclide_decay_probability(make_single):
    scheme = make_single(half_life=1.5)
    assert survival_probability(1.5, [1.5, 3.0]).tolist() == pytest.approx([0.5, 0.25], abs=1e-15)
    q = decay_probabilities(scheme, "Sa-400", [0.0, 1.5, 3.0])[:, 0]
    assert q == pytest.approx([0.0, 0.5, 0.75], abs=1e-14)


def test_decay_probabilities_tend_to_path_visit_probabilities(branching):
    q = decay_probabilities(branching, "Xa-400", 400.0)[0]
    visits = {"Xa-400": 1.0, "Xb-396": 0.64, "Xd-396": 0.64, "Xc-400m": 0.36, "Xc-400": 0.36, "Xe-396": 0.36}
    for name, probability in visits.items():
        assert q[branching.index[name]] == pytest.approx(probability, abs=1e-12)
    assert np.all(q[[branching.index[n] for n in ("Xf-392", "Xg-396", "Xh-395")]] == 0.0)


def test_decay_probability_equals_entries_minus_occupancy(branching):
    """q_j(t) = P(entered j by t) - P(in j at t); entered = initial or decayed into j."""
    times = [0.7, 3.0]
    q = decay_probabilities(branching, "Xa-400", times)
    p = state_probabilities(branching, "Xa-400", times)
    xb = branching.index["Xb-396"]
    entered = 0.64 * q[:, branching.index["Xa-400"]]
    assert np.allclose(q[:, xb], entered - p[:, xb], atol=1e-14, rtol=0)


def test_expected_populations_mix_initial_states(branching):
    mixed = expected_populations(branching, {"Xa-400": 10, "Xe-396": 5}, [2.0])
    separate = 10 * state_probabilities(branching, "Xa-400", 2.0) + 5 * state_probabilities(branching, "Xe-396", 2.0)
    assert np.allclose(mixed, separate, rtol=0, atol=1e-12)


def test_transient_equilibrium_is_reached():
    from rnglab.nuclear import Channel, DecayScheme, Mode, Nuclide
    scheme = DecayScheme("eq", [Nuclide("Ya-410", 160, 410, 10.0, (Channel("Yb-410", Mode.BETA_MINUS, 1.0),)),
                                Nuclide("Yb-410", 161, 410, 1.0, (Channel("Yc-406", Mode.ALPHA, 1.0),)),
                                Nuclide("Yc-406", 159, 406, math.inf)])
    p = state_probabilities(scheme, "Ya-410", [30.0])[0]
    ratio = (math.log(2) / 1.0 * p[1]) / (math.log(2) / 10.0 * p[0])
    assert ratio == pytest.approx(transient_equilibrium_ratio(10.0, 1.0), rel=1e-6)
    assert transient_equilibrium_ratio(10.0, 1.0) == pytest.approx(10 / 9, rel=1e-15)
    with pytest.raises(ValueError):
        transient_equilibrium_ratio(1.0, 10.0)


def test_truncated_exponential_cdf():
    cdf = truncated_exponential_cdf(0.5, 4.0)
    assert cdf(0.0) == 0.0 and cdf(4.0) == pytest.approx(1.0, abs=1e-15) and cdf(10.0) == pytest.approx(1.0)
    assert cdf(1.0) == pytest.approx((1 - math.exp(-0.5)) / (1 - math.exp(-2.0)), rel=1e-14)
    assert truncated_exponential_cdf(0.5, math.inf)(2.0) == pytest.approx(1 - math.exp(-1.0), rel=1e-14)
    with pytest.raises(ValueError):
        truncated_exponential_cdf(0.0, 1.0)


def test_times_are_validated(branching):
    with pytest.raises(ValueError):
        state_probabilities(branching, "Xa-400", [-1.0])
    with pytest.raises(ValueError):
        state_probabilities(branching, "Xa-400", [math.inf])


def _chain_reference(rates, t):
    """Matrix exponential of the linear-chain generator, on rates scaled to order one."""
    from scipy.linalg import expm
    scale = max(rates)
    matrix = np.zeros((len(rates), len(rates)))
    for i, rate in enumerate(rates):
        matrix[i, i] = -rate / scale
        if i + 1 < len(rates):
            matrix[i + 1, i] = rate / scale
    return expm(matrix * t * scale)[-1, 0]


def test_closed_form_refuses_ill_conditioned_or_non_finite_evaluations():
    """Close rates pass the pairwise guard but cancel catastrophically (an independent review's example)."""
    with pytest.raises(ValueError, match="ill-conditioned"):
        bateman_chain([1.0, 1.00001, 1.00002, 1.00003, 1.00004], [1.0])
    with pytest.raises(ValueError):
        bateman_chain([1.0, 2.0, 3.0], [1e308])


@pytest.mark.parametrize("scale", [1e-200, 1e-75, 1e-3, 1.0, 1e3, 1e75, 1e200])
@pytest.mark.parametrize("spacing", [2e-6, 1e-4, 1e-2, 0.3])
def test_closed_form_is_refused_or_within_the_rounding_limit(scale, spacing):
    """The contract: either ValueError, or within ROUNDING_LIMIT (absolute) of the matrix exponential."""
    from rnglab.nuclear.bateman import ROUNDING_LIMIT
    rates = [(1 + i * spacing) * scale for i in range(5)]
    for t in (0.1 / scale, 1.0 / scale, 45.0 / scale):
        try:
            value = bateman_chain(rates, [t])[0]
        except ValueError:
            continue
        assert abs(value - _chain_reference(rates, t)) <= ROUNDING_LIMIT


def test_well_separated_paths_stay_within_the_rounding_limit(branching, competing):
    for scheme in (branching, competing):
        for nuclide in scheme.nuclides:
            path_state_probabilities(scheme, nuclide.name, [0.0, 0.5, 3.0, 40.0])


def test_analytic_references_use_the_sampled_branching_law():
    from rnglab.nuclear import Channel, DecayScheme, Mode, Nuclide
    third = 1.0 / 3.0
    scheme = DecayScheme("thirds", [
        Nuclide("A", 150, 400, 1.0, (Channel("B", Mode.ALPHA, third), Channel("C", Mode.BETA_MINUS, third),
                                     Channel("D", Mode.BETA_PLUS, third))),
        Nuclide("B", 148, 396, math.inf), Nuclide("C", 151, 400, math.inf), Nuclide("D", 149, 400, math.inf)])
    probabilities = scheme.nuclides[0].channel_probabilities()
    matrix = generator_matrix(scheme)
    assert matrix[1:, 0].tolist() == [p * math.log(2) for p in probabilities]
    assert abs(sum(probabilities) - 1.0) <= 1e-15 and max(abs(p - third) for p in probabilities) <= 1e-15


def test_tiny_close_rates_no_longer_hit_subnormal_products():
    """Rates near 1e-75 with relative spacing 2e-6 (a confirmation review's example) gave -15.6 before scaling."""
    from rnglab.nuclear.bateman import ROUNDING_LIMIT
    rates = [(1 + i * 2e-6) * 1e-75 for i in range(5)]
    value = bateman_chain(rates, [45 / 1e-75])[0]
    assert abs(value - _chain_reference(rates, 45 / 1e-75)) <= ROUNDING_LIMIT


def test_closed_form_is_scale_invariant():
    rates = [0.7, 0.2, 1.3, 0.0]
    t = np.array([0.5, 3.0, 9.0])
    for scale in (1e-200, 1e-60, 1.0, 1e60, 1e200):
        assert np.allclose(bateman_chain([r * scale for r in rates], t / scale), bateman_chain(rates, t),
                           rtol=1e-12, atol=1e-15)
