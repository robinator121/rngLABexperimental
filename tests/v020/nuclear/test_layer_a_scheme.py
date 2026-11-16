"""Layer A decay-scheme data model: bookkeeping, normalisation and acyclicity checks."""
import math

import numpy as np
import pytest

from rnglab.nuclear import Channel, DecayScheme, Mode, Nuclide


def _scheme(*nuclides):
    return DecayScheme("t", list(nuclides))


def test_branching_scheme_structure(branching):
    assert branching.depth == 4 and branching.words_per_nucleus == 12 and branching.max_channels == 2
    assert branching.precursors() == pytest.approx({"Xe-396": 0.1}, abs=1e-15)
    rates = branching.decay_constants
    assert rates[branching.index["Xa-400"]] == math.log(2) / 2.0 and rates[branching.index["Xf-392"]] == 0.0
    assert branching.cum_thresholds[branching.index["Xa-400"], 0] == 0.64
    assert np.isinf(branching.cum_thresholds[branching.index["Xb-396"], 0])


def test_round_trip_and_digest_are_stable(branching, competing):
    for scheme in (branching, competing):
        again = DecayScheme.from_dict(scheme.to_dict())
        assert again.to_dict() == scheme.to_dict() and again.sha256() == scheme.sha256()
    assert branching.sha256() != competing.sha256()


def test_paths_include_parallel_channels_to_one_daughter(competing):
    prefixes = list(competing.paths("Pa-300"))
    to_b = [w for p, w in prefixes if p == (0, 1)]
    assert sorted(to_b) == [0.3, 0.7]
    total_at_depth_one = math.fsum(w for p, w in prefixes if len(p) == 2)
    assert total_at_depth_one == pytest.approx(1.0, abs=1e-15)


def test_arrays_are_read_only(branching):
    with pytest.raises(ValueError):
        branching.decay_constants[0] = 1.0


@pytest.mark.parametrize("make, message", [
    (lambda: Nuclide("A", 150, 400, 1.0, (Channel("B", Mode.ALPHA, 0.5), Channel("C", Mode.BETA_MINUS, 0.4))),
     "sum to"),
    (lambda: Nuclide("A", 150, 400, 1.0), "stable if and only if"),
    (lambda: Nuclide("A", 150, 400, math.inf, (Channel("B", Mode.ALPHA, 1.0),)), "stable if and only if"),
    (lambda: Nuclide("A", 150, 400, float("nan")), "NaN"),
    (lambda: Nuclide("A", 150, 400, -1.0), "positive"),
    (lambda: Nuclide("A", 401, 400, math.inf), "Z <= A"),
    (lambda: Nuclide("A", 150, 400, 1.0, (Channel("B", Mode.ALPHA, 1.0), Channel("B", Mode.ALPHA, 1e-13))),
     "duplicate"),
    (lambda: Channel("B", Mode.ALPHA, 0.0), "branching"),
    (lambda: Channel("B", Mode.ALPHA, 1.0, 1.5), "photon yield"),
])
def test_invalid_nuclides_and_channels_are_rejected(make, message):
    with pytest.raises((ValueError, TypeError), match=message):
        make()


def test_bool_and_wrong_types_are_rejected():
    with pytest.raises(TypeError):
        Channel("B", Mode.ALPHA, True)
    with pytest.raises(TypeError):
        Channel("B", "alpha", 1.0)
    with pytest.raises(TypeError):
        Nuclide("A", 150.0, 400, math.inf)


def test_normalisation_tolerance_accepts_rounding_but_does_not_renormalise():
    third = 1.0 / 3.0
    nuclide = Nuclide("A", 150, 400, 1.0, (Channel("B", Mode.ALPHA, third), Channel("C", Mode.BETA_MINUS, third),
                                           Channel("D", Mode.BETA_PLUS, third)))
    assert [c.branching for c in nuclide.channels] == [third] * 3


def test_mode_bookkeeping_is_enforced():
    with pytest.raises(ValueError, match="contradicts"):
        _scheme(Nuclide("A", 150, 400, 1.0, (Channel("B", Mode.ALPHA, 1.0),)), Nuclide("B", 149, 396, math.inf))
    with pytest.raises(ValueError, match="lower the level"):
        _scheme(Nuclide("A", 150, 400, 1.0, (Channel("B", Mode.ISOMERIC_TRANSITION, 1.0),)),
                Nuclide("B", 150, 400, math.inf, level=1))
    with pytest.raises(ValueError, match="unknown daughter"):
        _scheme(Nuclide("A", 150, 400, 1.0, (Channel("Z", Mode.ALPHA, 1.0),)))
    with pytest.raises(ValueError, match="share"):
        _scheme(Nuclide("A", 150, 400, math.inf), Nuclide("B", 150, 400, math.inf))


def test_cycles_are_rejected():
    with pytest.raises(ValueError, match="cycle"):
        _scheme(Nuclide("A", 150, 400, 1.0, (Channel("B", Mode.BETA_MINUS, 1.0),)),
                Nuclide("B", 151, 400, 1.0, (Channel("A", Mode.BETA_PLUS, 1.0),)))


def test_initial_states_follow_scheme_order(branching):
    states = branching.initial_states({"Xe-396": 2, "Xa-400": 3})
    assert states.tolist() == [0, 0, 0, 5, 5]
    with pytest.raises(ValueError):
        branching.initial_states({"nope": 1})
    with pytest.raises(ValueError):
        branching.initial_states({"Xa-400": -1})
    with pytest.raises(ValueError):
        branching.initial_states({"Xa-400": 1.0})


def test_all_stable_scheme_has_depth_zero():
    scheme = _scheme(Nuclide("A", 150, 400, math.inf))
    assert scheme.depth == 0 and scheme.words_per_nucleus == 0 and scheme.precursors() == {}


def test_half_lives_outside_the_supported_range_are_refused():
    for half_life in (1e-310, 1e-301, 1e301):
        with pytest.raises(ValueError, match="must lie in"):
            Nuclide("A", 150, 400, half_life, (Channel("B", Mode.ALPHA, 1.0),))
    Nuclide("A", 150, 400, 1e-300, (Channel("B", Mode.ALPHA, 1.0),))


def test_channels_rounded_away_by_the_thresholds_are_refused():
    """Declared [1.0, 5e-13] passes the sum tolerance, but the sampler could never choose the second channel."""
    with pytest.raises(ValueError, match="no probability left"):
        Nuclide("A", 150, 400, 1.0, (Channel("B", Mode.ALPHA, 1.0), Channel("C", Mode.BETA_MINUS, 5e-13)))


def test_channel_probabilities_follow_the_thresholds(branching):
    xa = branching.nuclides[branching.index["Xa-400"]]
    assert xa.channel_probabilities() == (0.64, 1.0 - 0.64)
    assert branching.cum_thresholds[branching.index["Xa-400"], 0] == xa.channel_probabilities()[0]


def test_stable_nuclides_have_no_channel_probabilities():
    assert Nuclide("S", 150, 400, math.inf).channel_probabilities() == ()
