"""Layer A engines: layout, replay, consumption, invariants and edge words.

Statistical assertions here use fixed seeds of the test namespace and are
regression checks, not evidence. Evidence about the physics comes from the
frozen study research/v020-nuclear-a-verify.
"""
import math

import numpy as np
import pytest

from rnglab.nuclear import (
    PARTICLE_CODES, NumpyPCG64Source, WordStream, labelled_source, simulate_decays,
    simulate_decays_reference, state_probabilities,
)
from rnglab.nuclear.decay import EVENT_COLUMNS, MODE_BY_CODE, MODE_SPECS, Mode
from rnglab14.sources.base import RandomSource, SourceManifest


def _run(scheme, ns, label, n, horizon, *, engine=simulate_decays, start="Xa-400", **kw):
    stream = WordStream(labelled_source(ns, "engine", label))
    return engine(scheme, {start: n}, horizon, stream, **kw), stream


class ConstantWords(RandomSource):
    """Every byte equals ``value``: all-zero words give u = 0, all-ones words give u = 1 - 2**-53."""

    def __init__(self, value: int):
        super().__init__(SourceManifest(source_id=f"const-{value}", family="test", is_deterministic=True,
                                        is_physical=False))
        self.value = value

    def _read_bytes(self, n):
        return bytes([self.value]) * n


@pytest.mark.parametrize("scheme_name, start", [("branching", "Xa-400"), ("competing", "Pa-300")])
def test_vectorised_engine_equals_reference_bit_for_bit(request, ns, scheme_name, start):
    scheme = request.getfixturevalue(scheme_name)
    fast, fast_stream = _run(scheme, ns, f"xref-{scheme_name}", 1500, 6.0, chunk=97, start=start)
    slow, slow_stream = _run(scheme, ns, f"xref-{scheme_name}", 1500, 6.0, engine=simulate_decays_reference,
                             start=start)
    for name in EVENT_COLUMNS:
        assert np.array_equal(fast.events[name], slow.events[name]), name
    assert np.array_equal(fast.final_state, slow.final_state)
    assert fast.digest() == slow.digest()
    assert fast_stream.source.stream_digest == slow_stream.source.stream_digest


def test_chunk_size_does_not_change_the_result(branching, ns):
    digests = {_run(branching, ns, "chunks", 2000, 9.0, chunk=c)[0].digest() for c in (1, 7, 64, 2000, 10**6)}
    assert len(digests) == 1


def test_replay_and_seed_sensitivity(branching, ns):
    a, _ = _run(branching, ns, "replay", 3000, 9.0)
    b, _ = _run(branching, ns, "replay", 3000, 9.0)
    c, _ = _run(branching, ns, "replay-other", 3000, 9.0)
    assert a.digest() == b.digest() != c.digest()


@pytest.mark.parametrize("horizon", [0.0, 0.3, 12.0, math.inf])
def test_consumption_is_fixed_by_the_layout(branching, ns, horizon):
    run, stream = _run(branching, ns, "consumption", 1234, horizon)
    assert run.words_per_nucleus == 12 and run.words_used == 1234 * 12
    assert stream.words_read == 1234 * 12 and stream.source.consumed_bytes == 8 * 1234 * 12


def test_horizon_zero_and_infinity(branching, ns):
    none, _ = _run(branching, ns, "h0", 500, 0.0)
    assert none.n_events == 0 and np.array_equal(none.final_state, none.initial_state)
    full, _ = _run(branching, ns, "hinf", 500, math.inf)
    stable = {i for i, n in enumerate(branching.nuclides) if n.stable}
    assert set(full.final_state.tolist()) <= stable
    steps = np.bincount(full.events["nucleus"], minlength=500)
    assert set(steps.tolist()) <= {3, 4}


def test_event_invariants(competing, ns):
    run, _ = _run(competing, ns, "invariants", 4000, 5.0, start="Pa-300")
    ev = run.events
    assert np.all(ev["time"] >= 0) and np.all(ev["time"] <= 5.0)
    assert np.array_equal(ev["event_id"], ev["nucleus"] * competing.depth + ev["step"])
    assert np.unique(ev["event_id"]).size == run.n_events
    later = ev["step"] > 0
    position = np.searchsorted(ev["event_id"], ev["previous_event"][later])
    assert np.array_equal(ev["event_id"][position], ev["previous_event"][later])
    assert np.array_equal(ev["daughter"][position], ev["parent"][later])
    assert np.all(ev["time"][position] <= ev["time"][later])
    first = ~later
    assert np.array_equal(ev["parent"][first], run.initial_state[ev["nucleus"][first]])
    for code in np.unique(ev["mode"]):
        spec = MODE_SPECS[MODE_BY_CODE[int(code)]]
        rows = ev["mode"] == code
        dz = [competing.nuclides[d].z - competing.nuclides[p].z for p, d in zip(ev["parent"][rows], ev["daughter"][rows])]
        da = [competing.nuclides[d].a - competing.nuclides[p].a for p, d in zip(ev["parent"][rows], ev["daughter"][rows])]
        assert set(dz) == {spec.delta_z} and set(da) == {spec.delta_a}
    assert np.array_equal(ev["neutron"], ev["mode"] == MODE_SPECS[Mode.BETA_DELAYED_NEUTRON].code)
    photon_yield = competing.photon_yield[ev["parent"], ev["channel"]]
    assert not np.any(ev["photon"][photon_yield == 0.0]) and np.all(ev["photon"][photon_yield == 1.0])
    last_rows = np.flatnonzero(np.r_[ev["nucleus"][1:] != ev["nucleus"][:-1], True])
    final = run.initial_state.copy()
    final[ev["nucleus"][last_rows]] = ev["daughter"][last_rows]
    assert np.array_equal(final, run.final_state)


def test_populations_match_a_per_nucleus_reconstruction(branching, ns):
    run, _ = _run(branching, ns, "populations", 800, 10.0)
    times = [0.0, 0.25, 1.0, 3.3, 10.0]
    expected = np.zeros((len(times), len(branching)), dtype=np.int64)
    for k, t in enumerate(times):
        for nucleus in range(run.n_nuclei):
            rows = (run.events["nucleus"] == nucleus) & (run.events["time"] <= t)
            state = run.events["daughter"][rows][-1] if rows.any() else run.initial_state[nucleus]
            expected[k, state] += 1
    assert np.array_equal(run.populations(times), expected)
    with pytest.raises(ValueError):
        run.populations([11.0])


def test_emissions_are_time_ordered_and_complete(competing, ns):
    run, _ = _run(competing, ns, "emissions", 3000, 8.0, start="Pa-300")
    em = run.emissions()
    assert np.all(np.diff(em["time"]) >= 0)
    mode = run.events["mode"]
    count = {name: int((em["particle"] == code).sum()) for name, code in PARTICLE_CODES.items()}
    code = {m: MODE_SPECS[m].code for m in Mode}
    assert count["alpha"] == int((mode == code[Mode.ALPHA]).sum())
    assert count["positron"] == int((mode == code[Mode.BETA_PLUS]).sum())
    assert count["neutron"] == int((mode == code[Mode.BETA_DELAYED_NEUTRON]).sum()) > 0
    assert count["electron"] == int(((mode == code[Mode.BETA_MINUS]) | (mode == code[Mode.BETA_DELAYED_NEUTRON])).sum())
    ec = mode == code[Mode.ELECTRON_CAPTURE]
    assert count["xray"] == int((run.events["photon"] & ec).sum())
    assert count["gamma"] == int((run.events["photon"] & ~ec).sum())


def test_only_beta_delayed_neutron_channels_emit_neutrons(branching, ns):
    run, _ = _run(branching, ns, "neutrons", 5000, 20.0)
    neutron_rows = run.events["neutron"]
    assert set(run.events["parent"][neutron_rows].tolist()) == {branching.index["Xe-396"]}
    assert int(neutron_rows.sum()) == int(run.channel_counts()[branching.index["Xe-396"], 1])


def test_all_zero_words_decay_instantly_along_first_channels(branching):
    run = simulate_decays(branching, {"Xa-400": 3}, 1.0, WordStream(ConstantWords(0)))
    assert np.all(run.events["time"] == 0.0) and np.all(run.events["channel"] == 0)
    assert [branching.nuclides[s].name for s in run.final_state] == ["Xf-392"] * 3
    photon_yield = branching.photon_yield[run.events["parent"], 0]
    assert np.array_equal(run.events["photon"], photon_yield > 0)


def test_all_one_words_give_the_longest_waits_and_last_channels(branching):
    run = simulate_decays(branching, {"Xa-400": 2}, math.inf, WordStream(ConstantWords(255)))
    longest = 53 * math.log(2) / (math.log(2) / 2.0)
    assert run.events["time"][0] == pytest.approx(longest, rel=1e-12)
    assert np.all(run.events["channel"] == np.array([1, 0, 0, 1] * 2))
    assert run.events["photon"][run.events["parent"] == branching.index["Xc-400m"]].all()
    assert not run.events["photon"][run.events["parent"] == branching.index["Xa-400"]].any()


def test_runs_share_a_stream_back_to_back_like_an_offset(branching, ns):
    stream = WordStream(labelled_source(ns, "engine", "offset"))
    first = simulate_decays(branching, {"Xa-400": 100}, 5.0, stream)
    second = simulate_decays(branching, {"Xa-400": 100}, 5.0, stream)
    assert second.first_word == first.words_used == 1200
    seed = labelled_source(ns, "engine", "offset").manifest.seed
    jumped = simulate_decays(branching, {"Xa-400": 100}, 5.0, WordStream(NumpyPCG64Source(seed, start_word=1200)))
    assert jumped.digest() == second.digest()


def test_large_sample_regression_against_exact_law(branching, ns):
    run, _ = _run(branching, ns, "regression", 50000, 4.0)
    p = state_probabilities(branching, "Xa-400", 4.0)[0]
    observed = run.populations([4.0])[0]
    sd = np.sqrt(50000 * p * (1 - p))
    z = (observed - 50000 * p)[p > 0] / sd[p > 0]
    assert np.all(np.abs(z) < 5)


def test_bad_arguments(branching, ns):
    stream = WordStream(labelled_source(ns, "engine", "bad"))
    for horizon in (-1.0, float("nan")):
        with pytest.raises(ValueError):
            simulate_decays(branching, {"Xa-400": 1}, horizon, stream)
    with pytest.raises(ValueError):
        simulate_decays(branching, {"Xa-400": 1}, 1.0, stream, chunk=0)
    with pytest.raises(TypeError):
        simulate_decays(branching, {"Xa-400": 1}, 1.0, labelled_source(ns, "engine", "bad"))
