"""B2 engine: bookkeeping, layout, replay, checkpoints, bounded memory, agreement with exact laws.

Statistical assertions use fixed seeds of the B2 test namespace; they are
regression checks, not evidence. Evidence comes from research/v020-nuclear-b2-verify.
"""
import json
import math

import numpy as np
import pytest

from rnglab.nuclear import NumpyPCG64Source, WordStream, labelled_source, simulate_decays, state_probabilities
from rnglab.nuclear.aggregate import (
    INT64_MAX, AggregatePlan, simulate_aggregate, simulate_aggregate_reference,
)
from rnglab.nuclear.checks import multinomial_check
from rnglab.nuclear.decay import Mode

TIMES = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]


def _run(scheme, ns, label, initial, times=TIMES, engine=simulate_aggregate, **kw):
    stream = WordStream(labelled_source(ns, "b2-engine", label))
    return engine(scheme, initial, times, stream, **kw), stream


@pytest.mark.parametrize("name, initial, part_max", [
    ("branching", {"Xa-400": 50000}, 2 ** 53), ("branching", {"Xa-400": 3000, "Xe-396": 70}, 997),
    ("competing", {"Pa-300": 20000}, 2 ** 53), ("precursors", {"G1": 4000, "G2": 4000, "G3": 50}, 1500)])
def test_vectorised_equals_reference(request, ns, name, initial, part_max):
    scheme = request.getfixturevalue(name)
    fast, fs = _run(scheme, ns, f"xref-{name}-{part_max}", initial, part_max=part_max)
    slow, ss = _run(scheme, ns, f"xref-{name}-{part_max}", initial, engine=simulate_aggregate_reference,
                    part_max=part_max)
    assert fast.digest == slow.digest
    for key in ("populations", "path_counts", "decays", "photons", "neutrons"):
        assert np.array_equal(getattr(fast, key), getattr(slow, key)), key
    assert fs.source.stream_digest == ss.source.stream_digest and fs.words_read == ss.words_read


def test_bookkeeping_invariants(branching, competing, precursors, ns):
    for scheme, initial in ((branching, {"Xa-400": 30000, "Xc-400": 500}), (competing, {"Pa-300": 30000}),
                            (precursors, {"G1": 9000, "G3": 900})):
        run, stream = _run(scheme, ns, f"inv-{scheme.name}", initial)
        total = sum(initial.values())
        assert np.all(run.populations >= 0) and set(run.populations.sum(axis=1).tolist()) == {total}
        inflow = np.zeros_like(run.populations[1:])
        for j, nuclide in enumerate(scheme.nuclides):
            for k, channel in enumerate(nuclide.channels):
                inflow[:, scheme.index[channel.daughter]] += run.decays[:, j, k]
        outflow = run.decays.sum(axis=2)
        assert np.array_equal(np.diff(run.populations, axis=0), inflow - outflow)
        assert np.all(outflow <= run.populations[:-1] + inflow)  # nuclei can enter and leave within one interval
        assert np.all((run.photons >= 0) & (run.photons <= run.decays))
        yields = np.zeros(run.decays.shape[1:])
        neutron = np.zeros(run.decays.shape[1:], dtype=bool)
        for j, nuclide in enumerate(scheme.nuclides):
            for k, channel in enumerate(nuclide.channels):
                yields[j, k] = channel.photon_yield
                neutron[j, k] = channel.mode is Mode.BETA_DELAYED_NEUTRON
        assert np.all(run.photons[:, yields == 0.0] == 0)
        assert np.array_equal(run.photons[:, yields == 1.0], run.decays[:, yields == 1.0])
        assert np.array_equal(run.neutrons, run.decays[:, neutron].sum(axis=1))
        assert stream.words_read == run.words_read == (len(TIMES) - 1) * run.plan.words_per_interval
        assert stream.source.consumed_bytes == 8 * stream.words_read


def test_consumption_does_not_depend_on_the_population(branching, ns):
    plans = [AggregatePlan(branching, TIMES, n) for n in (0, 1, 10 ** 6, 2 ** 53)]
    assert len({p.words_per_interval for p in plans}) == 1 and plans[0].parts == 1
    run, stream = _run(branching, ns, "empty", {"Xa-400": 0})
    assert stream.words_read == 6 * plans[0].words_per_interval and np.all(run.populations == 0)
    assert AggregatePlan(branching, TIMES, 2 ** 53 + 1).parts == 2


def test_replay_and_seed_sensitivity(branching, ns):
    a, _ = _run(branching, ns, "replay", {"Xa-400": 10 ** 6})
    b, _ = _run(branching, ns, "replay", {"Xa-400": 10 ** 6})
    c, _ = _run(branching, ns, "replay-other", {"Xa-400": 10 ** 6})
    assert a.digest == b.digest != c.digest


def test_checkpoint_resume_is_bit_identical(branching, ns):
    checkpoints = []
    full, _ = _run(branching, ns, "checkpoint", {"Xa-400": 10 ** 7}, checkpoint_every=1,
                   on_checkpoint=checkpoints.append)
    seed = labelled_source(ns, "b2-engine", "checkpoint").manifest.seed
    for cut in (1, 3, 5):
        checkpoint = json.loads(json.dumps(checkpoints[cut - 1]))
        stream = WordStream(NumpyPCG64Source(seed, start_word=checkpoint["words_consumed"]))
        rest = simulate_aggregate(branching, {"Xa-400": 10 ** 7}, TIMES, stream, resume_from=checkpoint)
        assert rest.digest == full.digest and rest.first_interval == cut
        assert np.array_equal(rest.populations, full.populations[cut:])
        assert np.array_equal(rest.decays, full.decays[cut:])


def test_malformed_checkpoints_are_refused_not_repaired(branching, ns):
    """Fractional or wrongly typed fields and a bad digest must be refused (an independent review's finding)."""
    checkpoints = []
    _run(branching, ns, "checkpoint-malformed", {"Xa-400": 1000}, times=[0.0, 1.0, 2.0], checkpoint_every=1,
         on_checkpoint=checkpoints.append)
    good = checkpoints[1]
    stream = WordStream(labelled_source(ns, "b2-engine", "unused"))
    pops = good["populations"]
    cases = [("next_interval", 2.9), ("next_interval", True), ("next_interval", "2"), ("words_consumed", 50.0),
             ("populations", [float(x) for x in pops]), ("populations", pops[:-1]),
             ("populations", [pops[0] + 0.5, pops[1] - 0.5] + pops[2:]), ("running_sha256", ""),
             ("running_sha256", "zz" * 32), ("running_sha256", good["running_sha256"].upper()), ("schema", "other")]
    for key, value in cases:
        with pytest.raises(ValueError):
            simulate_aggregate(branching, {"Xa-400": 1000}, [0.0, 1.0, 2.0], stream, resume_from=dict(good, **{key: value}))
    with pytest.raises(ValueError):
        simulate_aggregate(branching, {"Xa-400": 1000}, [0.0, 1.0, 2.0], stream, resume_from=dict(good, extra=1))
    missing = {k: v for k, v in good.items() if k != "running_sha256"}
    with pytest.raises(ValueError):
        simulate_aggregate(branching, {"Xa-400": 1000}, [0.0, 1.0, 2.0], stream, resume_from=missing)
    assert stream.words_read == 0
    with pytest.raises(ValueError, match="word position"):   # a valid checkpoint, but not this stream
        simulate_aggregate(branching, {"Xa-400": 1000}, [0.0, 1.0, 2.0], stream, resume_from=good)
    seed = labelled_source(ns, "b2-engine", "checkpoint-malformed").manifest.seed
    right = WordStream(NumpyPCG64Source(seed, start_word=good["stream"]["word_position"]))
    done = simulate_aggregate(branching, {"Xa-400": 1000}, [0.0, 1.0, 2.0], right, resume_from=good)
    assert done.first_interval == 2 and right.words_read == 0 and done.populations.shape[0] == 1


def test_a_sink_cannot_alter_the_run(branching, ns):
    """A sink gets detached copies in a dict of its own (an independent review's finding)."""
    reference, _ = _run(branching, ns, "sink-safe", {"Xa-400": 5000})

    def replace_entries(interval, record):
        record["populations"] = record["populations"] * 0
        record["decays"] = None

    replaced, _ = _run(branching, ns, "sink-safe", {"Xa-400": 5000}, sink=replace_entries)
    assert replaced.digest == reference.digest
    for key in ("populations", "path_counts", "decays", "photons", "neutrons"):
        assert np.array_equal(getattr(replaced, key), getattr(reference, key)), key
    for key in ("populations", "path_counts", "decays", "photons"):
        def write(interval, record, key=key):
            record[key].fill(0)                  # its own copy: allowed, and without effect on the run
        written, _ = _run(branching, ns, "sink-safe", {"Xa-400": 5000}, sink=write)
        assert written.digest == reference.digest and np.array_equal(written.populations, reference.populations)


def test_neutron_totals_that_could_overflow_are_refused_before_any_draw(ns):
    from rnglab.nuclear import Channel, DecayScheme, Nuclide
    scheme = DecayScheme("two-neutrons", [
        Nuclide("N0", 150, 400, 1.0, (Channel("N1", Mode.BETA_DELAYED_NEUTRON, 1.0),)),
        Nuclide("N1", 151, 399, 2.0, (Channel("N2", Mode.BETA_DELAYED_NEUTRON, 1.0),)),
        Nuclide("N2", 152, 398, math.inf)])
    stream = WordStream(labelled_source(ns, "b2-engine", "neutron-overflow"))
    with pytest.raises(ValueError, match="neutron count"):
        simulate_aggregate(scheme, {"N0": 2 ** 62}, [0.0, 200.0], stream)
    assert stream.words_read == 0
    run = simulate_aggregate(scheme, {"N0": 2 ** 61}, [0.0, 200.0], stream)
    assert int(run.neutrons[0]) == 2 ** 62 and run.plan.max_neutrons_per_nucleus == 2


def test_mismatched_checkpoints_are_refused(branching, ns):
    checkpoints = []
    _run(branching, ns, "checkpoint-bad", {"Xa-400": 1000}, checkpoint_every=2, on_checkpoint=checkpoints.append)
    good = checkpoints[0]
    stream = WordStream(labelled_source(ns, "b2-engine", "unused"))
    for key, value in (("times", ["0.0", "1.0"]), ("initial", {"Xa-400": 999}), ("words_consumed", 7),
                       ("populations", [1001] + [0] * 8), ("scheme_sha256", "0" * 64)):
        bad = dict(good, **{key: value})
        with pytest.raises(ValueError):
            simulate_aggregate(branching, {"Xa-400": 1000}, TIMES, stream, resume_from=bad)


def test_streaming_without_kept_outputs(branching, ns):
    records = []
    kept, _ = _run(branching, ns, "sink", {"Xa-400": 10 ** 5})
    lean, _ = _run(branching, ns, "sink", {"Xa-400": 10 ** 5}, keep_outputs=False,
                   sink=lambda k, record: records.append((k, record)))
    assert lean.digest == kept.digest and lean.populations.shape[0] == 0
    assert [k for k, _ in records] == list(range(6))
    assert np.array_equal(np.array([r["populations"] for _, r in records]), kept.populations[1:])
    assert np.array_equal(lean.final_populations, kept.populations[-1])


def test_bad_inputs(branching, ns):
    stream = WordStream(labelled_source(ns, "b2-engine", "bad"))
    with pytest.raises(ValueError):
        simulate_aggregate(branching, {"Xa-400": INT64_MAX, "Xb-396": 1}, TIMES, stream)
    with pytest.raises(ValueError):
        simulate_aggregate(branching, {"nope": 1}, TIMES, stream)
    for times in ([0.0], [0.0, 0.0, 1.0], [0.0, 2.0, 1.0], [-1.0, 1.0], [0.0, float("inf")]):
        with pytest.raises(ValueError):
            simulate_aggregate(branching, {"Xa-400": 1}, times, stream)
    with pytest.raises(TypeError):
        simulate_aggregate(branching, {"Xa-400": 1}, TIMES, labelled_source(ns, "b2-engine", "bad"))


def test_regression_against_the_exact_multinomial_law(branching, ns):
    run, _ = _run(branching, ns, "law", {"Xa-400": 10 ** 6})
    for t, counts in zip(TIMES[1:], run.populations[1:]):
        assert multinomial_check(counts, state_probabilities(branching, "Xa-400", t)[0])["p_value"] > 1e-4


def test_huge_populations_are_split_and_conserved(branching, ns):
    n0 = 2 ** 56 + 12345
    run, _ = _run(branching, ns, "huge", {"Xa-400": n0})
    assert run.plan.parts == 9 and set(run.populations.sum(axis=1).tolist()) == {n0}
    p = state_probabilities(branching, "Xa-400", 12.0)[0]
    z = (run.populations[-1] - n0 * p) / np.sqrt(n0 * p * (1 - p))
    assert np.all(np.abs(z[p > 0]) < 5)


def test_layer_a_and_b2_means_agree(branching, ns):
    """Ensemble means of both engines against each other (fixed-seed regression, not evidence)."""
    a = np.array([simulate_decays(branching, {"Xa-400": 400}, 12.0,
                                  WordStream(labelled_source(ns, "b2-engine", "ens-a", r))).populations(TIMES)
                  for r in range(60)])
    b = np.array([_run(branching, ns, f"ens-b2-{r}", {"Xa-400": 400})[0].populations for r in range(60)])
    sd = np.sqrt((a.var(axis=0, ddof=1) + b.var(axis=0, ddof=1)) / 60)
    z = (a.mean(axis=0) - b.mean(axis=0))[sd > 0] / sd[sd > 0]
    assert np.all(np.abs(z) < 5)
