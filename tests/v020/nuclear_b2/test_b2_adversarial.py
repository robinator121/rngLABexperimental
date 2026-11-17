"""Adversarial tests for the three defect classes an independent review confirmed in B2.

1. State aliasing: nothing a caller, a sink or a checkpoint callback does to the objects it holds
   may change a run, its outputs or a later run.
2. Checkpoint integrity: a corrupted, edited or misplaced checkpoint is refused, never repaired.
3. Floating-point boundaries: the binomial inversion agrees with exact midpoint inversion at the edges
   of its domain (extreme uniforms, regime switches, the part limit, vanishing probabilities). Next to a
   step of F, where the evaluated F may err by a few ulps, see test_b2_guards.

Fixed seeds of the B2 test namespace; exact references are rational arithmetic.
"""
import copy
import json
import math
from fractions import Fraction

import numpy as np
import pytest

from rnglab.nuclear import NumpyPCG64Source, WordStream, labelled_source
from rnglab.nuclear.aggregate import (
    AggregatePlan, _checkpoint_checksum, simulate_aggregate, simulate_aggregate_reference, stream_identity,
)
from rnglab.nuclear.binomial import (
    PART_MAX, TABLE_MEAN_MAX, UPPER_TAIL_MAX_N, UPPER_TAIL_SWITCH, binomial_draw, binomial_inverse,
)
from rnglab.nuclear.paths import exact_channel_probabilities
from rnglab14.sources import Xoshiro256PlusPlusSource

TIMES = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]
OUTPUTS = ("populations", "path_counts", "decays", "photons", "neutrons")
HALF_CELL = Fraction(1, 2 ** 54)


def _source(ns, label, start_word=0):
    seed = labelled_source(ns, "b2-adversarial", label).manifest.seed
    return NumpyPCG64Source(seed, start_word=start_word)


def _run(scheme, ns, label, initial, times=TIMES, **kw):
    return simulate_aggregate(scheme, initial, times, WordStream(_source(ns, label)), **kw)


def _same(a, b):
    return a.digest == b.digest and all(np.array_equal(getattr(a, k), getattr(b, k)) for k in OUTPUTS)


# ---------------------------------------------------------------------------------------------
# 1. State aliasing
# ---------------------------------------------------------------------------------------------

def test_a_sink_that_forces_writes_cannot_reach_the_engine(branching, ns):
    reference = _run(branching, ns, "sink", {"Xa-400": 20000})
    hoard = []

    def vandal(interval, record):
        for value in record.values():
            if isinstance(value, np.ndarray):
                value.setflags(write=True)      # a determined sink can re-enable writing on what it was given
                value[...] = -7
        hoard.append(record)
        for earlier in hoard:                   # and keep scribbling on earlier intervals' records
            for key in list(earlier):
                earlier[key] = None

    attacked = _run(branching, ns, "sink", {"Xa-400": 20000}, sink=vandal)
    assert _same(attacked, reference) and np.array_equal(attacked.final_populations, reference.final_populations)
    lean = _run(branching, ns, "sink", {"Xa-400": 20000}, sink=vandal, keep_outputs=False)
    assert lean.digest == reference.digest and np.array_equal(lean.final_populations, reference.final_populations)


def test_sink_records_are_detached_from_each_other_and_from_the_results(branching, ns):
    records = []
    run = _run(branching, ns, "detached", {"Xa-400": 5000}, sink=lambda k, record: records.append(record))
    arrays = [value for record in records for value in record.values() if isinstance(value, np.ndarray)]
    results = [getattr(run, key) for key in OUTPUTS] + [run.final_populations, run.plan.category_steps]
    for i, a in enumerate(arrays):
        assert not any(np.shares_memory(a, b) for b in arrays[i + 1:])
        assert not any(np.shares_memory(a, b) for b in results)
    assert np.array_equal(np.array([r["populations"] for r in records]), run.populations[1:])


def test_callers_inputs_are_not_aliased(branching, ns):
    times = np.array(TIMES)
    initial = {"Xa-400": 8000}
    reference = _run(branching, ns, "inputs", {"Xa-400": 8000}, times=list(TIMES))
    plan = AggregatePlan(branching, times, 8000)

    def meddle(interval, record):               # the caller's own objects change while the run is under way
        times[:] = np.arange(times.size) * 100.0
        initial["Xa-400"] = 1
        initial["Xb-396"] = 10 ** 6

    run = simulate_aggregate(branching, initial, TIMES, WordStream(_source(ns, "inputs")), plan=plan, sink=meddle)
    assert _same(run, reference) and run.initial == {n.name: (8000 if n.name == "Xa-400" else 0)
                                                     for n in branching.nuclides}
    assert plan.times.tolist() == TIMES and not np.shares_memory(plan.times, times)


def test_plans_and_results_are_read_only(branching, ns):
    run = _run(branching, ns, "read-only", {"Xa-400": 3000})
    plan = run.plan
    source = next(iter(plan.interval_plans[float(plan.deltas[0])].values()))
    for array in (plan.times, plan.deltas, plan.category_end, plan.category_steps, plan.first_category,
                  source.p_small, source.flip, run.populations, run.path_counts, run.decays, run.photons,
                  run.neutrons, run.final_populations):
        with pytest.raises(ValueError, match="read-only"):
            array[...] = 0


def test_forced_writes_to_returned_arrays_do_not_leak_into_later_runs(branching, ns):
    plan = AggregatePlan(branching, TIMES, 6000)
    first = simulate_aggregate(branching, {"Xa-400": 6000}, TIMES, WordStream(_source(ns, "leak")), plan=plan)
    digest, populations = first.digest, first.populations.copy()
    for key in OUTPUTS + ("final_populations",):
        array = getattr(first, key)
        array.setflags(write=True)
        array[...] = 123
    second = simulate_aggregate(branching, {"Xa-400": 6000}, TIMES, WordStream(_source(ns, "leak")), plan=plan)
    assert second.digest == digest and np.array_equal(second.populations, populations)


def test_checkpoints_are_detached_from_the_engine(branching, ns):
    reference = _run(branching, ns, "checkpoint-objects", {"Xa-400": 9000})
    seen = []

    def vandal(checkpoint):
        seen.append(copy.deepcopy(checkpoint))
        checkpoint["populations"][:] = [0] * len(checkpoint["populations"])
        checkpoint["layout"]["parts"] = 99
        checkpoint["stream"]["word_position"] = -1
        checkpoint["initial"].clear()
        checkpoint["next_interval"] = 0

    attacked = _run(branching, ns, "checkpoint-objects", {"Xa-400": 9000}, checkpoint_every=1, on_checkpoint=vandal)
    assert _same(attacked, reference) and len(seen) == 6
    for checkpoint in seen:                      # every checkpoint was valid when it was handed over
        body = {k: v for k, v in checkpoint.items() if k != "checkpoint_sha256"}
        assert checkpoint["checkpoint_sha256"] == _checkpoint_checksum(body)


def test_both_engines_return_independent_read_only_results(competing, ns):
    fast = _run(competing, ns, "engines", {"Pa-300": 4000})
    slow = simulate_aggregate_reference(competing, {"Pa-300": 4000}, TIMES, WordStream(_source(ns, "engines")))
    assert _same(fast, slow)
    for key in OUTPUTS:
        assert not np.shares_memory(getattr(fast, key), getattr(slow, key))
        with pytest.raises(ValueError, match="read-only"):
            getattr(slow, key)[...] = 0


# ---------------------------------------------------------------------------------------------
# 2. Checkpoint integrity
# ---------------------------------------------------------------------------------------------

def _checkpointed(scheme, ns, label, initial, times=TIMES, **kw):
    checkpoints = []
    run = _run(scheme, ns, label, initial, times, checkpoint_every=1, on_checkpoint=checkpoints.append, **kw)
    return run, [json.loads(json.dumps(c)) for c in checkpoints]


def _resume(scheme, ns, label, initial, checkpoint, times=TIMES, stream=None, **kw):
    stream = stream or WordStream(_source(ns, label, start_word=checkpoint["stream"]["word_position"]))
    return simulate_aggregate(scheme, initial, times, stream, resume_from=checkpoint, **kw), stream


def test_resume_is_bit_identical_through_json_at_every_interval_with_huge_counts(branching, ns):
    initial = {"Xa-400": 2 ** 54 + 3}           # 17 parts; counts beyond 2**53 must survive JSON exactly
    full, checkpoints = _checkpointed(branching, ns, "huge", initial)
    assert full.plan.parts == -(-(2 ** 54 + 3) // PART_MAX) == 17 and len(checkpoints) == 6
    for cut, checkpoint in enumerate(checkpoints, start=1):
        assert sum(checkpoint["populations"]) == 2 ** 54 + 3
        rest, stream = _resume(branching, ns, "huge", initial, checkpoint)
        assert rest.digest == full.digest and rest.first_interval == cut
        assert np.array_equal(rest.populations, full.populations[cut:])
        assert np.array_equal(rest.decays, full.decays[cut:])
        assert stream.words_read == (6 - cut) * full.plan.words_per_interval


def test_resume_at_the_final_interval_reads_nothing(branching, ns):
    full, checkpoints = _checkpointed(branching, ns, "final", {"Xa-400": 7000})
    rest, stream = _resume(branching, ns, "final", {"Xa-400": 7000}, checkpoints[-1])
    assert rest.digest == full.digest and stream.words_read == 0 and rest.decays.shape[0] == 0
    assert np.array_equal(rest.final_populations, full.final_populations)


def _edits(checkpoint, words_per_interval):
    """Type-valid edits of a checkpoint: (name, edited copy, still refused if its checksum is recomputed)."""
    pops = checkpoint["populations"]
    moved = list(pops)
    donor = max(range(len(moved)), key=lambda i: moved[i])
    receiver = (donor + 1) % len(moved)
    moved[donor] -= 1
    moved[receiver] += 1
    k = checkpoint["next_interval"]
    yield "one nucleus moved (total conserved)", dict(checkpoint, populations=moved), False
    yield "populations permuted (total conserved)", dict(checkpoint, populations=pops[1:] + pops[:1]), False
    yield "another interval, word count consistent", dict(checkpoint, next_interval=k + 1,
                                                          words_consumed=(k + 1) * words_per_interval), False
    yield "another interval only", dict(checkpoint, next_interval=k - 1), True
    yield "another word count only", dict(checkpoint, words_consumed=checkpoint["words_consumed"] + 1), True
    flipped = ("0" if checkpoint["running_sha256"][0] != "0" else "1") + checkpoint["running_sha256"][1:]
    yield "running digest changed", dict(checkpoint, running_sha256=flipped), False
    yield "initial changed", dict(checkpoint, initial=dict(checkpoint["initial"], **{"Xa-400": 1})), True
    yield "grid changed", dict(checkpoint, times=checkpoint["times"][:-1] + ["13.0"]), True
    yield "layout changed", dict(checkpoint, layout=dict(checkpoint["layout"], parts=2)), True
    yield "scheme changed", dict(checkpoint, scheme_sha256="0" * 64), True
    yield "seed changed", dict(checkpoint, stream=dict(checkpoint["stream"], seed=checkpoint["stream"]["seed"] + 1)), True
    yield "position changed", dict(checkpoint, stream=dict(checkpoint["stream"],
                                                           word_position=checkpoint["stream"]["word_position"] + 1)), True


def test_every_type_valid_edit_of_a_checkpoint_is_refused(branching, ns):
    full, checkpoints = _checkpointed(branching, ns, "edits", {"Xa-400": 10 ** 5})
    good = checkpoints[2]
    edits = list(_edits(good, full.plan.words_per_interval))
    edits.append(("checksum changed", dict(good, checkpoint_sha256="f" * 64), True))
    assert len(edits) == 13
    for name, edited, _ in edits:
        stream = WordStream(_source(ns, "edits", start_word=good["stream"]["word_position"]))
        with pytest.raises(ValueError):
            simulate_aggregate(branching, {"Xa-400": 10 ** 5}, TIMES, stream, resume_from=edited)
        assert stream.words_read == 0, name      # a refused resume reads no word
    rest, _ = _resume(branching, ns, "edits", {"Xa-400": 10 ** 5}, good)
    assert rest.digest == full.digest


def test_a_checkpoint_is_bound_to_its_generator_and_position(branching, ns):
    full, checkpoints = _checkpointed(branching, ns, "stream", {"Xa-400": 50000})
    good = checkpoints[3]
    position = good["stream"]["word_position"]
    assert position == 4 * full.plan.words_per_interval == good["words_consumed"]
    seed = _source(ns, "stream").manifest.seed
    wrong = {"fresh stream at the origin": NumpyPCG64Source(seed),
             "one word early": NumpyPCG64Source(seed, start_word=position - 1),
             "one word late": NumpyPCG64Source(seed, start_word=position + 1),
             "another seed": NumpyPCG64Source(seed + 1, start_word=position),
             "another generator": Xoshiro256PlusPlusSource(seed=seed)}
    for name, source in wrong.items():
        with pytest.raises(ValueError, match="word position"):
            simulate_aggregate(branching, {"Xa-400": 50000}, TIMES, WordStream(source), resume_from=good)
        assert source.consumed_bytes == 0, name
    advanced = NumpyPCG64Source(seed)            # reading up to the position is as good as jumping there
    advanced.read_bytes(8 * position)
    rest = simulate_aggregate(branching, {"Xa-400": 50000}, TIMES, WordStream(advanced), resume_from=good)
    assert rest.digest == full.digest
    misaligned = NumpyPCG64Source(seed, start_word=position - 1)
    misaligned.read_bytes(11)                    # 8 bytes to the position, 3 into the next word
    with pytest.raises(ValueError):
        simulate_aggregate(branching, {"Xa-400": 50000}, TIMES, WordStream(misaligned), resume_from=good)


def test_resume_with_a_generator_that_cannot_jump(branching, ns):
    seed = _source(ns, "xoshiro").manifest.seed
    checkpoints = []
    full = simulate_aggregate(branching, {"Xa-400": 3000}, TIMES, WordStream(Xoshiro256PlusPlusSource(seed=seed)),
                              checkpoint_every=2, on_checkpoint=checkpoints.append)
    checkpoint = json.loads(json.dumps(checkpoints[0]))
    replayed = Xoshiro256PlusPlusSource(seed=seed)
    replayed.read_bytes(8 * checkpoint["stream"]["word_position"])
    rest = simulate_aggregate(branching, {"Xa-400": 3000}, TIMES, WordStream(replayed), resume_from=checkpoint)
    assert rest.digest == full.digest and rest.first_interval == 2


def test_a_forged_checkpoint_with_a_recomputed_checksum_is_the_documented_limit(branching, ns):
    """The checksum is an integrity check, not an authentication: a consistent forgery resumes.
    The final digest then differs from the uninterrupted run, which is how it is found."""
    full, checkpoints = _checkpointed(branching, ns, "forgery", {"Xa-400": 10 ** 5})
    good = checkpoints[2]
    accepted = []
    for name, edited, refused_anyway in _edits(good, full.plan.words_per_interval):
        forged = dict(edited)
        forged["checkpoint_sha256"] = _checkpoint_checksum({k: v for k, v in forged.items() if k != "checkpoint_sha256"})
        stream = WordStream(_source(ns, "forgery", start_word=good["stream"]["word_position"]))
        if refused_anyway:                       # the other validations still catch it
            with pytest.raises(ValueError):
                simulate_aggregate(branching, {"Xa-400": 10 ** 5}, TIMES, stream, resume_from=forged)
            assert stream.words_read == 0, name
        else:                                    # it resumes, and the final digest gives it away
            rest = simulate_aggregate(branching, {"Xa-400": 10 ** 5}, TIMES, stream, resume_from=forged)
            assert rest.digest != full.digest, name
            accepted.append(name)
    assert len(accepted) == 4


def test_checkpoints_of_other_runs_are_refused(branching, competing, ns):
    _, checkpoints = _checkpointed(branching, ns, "other", {"Xa-400": 4000})
    good = checkpoints[0]

    def stream():
        return WordStream(_source(ns, "other", start_word=good["stream"]["word_position"]))

    with pytest.raises(ValueError):              # another population
        simulate_aggregate(branching, {"Xa-400": 4001}, TIMES, stream(), resume_from=good)
    with pytest.raises(ValueError):              # another grid
        simulate_aggregate(branching, {"Xa-400": 4000}, TIMES[:-1] + [13.0], stream(), resume_from=good)
    with pytest.raises(ValueError):              # another part size, hence another layout
        simulate_aggregate(branching, {"Xa-400": 4000}, TIMES, stream(), resume_from=good, part_max=1000)
    with pytest.raises(ValueError):              # another scheme
        simulate_aggregate(competing, {"Pa-300": 4000}, TIMES, stream(), resume_from=good)
    for junk in (object(), [], "checkpoint", 3):
        with pytest.raises(ValueError):
            simulate_aggregate(branching, {"Xa-400": 4000}, TIMES, stream(), resume_from=junk)
    nan = dict(good, populations=[float("nan")] * len(good["populations"]))
    with pytest.raises(ValueError):
        simulate_aggregate(branching, {"Xa-400": 4000}, TIMES, stream(), resume_from=nan)


def test_checkpoint_arguments_are_validated(branching, ns):
    stream = WordStream(_source(ns, "arguments"))
    for kw in (dict(checkpoint_every=1), dict(on_checkpoint=print), dict(checkpoint_every=-1, on_checkpoint=print),
               dict(checkpoint_every=True, on_checkpoint=print), dict(checkpoint_every=1.0, on_checkpoint=print)):
        with pytest.raises(ValueError):
            simulate_aggregate(branching, {"Xa-400": 10}, TIMES, stream, **kw)
    plan = AggregatePlan(branching, TIMES, 10, part_max=7)
    with pytest.raises(ValueError):
        simulate_aggregate(branching, {"Xa-400": 10}, TIMES, stream, plan=plan, part_max=8)
    with pytest.raises(TypeError):
        simulate_aggregate(branching, {"Xa-400": 10}, TIMES, stream, plan="plan")
    assert stream.words_read == 0
    assert simulate_aggregate(branching, {"Xa-400": 10}, TIMES, stream, plan=plan).plan.parts == 2


def test_stream_identity_tracks_the_absolute_position(ns):
    source = _source(ns, "identity", start_word=40)
    stream = WordStream(source)
    assert stream_identity(stream)["word_position"] == 40 and stream_identity(stream)["byte_aligned"]
    stream.words(3)
    assert stream_identity(stream)["word_position"] == 43
    source.read_bytes(1)
    assert not stream_identity(stream)["byte_aligned"]


# ---------------------------------------------------------------------------------------------
# 3. Floating-point boundaries of the binomial inversion
# ---------------------------------------------------------------------------------------------

def _exact(n, p, us):
    """min{x : F(x) > u + 2**-54}, F in exact integer arithmetic (p = m / d exactly, everything scaled by d**n)."""
    m, d = Fraction(p).numerator, Fraction(p).denominator
    q = d - m
    cdf, total, a, b, c = [], 0, 1, q ** n, 1          # a = m**x, b = q**(n - x), c = C(n, x)
    for x in range(n + 1):
        total += c * a * b
        cdf.append(total)
        if x < n:
            a, b, c = a * m, b // q, c * (n - x) // (x + 1)
    scale = d ** n
    out = []
    for u in us:
        v = Fraction(u) + HALF_CELL
        out.append(next(x for x in range(n + 1) if cdf[x] * v.denominator > v.numerator * scale))
    return out


def _draws(n, p, us):
    return binomial_inverse(np.full(len(us), n), np.full(len(us), p), np.array(us)).tolist()


ULP = 2.0 ** -53
EXTREMES = [0.0, ULP, 2 * ULP, 0.5 - ULP, 0.5, 0.5 + 2 * ULP, 1.0 - 2 * ULP, 1.0 - ULP]


def test_dyadic_steps_fall_on_the_right_side():
    """n = 2, p = 1/2: F = 1/4, 3/4. A uniform one ulp below a step stays below it; at the step it moves up."""
    us = [0.25 - ULP, 0.25, 0.75 - ULP, 0.75, 0.0, 1.0 - ULP]
    assert _draws(2, 0.5, us) == [0, 1, 1, 2, 0, 2] == _exact(2, 0.5, us)
    us = [0.5 - ULP, 0.5]
    assert _draws(1, 0.5, us) == [0, 1] == _exact(1, 0.5, us)


@pytest.mark.parametrize("n, p", [(512, 0.5), (513, 0.5), (1024, 0.25), (1025, 0.25), (2000, 0.128), (2001, 0.128)])
def test_both_sides_of_the_regime_switch_are_exact_at_the_extremes(n, p):
    """n*p = 256 is the last table-regime mean; one more trial switches to the incomplete beta."""
    us = EXTREMES + np.linspace(0.0, 1.0 - ULP, 41).tolist()
    assert _draws(n, p, us) == _exact(n, p, us)
    assert (n * p <= TABLE_MEAN_MAX) == (n in (512, 1024, 2000))


@pytest.mark.parametrize("n, p", [(600, 0.45), (3000, 0.1)])
def test_the_upper_tail_switch_is_seamless(n, p):
    """In the incomplete-beta regime the upper tail is consulted for 1 - u < 2**-10: check both sides."""
    assert n * p > TABLE_MEAN_MAX
    edge = 1.0 - UPPER_TAIL_SWITCH
    us = [edge - 2 * ULP, edge - ULP, edge, edge + ULP, edge + 2 * ULP, 1.0 - 2.0 ** -20, 1.0 - 2.0 ** -40,
          1.0 - 2 * ULP, 1.0 - ULP]
    assert _draws(n, p, us) == _exact(n, p, us)


def test_the_population_limit_of_the_upper_tail_is_continuous():
    """Above n = 1e9 only the lower tail is used. Draws just below and above the limit agree to a
    fraction of a standard deviation, and the last cell stays in the far tail (documented limit)."""
    p = 0.3
    us = np.array([1.0 - 2.0 ** -11, 1.0 - 2.0 ** -20, 1.0 - 2.0 ** -30, 1.0 - ULP])
    z = []
    for n in (UPPER_TAIL_MAX_N, UPPER_TAIL_MAX_N + 1):
        x = binomial_inverse(np.full(us.size, n), np.full(us.size, p), us)
        z.append((x - n * p) / math.sqrt(n * p * (1 - p)))
    assert np.all(np.abs(z[0][:3] - z[1][:3]) < 0.02)
    assert 7.0 < z[0][3] < 9.5 and 7.0 < z[1][3] < 9.5


def test_counts_at_and_around_the_part_limit():
    us = np.array([0.0, 0.25, 0.5, 0.75, 1.0 - ULP])
    for n in (PART_MAX - 1, PART_MAX):
        x = binomial_inverse(np.full(us.size, n), np.full(us.size, 0.25), us)
        assert np.all(np.diff(x) > 0) and x[0] > 0 and x[-1] < n
    with pytest.raises(ValueError):
        binomial_inverse(PART_MAX + 1, 0.25, 0.5)
    for u1, u2 in ((0.3, 0.1), (0.3, 0.9), (1.0 - ULP, 0.0)):
        whole = int(binomial_inverse(PART_MAX, 0.25, u1)) + int(binomial_inverse(1, 0.25, u2))
        assert binomial_draw(PART_MAX + 1, 0.25, False, np.array([u1, u2])) == whole
        assert binomial_draw(PART_MAX + 1, 0.25, True, np.array([u1, u2])) == PART_MAX + 1 - whole


def test_vanishing_probabilities():
    us = [0.0, 0.5, 1.0 - ULP]
    for p in (5e-324, 1e-320, 1e-300, 2.0 ** -108):
        assert _draws(PART_MAX, p, us) == [0, 0, 0]
    # P(X >= 1) = 1 - (1 - p)**n is about 2**-53 > 2**-54 for n = 2**50, p = 2**-103: only the last cell gives 1.
    assert _draws(PART_MAX, 2.0 ** -103, us) == [0, 0, 1]
    assert _draws(PART_MAX, 2.0 ** -106, us) == [0, 0, 0]   # about 2**-56 < 2**-54: below the grid's resolution
    assert _draws(1, 2.0 ** -53, [1.0 - 2 * ULP, 1.0 - ULP]) == [0, 1]
    assert _draws(7, 0.0, us) == [0, 0, 0] and _draws(0, 0.3, us) == [0, 0, 0]


def test_probability_one_half_and_its_neighbours_in_the_engine(ns):
    """A stay probability of exactly 1/2 is sampled directly; one ulp above it, through the complement."""
    from rnglab.nuclear import Channel, DecayScheme, Mode, Nuclide
    scheme = DecayScheme("half", [Nuclide("Ha-400", 150, 400, 1.0, (Channel("Hb-400", Mode.BETA_MINUS, 1.0),)),
                                  Nuclide("Hb-400", 151, 400, math.inf)])
    flips = []
    deltas = [1.1, math.nextafter(1.0, 2.0), 1.0, math.nextafter(1.0, 0.0), 1.0 - 2.0 ** -52, 0.9]
    for delta in deltas:
        plan = AggregatePlan(scheme, [0.0, delta], 10 ** 6)
        source = plan.interval_plans[delta][0]
        flips.append(bool(source.flip[0]))
        fast = simulate_aggregate(scheme, {"Ha-400": 10 ** 6}, [0.0, delta], WordStream(_source(ns, "half")), plan=plan)
        slow = simulate_aggregate_reference(scheme, {"Ha-400": 10 ** 6}, [0.0, delta], WordStream(_source(ns, "half")))
        assert fast.digest == slow.digest and int(fast.populations[1].sum()) == 10 ** 6
        stay = 2.0 ** -delta
        assert abs(int(fast.populations[1][0]) - 10 ** 6 * stay) < 5 * math.sqrt(10 ** 6 * stay * (1 - stay))
        assert 0.0 < float(source.p_small[0]) <= 0.5
    # exactly 1/2 at delta = 1 is sampled directly; shorter intervals stay with higher probability and flip
    assert flips[0] is False and flips[2] is False and flips[-1] is True and flips == sorted(flips)


def test_uniform_domain():
    assert int(binomial_inverse(10, 0.5, math.nextafter(1.0, 0.0))) == 10
    assert int(binomial_inverse(10, 0.5, -0.0)) == 0
    for bad in (1.0, math.nextafter(1.0, 2.0), -ULP, float("nan"), float("inf")):
        with pytest.raises(ValueError):
            binomial_inverse(10, 0.5, bad)


def test_branching_laws_of_both_engines_agree_at_float_precision(branching, competing):
    """B2 uses the exact differences of the thresholds Layer A compares uniforms with."""
    for scheme in (branching, competing):
        for i, nuclide in enumerate(scheme.nuclides):
            exact = [Fraction(x) for x in exact_channel_probabilities(nuclide)]
            if not nuclide.channels:
                assert exact == [1]              # no channel: the single edge pair (0, 1)
                continue
            assert sum(exact) == 1
            edges = [0.0] + [float(x) for x in scheme.cum_thresholds[i][:len(nuclide.channels) - 1]] + [1.0]
            for k, probability in enumerate(exact):
                assert probability == Fraction(edges[k + 1]) - Fraction(edges[k])
                assert float(probability) == pytest.approx(nuclide.channel_probabilities()[k], abs=2.0 ** -52)
