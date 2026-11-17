"""Regression tests for the third review of B2 (architecture and performance roles).

Each test reproduces a confirmed finding and checks its correction:
- a plan, a scheme or a result cannot be changed under a run: plans and results are immutable, a
  scheme that changed since planning is refused, and conservation is checked after every interval;
- a checkpoint identifies its generator by every parameter that fixes the sequence (another PCG64
  stream id is refused); sources whose sequence a checkpoint cannot identify are refused before
  any draw;
- a run owns its source: a read by a callback stops it; manifests come from the run's own record;
- both engines leave the stream at an interval boundary when a draw fails;
- comparison tests count in Python integers, refuse non-integer counts, and never read a
  degenerate degree of freedom or a NaN tail as p = 1;
- once-only record runs reserve an attempt atomically, record how it ended and never overwrite;
- every namespace a protocol names (at any depth) is reserved, statically and at run time;
- checkpoints do not grow with the grid; category counts are bounded before enumeration;
- benchmark timings are medians with their samples; calibration forces one thread per worker.
Draws only from the B2 test namespace (or none at all).
"""
import dataclasses
import importlib.util
import json
import math
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pytest

import rnglab.nuclear.aggregate as aggregate_module
from rnglab.nuclear import Channel, DecayScheme, Mode, Nuclide, NumpyPCG64Source, WordStream, labelled_source
from rnglab.nuclear import binomial as binomial_module
from rnglab.nuclear.aggregate import (
    REPLAYABLE, AggregatePlan, aggregate_manifest, category_counts, replay_identity, simulate_aggregate,
    simulate_aggregate_reference,
)
from rnglab.nuclear.comparison import (
    allocation_check, fixed_cells_chisquare, homogeneity_check, scaled_chisquare_variance_check,
    variance_ratio_f_check,
)
from rnglab.nuclear.paths import path_categories
from rnglab.nuclear.precursors import PrecursorGroup, precursor_scheme
from rnglab.seeds import namespace_seed
from rnglab14.sources import MersenneTwisterSource, PCG64Source, Xoshiro256PlusPlusSource

from .conftest import B2_TEST_NAMESPACE, ReservedNamespaceError, reserved_namespaces

REPO = Path(__file__).resolve().parents[3]
TIMES = [0.0, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0]


def _load(name: str, relative: str):
    if name not in sys.modules:
        spec = importlib.util.spec_from_file_location(name, REPO / relative)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
    return sys.modules[name]


def _verification():
    return _load("b2_verify_run_verification", "research/v020-nuclear-b2-verify/run_verification.py")


def _calibration():
    return _load("b2_integrity_calibration", "research/v020-nuclear-b2-verify/calibration.py")


def _stream(ns, *label):
    return WordStream(labelled_source(ns, "integrity", *label))


# --- 1. Plans, schemes and results cannot change under a run ----------------------------------------------

def test_a_sink_cannot_change_the_plan_it_runs_with(branching, ns):
    """The review's probe: p.sources.clear() in a sink emptied a nuclide and broke conservation."""
    reference = simulate_aggregate(branching, {"Xa-400": 10000}, TIMES, _stream(ns, "plan"))
    plan = AggregatePlan(branching, TIMES, 10000)
    delta = float(plan.deltas[0])
    attempts = {
        "sources.clear()": lambda: plan.sources.clear(),
        "slot_offset[0] = 5": lambda: plan.slot_offset.__setitem__(0, 5),
        "interval_plans[delta][0] = None": lambda: plan.interval_plans[delta].__setitem__(0, None),
        "interval_plans.clear()": lambda: plan.interval_plans.clear(),
        "diagnostics[...] = 0": lambda: plan.diagnostics.__setitem__("tree_comparisons", 0),
        "categories[0].append": lambda: plan.categories[0].append(None),
        "total = 1": lambda: setattr(plan, "total", 1),
        "sources rebound": lambda: setattr(plan, "sources", ()),
        "del categories": lambda: delattr(plan, "categories"),
        "p_small made writable": lambda: plan.interval_plans[delta][0].p_small.setflags(write=True),
        "law replaced": lambda: setattr(plan.interval_plans[delta][0], "law", None),
    }
    succeeded = []

    def vandal(interval, record):
        for name, attempt in attempts.items():
            try:
                attempt()
            except (AttributeError, TypeError, ValueError, dataclasses.FrozenInstanceError):
                continue
            succeeded.append(name)

    run = simulate_aggregate(branching, {"Xa-400": 10000}, TIMES, _stream(ns, "plan"), plan=plan, sink=vandal)
    assert succeeded == []
    assert run.digest == reference.digest
    again = simulate_aggregate(branching, {"Xa-400": 10000}, TIMES, _stream(ns, "plan"), plan=plan)
    assert again.digest == reference.digest       # the reused plan is unchanged


def test_a_scheme_changed_after_planning_is_refused(branching, ns):
    plan = AggregatePlan(branching, TIMES, 4000)
    run = simulate_aggregate(branching, {"Xa-400": 4000}, TIMES, _stream(ns, "scheme"), plan=plan)
    first = branching.nuclides[0]
    branching.nuclides = (dataclasses.replace(first, half_life=2 * first.half_life),) + branching.nuclides[1:]
    stream = _stream(ns, "scheme")
    with pytest.raises(ValueError, match="changed since the plan"):
        simulate_aggregate(branching, {"Xa-400": 4000}, TIMES, stream, plan=plan)
    assert stream.words_read == 0
    with pytest.raises(ValueError, match="changed since the run"):
        aggregate_manifest(experiment="e", run_id="r", run=run, code={})
    assert plan.nuclides[0] is first               # the plan kept its own snapshot


def test_results_are_immutable(branching, ns):
    run = simulate_aggregate(branching, {"Xa-400": 3000}, TIMES, _stream(ns, "result"))
    with pytest.raises(TypeError):
        run.initial["Xa-400"] = 5
    with pytest.raises(dataclasses.FrozenInstanceError):
        run.digest = "0" * 64
    provenance = run.stream_provenance()
    provenance["words_read"] = -1
    assert run.stream_provenance()["words_read"] == run.words_read


def _corrupt_vector(monkeypatch, how):
    real = aggregate_module._vector_binomials

    def corrupted(n, p_small, flip, uniforms, parts, part_max):
        out = real(n, p_small, flip, uniforms, parts, part_max)
        return n + 1 if how == "too many" else out
    monkeypatch.setattr(aggregate_module, "_vector_binomials", corrupted)


def test_conservation_is_checked_after_every_interval(branching, ns, monkeypatch):
    _corrupt_vector(monkeypatch, "too many")     # a draw above its trials makes the remainder negative
    with pytest.raises(RuntimeError, match="conservation"):
        simulate_aggregate(branching, {"Xa-400": 3000}, TIMES, _stream(ns, "conservation"))
    monkeypatch.undo()
    real = aggregate_module.binomial_draw
    monkeypatch.setattr(aggregate_module, "binomial_draw", lambda n, *args, **kw: n + 1 if n else real(n, *args, **kw))
    with pytest.raises(RuntimeError, match="conservation"):
        simulate_aggregate_reference(branching, {"Xa-400": 3000}, TIMES, _stream(ns, "conservation"))


# --- 2. Checkpoints identify the whole sequence; other sources are refused ---------------------------------

def test_a_checkpoint_refuses_another_pcg64_stream_id(branching, ns):
    """The review's probe: stream_id 2 resumed a checkpoint of stream_id 1 (same seed, same position)."""
    seed = namespace_seed(ns, "integrity", "pcg64")
    checkpoints = []
    full = simulate_aggregate(branching, {"Xa-400": 10000}, TIMES, WordStream(PCG64Source(seed, stream_id=1)),
                              checkpoint_every=1, on_checkpoint=checkpoints.append)
    checkpoint = json.loads(json.dumps(checkpoints[2]))
    assert checkpoint["stream"]["parameters"] == {"seed": seed, "stream_id": 1}
    position = checkpoint["stream"]["word_position"]
    wrong = PCG64Source(seed, stream_id=2)
    wrong.read_bytes(8 * position)
    with pytest.raises(ValueError, match="word position"):
        simulate_aggregate(branching, {"Xa-400": 10000}, TIMES, WordStream(wrong), resume_from=checkpoint)
    right = PCG64Source(seed, stream_id=1)
    right.read_bytes(8 * position)
    rest = simulate_aggregate(branching, {"Xa-400": 10000}, TIMES, WordStream(right), resume_from=checkpoint)
    assert rest.digest == full.digest


class _Subclass(NumpyPCG64Source):
    """Same constructor, another class: a subclass could change the sequence, so it is not replayable."""


def test_sources_a_checkpoint_cannot_identify_are_refused_before_any_draw(branching, ns):
    seed = namespace_seed(ns, "integrity", "other-sources")
    assert replay_identity(NumpyPCG64Source(seed, start_word=7)) == {
        "generator": REPLAYABLE[0], "parameters": {"seed": seed}, "start_word": 7}
    assert replay_identity(Xoshiro256PlusPlusSource(seed=seed))["generator"] == REPLAYABLE[2]
    for source in (MersenneTwisterSource(seed=seed), _Subclass(seed)):
        assert replay_identity(source) is None
        with pytest.raises(ValueError, match="not a replayable generator"):
            simulate_aggregate(branching, {"Xa-400": 2000}, TIMES, WordStream(source), checkpoint_every=1,
                               on_checkpoint=lambda checkpoint: None)
        assert source.consumed_bytes == 0
    plain = simulate_aggregate(branching, {"Xa-400": 2000}, TIMES, WordStream(MersenneTwisterSource(seed=seed)))
    assert plain.stream_provenance()["replay_identity"] is None   # runs without checkpoints are unaffected


# --- 3. A run owns its source; manifests come from the run's record ----------------------------------------

@pytest.mark.parametrize("thief", ["another stream on the source", "the run's own stream", "the source itself"])
def test_a_callback_that_reads_the_run_source_stops_the_run(branching, ns, thief):
    stream = _stream(ns, "exclusive")

    def read():
        if thief == "another stream on the source":
            WordStream(stream.source).raw(1)
        elif thief == "the run's own stream":
            stream.raw(1)
        else:
            stream.source.read_bytes(8)

    with pytest.raises(RuntimeError, match="outside its layout"):
        simulate_aggregate(branching, {"Xa-400": 5000}, TIMES, stream, sink=lambda interval, record: read())
    stream = _stream(ns, "exclusive")
    with pytest.raises(RuntimeError, match="outside its layout"):
        simulate_aggregate(branching, {"Xa-400": 5000}, TIMES, stream, checkpoint_every=2,
                           on_checkpoint=lambda checkpoint: read())


def test_a_callback_may_use_a_source_of_its_own(branching, ns):
    reference = simulate_aggregate(branching, {"Xa-400": 5000}, TIMES, _stream(ns, "own"))
    other = _stream(ns, "callback's own")
    run = simulate_aggregate(branching, {"Xa-400": 5000}, TIMES, _stream(ns, "own"),
                             sink=lambda interval, record: other.raw(3))
    assert run.digest == reference.digest and other.words_read == 3 * 6


def test_manifests_are_built_from_the_run_record(branching, ns):
    stream = _stream(ns, "manifest")
    run = simulate_aggregate(branching, {"Xa-400": 5000}, TIMES, stream)
    at_end = stream.provenance()
    stream.raw(3)                                # read after the run: not the run's
    rng = aggregate_manifest(experiment="e", run_id="r", run=run, code={})["rng"]
    assert rng["source_consumed_bytes"] == at_end["source_consumed_bytes"] == 8 * run.words_read
    assert rng["source_stream_sha256"] == at_end["source_stream_sha256"]
    assert rng["words_read_by_run"] == run.words_read == 6 * run.plan.words_per_interval
    assert rng["replay_identity"]["generator"] == REPLAYABLE[0]
    with pytest.raises(TypeError):               # no live stream can be attributed to the run
        aggregate_manifest(experiment="e", run_id="r", run=run, stream=stream, code={})


# --- 4. A failed draw leaves the stream at an interval boundary in both engines ----------------------------

def test_both_engines_consume_whole_intervals_when_a_draw_fails(ns, monkeypatch):
    scheme = precursor_scheme([PrecursorGroup("G", 1.0, 0.5)])
    plan = AggregatePlan(scheme, [0.0, 1.0], 10000)
    monkeypatch.setattr(binomial_module.special, "betaincc",
                        lambda a, b, x: np.full(np.broadcast(a, b, x).shape, np.nan))
    for engine in (simulate_aggregate, simulate_aggregate_reference):
        stream = _stream(ns, "failure", engine.__name__)
        with pytest.raises(FloatingPointError):
            engine(scheme, {"G": 10000}, [0.0, 1.0], stream)
        assert stream.words_read == plan.words_per_interval == 2, engine.__name__


# --- 5. Comparison tests at their boundaries ---------------------------------------------------------------

def test_allocation_sums_are_exact_and_inputs_are_integers():
    """The review's probe: cell sums of [2**62, 2**62] overflowed int64 and the cell was skipped (p = 1)."""
    result = allocation_check([0, 0, 0], [2 ** 62, 2 ** 62, 0], 0.5, [[0, 1], [2]])
    assert result["cell_trials"] == [2 ** 63, 0] and result["p_value"] < 1e-300
    for successes, trials in (([0.9, 0.9], [1.9, 1.9]), ([True, 0], [1, 1]), ([-1, 0], [1, 1]), ([2, 0], [1, 1])):
        with pytest.raises(ValueError):
            allocation_check(successes, trials, 0.5, [[0], [1]])


def test_homogeneity_and_fixed_cells_refuse_non_integer_counts_and_do_not_overflow():
    with pytest.raises(ValueError):
        homogeneity_check([1.5, 2], [1, 2])
    with pytest.raises(ValueError):
        fixed_cells_chisquare([0.9, 1], [0.5, 0.5], [[0], [1]])
    huge = fixed_cells_chisquare([2 ** 62, 2 ** 62], [0.5, 0.5], [[0], [1]])
    assert huge["total"] == 2 ** 63 and huge["statistic"] == 0.0


def test_degenerate_degrees_of_freedom_raise_instead_of_reading_as_p_one():
    """The review's probe: nu underflowed to 0, SciPy returned NaN and min(1, NaN) reported p = 1."""
    with pytest.raises(ValueError, match="degrees of freedom"):
        scaled_chisquare_variance_check(1.0, 1e-200, 1e-200)
    with pytest.raises(ValueError, match="degrees of freedom"):
        variance_ratio_f_check(1.0, 2.0, 1e-200, 1e-200)
    with pytest.raises(ValueError, match="degrees of freedom"):
        scaled_chisquare_variance_check(1.0, 1e200, 1e-200)


def test_a_nan_tail_is_never_read_as_p_one(monkeypatch):
    from rnglab.nuclear import comparison
    monkeypatch.setattr(comparison.stats.chi2, "cdf", lambda x, nu: float("nan"))
    with pytest.raises(ValueError):
        scaled_chisquare_variance_check(1.0, 1.0, 0.02)
    monkeypatch.setattr(comparison.stats.f, "sf", lambda x, a, b: float("nan"))
    with pytest.raises(ValueError):
        variance_ratio_f_check(1.0, 1.1, 1.0, 0.02)


# --- 6. Once-only record runs: attempts, atomic publication ------------------------------------------------

def test_attempts_are_reserved_once_and_record_how_they_ended(tmp_path):
    rv = _verification()
    code = {"git_commit": "c1", "git_dirty": False}
    first = rv.begin_attempt(tmp_path, "x", code)
    assert first["path"].name == "x-1.json" and first["dir"].is_dir()
    assert json.loads(first["path"].read_text())["status"] == "started"
    with pytest.raises(SystemExit, match="no end record"):
        rv.begin_attempt(tmp_path, "x", code)                     # still marked as running
    second = rv.begin_attempt(tmp_path, "x", code, previous_ended=True)
    assert second["record"]["previous_attempts"] == [{"attempt": 1, "status": "started", "git_commit": "c1"}]

    def broken():
        raise KeyError("boom")
    with pytest.raises(KeyError):
        rv.run_attempt(second, broken, tmp_path / "result.json")
    assert json.loads(second["path"].read_text())["status"] == "failed"
    with pytest.raises(SystemExit, match="deviation"):
        rv.begin_attempt(tmp_path, "x", {"git_commit": "c2"}, previous_ended=True)
    third = rv.begin_attempt(tmp_path, "x", code, previous_ended=True)

    def works():
        rv.write_new(tmp_path / "result.json", "{}\n")
        return 0
    assert rv.run_attempt(third, works, tmp_path / "result.json") == 0
    record = json.loads(third["path"].read_text())
    assert record["status"] == "completed" and record["result_sha256"] == rv.file_sha256(tmp_path / "result.json")
    with pytest.raises(SystemExit, match="completed attempt"):
        rv.begin_attempt(tmp_path, "x", code, previous_ended=True, deviation="any")


def test_results_are_published_atomically_and_never_overwritten(tmp_path):
    rv = _verification()
    target = tmp_path / "verification.json"
    rv.write_new(target, "first\n")
    with pytest.raises(FileExistsError):
        rv.write_new(target, "second\n")
    assert target.read_text() == "first\n"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["verification.json"]   # no temporary file left


def test_sigterm_ends_the_attempt_as_terminated(tmp_path):
    rv = _verification()
    attempt = rv.begin_attempt(tmp_path, "x", {"git_commit": "c"})
    before = signal.getsignal(signal.SIGTERM)

    def killed():
        os.kill(os.getpid(), signal.SIGTERM)
        time.sleep(5)                             # the handler raises before this returns
        return 0
    with pytest.raises(rv.Terminated):
        rv.run_attempt(attempt, killed, tmp_path / "result.json")
    assert json.loads(attempt["path"].read_text())["status"] == "terminated"
    assert signal.getsignal(signal.SIGTERM) is before


def test_the_replay_check_compares_chunks_of_earlier_attempts(tmp_path):
    cal = _calibration()
    folder = tmp_path / "attempts"
    for k in (1, 2):
        (folder / f"calibration-{k}").mkdir(parents=True)
    line = lambda first, digest: json.dumps({"first_battery": first, "batteries": 10, "sha256": digest}) + "\n"  # noqa: E731
    (folder / "calibration-1" / "chunks.jsonl").write_text(line(0, "a") + line(10, "b") + '{"first_battery": 20, "ba')
    (folder / "calibration-2" / "chunks.jsonl").write_text(line(0, "a") + line(10, "c") + line(20, "d"))
    attempt = {"dir": folder / "calibration-2",
               "record": {"code": {"git_commit": "c"},
                          "previous_attempts": [{"attempt": 1, "status": "terminated", "git_commit": "c"}]}}
    assert cal.replay_check(attempt) == {"chunks_compared": 2,
                                         "mismatches": [{"attempt": 1, "first_battery": 10, "batteries": 10}]}
    assert cal.chunk_digest({"b": 1, "a": [1, 2]}) == cal.chunk_digest({"a": [1, 2], "b": 1})


def test_calibration_sets_one_thread_per_worker_whatever_was_inherited(monkeypatch):
    """The review's probe: setdefault kept an inherited OPENBLAS_NUM_THREADS=4."""
    cal = _calibration()
    for name in cal.THREAD_VARIABLES:
        monkeypatch.setenv(name, "4")
    inherited, applied = cal.single_thread_environment()
    assert inherited == {name: "4" for name in cal.THREAD_VARIABLES}
    assert applied == {name: "1" for name in cal.THREAD_VARIABLES} == {n: os.environ[n] for n in cal.THREAD_VARIABLES}
    assert 1 <= cal.allocated_cpus() <= (os.cpu_count() or 1)


# --- 7. Reserved namespaces, at any depth and at run time --------------------------------------------------

def test_every_namespace_a_protocol_names_is_reserved_and_absent_from_the_tests():
    protocol = json.loads((REPO / "research/v020-nuclear-b2-verify/protocol.json").read_text())
    reserved = reserved_namespaces()
    nested = {protocol["calibration"]["namespace"], protocol["calibration"]["development_namespace"]}
    assert nested <= reserved and {protocol["namespace"], protocol["smoke_namespace"]} <= reserved
    assert all(not B2_TEST_NAMESPACE.startswith(r) and not r.startswith(B2_TEST_NAMESPACE) for r in reserved)
    frozen_guard = REPO / "tests/v020/nuclear/test_layer_a_namespaces.py"   # names Layer A's namespace to check it
    for path in (REPO / "tests").rglob("*.py"):
        if path == frozen_guard:
            continue
        text = path.read_text()
        assert not any(name in text for name in reserved), path


def test_namespace_seed_refuses_reserved_namespaces_while_a_test_runs():
    protocol = json.loads((REPO / "research/v020-nuclear-b2-verify/protocol.json").read_text())
    record = protocol["calibration"]["namespace"]           # reached indirectly, as a script would
    import rnglab.seeds
    cal = _calibration()
    for seed_function in (rnglab.seeds.namespace_seed, namespace_seed, cal.namespace_seed):
        with pytest.raises(ReservedNamespaceError):
            seed_function(record, "battery", 0)
        with pytest.raises(ReservedNamespaceError):
            seed_function(record + "/battery", 0)                # the same label, spelled differently
    with pytest.raises(ReservedNamespaceError):
        labelled_source(protocol["namespace"], "anything")
    assert isinstance(namespace_seed(B2_TEST_NAMESPACE, "x"), int)   # other namespaces are unaffected


# --- 8. Costs: checkpoint size, category counts --------------------------------------------------------------

def test_checkpoint_size_does_not_grow_with_the_grid(precursors, ns):
    """The review's probe: every checkpoint serialised the whole grid (O(I**2 / k) in all)."""
    sizes = []
    for intervals in (10, 2000):
        seen = []
        simulate_aggregate(precursors, {"G1": 1000}, [float(k) for k in range(intervals + 1)],
                           _stream(ns, "size", intervals), checkpoint_every=intervals, on_checkpoint=seen.append,
                           keep_outputs=False)
        sizes.append(len(json.dumps(seen[-1])))
    assert abs(sizes[1] - sizes[0]) < 40


def _ladder(depth: int) -> DecayScheme:
    """A chain whose every nuclide decays to the next by two channels: 2**(depth + 1) - 1 categories."""
    nuclides = []
    for d in range(depth + 1):
        channels = () if d == depth else (Channel(f"L{d + 1}", Mode.BETA_PLUS, 0.5),
                                          Channel(f"L{d + 1}", Mode.ELECTRON_CAPTURE, 0.5))
        nuclides.append(Nuclide(f"L{d}", 150 - d, 300, math.inf if d == depth else 1.0 + 0.37 * d, channels))
    return DecayScheme(f"ladder-{depth}", nuclides)


def test_categories_are_counted_exactly_and_bounded_before_enumeration():
    small = _ladder(4)
    assert category_counts(small) == [len(path_categories(small, i)) for i in range(len(small))]
    assert category_counts(small)[0] == 31
    start = time.perf_counter()
    with pytest.raises(ValueError, match="max_categories"):
        AggregatePlan(_ladder(40), [0.0, 1.0], 100)            # 2**41 - 1 categories: refused, not enumerated
    assert time.perf_counter() - start < 1.0
    with pytest.raises(ValueError, match="max_categories"):
        AggregatePlan(_ladder(10), [0.0, 1.0], 100)            # 2047 > 1024
    assert AggregatePlan(small, [0.0, 1.0], 100, max_categories=31).n_categories == sum(category_counts(small))


# --- 9. Benchmark summaries --------------------------------------------------------------------------------

def test_benchmark_timings_are_medians_with_their_samples(monkeypatch):
    """The review's probe: plan_seconds was copied from the first repeat instead of the median."""
    bench = _load("b2_integrity_benchmark", "research/v020-nuclear-b2-bench/benchmark.py")
    base = {"words": 6, "parts": 1, "traced_peak_bytes": None}
    outputs = iter([dict(base, wall_seconds=3.0, cpu_seconds=2.5, plan_seconds=0.3, max_rss_bytes=30),
                    dict(base, wall_seconds=1.0, cpu_seconds=0.5, plan_seconds=0.1, max_rss_bytes=10),
                    dict(base, wall_seconds=2.0, cpu_seconds=1.5, plan_seconds=0.2, max_rss_bytes=20),
                    dict(base, wall_seconds=9.0, cpu_seconds=9.0, plan_seconds=0.9, max_rss_bytes=90,
                         traced_peak_bytes=99)])
    monkeypatch.setattr(bench.subprocess, "run", lambda *args, **kw: subprocess.CompletedProcess(
        args, 0, stdout=json.dumps(next(outputs)) + "\n", stderr=""))
    summary = bench.measure({"engine": "b2", "n": 10}, repeats=3, threads="1")
    assert (summary["plan_seconds"], summary["wall_seconds"], summary["cpu_seconds"]) == (0.2, 2.0, 1.5)
    assert summary["samples"]["plan_seconds"] == [0.3, 0.1, 0.2] and summary["traced_peak_bytes"] == 99
    outputs = iter([dict(base, wall_seconds=1.0, cpu_seconds=1.0), dict(base, words=7, wall_seconds=1.0,
                                                                       cpu_seconds=1.0)])
    with pytest.raises(RuntimeError, match="disagree"):
        bench.measure({"engine": "b2", "n": 10}, repeats=2, threads="1")


class _InlineExecutor:
    """Stands in for the process pool: the same initializer and map, in this process."""

    def __init__(self, max_workers, initializer):
        initializer()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False

    def map(self, function, *iterables):
        return map(function, *iterables)


def test_the_record_calibration_path_logs_chunks_replays_and_publishes_once(tmp_path, monkeypatch):
    """calibrate() as the record run uses it (an attempt), in the B2 test namespace, without a process pool."""
    from types import SimpleNamespace
    cal, rv = _calibration(), _verification()
    monkeypatch.setattr(cal, "ProcessPoolExecutor", _InlineExecutor)
    for name in cal.THREAD_VARIABLES:
        monkeypatch.setenv(name, "1")
    namespace = f"{B2_TEST_NAMESPACE}-calibration-path"
    args = SimpleNamespace(chunk=2, development=False)
    code = {"git_commit": "c", "git_dirty": False}
    results = tmp_path / "results"
    first = rv.begin_attempt(results, "calibration", code)
    assert cal.calibrate(rv, namespace, 4, results / "a.json", {}, args, 1, first) == 0
    a = json.loads((results / "a.json").read_text())
    assert a["batteries"] == 4 and a["attempt"]["number"] == 1
    assert a["replay_check"] == {"chunks_compared": 0, "mismatches": []}
    assert len((first["dir"] / "chunks.jsonl").read_text().splitlines()) == 2
    second = rv.begin_attempt(results, "calibration", code, previous_ended=True)   # the first never ended
    assert cal.calibrate(rv, namespace, 4, results / "b.json", {}, args, 1, second) == 0
    b = json.loads((results / "b.json").read_text())
    assert b["replay_check"] == {"chunks_compared": 2, "mismatches": []}            # the same draws again
    assert a["families"] == b["families"]
    with pytest.raises(FileExistsError):
        cal.calibrate(rv, namespace, 4, results / "a.json", {}, args, 1,
                      rv.begin_attempt(results, "calibration", code, previous_ended=True))
