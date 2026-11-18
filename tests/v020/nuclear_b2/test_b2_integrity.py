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
        "_source_plans[...] = another": lambda: plan._source_plans.__setitem__((delta, 0), None),
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
    assert again.digest == reference.digest                        # the next run with the same plan too
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
    with pytest.raises(AttributeError):
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
    with pytest.raises(SystemExit, match="another x attempt is running"):
        rv.begin_attempt(tmp_path, "x", code)                     # the first holds the lock
    rv.release_attempt(first)                    # as if its process died: the lock is gone, no end record
    second = rv.begin_attempt(tmp_path, "x", code)
    assert second["record"]["previous_attempts"] == [{"attempt": 1, "status": "started", "git_commit": "c1"}]
    assert second["record"]["interrupted_attempts"] == [1]

    def broken():
        raise KeyError("boom")
    with pytest.raises(KeyError):
        rv.run_attempt(second, broken, tmp_path / "result.json")
    assert json.loads(second["path"].read_text())["status"] == "failed"
    with pytest.raises(SystemExit, match="deviation"):
        rv.begin_attempt(tmp_path, "x", {"git_commit": "c2"})
    third = rv.begin_attempt(tmp_path, "x", code)

    def works():
        rv.write_new(tmp_path / "result.json", "{}\n")
        return 0
    assert rv.run_attempt(third, works, tmp_path / "result.json") == 0
    record = json.loads(third["path"].read_text())
    assert record["status"] == "completed" and record["result_sha256"] == rv.file_sha256(tmp_path / "result.json")
    with pytest.raises(SystemExit, match="completed attempt"):
        rv.begin_attempt(tmp_path, "x", code, deviation="any")


def test_overlapping_launches_are_refused_and_records_are_complete_or_absent(tmp_path):
    """The review's probes: a second recovery launch treated the running attempt as dead (2 and 3 both
    'started'); a reservation cut off mid-write ({"schema":) raised JSONDecodeError."""
    rv = _verification()
    code = {"git_commit": "c1"}
    running = rv.begin_attempt(tmp_path, "y", code)
    for _ in range(2):                           # however often it is tried while the first runs
        with pytest.raises(SystemExit, match="another y attempt is running"):
            rv.begin_attempt(tmp_path, "y", code)
    assert sorted(p.name for p in (tmp_path / "attempts").glob("y-*.json")) == ["y-1.json"]
    rv.release_attempt(running)
    leftover = tmp_path / "attempts" / ".y-2.json.999.tmp"
    leftover.write_text('{"schema":')            # a crash between the temporary file and its link
    nxt = rv.begin_attempt(tmp_path, "y", code)
    assert nxt["record"]["attempt"] == 2 and nxt["record"]["leftover_temporary_files"] == [leftover.name]
    rv.release_attempt(nxt)
    (tmp_path / "attempts" / "y-2.json").write_text('{"schema":')   # edited by hand
    with pytest.raises(SystemExit, match="cannot be read"):
        rv.begin_attempt(tmp_path, "y", code)


def _git(repo, *args):
    subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@example.invalid", *args],
                   check=True, capture_output=True)


def test_a_record_run_may_start_with_files_of_earlier_attempts_and_nothing_else(tmp_path):
    """The review's probe: an interrupted attempt's untracked files made the tree dirty, so the documented
    rerun on the same commit was refused by the clean-tree guard."""
    rv = _verification()
    _git(tmp_path, "init", "-q")
    (tmp_path / "tracked.txt").write_text("a\n")
    _git(tmp_path, "add", "tracked.txt")
    _git(tmp_path, "commit", "-q", "-m", "c")
    attempts = tmp_path / rv.ATTEMPTS
    (attempts / "calibration-1").mkdir(parents=True)
    (attempts / "calibration-1.json").write_text("{}\n")
    (attempts / "calibration-1" / "chunks.jsonl").write_text("")
    tree = rv.require_record_tree(tmp_path, "calibration")
    assert len(tree["attempt_files"]) == 2 and tree["other_changes"] == []
    (tmp_path / "tracked.txt").write_text("b\n")
    with pytest.raises(SystemExit, match="tracked.txt"):
        rv.require_record_tree(tmp_path, "calibration")
    (tmp_path / "tracked.txt").write_text("a\n")
    (attempts.parent / "calibration.json").write_text("{}\n")   # results outside attempts/ must be committed
    with pytest.raises(SystemExit, match="calibration.json"):
        rv.require_record_tree(tmp_path, "verification")


def test_r0_never_reads_an_upper_bound_above_the_limit_as_an_excess():
    """The review's probe: 950 of 50,000 (1.9%, interval 1.78% to 2.02%) was to be reported as an excess.
    Fifth review: R0 trusted the stored interval (950 failures with a stored [0.001, 0.005] read as C1 holding),
    and R1's own false-FAIL rate (a FAIL in either family) was never reported."""
    rv = _verification()

    def family(k):
        return {"failed_batteries": k, "false_fail_ci95": [0.001, 0.005]}     # a stored interval is not read
    reading = rv.r0_reading({"inconclusive": family(950), "holds": family(500), "excess": family(1300)}, 50000,
                            {"failed_batteries": 1800})
    assert reading["inconclusive"]["false_fail_ci95"] == rv.clopper_pearson(950, 50000)
    assert reading["R1_event"]["false_fail_ci95"][0] > 0.03 and "NC1 or NC2" in reading["R1_event"]["statement"]
    assert not reading["inconclusive"]["C1_holds"] and not reading["inconclusive"]["excess_established"]
    assert "did not establish" in reading["inconclusive"]["statement"]
    assert reading["holds"]["C1_holds"] and "C1 holds" in reading["holds"]["statement"]
    assert reading["excess"]["excess_established"] and "exceeds" in reading["excess"]["statement"]


def _calibration_record(rv, repo, results, protocol, frozen, **change):
    """A record calibration and its completed attempt record, as the record procedure leaves them."""
    spec, batteries = protocol["calibration"], protocol["calibration"]["batteries"]
    code = rv.record_code(repo)
    planned = json.loads((REPO / rv.STUDY / "design.json").read_text())["planned_tests"]

    def family(name, k):
        interval = rv.clopper_pearson(k, batteries)
        return {"tests": len(planned[name]), "failed_batteries": k, "false_fail_rate": k / batteries,
                "false_fail_ci95": interval, "C1_holds": interval[1] <= rv.C1_LIMIT}
    data = {"development": False, "namespace": spec["namespace"], "batteries": batteries, "protocol_files": frozen,
            "code": code, "environment": rv.environment(),
            "attempt": {"number": 1, "record": "attempts/calibration-1.json", "deviation": None},
            "families": {"NC1": family("NC1", 500), "NC2": family("NC2", 400)}, "C1": True,
            "R1_event": {"failed_batteries": 900, "batteries": batteries, "false_fail_rate": 900 / batteries,
                         "false_fail_ci95": rv.clopper_pearson(900, batteries)},
            "replay_check": {"chunks_compared": 0, "mismatches": []}}
    data.update(change)
    (results / "attempts").mkdir(parents=True, exist_ok=True)
    (results / "calibration.json").write_text(json.dumps(data))
    record = {"status": "completed", "result_sha256": rv.file_sha256(results / "calibration.json")}
    (results / "attempts" / "calibration-1.json").write_text(json.dumps(record))
    return data


def test_the_verification_needs_a_matching_completed_record_calibration(tmp_path):
    """The review's probe: the verification reached the study without results/calibration.json. Fifth review:
    an uncommitted or ignored calibration, an uncommitted attempt record, other rnglab14 or seed code, a stored
    interval that its counts contradict and a failed replay check were all accepted."""
    rv = _verification()
    protocol = json.loads((REPO / "research/v020-nuclear-b2-verify/protocol.json").read_text())
    frozen = {"frozen": True, "files": {"protocol.json": "0" * 64}}
    repo, results = tmp_path / "repo", tmp_path / "repo" / rv.RESULTS
    (repo / "src").mkdir(parents=True)
    (repo / "src" / "seeds.py").write_text("x = 1\n")
    (repo / "README").write_text("r\n")
    _git(repo, "init", "-q")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "frozen")
    code = rv.record_code(repo)
    with pytest.raises(SystemExit, match="missing"):
        rv.record_calibration(protocol, frozen, code, results, repo)
    data = _calibration_record(rv, repo, results, protocol, frozen)
    with pytest.raises(SystemExit, match="not committed"):
        rv.record_calibration(protocol, frozen, code, results, repo)
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "calibration")
    info = rv.record_calibration(protocol, frozen, code, results, repo)
    assert info["R0"]["NC1"]["C1_holds"] and info["attempt"] == 1 and "R1_event" in info["R0"]
    design = json.loads((REPO / rv.STUDY / "design.json").read_text())
    assert rv.record_calibration(protocol, frozen, code, results, repo, design)["attempt"] == 1
    short = {"planned_tests": dict(design["planned_tests"], NC1=design["planned_tests"]["NC1"][:3])}
    with pytest.raises(SystemExit, match="number of tests for NC1"):
        rv.record_calibration(protocol, frozen, code, results, repo, short)
    nc1, r1 = data["families"]["NC1"], data["R1_event"]
    big = {"NC1": dict(nc1, failed_batteries=30000), "NC2": dict(nc1, failed_batteries=30000)}
    for change, message in ((dict(data, development=True), "development"),
                             (dict(data, protocol_files={"frozen": True, "files": {}}), "frozen files"),
                             (dict(data, code=dict(data["code"], module_files_sha256={})), "module files"),
                             (dict(data, code=dict(data["code"], source_files_sha256="0")), "source files"),
                             (dict(data, code=dict(data["code"], git_commit="0" * 40)), "not an ancestor"),
                             (dict(data, environment=dict(data["environment"], numpy="0.0")), "numpy"),
                             (dict(data, replay_check={"mismatches": [{"attempt": 1}]}), "replay check"),
                             (dict(data, families={"NC1": {"failed_batteries": 950,
                                                            "false_fail_ci95": [0.001, 0.005]}}), "counts"),
                             (dict(data, R1_event=dict(data["R1_event"], failed_batteries=5)), "R1 event"),
                             (dict(data, families={"NC1": data["families"]["NC1"]}), "not NC1 and NC2"),
                             (dict(data, families=dict(data["families"], NC2=dict(data["families"]["NC2"],
                                                                                 C1_holds=False))), "rate or C1"),
                             (dict(data, C1=False), "its C1"),
                             (dict(data, families=dict(data["families"], NC1=dict(nc1, false_fail_rate=0.5))),
                              "rate or C1"),
                             (dict(data, families=dict(data["families"], NC1=dict(nc1, failed_batteries=-1))),
                              "missing or impossible"),
                             (dict(data, families=dict(data["families"], NC1=dict(nc1, failed_batteries=True))),
                              "missing or impossible"),
                             (dict(data, families=big, R1_event=dict(r1, failed_batteries=55000)), "R1 event"),
                             (dict(data, R1_event=dict(r1, failed_batteries=True)), "R1 event"),
                             (dict(data, R1_event=dict(r1, false_fail_rate=0.5)), "R1 event"),
                             (dict(data, R1_event=[1, 2]), "R1 event"),
                             (dict(data, families=["NC1", "NC2"]), "not a mapping"),
                             (dict(data, families={"NC1": 3, "NC2": None}), "not a mapping"),
                             (dict(data, batteries="50000"), "positive integer"),
                             (dict(data, code=[1]), "module files"),
                             (dict(data, attempt=dict(data["attempt"], record="../calibration-1.json")),
                              "named as")):
        (results / "calibration.json").write_text(json.dumps(change))   # uncommitted, and no longer the record's
        with pytest.raises(SystemExit, match=message):
            rv.record_calibration(protocol, frozen, code, results, repo)
    _git(repo, "checkout", "-q", "--", ".")
    _git(repo, "mv", "README", rv.RESULTS + "README")      # the review's probe: rename detection hid the deletion
    _git(repo, "commit", "-q", "-m", "moved")
    with pytest.raises(SystemExit, match="outside results.*README"):
        rv.record_calibration(protocol, frozen, rv.record_code(repo), results, repo)
    _git(repo, "mv", rv.RESULTS + "README", "README")
    _git(repo, "commit", "-q", "-m", "moved back")
    (repo / "src" / "seeds.py").write_text("x = 2\n")                   # committed after the calibration
    _git(repo, "commit", "-q", "-am", "another seed derivation")
    with pytest.raises(SystemExit, match="outside results"):
        rv.record_calibration(protocol, frozen, rv.record_code(repo), results, repo)


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
    rv.release_attempt(first)                    # as if its process died: no end record
    second = rv.begin_attempt(results, "calibration", code)
    assert cal.calibrate(rv, namespace, 4, results / "b.json", {}, args, 1, second) == 0
    b = json.loads((results / "b.json").read_text())
    assert b["replay_check"] == {"chunks_compared": 2, "mismatches": []}            # the same draws again
    assert a["families"] == b["families"]
    rv.release_attempt(second)
    with pytest.raises(FileExistsError):
        cal.calibrate(rv, namespace, 4, results / "a.json", {}, args, 1,
                      rv.begin_attempt(results, "calibration", code))


# --- 10. Fourth review: corrections to the corrections -------------------------------------------------------

def _simple() -> DecayScheme:
    return DecayScheme("simple", [Nuclide("P", 1, 2, 1.0, (Channel("D", Mode.BETA_MINUS, 1.0),)),
                                  Nuclide("D", 2, 2, math.inf)])


def test_reassigning_dtype_or_shape_of_a_returned_array_changes_nothing(ns):
    """The review's probe: p_small.dtype = '>f8' in a sink made the second interval decay nothing,
    and final_populations.dtype = '>i8' changed the reported populations."""
    scheme, t = _simple(), [0.0, 1.0, 2.0]
    reference = simulate_aggregate(scheme, {"P": 10000}, t, _stream(ns, "dtype"))
    plan = AggregatePlan(scheme, t, 10000)

    def vandal(interval, record):
        source = plan.interval_plans[1.0][0]
        for view in (source.p_small, source.flip, plan.category_steps, plan.first_category, plan.deltas):
            if view.ndim > 1:
                with pytest.warns(DeprecationWarning):
                    view.shape = (-1,)
            with pytest.warns(DeprecationWarning):
                view.dtype = ">f8" if view.dtype.kind == "f" else ">i8" if view.dtype.kind == "i" else "u1"

    run = simulate_aggregate(scheme, {"P": 10000}, t, _stream(ns, "dtype"), plan=plan, sink=vandal)
    assert run.digest == reference.digest and run.populations.tolist() == reference.populations.tolist()
    final = run.final_populations.tolist()
    view = run.final_populations
    with pytest.warns(DeprecationWarning):
        view.dtype = ">i8"
    assert run.final_populations.tolist() == final and run.summary()["final_populations"] == dict(zip("PD", final))
    assert run.final_populations is not run.final_populations      # every read is a fresh view


def test_a_stream_given_another_source_during_the_run_is_refused(ns):
    """The review's probe: w.source = other in a sink; the run used two sources and named only the second."""
    scheme, t = _simple(), [0.0, 1.0, 2.0]
    stream = _stream(ns, "swap-a")
    other = labelled_source(ns, "integrity", "swap-b")
    words = AggregatePlan(scheme, t, 10000).words_per_interval
    other.read_bytes(8 * words)        # as in the review: the other source's counters match the layout's

    def swap(interval, record):
        stream.source = other
    with pytest.raises(RuntimeError, match="another source"):
        simulate_aggregate(scheme, {"P": 10000}, t, stream, sink=swap)
    assert other.consumed_bytes == 8 * words


def test_unaligned_sources_cannot_publish_checkpoints(ns):
    """The review's probe: a source one byte into a word published checkpoints no resume could use."""
    source = labelled_source(ns, "integrity", "unaligned")
    source.read_bytes(1)
    with pytest.raises(ValueError, match="word boundary"):
        simulate_aggregate(_simple(), {"P": 100}, [0.0, 1.0, 2.0], WordStream(source), checkpoint_every=1,
                           on_checkpoint=lambda checkpoint: None)
    assert source.consumed_bytes == 1


def test_the_allocation_statistic_is_exact_at_probabilities_near_one():
    """The review's probe: at p = 1 - 2**-53 and T = 2**62 - 1 the float residual was 0 (p = 1);
    the exact statistic is 30.25 (p = 2.7e-7)."""
    from fractions import Fraction
    n, p = 2 ** 62 - 1, 1 - 2 ** -53
    result = allocation_check([n - 600, n - 600], [n, n], p, [[0], [1]])
    exact = 2 * (Fraction(n - 600) - Fraction(p) * n) ** 2 / (Fraction(p) * (1 - Fraction(p)) * n)
    assert result["statistic"] == float(exact) and 2.6e-7 < result["p_value"] < 2.8e-7
    small = allocation_check([3, 5, 9], [10, 10, 20], 0.37, [[0, 1], [2]])
    terms = [(Fraction(8) - Fraction(0.37) * 20) ** 2 / (Fraction(0.37) * (1 - Fraction(0.37)) * 20),
             (Fraction(9) - Fraction(0.37) * 20) ** 2 / (Fraction(0.37) * (1 - Fraction(0.37)) * 20)]
    assert small["statistic"] == float(sum(terms))


def test_sampler_tails_refuse_a_non_finite_value(monkeypatch):
    monkeypatch.setattr(binomial_module.special, "betaincc", lambda a, b, x: float("nan"))
    with pytest.raises(FloatingPointError):
        binomial_module.sampler_tails(5000, 10000, 0.5)


def test_chernoff_tails_do_not_underflow():
    """The review's probe: float64 gave 0 for the bound 3.07e-783 at n = 10**6, p = 1/2, k = 469,994."""
    from decimal import Decimal, localcontext
    numerics = _load("b2_integrity_numerics", "research/v020-nuclear-b2-numerics/numerics.py")
    value = numerics._chernoff_tail(10 ** 6, 0.5, 469994, True)
    with localcontext() as context:
        context.prec = 80
        x = Decimal(469994) / 10 ** 6
        kl = x * (x / Decimal("0.5")).ln() + (1 - x) * ((1 - x) / Decimal("0.5")).ln()
        reference = (-10 ** 6 * kl).exp()
        assert value > 0 and abs(value / reference - 1) < Decimal("1e-50")


def test_benchmark_keeps_stochastic_fields_of_real_layer_a_repeats(ns, monkeypatch):
    """The review's probe: Layer A event counts differ between repeats (own seeds), and the consistency
    check stopped the benchmark on its first Layer A configuration."""
    bench = _load("b2_integrity_benchmark", "research/v020-nuclear-b2-bench/benchmark.py")
    monkeypatch.setattr(bench, "NAMESPACE", f"{ns}-benchmark")
    samples = [bench.worker({"engine": "layer_a", "n": 1000, "repeat": k}) for k in range(3)]
    assert len({s["events"] for s in samples}) > 1 and len({s["words"] for s in samples}) == 1
    traced = dict(samples[0], traced_peak_bytes=1)
    outputs = iter(samples + [traced])
    monkeypatch.setattr(bench.subprocess, "run", lambda *args, **kw: subprocess.CompletedProcess(
        args, 0, stdout=json.dumps(next(outputs)) + "\n", stderr=""))
    summary = bench.measure({"engine": "layer_a", "n": 1000}, repeats=3, threads="1")
    assert summary["samples"]["events"] == [s["events"] for s in samples]
    assert summary["events"] == sorted(s["events"] for s in samples)[1] and summary["words"] == samples[0]["words"]


def test_the_design_bounds_the_calibration_null_lattice():
    """numpy's binomial at counts above 2**53 returns integers on a lattice; design.json bounds its spacing
    against the standard deviation of the BIG state counts."""
    design = json.loads((REPO / "research/v020-nuclear-b2-verify/design.json").read_text())
    lattice = design["calibration_null_lattice"]
    assert lattice["cells_above_2_53"] > 0 and lattice["worst"]["spacing_over_sd"] < 1e-6


# --- 11. Fifth review ------------------------------------------------------------------------------------------

class _Grab:
    """Receives the dict a mapping proxy hands to the other operand of == or |."""

    __hash__ = None

    def __eq__(self, other):
        self.got = other
        return NotImplemented

    def __ror__(self, other):
        self.got = other
        return NotImplemented


def _leak(mapping):
    grabs = []
    for take in (lambda g: mapping == g, lambda g: mapping | g):
        grab = _Grab()
        try:
            take(grab)
        except TypeError:
            pass
        grabs.append(grab.got)
    return grabs


def test_a_source_given_another_manifest_or_first_byte_during_the_run_is_refused(ns):
    """The review's probe: st.source.manifest = other.manifest in a sink; the run completed and its provenance
    named the other source (seed and start word included)."""
    scheme, t = _simple(), [0.0, 1.0, 2.0]
    other = labelled_source(ns, "integrity", "someone-else")
    vandals = {"manifest": lambda stream: setattr(stream.source, "manifest", other.manifest),
               "first byte": lambda stream: setattr(stream, "first_source_byte", 4242)}
    for name, vandal in vandals.items():
        stream = _stream(ns, "manifest", name)
        with pytest.raises(RuntimeError, match="another manifest"):
            simulate_aggregate(scheme, {"P": 10000}, t, stream, sink=lambda i, r, s=stream, v=vandal: v(s))


def test_a_manifest_edited_in_place_fails_the_run_and_checkpoints_keep_the_recorded_identity(ns):
    """The review's probe: metadata['start_word'] = 999 in a sink put a false start word into the run's
    provenance; checkpoints read the live manifest too."""
    scheme, t = _simple(), [0.0, 1.0, 2.0, 3.0]
    plan = AggregatePlan(scheme, t, 10000)
    stream = _stream(ns, "manifest", "in-place")
    seed, published = stream.source.manifest.seed, []

    def vandal(interval, record):
        if interval == 0:
            stream.source.manifest.metadata["start_word"] = 999
            stream.source.manifest.metadata["seed_label"] = "forged"
    with pytest.raises(RuntimeError, match="edited"):
        simulate_aggregate(scheme, {"P": 10000}, t, stream, plan=plan, sink=vandal, checkpoint_every=1,
                           on_checkpoint=published.append)
    assert [c["stream"] for c in published] == [
        {"generator": REPLAYABLE[0], "parameters": {"seed": seed}, "word_position": (k + 1) * plan.words_per_interval,
         "byte_aligned": True} for k in range(3)]


def test_a_mapping_leaked_through_comparison_or_union_is_a_throwaway_copy(ns):
    """The review's probe: a mapping proxy hands its dict to the other operand of == and |; editing it in a
    sink changed the slots the engine read, and editing run.initial's changed the manifest's populations."""
    scheme, t = _simple(), [0.0, 1.0, 2.0]
    reference = simulate_aggregate(scheme, {"P": 10000}, t, _stream(ns, "leak"))
    plan = AggregatePlan(scheme, t, 10000)

    def vandal(interval, record):
        for mapping in (plan.slot_offset, plan.diagnostics, plan.interval_plans, plan.interval_plans[1.0]):
            for leaked in _leak(mapping):
                key = next(iter(leaked))
                leaked[key] = 5
    run = simulate_aggregate(scheme, {"P": 10000}, t, _stream(ns, "leak"), plan=plan, sink=vandal)
    assert run.digest == reference.digest and dict(plan.slot_offset) == {0: 0}
    assert dict(plan.diagnostics) == dict(reference.plan.diagnostics)
    for leaked in _leak(run.initial):
        leaked["P"] = 5
    assert dict(run.initial) == {"P": 10000, "D": 0}
    manifest = aggregate_manifest(experiment="integrity", run_id="leak", run=run, code={"commit": "test"})
    assert json.dumps(manifest).count('"P": 5') == 0


def test_source_plans_of_every_interval_length_are_taken_when_the_run_begins(ns):
    """The review's probe on the grid [0, 1, 3]: replacing the second length's source plan in a sink after
    interval 0 changed the second interval's draws, because views were taken lazily, per length."""
    scheme, t = _simple(), [0.0, 1.0, 3.0]
    reference = simulate_aggregate(scheme, {"P": 10000}, t, _stream(ns, "lengths"))
    plan = AggregatePlan(scheme, t, 10000)

    def vandal(interval, record):
        if interval == 0:          # a deliberate bypass of the plan's protections, after the run has begun
            later = plan.interval_plans[2.0][0]
            later.__dict__["_p_small_array"] = (np.array([0.5]).tobytes(), "<f8", (1,))
    run = simulate_aggregate(scheme, {"P": 10000}, t, _stream(ns, "lengths"), plan=plan, sink=vandal)
    assert run.digest == reference.digest


_CHILD = """
import importlib.util, sys, time
from pathlib import Path
spec = importlib.util.spec_from_file_location("child_run_verification", sys.argv[1])
rv = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rv)
attempt = rv.begin_attempt(Path(sys.argv[2]), "calibration", {"git_commit": "c"})   # holds the lock
print("ready", flush=True)
time.sleep(120)
"""


def test_a_deleted_or_replaced_lock_file_does_not_let_a_second_attempt_start(tmp_path):
    """The review's probe: with attempt 1 alive, rm (or cp + mv) of the lock file let attempt 2 start and
    list attempt 1 as interrupted."""
    rv = _verification()
    out, script = tmp_path / "results", tmp_path / "calibration.py"
    script.write_text(_CHILD)
    child = subprocess.Popen([sys.executable, str(script), str(REPO / rv.STUDY / "run_verification.py"), str(out)],
                             stdout=subprocess.PIPE, text=True, env=dict(os.environ, PYTHONPATH=str(REPO / "src")))
    try:
        assert child.stdout.readline().strip() == "ready"
        lock = out / "attempts" / "calibration.lock"
        with pytest.raises(SystemExit, match="another calibration attempt is running"):
            rv.begin_attempt(out, "calibration", {"git_commit": "c"})
        lock.unlink()
        with pytest.raises(SystemExit, match="still runs on this host"):
            rv.begin_attempt(out, "calibration", {"git_commit": "c"})
        copy = lock.with_name("copy")
        copy.write_text(lock.read_text())
        os.replace(copy, lock)
        with pytest.raises(SystemExit, match="still runs on this host"):
            rv.begin_attempt(out, "calibration", {"git_commit": "c"})
        assert sorted(p.name for p in (out / "attempts").glob("calibration-*.json")) == ["calibration-1.json"]
    finally:
        child.kill()
        child.wait()


def test_a_lock_taken_on_a_file_replaced_meanwhile_is_dropped_and_taken_again(tmp_path, monkeypatch):
    rv = _verification()
    real, answers = rv._same_file, iter([False])            # the first lock is on a file replaced meanwhile
    monkeypatch.setattr(rv, "_same_file", lambda handle, path: next(answers, None) is None and real(handle, path))
    attempt = rv.begin_attempt(tmp_path, "calibration", {"git_commit": "c"})
    assert attempt["record"]["attempt"] == 1
    rv.release_attempt(attempt)
    monkeypatch.setattr(rv, "_same_file", lambda handle, path: False)
    with pytest.raises(SystemExit, match="keeps being deleted or replaced"):
        rv.begin_attempt(tmp_path, "calibration", {"git_commit": "c"})


def test_a_result_published_before_its_record_was_finished_is_completed_on_the_next_launch(tmp_path):
    """The review's probe: a kill between publishing calibration.json and recording 'completed' left a
    published result that no later step could use."""
    rv = _verification()
    code = {"git_commit": "c"}
    attempt = rv.begin_attempt(tmp_path, "calibration", code)
    result = tmp_path / "calibration.json"
    rv.write_new(result, json.dumps({"attempt": {"number": 1}, "code": code}))
    rv.release_attempt(attempt)                                   # killed before end_attempt
    assert rv.complete_published_attempt(tmp_path, "calibration", result)
    record = json.loads(attempt["path"].read_text())
    assert record["status"] == "completed" and record["result_sha256"] == rv.file_sha256(result)
    assert record["completed_after_interruption"]["status_before"] == "started"
    assert not rv.complete_published_attempt(tmp_path, "calibration", result)       # nothing left to finish
    other = tmp_path / "other"
    attempt = rv.begin_attempt(other, "calibration", code)
    rv.end_attempt(attempt, "failed", error="OSError: [Errno 5] fsync after the link")  # a late publication error
    rv.release_attempt(attempt)
    rv.write_new(other / "calibration.json", json.dumps({"attempt": {"number": 1}, "code": code}))
    assert rv.complete_published_attempt(other, "calibration", other / "calibration.json")
    assert json.loads(attempt["path"].read_text())["completed_after_interruption"]["status_before"] == "failed"
    third = tmp_path / "third"
    attempt = rv.begin_attempt(third, "calibration", code)
    rv.release_attempt(attempt)
    rv.write_new(third / "calibration.json", json.dumps({"attempt": {"number": 1}, "code": {"git_commit": "d"}}))
    with pytest.raises(SystemExit, match="does not describe"):
        rv.complete_published_attempt(third, "calibration", third / "calibration.json")


def test_a_record_calibration_whose_code_changed_while_it_ran_publishes_nothing(tmp_path, monkeypatch):
    """The review's probe: comparison.py edited and committed 5 s into a record calibration; the result carried
    the end identity and was accepted."""
    from types import SimpleNamespace
    cal, rv = _calibration(), _verification()
    monkeypatch.setattr(cal, "ProcessPoolExecutor", _InlineExecutor)
    for name in cal.THREAD_VARIABLES:
        monkeypatch.setenv(name, "1")
    start = rv.record_code(REPO)
    monkeypatch.setattr(rv, "record_code", lambda repo: dict(start, source_files_sha256="changed"))
    attempt = rv.begin_attempt(tmp_path, "calibration", start)
    with pytest.raises(RuntimeError, match="changed while it ran"):
        cal.calibrate(rv, f"{B2_TEST_NAMESPACE}-code-change", 2, tmp_path / "c.json", {},
                      SimpleNamespace(chunk=2, development=False), 1, attempt, code=start)
    assert not (tmp_path / "c.json").exists()
    rv.release_attempt(attempt)
    with pytest.raises(RuntimeError, match="source_files_sha256"):
        rv.require_same_code(start, REPO, "verification")


def test_record_runs_refuse_ignored_modules_beside_the_code_and_bytecode_from_the_tree(tmp_path, monkeypatch):
    rv = _verification()
    _git(tmp_path, "init", "-q")
    (tmp_path / ".gitignore").write_text("__pycache__/\nsitecustomize.py\nsmoke/\n*.egg-info/\n")
    (tmp_path / "src" / "pkg.egg-info").mkdir(parents=True)                  # left by pip install -e .
    (tmp_path / "src" / "pkg.egg-info" / "PKG-INFO").write_text("Name: pkg\n")
    (tmp_path / "src" / "pkg" / "__pycache__").mkdir(parents=True)
    (tmp_path / "src" / "pkg" / "__pycache__" / "m.cpython-314.pyc").write_bytes(b"x")
    (tmp_path / rv.STUDY / "smoke").mkdir(parents=True)
    (tmp_path / rv.STUDY / "smoke" / "verification.json").write_text("{}")
    _git(tmp_path, "add", "-A")
    _git(tmp_path, "commit", "-q", "-m", "c")
    tree = rv.require_record_tree(tmp_path, "verification")
    assert tree["other_changes"] == [] and len(tree["ignored_allowed"]) == 3      # recorded in the result
    assert all(any(kind in entry for entry in tree["ignored_allowed"]) for kind in ("egg-info", "__pycache__", "smoke/"))
    (tmp_path / "src" / "sitecustomize.py").write_text("import os\n")
    with pytest.raises(SystemExit, match="sitecustomize.py"):
        rv.require_record_tree(tmp_path, "verification")
    for prefix in (None, str(tmp_path / "cache")):
        monkeypatch.setattr(sys, "pycache_prefix", prefix)
        with pytest.raises(SystemExit, match="PYTHONPYCACHEPREFIX"):
            rv.require_pycache_outside(tmp_path)
    monkeypatch.setattr(sys, "pycache_prefix", str(tmp_path.parent / "elsewhere"))
    rv.require_pycache_outside(tmp_path)


class _ScriptedContext:
    """Stands in for the calibration's context: NC1 fails in batteries 0 and 1, NC2 in 1 and 2."""

    alpha = 0.01
    planned = {"NC1": {"BT-STATE a", "BT-STATE b"}, "NC2": {"BT-STATE c"}}

    def __init__(self):
        self.calls = 0

    def battery(self, rng):
        k, self.calls = self.calls, self.calls + 1
        fail = {"NC1": k in (0, 1), "NC2": k in (1, 2)}
        return {f: {n: {"p_value": 1e-12 if fail[f] else 0.5} for n in sorted(self.planned[f])}
                for f in ("NC1", "NC2")}


def test_the_calibration_counts_the_event_r1_fails_on(monkeypatch):
    """The review's point: C1 is per family, R1 fails when either family fails, and the union was never counted."""
    cal = _calibration()
    monkeypatch.setattr(cal, "_CONTEXT", _ScriptedContext())
    raw = cal.merge(cal.run_chunk(f"{B2_TEST_NAMESPACE}-union", [0, 1]), cal.run_chunk(f"{B2_TEST_NAMESPACE}-union",
                                                                                       [2, 3]))
    assert [raw["families"][f]["failed_batteries"] for f in ("NC1", "NC2")] == [2, 2]
    event = cal.r1_event(raw)
    assert event["failed_batteries"] == 3 and event["false_fail_ci95"] == cal.clopper_pearson(3, 4)


def _record_repo(tmp_path, rv, monkeypatch):
    """A Git repository holding the study's protocol.json and design.json, with the scripts pointed at it."""
    repo = tmp_path / "repo"
    study = repo / rv.STUDY
    study.mkdir(parents=True)
    for name in ("protocol.json", "design.json"):
        (study / name).write_bytes((REPO / rv.STUDY / name).read_bytes())
    _git(repo, "init", "-q")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "frozen")
    monkeypatch.setattr(rv, "HERE", study)
    monkeypatch.setattr(rv, "REPO", repo)
    monkeypatch.setattr(rv, "frozen_hashes", lambda smoke=False: {"frozen": True, "files": {}})
    monkeypatch.setattr(sys, "pycache_prefix", str(tmp_path / "pycache"))
    return repo, study


def test_the_verification_main_needs_the_committed_record_calibration_and_passes_it_on(tmp_path, monkeypatch):
    """The review's mutation: with main() back on the old guards and no calibration passed on, the suite passed."""
    rv = _verification()
    repo, study = _record_repo(tmp_path, rv, monkeypatch)
    protocol, results, calls = json.loads((study / "protocol.json").read_text()), study / "results", []

    def fake_study(args, frozen, protocol, design, code, out_dir, manifests, attempt, calibration, tree):
        calls.append({"code": code, "attempt": attempt, "calibration": calibration, "tree": tree})
        return {"verdicts": {"R1": "PASS"}, "code": code}, {"m.json": "{}"}      # main() publishes it
    monkeypatch.setattr(rv, "study", fake_study)
    printed = []                                   # the summary comes after the completed record (item 63)
    monkeypatch.setattr(rv, "print_summary", lambda result: printed.append(json.loads(
        (results / "attempts" / "verification-1.json").read_text())["status"]))
    monkeypatch.setattr(sys, "argv", ["run_verification.py"])
    with pytest.raises(SystemExit, match="missing"):
        rv.main()
    _calibration_record(rv, repo, results, protocol, {"frozen": True, "files": {}}, development=True)
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "a development calibration")
    with pytest.raises(SystemExit, match="development"):
        rv.main()
    assert not calls and not (results / "attempts" / "verification-1.json").exists()
    _calibration_record(rv, repo, results, protocol, {"frozen": True, "files": {}})
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "the record calibration")
    assert rv.main() == 0 and printed == ["completed"]
    (call,) = calls
    assert call["calibration"]["R0"]["NC1"]["C1_holds"] and call["attempt"]["record"]["attempt"] == 1
    assert call["tree"] == {"attempt_files": [], "other_changes": [], "ignored_allowed": []}
    assert call["code"]["git_dirty"] is False
    assert json.loads((results / "attempts" / "verification-1.json").read_text())["status"] == "completed"
    published = json.loads((results / "verification.json").read_text())         # through publish_verification
    assert published["attempt"]["number"] == 1 and (results / "attempts" / "verification-1" / "manifests" /
                                                    "m.json").exists()


def test_the_calibration_main_needs_a_committed_tree_and_takes_its_code_when_it_starts(tmp_path, monkeypatch):
    rv, cal = _verification(), _calibration()
    repo, study = _record_repo(tmp_path, rv, monkeypatch)
    monkeypatch.setattr(cal, "HERE", study)
    monkeypatch.setattr(cal, "REPO", repo)
    seen = []

    def fake_calibrate(verify, namespace, batteries, out, frozen, args, workers, attempt, tree=None, *, code=None):
        seen.append({"attempt": attempt, "tree": tree, "code": code})
        verify.write_new(out, json.dumps({"attempt": {"number": attempt["record"]["attempt"]}, "code": code}))
        return 0
    monkeypatch.setattr(cal, "calibrate", fake_calibrate)
    printed = []                                   # the summary comes after the completed record (item 63)
    monkeypatch.setattr(cal, "print_summary", lambda verify, result, out: printed.append(json.loads(
        (out.parent / "attempts" / "calibration-1.json").read_text())["status"]))
    monkeypatch.setattr(sys, "argv", ["calibration.py", "--workers", "1"])
    (repo / "notes.txt").write_text("x")
    with pytest.raises(SystemExit, match="committed tree"):
        cal.main()
    (repo / "notes.txt").unlink()
    assert cal.main() == 0 and printed == ["completed"]
    (call,) = seen
    assert call["code"]["git_dirty"] is False and call["code"] == call["attempt"]["record"]["code"]
    assert call["tree"] == {"attempt_files": [], "other_changes": [], "ignored_allowed": []}
    with pytest.raises(SystemExit, match="committed tree"):           # its result is not committed yet
        cal.main()


def test_benchmark_repeats_must_agree_on_every_field_their_configuration_fixes(monkeypatch):
    """The review's probe: B2 output_bytes [6072, 9999, 6072] was summarised as 6072 without an error."""
    bench = _load("b2_integrity_benchmark", "research/v020-nuclear-b2-bench/benchmark.py")
    sample = {"wall_seconds": 1.0, "cpu_seconds": 1.0, "words": 10, "output_bytes": 6072}
    outputs = iter([sample, dict(sample, output_bytes=9999), sample, dict(sample, traced_peak_bytes=1)])
    monkeypatch.setattr(bench.subprocess, "run", lambda *args, **kw: subprocess.CompletedProcess(
        args, 0, stdout=json.dumps(next(outputs)) + "\n", stderr=""))
    with pytest.raises(RuntimeError, match="output_bytes"):
        bench.measure({"engine": "b2", "n": 10}, repeats=3, threads="1")
    assert "populations_sum" in bench.fixed_fields({"engine": "layer_a"})
    assert "output_bytes" not in bench.fixed_fields({"engine": "layer_a"})


def test_the_binomial_z_test_subtracts_exactly():
    """The review's probe: z_check rounds S and T p before subtracting; at T = 2**60 + 12345, p = 1 - 2**-53,
    S = T - 200 it gave z = 0 (p = 1) where the exact z is -6.36 (BIG-DECAY now uses binomial_z_check)."""
    from decimal import Decimal, localcontext
    from rnglab.nuclear.comparison import binomial_z_check, z_check
    size, p = 2 ** 60 + 12345, 1.0 - 2.0 ** -53
    a, b = p.as_integer_ratio()
    with localcontext() as context:
        context.prec = 60
        exact = Decimal((size - 200) * b - size * a) / Decimal(size * a * (b - a)).sqrt()
    assert binomial_z_check(size - 200, size, p)["z"] == float(exact) and abs(float(exact) + 6.364) < 1e-3
    assert z_check(size - 200, size * p, size * p * (1 - p))["z"] == 0.0           # the cancellation it replaces
    assert binomial_z_check(520, 1000, 0.5)["z"] == z_check(520, 500.0, 250.0)["z"]
    for bad in ((5, 4, 0.5), (1, 2, 1.0), (1, 2, float("nan")), (-1, 2, 0.5), (0, 0, 0.5)):
        with pytest.raises(ValueError):
            binomial_z_check(*bad)


def test_traced_draws_reach_both_tail_paths_and_count_refusals():
    """The review's probe: every traced draw had n >= 2**30, above UPPER_TAIL_MAX_N, so S was never evaluated,
    and a non-finite value would have aborted the section instead of being counted."""
    numerics = _load("b2_integrity_numerics", "research/v020-nuclear-b2-numerics/numerics.py")
    traced = numerics._traced_draws()
    assert traced["evaluations"]["F"] > 0 and traced["evaluations"]["S"] > 0
    assert traced["nonfinite"] == {"F": 0, "S": 0} and traced["refused_draws"] == 0


def test_numpy_binomial_is_measured_where_it_fails_and_where_the_big_null_calls_it(ns, monkeypatch):
    """The review's probe: numpy's binomial at n = 2**60 + 12345 has a variance 1.22 times too large at mean 115;
    the design's note said the BIG null was exact up to the lattice."""
    numerics = _load("b2_integrity_numerics", "research/v020-nuclear-b2-numerics/numerics.py")
    monkeypatch.setattr(numerics, "NAMESPACE", f"{ns}-numpy-binomial")
    rows = {(r["n"], round(r["mean"])): r for r in numerics.numpy_binomial_section(draws=100_000)["rows"]}
    assert rows[(2 ** 60 + 12345, 115)]["variance_ratio"] > 1.1
    assert rows[(2 ** 53 - 1, 115)]["variance_ratio"] > 1.02          # the fifth review's check: below 2**53 too
    design = json.loads((REPO / "research/v020-nuclear-b2-verify/design.json").read_text())
    calls = design["calibration_null_lattice"]["numpy_calls_large_n"]
    assert calls["n_above"] == 1e13 and calls["calls"] > 0 and calls["smallest_expected_count"] >= 1e12
    last = numerics.numpy_binomial_section(draws=1000)["rows"][-1]
    assert last["n"] == round(calls["smallest"]["n"]) and last["p"] == calls["smallest"]["p"]
    assert last["expected_count"] == pytest.approx(calls["smallest_expected_count"], rel=1e-9)


def test_withdrawn_claims_are_gone_from_the_docstrings():
    """The review's probe: three docstrings still claimed what items 34-35 withdrew."""
    for path in ("research/v020-nuclear-b2-numerics/numerics.py", "src/rnglab/nuclear/binomial.py"):
        text = (REPO / path).read_text()
        for claim in ("every point the bisection can reach", "evaluates it exactly", "tails (exact)"):
            assert claim not in text, (path, claim)


def test_the_tail_probe_draws_the_study_statistics(ns):
    """The review's request: per-test null tail rates from a script, with each test's own statistic."""
    probe = _load("b2_integrity_tail_probe", "research/v020-nuclear-b2-verify/tail_probe.py")
    probe._init_worker()
    neutron = probe.neutron_chunk(f"{ns}-tail-probe", "G1", 0, 500)        # checks 100 draws against chi2
    assert neutron["draws"] == 500 and neutron["checked"] == 100 and neutron["missing"] == 0
    photons = probe.photon_chunk(f"{ns}-tail-probe", "S1", 0, 50)
    assert photons and all(c["draws"] == 50 and c["missing"] == 0 for c in photons.values())


# --- 12. Sixth review ------------------------------------------------------------------------------------------

def test_source_plans_are_read_without_copying_the_nested_mapping():
    """The review's probe: every read of plan.interval_plans copies the whole nested mapping, so the loops that
    indexed it per interval and nuclide were quadratic in the number of distinct lengths."""
    plan = AggregatePlan(_simple(), [0.0, 1.0, 3.0, 3.5], 100)
    pairs = [(delta, i, source) for delta, plans in plan.interval_plans.items() for i, source in plans.items()]
    assert len(pairs) == 3 and all(plan.source_plan(delta, i) is source for delta, i, source in pairs)
    with pytest.raises(KeyError):
        plan.source_plan(0.25, 0)


def test_binomial_z_check_records_its_exact_inputs():
    """The review's probe: the record held estimate, expected and variance as floats, from which the z of a
    BIG-DECAY result could not be recomputed (0.0 against the recorded -6.36)."""
    from rnglab.nuclear.comparison import binomial_z_check
    size, p = 2 ** 60 + 12345, 1.0 - 2.0 ** -53
    result = json.loads(json.dumps(binomial_z_check(size - 200, size, p)))
    assert (result["successes"], result["trials"], result["probability"]) == (size - 200, size, p)
    assert binomial_z_check(result["successes"], result["trials"], result["probability"])["z"] == result["z"]


def test_a_live_process_is_an_earlier_attempt_only_if_its_start_time_matches(tmp_path):
    """The review's probe: a dead attempt's record carrying the pid of an unrelated live process whose command
    line contained the script's name refused the launch, with a message about a deleted lock file."""
    rv = _verification()
    out = tmp_path / "results"
    attempt = rv.begin_attempt(out, "calibration", {"git_commit": "c"})
    rv.release_attempt(attempt)                                         # as if its process had died
    other = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(60)", "calibration.py"])
    try:
        record = dict(json.loads(attempt["path"].read_text()), pid=other.pid)
        attempt["path"].write_text(json.dumps(record))
        assert record["process_start"] is not None and not rv._alive_here(record)
        second = rv.begin_attempt(out, "calibration", {"git_commit": "c"})
        assert second["record"]["interrupted_attempts"] == [1]
        rv.release_attempt(second)
        assert rv._alive_here(dict(record, process_start=rv._process_start(other.pid)))
    finally:
        other.kill()
        other.wait()


def test_summaries_follow_completion_and_a_closed_output_cannot_fail_an_attempt(tmp_path, monkeypatch):
    """The sixth review's probe: a broken pipe while printing the summary, after the result had been published,
    left the attempt 'failed' beside a published result; the seventh: printing the summary before publication
    showed statistics before completion. calibrate() now prints no statistic, main() prints the summary after
    the completed record, and a closed output is ignored."""
    import builtins
    from types import SimpleNamespace
    cal, rv = _calibration(), _verification()
    monkeypatch.setattr(cal, "ProcessPoolExecutor", _InlineExecutor)
    for name in cal.THREAD_VARIABLES:
        monkeypatch.setenv(name, "1")
    real_print = builtins.print

    def broken(*args, **kwargs):
        if args and str(args[0]).startswith("NC1:"):
            raise BrokenPipeError(32, "Broken pipe")
        real_print(*args, **kwargs)
    monkeypatch.setattr(builtins, "print", broken)
    attempt = rv.begin_attempt(tmp_path, "calibration", {"git_commit": "c"})
    result = tmp_path / "c.json"
    assert rv.run_attempt(attempt, lambda: cal.calibrate(rv, f"{B2_TEST_NAMESPACE}-publish-last", 2, result, {},
                                                         SimpleNamespace(chunk=2, development=False), 1, attempt),
                          result) == 0
    assert result.exists() and json.loads(attempt["path"].read_text())["status"] == "completed"
    output = open(tmp_path / "stdout.txt", "w")      # a real file, whatever pytest's capture mode
    monkeypatch.setattr(sys, "stdout", output)
    redirected = []
    monkeypatch.setattr(rv.os, "dup2", lambda fd, target: redirected.append(target))
    cal.print_summary(rv, json.loads(result.read_text()), result)     # the "NC1:" line hits a broken pipe
    import errno
    monkeypatch.setattr(builtins, "print", lambda *a, **k: (_ for _ in ()).throw(OSError(errno.ENOSPC, "full")))
    rv.quiet_print("a full disk")                                     # not only a broken pipe
    assert redirected == [output.fileno(), output.fileno()]
    output.close()


def test_the_verification_publishes_last_and_only_with_the_code_it_started_with(tmp_path, monkeypatch):
    """The review's mutation: removing require_same_code from the study went unnoticed by the suite."""
    rv = _verification()
    start = rv.record_code(REPO)
    attempt = rv.begin_attempt(tmp_path, "verification", start)
    manifests = tmp_path / "manifests"
    manifests.mkdir()
    monkeypatch.setattr(rv, "record_code", lambda repo: dict(start, git_commit="changed"))
    with pytest.raises(RuntimeError, match="changed while it ran"):
        rv.publish_verification({"verdicts": {}}, {"m.json": "{}"}, tmp_path, manifests, attempt, start)
    assert not (tmp_path / "verification.json").exists() and not (manifests / "m.json").exists()
    monkeypatch.setattr(rv, "record_code", lambda repo: start)
    rv.publish_verification({"verdicts": {}}, {"m.json": "{}"}, tmp_path, manifests, attempt, start)
    assert json.loads((tmp_path / "verification.json").read_text())["attempt"]["number"] == 1
    rv.release_attempt(attempt)


def test_both_main_functions_finish_a_published_but_unrecorded_attempt(tmp_path, monkeypatch):
    """The review's mutation: disabling complete_published_attempt in either main() went unnoticed."""
    rv, cal = _verification(), _calibration()
    repo, study = _record_repo(tmp_path, rv, monkeypatch)
    monkeypatch.setattr(cal, "HERE", study)
    monkeypatch.setattr(cal, "REPO", repo)
    for script, main, argv in (("calibration", cal.main, ["calibration.py", "--workers", "1"]),
                               ("verification", rv.main, ["run_verification.py"])):
        results, code = study / "results", rv.record_code(repo)
        attempt = rv.begin_attempt(results, script, code)
        rv.write_new(results / f"{script}.json",
                     json.dumps({"attempt": {"number": attempt["record"]["attempt"]}, "code": code}))
        rv.release_attempt(attempt)                                      # killed before its record said so
        monkeypatch.setattr(sys, "argv", argv)
        with pytest.raises(SystemExit, match="now shows it completed"):
            main()
        assert json.loads(attempt["path"].read_text())["status"] == "completed"


def test_the_verification_main_refuses_an_untracked_module_beside_the_code(tmp_path, monkeypatch):
    """The review's mutation: require_record_tree replaced by an empty tree went unnoticed."""
    rv = _verification()
    repo, study = _record_repo(tmp_path, rv, monkeypatch)
    (repo / "src").mkdir()
    (repo / "src" / "sitecustomize.py").write_text("import os\n")
    monkeypatch.setattr(sys, "argv", ["run_verification.py"])
    with pytest.raises(SystemExit, match="committed tree.*sitecustomize.py"):
        rv.main()


def test_the_report_refuses_a_calibration_other_than_the_one_the_verification_read(tmp_path, monkeypatch):
    """The review's mutations: disabling the report's sha256 check, or skipping a missing calibration file,
    went unnoticed."""
    rt = _load("b2_integrity_report_tables", "research/v020-nuclear-b2-verify/report_tables.py")
    (tmp_path / "results").mkdir()
    (tmp_path / "results" / "verification.json").write_text(json.dumps({"calibration": {"sha256": "0" * 64}}))
    monkeypatch.setattr(rt, "HERE", tmp_path)
    monkeypatch.setattr(sys, "argv", ["report_tables.py"])
    with pytest.raises(SystemExit, match="not the file the verification read"):
        rt.main()                                                     # calibration.json missing
    (tmp_path / "results" / "calibration.json").write_text("{}\n")
    with pytest.raises(SystemExit, match="not the file the verification read"):
        rt.main()                                                     # another calibration.json


# --- 13. Eighth review -----------------------------------------------------------------------------------------

def test_the_start_time_of_a_process_whose_name_is_not_utf8_is_read():
    """The review's probe: a process named b"\\xff\\xfe) z" gave no start time, so the live-attempt check fell
    back to its cautious answer."""
    rv = _verification()
    code = ("import ctypes, time; ctypes.CDLL(None).prctl(15, bytes([255, 254, 41, 32, 122]), 0, 0, 0); "
            "print('ready', flush=True); time.sleep(60)")
    child = subprocess.Popen([sys.executable, "-c", code], stdout=subprocess.PIPE, text=True)
    try:
        assert child.stdout.readline().strip() == "ready"
        assert bytes([255, 254, 41]) in Path(f"/proc/{child.pid}/stat").read_bytes()
        assert isinstance(rv._process_start(child.pid), int)
    finally:
        child.kill()
        child.wait()
