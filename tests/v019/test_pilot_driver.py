"""Stage D pilot driver: plan, storage, resume, corruption, manifests, authorization.

Synthetic fixtures only (tests/v019/_pilot_fixture.py). No sampler run in a pilot or smoke
namespace, and nothing is written under research/. Tests that need the measurement layer are
skipped when research/v019-d-pilot/measure.py is absent (set RNGLAB_DPILOT_MEASURE to supply it).
"""
import dataclasses
import importlib.util as _ilu
import json
import os
import signal
import subprocess
import sys as _sys
from pathlib import Path as _Path

import numpy as np
import pytest

_spec = _ilu.spec_from_file_location("_pilot_fixture", _Path(__file__).with_name("_pilot_fixture.py"))
fx = _sys.modules.get("_pilot_fixture") or _ilu.module_from_spec(_spec)
if "_pilot_fixture" not in _sys.modules:
    _sys.modules["_pilot_fixture"] = fx
    _spec.loader.exec_module(fx)

au, pc, pg = fx.au, fx.pc, fx.pg
PROTOCOL = json.loads((au.REPO / "research/v019-d-pilot/protocol.json").read_text())
# The independent side: it never imports the pilot's modules.
RC = au.recount()
IG = au.load("v019_audit_integrity", au.REPO / "research/v019-audit/integrity.py")
VERIFY = au.load("v019_audit_verify", au.REPO / "research/v019-audit/verify.py")
PF = au.load("v019_audit_preflight", au.REPO / "research/v019-audit/preflight.py")


def light(root, plan=fx.SMALL_PLAN):
    """A context that needs no measurement layer: enough for failure records and storage states."""
    return fx.context(plan, root, need_measure=False)


# --------------------------------------------------------------------------- plan and constants
def test_frozen_plan_and_constants_match_the_protocol():
    plan = pc.FROZEN_PLAN
    assert plan.total() == 5999 == PROTOCOL["design_number_snapshot"]["total_runs"] and plan.namespace == PROTOCOL["data"]["namespace"]
    for (key, family, count), declared in zip(plan.conditions, PROTOCOL["data"]["conditions"]):
        assert key == declared["key"] and family == declared["family"]
        assert declared["replicas"] == f"range(0,{count})"
    assert PROTOCOL["data"]["conditions"][2]["fitting"] == "range(0,600)" and list(plan.fitting) == list(range(600))
    assert PROTOCOL["data"]["conditions"][2]["evaluation"] == "range(600,1000)" and list(plan.evaluation) == list(range(600, 1000))
    assert plan.split("fault/DUP", 599) == "fitting" and plan.split("fault/DUP", 600) == "evaluation"
    assert plan.split("calibration", 5) == "calibration" and plan.split("healthy", 5) == "validation"
    assert pc.SMOKE_PLAN.namespace == PROTOCOL["smoke"]["namespace"] and pc.SMOKE_PLAN.total() == 10
    criteria = PROTOCOL["stop_criteria"]
    assert pc.MINIMUM_ELIGIBLE == criteria["2"]["minimum_eligible_per_family"] == criteria["3"]["minimum_eligible_per_family"]
    assert pc.MINIMUM_ELIGIBLE_EVALUATION == criteria["4"]["minimum_eligible_evaluation"]
    assert pc.CRITERION_LENGTH == criteria["2"]["L"] == criteria["3"]["L"] == criteria["4"]["L"]
    assert pc.CRITERION_4_LEVEL == PROTOCOL["design_number_snapshot"]["criterion_4"]["interval_level"]
    assert pc.OTHER_LEVEL == PROTOCOL["scores"]["bootstrap"]["other_level"]
    assert list(pc.LENGTHS) == PROTOCOL["application"]["L"] and pc.H == PROTOCOL["application"]["H"]
    assert "8 processes" in PROTOCOL["storage"]["parallelism"] and pc.PROTOCOL_WORKERS == 8
    assert all(name in PROTOCOL["storage"]["parallelism"] for name in pc.THREAD_VARIABLES)
    assert "2026-12-22" in criteria["5"]["calculation"] and pc.ENGINEERING_DEADLINE == "2026-12-22"
    assert pc.TARGETS == {"SHIFT": "frequency_signed_z", "LAG1": "lag1_feature", "DUP": "duplicate_rate"}


def test_plan_rejects_overlapping_or_misplaced_splits():
    conditions = fx.SMALL_PLAN.conditions
    with pytest.raises(ValueError, match="overlap"):
        pc.Plan(fx.NAMESPACE, conditions, range(0, 2), range(1, 2), 2, 2, "test")   # replica 1 on both sides
    with pytest.raises(ValueError, match="outside"):
        pc.Plan(fx.NAMESPACE, conditions, range(0, 1), range(1, 3), 2, 2, "test")
    with pytest.raises(ValueError, match="kind"):
        pc.Plan(fx.NAMESPACE, conditions, range(0, 1), range(1, 2), 2, 2, "production")
    assert not set(pc.FROZEN_PLAN.fitting) & set(pc.FROZEN_PLAN.evaluation)


def test_canonical_encoding_and_atomic_write(tmp_path):
    assert pc.canonical_json({"b": 1, "a": [1.5, None]}) == b'{\n  "a": [\n    1.5,\n    null\n  ],\n  "b": 1\n}\n'
    for bad in (float("nan"), float("inf"), -float("inf")):
        with pytest.raises(ValueError):
            pc.canonical_json({"x": bad})       # never JSON NaN or Infinity
        with pytest.raises(ValueError):
            pc.canonical_jsonl([{"x": bad}])
    target = tmp_path / "a" / "file.json"
    digest = pc.atomic_write(target, b"one")
    assert target.read_bytes() == b"one" and digest == pc.sha256_bytes(b"one") == pc.sha256_file(target)
    with pytest.raises(pc.PilotError, match="not overwritten"):
        pc.atomic_write(target, b"two")
    assert target.read_bytes() == b"one" and not list(target.parent.glob("*.tmp"))
    pc.atomic_write(target, b"two", overwrite=True)
    assert target.read_bytes() == b"two"
    assert pc.finite_or_none(float("nan")) is None and pc.finite_or_none(None) is None and pc.finite_or_none(2) == 2.0
    assert pc.finite_or_none("3") is None


# --------------------------------------------------------------------------- failure records and resume
def test_failure_records_are_complete_and_never_rerun(tmp_path):
    ctx = light(tmp_path)
    out = pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    assert out["counts"] == {"COMPLETE": 0, "FAILED_RECORD": 10, "MISSING": 0, "PARTIAL": 0, "CORRUPT": 0}
    record = pc.read_json(ctx.record_path("fault/LAG1", 1))
    assert record["status"] == "FAILED" and record["trace_path"] is None and record["split"] == "evaluation"
    assert record["failure"] == {"phase": "generation", "exception": {"class": "ValueError", "message": "synthetic failure"},
                                 "completed_sweeps": 7, "last_completed_sweep": 6, "local_bytes": 1024}
    assert record["seeds"] == ctx.dlib.seeds_for(fx.NAMESPACE, "fault/LAG1", 1)
    events = pg.journal_events(ctx)
    assert sum(e["event"] == "start" for e in events) == 10 == sum(e["event"] == "commit" for e in events)
    again = pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)     # a recorded failure is not run again
    assert again["counts"]["FAILED_RECORD"] == 10 and sum(e["event"] == "start" for e in pg.journal_events(ctx)) == 10
    manifest = pc.read_json(ctx.result("run_manifest.json"))
    assert manifest["planned"] == 10 and all(e["state"] == "FAILED_RECORD" and e["attempts"] == 1 for e in manifest["runs"])
    assert [(e["condition"], e["replica"]) for e in manifest["runs"]] == [(c, r) for c, _, r in ctx.plan.runs()]
    failures = [json.loads(line) for line in ctx.result("failures.jsonl").read_text().splitlines()]
    assert len(failures) == 10 and all(f["disposition"] == "recorded, not replaced" for f in failures)
    with pytest.raises(pc.PilotError, match="not replaced"):
        pg.commit_run(ctx, record, None)
    assert pg.require_complete(ctx) == {"complete": 0, "failed_records": 10}


# --------------------------------------------------------------------------- interrupted runs: regenerate_same_id
def events_of(ctx, condition, replica):
    return [e for e in pg.journal_events(ctx) if e.get("condition") == condition and e.get("replica") == replica]


def manifest_entry(ctx, condition, replica):
    return next(e for e in pc.read_json(ctx.result("run_manifest.json"))["runs"]
                if (e["condition"], e["replica"]) == (condition, replica))


def attempt_log(ctx):
    """The independent check of the audit log, on the same tree."""
    return VERIFY.check_attempt_log(ctx.root, ctx.plan.namespace, fx.audit_plan(ctx.plan))


def log_start(ctx, condition, replica, attempt=1):
    """What the driver writes before it draws; the test then behaves as if the process had died."""
    seeds = ctx.dlib.seeds_for(ctx.plan.namespace, condition, replica)
    return pg.journal(ctx, "start", condition=condition, replica=replica, attempt=attempt, seeds=seeds,
                      rng_input_sha256=pg.rng_input(ctx.dlib, seeds)["sha256"], environment_sha256=pg.environment_sha256())


def test_the_decided_policy_is_the_only_one():
    rules = pg.INTERRUPTED_RUN_RULES
    assert pg.POLICY == rules["policy"] == "regenerate_same_id" == RC.INTERRUPTED_RUN_POLICY
    assert pg.MAX_INFRASTRUCTURE_RETRIES == rules["max_infrastructure_retries"] == 3 == RC.MAX_INFRASTRUCTURE_RETRIES
    assert pg.MAX_ATTEMPTS == 4 == VERIFY.MAX_ATTEMPTS and rules["decided_by"] == "author"
    assert tuple(rules["infrastructure_interruptions"]) == RC.INFRASTRUCTURE_REASONS
    assert pg.INFRASTRUCTURE_PHASE == RC.INFRASTRUCTURE_PHASE and pg.UNRESOLVED == RC.UNRESOLVED_INFRASTRUCTURE
    assert not hasattr(pg, "POLICIES")
    assert PROTOCOL["failed_runs"]["recording"].startswith("Every planned ID gets a success or failure record")


def test_a_run_that_was_started_and_left_no_record_is_run_again_with_the_same_identity(tmp_path):
    ctx = light(tmp_path)
    first = log_start(ctx, "healthy", 1)                           # the process died before any record
    out = pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    assert out["counts"]["FAILED_RECORD"] == 10 and out["unresolved_infrastructure"] == [] and not out["stopped_on_request"]
    events = events_of(ctx, "healthy", 1)
    assert [e["event"] for e in events] == ["start", "interrupted", "restart", "start", "commit"]
    assert events[1]["attempt"] == 1 and events[1]["detected_by"] == "resume"
    assert events[1]["reason"] == {"kind": "process_ended_without_record", "detail": "no record and no trace"}
    assert events[3]["attempt"] == 2 and events[3]["seeds"] == first["seeds"] == ctx.dlib.seeds_for(fx.NAMESPACE, "healthy", 1)
    assert events[3]["rng_input_sha256"] == first["rng_input_sha256"] and events[3]["environment_sha256"] == pg.environment_sha256()
    assert events[4]["record_sha256"] == pc.sha256_file(ctx.record_path("healthy", 1)) and events[4]["status"] == "FAILED"
    entry = manifest_entry(ctx, "healthy", 1)
    assert entry["attempts"] == 2 and entry["infrastructure_status"] is None and len(entry["interruptions"]) == 1
    assert entry["attempt_environments"] == [pg.environment_sha256()]
    assert manifest_entry(ctx, "healthy", 0)["attempts"] == 1 and manifest_entry(ctx, "healthy", 0)["interruptions"] == []
    record = pc.read_json(ctx.record_path("healthy", 1))
    assert record["failure"]["phase"] == "generation" and record["seeds"] == first["seeds"]    # the run's own failure stands
    infrastructure = pc.read_json(ctx.result("run_manifest.json"))["infrastructure"]
    assert infrastructure == {"policy": "regenerate_same_id", "max_infrastructure_retries": 3, "runs_with_more_than_one_attempt": 1,
                              "interruptions": 1, "unresolved": [], "environments": [pg.environment_sha256()]}
    check = attempt_log(ctx)
    assert check["agrees"] and check["runs_with_more_than_one_attempt"] == 1 and check["interruptions"] == 1
    assert check["unresolved"] == [] and check["status"] == "CHECKED"


@pytest.mark.parametrize("kind,reason", [("OSError", "operating_system_error"), ("MemoryError", "memory_error"),
                                         ("wrapped", "operating_system_error")])
def test_an_error_of_the_system_leaves_no_record_and_is_retried(tmp_path, kind, reason):
    ctx = light(tmp_path)
    fx.set_interruptions(tmp_path, {"healthy/0": [kind]})
    with pytest.raises(pg.InfrastructureInterruption, match="no record was written") as caught:
        pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert caught.value.reason["kind"] == reason
    assert pg.classify(ctx, "healthy", 0)["state"] == "MISSING"                 # no failure record for an interruption
    assert pg.classify(ctx, "calibration", 1)["state"] == "FAILED_RECORD"       # what came before is stored
    assert pg.classify(ctx, "healthy", 1)["state"] == "MISSING"                 # the process stopped: nothing after it was started
    events = events_of(ctx, "healthy", 0)
    assert [e["event"] for e in events] == ["start", "interrupted"]
    assert events[1]["reason"]["kind"] == reason and events[1]["detected_by"] == "attempt" and events[1]["attempt"] == 1
    assert not events_of(ctx, "healthy", 1)
    assert pg.journal_events(ctx)[-1]["event"] == "generate_end"                # the invocation still closed its log
    out = pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert out["counts"]["FAILED_RECORD"] == 10
    assert [e["event"] for e in events_of(ctx, "healthy", 0)] == ["start", "interrupted", "restart", "start", "commit"]
    assert manifest_entry(ctx, "healthy", 0)["attempts"] == 2
    assert manifest_entry(ctx, "healthy", 0)["interruptions"][0]["reason"]["kind"] == reason
    assert attempt_log(ctx)["agrees"]


def test_the_runs_own_failure_is_evidence_and_is_never_retried(tmp_path):
    ctx = light(tmp_path)
    fx.set_interruptions(tmp_path, {"fault/DUP/0": ["ValueError"]})             # an exception the engine did not classify
    out = pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert out["counts"]["FAILED_RECORD"] == 10 and out["unresolved_infrastructure"] == []
    record = pc.read_json(ctx.record_path("fault/DUP", 0))
    assert record["failure"] == {"phase": "execution", "exception": {"class": "ValueError", "message": "a numerical failure"}}
    assert [e["event"] for e in events_of(ctx, "fault/DUP", 0)] == ["start", "commit"]
    pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert [e["event"] for e in events_of(ctx, "fault/DUP", 0)] == ["start", "commit"]     # not run again
    assert all(e["attempts"] == 1 for e in pc.read_json(ctx.result("run_manifest.json"))["runs"])
    for error in (ValueError("x"), ArithmeticError("x"), FloatingPointError("x"), RuntimeError("serialization: nonfinite array")):
        assert pg.infrastructure_cause(error) is None
    chained = RuntimeError("outer")
    chained.__context__ = MemoryError()
    assert isinstance(pg.infrastructure_cause(chained), MemoryError)
    assert isinstance(pg.infrastructure_cause(FileNotFoundError("gone")), OSError)


def test_a_memory_error_recorded_inside_the_measurement_is_an_interruption(tmp_path):
    ctx = light(tmp_path)
    fx.set_interruptions(tmp_path, {"calibration/1": ["absorbed"]})
    with pytest.raises(pg.InfrastructureInterruption, match="memory_error"):
        pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert pg.classify(ctx, "calibration", 1)["state"] == "MISSING"            # the record that held it was not stored
    events = events_of(ctx, "calibration", 1)
    assert [e["event"] for e in events] == ["start", "interrupted"] and events[1]["reason"]["kind"] == "memory_error"
    assert "Unable to allocate" in events[1]["reason"]["detail"]
    pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert pg.classify(ctx, "calibration", 1)["state"] == "FAILED_RECORD" and manifest_entry(ctx, "calibration", 1)["attempts"] == 2
    assert pg.recorded_memory_error(pc.read_json(ctx.record_path("calibration", 1))) is None
    assert pg.recorded_memory_error({"conventional": [{"failures": [{"exception": {"class": "MemoryError", "message": "m"}}]}]}) == {
        "class": "MemoryError", "message": "m"}
    for harmless in ({"exception": {"class": "ValueError", "message": "MemoryError"}}, {"class": "ZeroDivisionError"}, [], None, "MemoryError"):
        assert pg.recorded_memory_error(harmless) is None
    assert attempt_log(ctx)["agrees"]


def test_three_retries_then_the_run_is_unresolved_and_never_run_again(tmp_path):
    ctx = light(tmp_path)
    fx.set_interruptions(tmp_path, {"fault/LAG1/1": ["OSError", "MemoryError", "OSError", "OSError", "OSError"]})
    for attempt in range(1, 5):
        with pytest.raises(pg.InfrastructureInterruption, match=f"attempt {attempt}"):
            pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
        assert pg.classify(ctx, "fault/LAG1", 1)["state"] == "MISSING"
        assert pg.classify(ctx, "fault/DUP", 0)["state"] == "MISSING"          # never reached while the cause persists
    out = pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)   # the fifth invocation does not try again
    assert out["unresolved_infrastructure"] == ["fault/LAG1/00001"] and out["counts"]["FAILED_RECORD"] == 10
    record = pc.read_json(ctx.record_path("fault/LAG1", 1))
    failure = record["failure"]
    assert record["status"] == "FAILED" and record["trace_path"] is None and record["seeds"] == ctx.dlib.seeds_for(fx.NAMESPACE, "fault/LAG1", 1)
    assert failure["phase"] == "infrastructure" and failure["infrastructure_status"] == "UNRESOLVED_INFRASTRUCTURE"
    assert failure["exception"] is None and failure["attempts"] == 4 and failure["max_infrastructure_retries"] == 3
    assert [(i["attempt"], i["reason"]["kind"], i["detected_by"]) for i in failure["interruptions"]] == [
        (1, "operating_system_error", "attempt"), (2, "memory_error", "attempt"), (3, "operating_system_error", "attempt"),
        (4, "operating_system_error", "attempt")]
    kinds = [e["event"] for e in events_of(ctx, "fault/LAG1", 1)]
    assert kinds == ["start", "interrupted", "restart"] * 3 + ["start", "interrupted", "unresolved"]
    pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert [e["event"] for e in events_of(ctx, "fault/LAG1", 1)] == kinds        # final: no fifth attempt, ever
    assert json.loads((tmp_path / fx.CONTROL_FILE).read_text())["fault/LAG1/1"] == ["OSError"]   # the fifth instruction was never used
    entry = manifest_entry(ctx, "fault/LAG1", 1)
    assert entry["state"] == "FAILED_RECORD" and entry["attempts"] == 4 and entry["infrastructure_status"] == "UNRESOLVED_INFRASTRUCTURE"
    assert pc.read_json(ctx.result("run_manifest.json"))["infrastructure"]["unresolved"] == ["fault/LAG1/00001"]
    failures = [json.loads(line) for line in ctx.result("failures.jsonl").read_text().splitlines()]
    mine = [f for f in failures if (f["condition"], f["replica"]) == ("fault/LAG1", 1)]
    assert [f["disposition"] for f in mine] == ["infrastructure unresolved, not an experimental failure"]
    assert sum(f["disposition"] == "recorded, not replaced" for f in failures) == 9
    assert pg.require_complete(ctx) == {"complete": 0, "failed_records": 10}       # every planned ID has its record
    unresolved = fx.pa.unresolved_infrastructure(ctx)
    assert unresolved["status"] == "UNRESOLVED_INFRASTRUCTURE" and unresolved["runs"] == [{"condition": "fault/LAG1", "replica": 1}]
    assert unresolved["criteria_that_cannot_pass"] == ["2", "3", "4"] and unresolved["per_condition"]["fault/LAG1"] == 1
    check = attempt_log(ctx)
    assert check["agrees"], check["discrepancies"]
    assert check["unresolved"] == [{"condition": "fault/LAG1", "replica": 1}] and check["criteria_that_cannot_pass"] == [2, 3, 4]
    assert check["interruptions"] == 4 and check["runs_with_more_than_one_attempt"] == 1


def test_seeds_and_random_input_are_checked_against_the_manifest_before_anything_is_drawn(tmp_path):
    ctx = light(tmp_path)
    engine = fx.FailingEngine(ctx)
    seeds = ctx.dlib.seeds_for(fx.NAMESPACE, "healthy", 0)
    fingerprint = pg.rng_input(ctx.dlib, seeds)
    assert set(fingerprint["roles"]) == set(pg.RNG_ROLES) == set(seeds) and len(set(fingerprint["roles"].values())) == 5
    good = {("healthy", 0): {"seeds": seeds, "rng_input_sha256": fingerprint["sha256"]}}
    for planned in ({("healthy", 0): {"seeds": {**seeds, "fault": seeds["fault"] + 1}, "rng_input_sha256": fingerprint["sha256"]}},
                    {("healthy", 0): {"seeds": seeds, "rng_input_sha256": "0" * 64}}, {}):
        with pytest.raises(pc.PilotError, match="nothing was drawn"):
            pg.run_one(ctx, engine, "healthy", "healthy", 0, 1, planned)
    assert not ctx.result("journal.jsonl").exists() and pg.classify(ctx, "healthy", 0)["state"] == "MISSING"
    pg.run_one(ctx, engine, "healthy", "healthy", 0, 1, good)
    start = events_of(ctx, "healthy", 0)[0]
    assert start["seeds"] == seeds and start["rng_input_sha256"] == fingerprint["sha256"]
    manifest = {"planned_runs": [{"condition": "healthy", "replica": 0, "seeds": seeds, "rng_input_sha256": fingerprint["sha256"]}]}
    assert pg.planned_inputs(manifest) == good and pg.planned_inputs(None) is None


def test_the_random_input_hash_is_reproduced_by_the_independent_generators():
    dlib = au.dlib()
    for condition, replica in (("calibration", 0), ("healthy", 7), ("fault/SHIFT", 3), ("fault/LAG1", 15), ("fault/DUP", 1)):
        seeds = dlib.seeds_for(fx.NAMESPACE, condition, replica)
        assert seeds == RC.run_seeds(fx.NAMESPACE, condition, replica)
        mine = pg.rng_input(dlib, seeds)["sha256"]
        assert mine == RC.rng_input_sha256(seeds) and len(mine) == 64
        for role in pg.RNG_ROLES:                                   # every one of the five generators is in the hash
            assert RC.rng_input_sha256({**seeds, role: seeds[role] + 1}) != mine
            assert pg.rng_input(dlib, {**seeds, role: seeds[role] + 1})["sha256"] != mine
    # A fixed value of the audit namespace: the definition does not change unnoticed.
    seeds = dlib.seeds_for(fx.NAMESPACE, "fault/LAG1", 15)
    words = dlib.PCG64Source(seeds["local"]).read_bytes(8 * 128)
    assert pg.rng_input(dlib, seeds)["roles"]["local"] == __import__("hashlib").sha256(words).hexdigest()
    assert words == RC.words_to_bytes(RC.pcg64_words(seeds["local"], 128))
    assert pg.environment_sha256() == pg.environment_sha256() and len(pg.environment_sha256()) == 64


def test_a_record_stored_before_the_log_said_so_is_reconciled_not_rerun(tmp_path):
    ctx = light(tmp_path)
    log_start(ctx, "fault/SHIFT", 1)
    record, _ = fx.FailingEngine(ctx).execute("fault/SHIFT", "SHIFT", 1)
    pg.commit_run(ctx, record, None)                                # stored; the process ended before the commit entry
    pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    events = events_of(ctx, "fault/SHIFT", 1)
    assert [e["event"] for e in events] == ["start", "commit"] and events[1]["detected_by"] == "resume"
    assert events[1]["record_sha256"] == pc.sha256_file(ctx.record_path("fault/SHIFT", 1)) and events[1]["attempt"] == 1
    assert manifest_entry(ctx, "fault/SHIFT", 1)["attempts"] == 1 and attempt_log(ctx)["agrees"]


def test_the_audit_log_is_append_only_and_chained(tmp_path):
    ctx = light(tmp_path)
    pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    path = ctx.result("journal.jsonl")
    original = path.read_bytes()
    events = pg.journal_events(ctx)
    assert [e["seq"] for e in events] == list(range(1, len(events) + 1)) and events[0]["previous_sha256"] is None
    lines = original.splitlines()
    assert all(e["previous_sha256"] == pc.sha256_bytes(line) for e, line in zip(events[1:], lines))
    assert events[0]["event"] == "generate_begin" and events[0]["policy"] == "regenerate_same_id" and events[-1]["event"] == "generate_end"
    assert IG.read_journal(tmp_path)["problems"] == [] and IG.read_journal(tmp_path)["events"] == events
    edited = json.loads(lines[4])
    edited["attempt"] = 3
    tampered = {"an edited entry": lines[:4] + [json.dumps(edited, sort_keys=True).encode()] + lines[5:],
                "a removed entry": lines[:4] + lines[5:], "two entries exchanged": lines[:4] + [lines[5], lines[4]] + lines[6:],
                "a removed first entry": lines[1:]}
    for name, changed in tampered.items():
        path.write_bytes(b"\n".join(changed) + b"\n")
        with pytest.raises(pc.PilotError, match="does not verify"):
            pg.journal_events(ctx)
        with pytest.raises(pc.PilotError, match="does not verify"):
            pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)       # nothing is run on a log that does not verify
        assert IG.read_journal(tmp_path)["problems"], name
        check = attempt_log(ctx)
        assert not check["agrees"] and any(d["field"] == "chain" for d in check["discrepancies"]), name
    path.write_bytes(original + b"not an entry\n")
    with pytest.raises(pc.PilotError, match="is not a journal entry"):
        pg.journal(ctx, "start", condition="healthy", replica=0)
    assert IG.read_journal(tmp_path)["problems"] == ["the journal ends with a line that is not an entry"]
    path.write_bytes(original)
    assert attempt_log(ctx)["agrees"]


def test_a_write_cut_short_is_left_in_place_and_declared(tmp_path):
    ctx = light(tmp_path)
    pg.journal(ctx, "generate_begin", workers=0)
    path = ctx.result("journal.jsonl")
    with open(path, "ab") as handle:
        handle.write(b'{"seq": 2, "event": "sta')                    # power was lost in the middle of a line
    entry = pg.journal(ctx, "generate_begin", workers=0)
    assert entry["seq"] == 2 and entry["after_incomplete_line"] is True
    assert b'{"seq": 2, "event": "sta\n' in path.read_bytes()        # the fragment is kept
    assert [e["seq"] for e in pg.journal_events(ctx)] == [1, 2]
    read = IG.read_journal(tmp_path)
    assert read["problems"] == [] and read["incomplete_lines"] == 1 and len(read["events"]) == 2
    third = pg.journal(ctx, "generate_begin", workers=0)
    assert third["seq"] == 3 and "after_incomplete_line" not in third and len(pg.journal_events(ctx)) == 3
    big = pg.journal(ctx, "generate_end", unresolved=["x" * 40] * 5000)   # one line longer than the block read from the end
    after = pg.journal(ctx, "generate_begin", workers=0)
    assert after["seq"] == big["seq"] + 1 and len(pg.journal_events(ctx)) == 5 and IG.read_journal(tmp_path)["problems"] == []


def test_only_one_generate_works_on_an_output_tree(tmp_path):
    import fcntl
    ctx = light(tmp_path)
    tmp_path.mkdir(exist_ok=True)
    holder = os.open(tmp_path / ".generate.lock", os.O_RDWR | os.O_CREAT, 0o644)
    try:
        fcntl.flock(holder, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(pc.PilotError, match="another generate"):
            pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
        assert not ctx.result("journal.jsonl").exists()
    finally:
        os.close(holder)
    assert pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)["counts"]["FAILED_RECORD"] == 10


def test_a_stop_request_finishes_the_run_in_progress_and_costs_no_attempt(tmp_path):
    ctx = light(tmp_path)
    fx.set_interruptions(tmp_path, {"healthy/1": ["SIGINT"]})
    before = signal.getsignal(signal.SIGINT)
    out = pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert signal.getsignal(signal.SIGINT) is before                  # the usual handler is back
    assert out["stopped_on_request"] is True and out["counts"] == {"COMPLETE": 0, "FAILED_RECORD": 4, "MISSING": 6,
                                                                   "PARTIAL": 0, "CORRUPT": 0}
    assert [e["event"] for e in events_of(ctx, "healthy", 1)] == ["start", "commit"]      # the run in progress was finished
    events = pg.journal_events(ctx)
    stop = next(e for e in events if e["event"] == "stop_requested")
    assert stop["signal"] == int(signal.SIGINT) and stop["runs_not_started"] == 6
    assert not any(e["event"] == "interrupted" for e in events)
    out = pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert out["stopped_on_request"] is False and out["counts"]["FAILED_RECORD"] == 10
    assert all(e["attempts"] == 1 for e in pc.read_json(ctx.result("run_manifest.json"))["runs"])
    assert attempt_log(ctx)["agrees"]


def test_the_independent_check_of_the_attempt_log_finds_what_was_changed(tmp_path):
    ctx = light(tmp_path)
    fx.set_interruptions(tmp_path, {"healthy/0": ["OSError"]})
    with pytest.raises(pg.InfrastructureInterruption):
        pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    pg.generate(ctx, engine_spec=fx.INTERRUPTED_FAILING_SPEC, workers=0)
    assert attempt_log(ctx)["agrees"]

    def fields(mutate):
        """Discrepancy fields after a change to a copy of the tree."""
        import shutil
        copy = tmp_path.parent / (tmp_path.name + "-copy")
        shutil.rmtree(copy, ignore_errors=True)
        shutil.copytree(tmp_path, copy)
        mutate(fx.context(ctx.plan, copy, need_measure=False))
        out = VERIFY.check_attempt_log(copy, fx.NAMESPACE, fx.audit_plan(ctx.plan))
        return {(d["id"], d["field"]) for d in out["discrepancies"]}

    def rewrite_record(c, condition, replica, change):
        record = pc.read_json(c.record_path(condition, replica))
        change(record)
        c.record_path(condition, replica).write_bytes(pc.canonical_json(record))

    def other_seed(record):
        record["seeds"]["local"] += 1

    assert ("fault/DUP/00001", "commit record_sha256") in fields(lambda c: rewrite_record(c, "fault/DUP", 1, lambda r: r["failure"].update(local_bytes=8)))
    found = fields(lambda c: rewrite_record(c, "fault/DUP", 1, other_seed))
    assert ("fault/DUP/00001", "record seeds") in found and ("fault/DUP/00001", "commit record_sha256") in found
    found = fields(lambda c: c.record_path("fault/DUP", 1).unlink())
    assert ("fault/DUP/00001", "commit entries of a run without a record") in found
    # A run started again although it had its record: a result replaced by another draw.
    found = fields(lambda c: log_start(c, "fault/SHIFT", 0, 2))
    assert ("fault/SHIFT/00000", "order") in found and ("fault/SHIFT/00000", "run manifest attempts") in found
    # An attempt under other seeds, and one whose random input is not the planned one.
    found = fields(lambda c: pg.journal(c, "start", condition="fault/SHIFT", replica=0, attempt=2,
                                        seeds={**c.dlib.seeds_for(fx.NAMESPACE, "fault/SHIFT", 0), "profile": 5}, rng_input_sha256="0" * 64))
    assert {("fault/SHIFT/00000", "attempt 2 seeds"), ("fault/SHIFT/00000", "attempt 2 rng_input_sha256")} <= found
    # An interruption for a reason that is not of the infrastructure, and one that names no started attempt.
    found = fields(lambda c: pg.journal(c, "interrupted", condition="healthy", replica=1, attempt=1, reason={"kind": "bad statistic"}))
    assert ("healthy/00001", "interruption reason") in found and ("healthy/00001", "order") in found
    found = fields(lambda c: pg.journal(c, "start", condition="healthy", replica=9, attempt=1))
    assert ("healthy/9", "planned") in found
    # The run manifest must tell the same history.
    def fewer_attempts(c):
        manifest = pc.read_json(c.result("run_manifest.json"))
        next(e for e in manifest["runs"] if (e["condition"], e["replica"]) == ("healthy", 0))["attempts"] = 1
        c.result("run_manifest.json").write_bytes(pc.canonical_json(manifest))

    assert fields(fewer_attempts) == {("healthy/00000", "run manifest attempts")}
    assert fields(lambda c: c.result("journal.jsonl").unlink()) == set()        # a test namespace needs no log
    absent = VERIFY.check_attempt_log(tmp_path / "nowhere", pc.SMOKE_NAMESPACE, fx.audit_plan(pc.SMOKE_PLAN))
    assert absent["status"] == "ABSENT" and not absent["agrees"] and absent["discrepancies"][0]["field"] == "present"


def damage(ctx, kind):
    record_path, trace_path = ctx.record_path("fault/SHIFT", 0), ctx.trace_path("fault/SHIFT", 0)
    record = pc.read_json(record_path)
    if kind == "truncated":
        record_path.write_bytes(record_path.read_bytes()[:-9])
    elif kind == "empty":
        record_path.write_bytes(b"")
    elif kind == "list":
        record_path.write_text("[]\n")
    elif kind == "other replica":
        record_path.write_bytes(pc.canonical_json({**record, "replica": 1}))
    elif kind == "other namespace":
        record_path.write_bytes(pc.canonical_json({**record, "namespace": pc.SMOKE_NAMESPACE}))
    elif kind == "claims a trace":
        record_path.write_bytes(pc.canonical_json({**record, "trace_path": "traces/fault/SHIFT/00000.npz", "trace_sha256": "0" * 64}))
    elif kind == "success without trace":
        record_path.write_bytes(pc.canonical_json({**record, "status": "SUCCESS", "failure": None}))
    elif kind == "trace beside failure record":
        trace_path.parent.mkdir(parents=True, exist_ok=True)
        trace_path.write_bytes(b"not a trace")
    return record_path, trace_path


@pytest.mark.parametrize("kind", ["truncated", "empty", "list", "other replica", "other namespace", "claims a trace",
                                  "success without trace", "trace beside failure record"])
def test_a_store_that_does_not_verify_stops_everything_and_is_left_as_found(tmp_path, kind):
    ctx = light(tmp_path)
    pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    record_path, trace_path = damage(ctx, kind)
    before = (record_path.read_bytes(), trace_path.read_bytes() if trace_path.exists() else None)
    state = pg.classify(ctx, "fault/SHIFT", 0)
    assert state["state"] == "CORRUPT" and state["reason"]
    with pytest.raises(pc.PilotError, match="do not verify and are not repaired"):
        pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    with pytest.raises(pc.PilotError, match="no valid record"):
        pg.require_complete(ctx)
    assert (record_path.read_bytes(), trace_path.read_bytes() if trace_path.exists() else None) == before
    assert not (tmp_path / "interrupted").exists()


def test_a_record_changed_after_generation_stops_the_analysis(tmp_path):
    ctx = light(tmp_path)
    pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    path = ctx.record_path("healthy", 0)
    record = pc.read_json(path)
    record["failure"]["local_bytes"] += 8          # still a well-formed failure record
    path.write_bytes(pc.canonical_json(record))
    assert pg.classify(ctx, "healthy", 0)["state"] == "FAILED_RECORD"
    with pytest.raises(pc.PilotError, match="differ from the run manifest"):
        pg.require_complete(ctx)
    ctx.result("run_manifest.json").unlink()
    with pytest.raises(pc.PilotError, match="not been finalised"):
        pg.require_complete(ctx)


def test_a_trace_without_a_record_is_set_aside_not_deleted(tmp_path):
    ctx = light(tmp_path)
    orphan = ctx.trace_path("fault/DUP", 0)
    orphan.parent.mkdir(parents=True)
    orphan.write_bytes(b"half a trace")
    assert pg.classify(ctx, "fault/DUP", 0)["state"] == "PARTIAL"
    with pytest.raises(pc.PilotError, match="origin is unknown"):     # no attempt was ever logged for it
        pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    assert orphan.read_bytes() == b"half a trace" and not (tmp_path / "interrupted").exists()
    assert not list((tmp_path / "results" / "runs").rglob("*.json"))
    log_start(ctx, "fault/DUP", 0)                                    # an attempt that died after writing its trace
    pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    kept = list((tmp_path / "interrupted").rglob("00000.npz.*"))
    assert len(kept) == 1 and kept[0].read_bytes() == b"half a trace" and not orphan.exists()
    events = events_of(ctx, "fault/DUP", 0)
    assert [e["event"] for e in events] == ["start", "interrupted", "restart", "start", "commit"]
    assert events[1]["reason"] == {"kind": "process_ended_without_record", "detail": "a trace without a record"}
    assert events[1]["trace_set_aside"].startswith("interrupted/traces/fault/DUP/")
    assert pc.read_json(ctx.record_path("fault/DUP", 0))["failure"]["phase"] == "generation"     # the second attempt's own result
    assert (tmp_path / ".staging").exists() is False and attempt_log(ctx)["agrees"]


def test_shards_partition_the_pending_runs():
    pending = list(pc.FROZEN_PLAN.runs())[:997]
    for workers in (1, 6, 8):
        parts = [pg.shard(pending, f"{k}/{workers}") for k in range(workers)]
        assert sorted(sum(parts, [])) == sorted(pending) and sum(map(len, parts)) == len(pending)
        assert max(map(len, parts)) - min(map(len, parts)) <= 1


def test_reserved_namespaces_and_engines_are_refused(tmp_path):
    reserved = dataclasses.replace(fx.SMALL_PLAN, namespace=pc.PILOT_NAMESPACE)
    with pytest.raises(pc.PilotError, match="must not use the pilot or smoke namespace"):
        pg.generate(light(tmp_path, reserved), engine_spec=fx.FAILING_SPEC, workers=0)
    smoke = pc.Context(pc.SMOKE_PLAN, tmp_path, None, au.dlib(), None)
    with pytest.raises(pc.PilotError, match="only the production engine"):
        pg.load_engine(fx.FAILING_SPEC, smoke)
    with pytest.raises(pc.PilotError, match="manifest.json is missing"):
        pg.generate(smoke, engine_spec="production", workers=0)
    assert not list(tmp_path.iterdir())            # nothing was created


# --------------------------------------------------------------------------- manifest and authorization
def test_pre_generation_manifest_contents():
    ctx = pc.Context(pc.FROZEN_PLAN, _Path("/nonexistent"), None, au.dlib(), None)
    for other in (None, "retry", "record_failure"):
        with pytest.raises(pc.PilotError, match="is not implemented"):
            pg.build_pre_generation_manifest(ctx, interrupted_run_policy=other, tdk_deadline=None)
    manifest = pg.build_pre_generation_manifest(ctx, tdk_deadline=None, workers=6)
    assert manifest["interrupted_run_policy"] == "regenerate_same_id" and manifest["interrupted_run_rules"] == pg.INTERRUPTED_RUN_RULES
    assert manifest["interrupted_run_rules"]["max_infrastructure_retries"] == 3
    assert manifest["environment_sha256"] == pg.environment_sha256()
    for index in (0, 998, 999, 2999, 4006, 5998):        # the hash of every run's random input, as the independent code has it
        run = manifest["planned_runs"][index]
        assert run["rng_input_sha256"] == RC.rng_input_sha256(RC.run_seeds(pc.PILOT_NAMESPACE, run["condition"], run["replica"]))
    assert len({run["rng_input_sha256"] for run in manifest["planned_runs"]}) == 5999
    assert PF.planned_run_problems(manifest, pc.PILOT_NAMESPACE, PF.PLAN) == []
    assert manifest["protocol_sha256"] == (au.REPO / "research/v019-d-pilot/protocol.sha256").read_text().split()[0]
    assert len(manifest["planned_runs"]) == 5999 and manifest["tdk_deadline"] is None
    run = manifest["planned_runs"][999 + 2000 + 1000 + 7]
    assert (run["condition"], run["replica"], run["split"]) == ("fault/LAG1", 7, "fitting")
    assert run["seeds"] == ctx.dlib.seeds_for(pc.PILOT_NAMESPACE, "fault/LAG1", 7)
    assert run["seed_labels"]["profile"] == "rnglab-v019-d-pilot/fault/LAG1/profile/00007"
    seeds = [s for r in manifest["planned_runs"] for s in r["seeds"].values()]
    assert len(set(seeds)) == len(seeds) == 5 * 5999          # no seed is used twice
    verify = au.load("v019_audit_verify", au.REPO / "research/v019-audit/verify.py")
    assert manifest["i4_replay_selection"] == verify.deep_sample(pc.PILOT_NAMESPACE)   # same 25 runs as the independent audit
    assert sum(map(len, manifest["i4_replay_selection"].values())) == 25
    for path in ("research/v019-d-pilot/dlib.py", "research/v019-d-pilot/pilot_generate.py", "research/v019-d-pilot/protocol.json",
                 "research/v018-h1-detection/h1lib.py", "src/rnglab/seeds.py", "src/rnglab14/models/pimc_harmonic.py",
                 "src/rnglab12/uncertainty.py"):
        assert manifest["code_sha256"][path] == pc.sha256_file(au.REPO / path)
    assert manifest["workers"] == 6 and manifest["workers_deviate_from_protocol"] is True
    assert set(manifest["thread_variables"]) == set(pc.THREAD_VARIABLES)
    with pytest.raises(pc.PilotError, match="test plan"):
        pg.build_pre_generation_manifest(light(_Path("/nonexistent")), tdk_deadline=None)


def test_generation_is_refused_without_the_authors_authorization(tmp_path, monkeypatch):
    for name in pc.THREAD_VARIABLES:
        monkeypatch.setenv(name, "1")
    ctx = pc.Context(pc.SMOKE_PLAN, tmp_path, None, au.dlib(), None)
    with pytest.raises(pc.PilotError, match="manifest.json is missing"):
        pg.require_generation_allowed(ctx)
    manifest = pg.build_pre_generation_manifest(ctx, tdk_deadline="2027-04-30")
    assert len(manifest["planned_runs"]) == 10 and PF.planned_run_problems(manifest, pc.SMOKE_NAMESPACE, PF.NAMESPACES[pc.SMOKE_NAMESPACE][0]) == []
    manifest_path, authorization_path = tmp_path / "pre_generation_manifest.json", tmp_path / "authorization.json"
    manifest_path.write_bytes(pc.canonical_json(manifest))
    with pytest.raises(pc.PilotError, match="authorization.json is missing"):
        pg.require_generation_allowed(ctx)
    good = {"authorized_by": "the author", "authorized_on": "2026-12-01", "namespace": pc.SMOKE_NAMESPACE,
            "pre_generation_manifest_sha256": pc.sha256_file(manifest_path)}
    for broken in ({**good, "pre_generation_manifest_sha256": "0" * 64}, {**good, "namespace": pc.PILOT_NAMESPACE},
                   {**good, "authorized_by": ""}, {k: v for k, v in good.items() if k != "authorized_on"}):
        authorization_path.write_bytes(pc.canonical_json(broken))
        with pytest.raises(pc.PilotError, match="does not authorise"):
            pg.require_generation_allowed(ctx)
    authorization_path.write_bytes(pc.canonical_json(good))
    assert pg.require_generation_allowed(ctx)["interrupted_run_policy"] == "regenerate_same_id"
    assert pg.planned_inputs(pg.require_generation_allowed(ctx))[("fault/DUP", 1)]["rng_input_sha256"] == RC.rng_input_sha256(
        RC.run_seeds(pc.SMOKE_NAMESPACE, "fault/DUP", 1))
    monkeypatch.setenv("OMP_NUM_THREADS", "4")
    with pytest.raises(pc.PilotError, match="OMP_NUM_THREADS must be 1"):
        pg.require_generation_allowed(ctx)
    monkeypatch.setenv("OMP_NUM_THREADS", "1")

    def reauthorised(changed):
        manifest_path.write_bytes(pc.canonical_json(changed))
        authorization_path.write_bytes(pc.canonical_json({**good, "pre_generation_manifest_sha256": pc.sha256_file(manifest_path)}))

    reauthorised({**manifest, "code_sha256": {**manifest["code_sha256"], "research/v019-d-pilot/dlib.py": "0" * 64}})
    with pytest.raises(pc.PilotError, match="code differs from the pre-generation manifest"):
        pg.require_generation_allowed(ctx)       # the code changed after the manifest was written
    reauthorised({**manifest, "code_sha256": {k: v for k, v in manifest["code_sha256"].items() if not k.endswith("pilot_common.py")}})
    with pytest.raises(pc.PilotError, match="not listed"):
        pg.require_generation_allowed(ctx)
    reauthorised({**manifest, "namespace": pc.PILOT_NAMESPACE})
    with pytest.raises(pc.PilotError, match="another namespace or protocol"):
        pg.require_generation_allowed(ctx)
    for changed in ({"interrupted_run_policy": None}, {"interrupted_run_policy": "record_failure"},
                    {"interrupted_run_rules": {**manifest["interrupted_run_rules"], "max_infrastructure_retries": 10}},
                    {"interrupted_run_rules": None}):
        reauthorised({**manifest, **changed})
        with pytest.raises(pc.PilotError, match="does not carry the interrupted-run policy the author decided"):
            pg.require_generation_allowed(ctx)
    assert sorted(p.name for p in tmp_path.iterdir()) == ["authorization.json", "pre_generation_manifest.json"]  # no draw happened


def test_amendments_need_a_recorded_approval_and_an_unaltered_document(tmp_path):
    assert pc.approved_amendments(tmp_path) == {} and pc.approved_amendments() == {}      # none exists today
    document = tmp_path / "AMENDMENT_01.md"
    document.write_text("text of the amendment\n")
    entry = {"document": "AMENDMENT_01.md", "sha256": pc.sha256_file(document), "approved_by": "the author", "approved_on": "2026-12-01"}
    (tmp_path / "amendments.json").write_bytes(pc.canonical_json({"A01": entry}))
    assert pc.approved_amendments(tmp_path) == {"A01": entry}
    (tmp_path / "amendments.json").write_bytes(pc.canonical_json({"A01": {**entry, "approved_by": ""}}))
    with pytest.raises(pc.PilotError, match="approval is not recorded"):
        pc.approved_amendments(tmp_path)
    (tmp_path / "amendments.json").write_bytes(pc.canonical_json({"A01": entry}))
    document.write_text("text of the amendment, edited\n")
    with pytest.raises(pc.PilotError, match="missing or altered"):
        pc.approved_amendments(tmp_path)


# --------------------------------------------------------------------------- engine pieces that need no measurement layer
def test_step_check_agrees_with_stage_c_and_detects_a_wrong_stream():
    dlib = au.dlib()
    for family in ("SHIFT", "LAG1", "DUP"):
        cap = np.float64(dlib.CAPS[family])
        profile = {"t0": 3, "T_cap": np.float64(7.5), "r": np.float64(cap / np.float64(7.5)), "s_max": cap}
        run = dlib.generate_run(fx.NAMESPACE, "fault/" + family, family, 4, sweeps=64, profile=profile)
        assert pg.step_check(dlib, run) == {"check": "I-3-step", "passed": True, "mismatches": [], "bytes": run["local_bytes"]}
        wrong = bytearray(run["raw_base"])
        for word in range(1, 65):          # the low bits of 64 words: a duplicated word hides its own raw value
            wrong[-8 * word + 7] ^= 1
        assert pg.step_check(dlib, {**run, "raw_base": bytes(wrong)})["mismatches"] == [{"field": "independent"}]
    healthy = dlib.generate_run(fx.NAMESPACE, "healthy", "healthy", 4, sweeps=64)
    assert pg.step_check(dlib, {**healthy, "family": "DUP", "profile": profile})["passed"] is False   # no onset: no B0


def test_a_sampler_exception_becomes_a_failure_record_with_the_progress_made(tmp_path, monkeypatch):
    dlib = au.dlib()

    class Broken:
        def __init__(self, **_):
            pass

        def run_simulation(self, local, **options):
            for _ in range(20):
                local.read_bytes(8)
            options["centroid_source"].read_bytes(16)        # one sweep boundary
            local.read_bytes(8)
            raise ArithmeticError("sampler broke in sweep 1")

    monkeypatch.setattr(dlib, "PathIntegralHarmonicOscillator", Broken)
    ctx = pc.Context(fx.SMALL_PLAN, tmp_path, None, dlib, None)
    record, arrays = pg.PilotEngine(ctx).execute("fault/SHIFT", "SHIFT", 1)
    assert arrays is None and record["status"] == "FAILED"
    assert record["failure"] == {"phase": "generation", "exception": {"class": "ArithmeticError", "message": "sampler broke in sweep 1"},
                                 "completed_sweeps": 1, "last_completed_sweep": 0, "local_bytes": 168}
    pg.commit_run(ctx, record, None)
    assert pg.classify(ctx, "fault/SHIFT", 1)["state"] == "FAILED_RECORD"


def test_command_line_stops_without_the_measurement_layer(capsys, tmp_path, monkeypatch):
    run_pilot = au.load("d_pilot_run_pilot_test", au.REPO / "research/v019-d-pilot/run_pilot.py")
    with pytest.raises(SystemExit):
        run_pilot.main(["calibrate", "--smoke"])     # no analysis of a smoke
    monkeypatch.setattr(run_pilot.pc, "HERE", tmp_path)          # a directory that does not hold the file
    assert run_pilot.main(["status"]) == 2
    assert "measure.py is not in this tree" in capsys.readouterr().err


def test_manifest_preview_shows_the_decision_and_writes_nothing(capsys):
    if au.measure_path() != au.REPO / "research/v019-d-pilot/measure.py":
        pytest.skip("the command line needs the measurement layer in the tree")
    run_pilot = au.load("d_pilot_run_pilot_test", au.REPO / "research/v019-d-pilot/run_pilot.py")
    directory = au.REPO / "research/v019-d-pilot"
    before = sorted(p.relative_to(directory).as_posix() for p in directory.rglob("*") if "__pycache__" not in p.parts)
    assert run_pilot.main(["manifest", "--smoke", "--preview", "--tdk-deadline", "2027-04-30"]) == 0
    out = capsys.readouterr().out
    assert "PREVIEW: nothing was written." in out and '"interrupted_run_policy": "regenerate_same_id"' in out
    assert '"max_infrastructure_retries": 3' in out and "10 runs, each with seed labels, seeds and rng_input_sha256" in out
    assert '"namespace": "rnglab-v019-d-pilot-smoke"' in out and '"tdk_deadline": "2027-04-30"' in out
    assert sorted(p.relative_to(directory).as_posix() for p in directory.rglob("*") if "__pycache__" not in p.parts) == before
    with pytest.raises(SystemExit):
        run_pilot.main(["manifest", "--interrupted-run-policy", "record_failure"])      # the choice is no longer offered


# --------------------------------------------------------------------------- with the measurement layer
def test_synthetic_runs_are_stored_whole_and_load_back(tmp_path):
    ctx = fx.context(fx.SMALL_PLAN, tmp_path)
    out = pg.generate(ctx, engine_spec=fx.ENGINE_SPEC, workers=0)
    assert out["counts"]["COMPLETE"] == 10 and not (tmp_path / ".staging").exists()
    engine = fx.SyntheticEngine(ctx)
    for condition, family, replica in ctx.plan.runs():
        record, arrays = ctx.measure.load_run(tmp_path, condition, replica)   # verifies file and array hashes
        _, expected = engine.execute(condition, family, replica)
        assert all(arrays[k].tobytes() == expected[k].tobytes() and arrays[k].dtype == expected[k].dtype for k in expected)
        assert record["trace_sha256"] == pc.sha256_file(ctx.trace_path(condition, replica))
        assert ctx.record_path(condition, replica).read_bytes() == pc.canonical_json(record)
    replay = pg.replay_check(ctx, engine, {"healthy": [1], "fault/LAG1": [0, 1]})
    assert replay["passed"] and len(replay["runs"]) == 3 and replay["ignored_fields"] == list(pg.TIMING_FIELDS)
    path = ctx.trace_path("fault/LAG1", 1)
    data = bytearray(path.read_bytes())
    data[len(data) // 2] ^= 1
    path.write_bytes(bytes(data))
    assert pg.classify(ctx, "fault/LAG1", 1) == {**pg.classify(ctx, "fault/LAG1", 1), "state": "CORRUPT",
                                                 "reason": "trace hash differs from the record"}
    with pytest.raises(ValueError, match="hash mismatch"):
        ctx.measure.load_run(tmp_path, "fault/LAG1", 1)
    assert pg.replay_check(ctx, engine, {"fault/LAG1": [1]})["passed"] is False


def test_replay_detects_a_changed_value(tmp_path):
    ctx = fx.context(fx.SMALL_PLAN, tmp_path)
    pg.generate(ctx, engine_spec=fx.ENGINE_SPEC, workers=0)

    class Drifting(fx.SyntheticEngine):
        def execute(self, condition, family, replica):
            record, arrays = super().execute(condition, family, replica)
            arrays["TD"][40959] = np.nextafter(arrays["TD"][40959], 9.0)        # one unit in the last place
            record["conventional"][3]["S"] += 1e-12
            record["measurement_seconds"] = 123.0                                # timing only: ignored
            return record, arrays

    result = pg.replay_check(ctx, Drifting(ctx), {"fault/DUP": [0]})
    assert result["passed"] is False
    assert {m["field"] for m in result["runs"][0]["mismatches"]} == {"array:TD", "conventional"}
    ctx.record_path("healthy", 0).unlink()
    ctx.trace_path("healthy", 0).unlink()
    assert pg.replay_check(ctx, fx.SyntheticEngine(ctx), {"healthy": [0]})["runs"][0]["reason"] == "no complete stored run to compare with"


def test_a_non_finite_trace_is_refused_by_the_storage_layer_and_recorded(tmp_path):
    ctx = fx.context(fx.SMALL_PLAN, tmp_path)
    out = pg.generate(ctx, engine_spec=fx.NON_FINITE_SPEC, workers=0)
    assert out["counts"]["COMPLETE"] == 9 and out["counts"]["FAILED_RECORD"] == 1
    record = pc.read_json(ctx.record_path("healthy", 0))
    assert record["failure"]["phase"] == "serialization" and "nonfinite" in record["failure"]["exception"]["message"]
    assert not ctx.trace_path("healthy", 0).exists() and not (tmp_path / ".staging").exists()


WORKER = r"""
import sys
sys.path.insert(0, sys.argv[1])
import _pilot_fixture as fx
plan = fx.pc.Plan(fx.NAMESPACE, (("calibration", "healthy", 2), ("healthy", "healthy", 2), ("fault/SHIFT", "SHIFT", 8),
                                 ("fault/LAG1", "LAG1", 8), ("fault/DUP", "DUP", 8)), range(0, 4), range(4, 8), 2, 2, "test")
ctx = fx.context(plan, sys.argv[2])
fx.pg.generate(ctx, engine_spec=fx.ENGINE_SPEC, workers=0)
print("done")
"""


def test_killed_generation_resumes_to_the_same_store(tmp_path):
    au.measure()                                    # skip without the measurement layer
    script = tmp_path / "worker.py"
    script.write_text(WORKER)
    environment = {**os.environ, "PYTHONDONTWRITEBYTECODE": "1", "PYTHONPATH": str(au.REPO / "src")}
    command = [_sys.executable, str(script), str(_Path(__file__).parent), str(tmp_path / "store")]
    plan = pc.Plan(fx.NAMESPACE, (("calibration", "healthy", 2), ("healthy", "healthy", 2), ("fault/SHIFT", "SHIFT", 8),
                                  ("fault/LAG1", "LAG1", 8), ("fault/DUP", "DUP", 8)), range(0, 4), range(4, 8), 2, 2, "test")
    ctx = fx.context(plan, tmp_path / "store")
    kills = 0
    for attempt in range(60):
        process = subprocess.Popen(command, env=environment, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            process.wait(timeout=2.0 + 0.37 * (attempt % 5))
        except subprocess.TimeoutExpired:
            process.send_signal(signal.SIGKILL)
            process.wait()
            kills += 1
        assert process.returncode in (0, -signal.SIGKILL), process.stderr.read().decode()[-2000:]
        states = pg.scan(ctx)
        assert not [s for s in states if s["state"] == "CORRUPT"]          # never a half-written run that claims to be whole
        if all(s["state"] == "COMPLETE" for s in states):
            break
    assert all(s["state"] == "COMPLETE" for s in states) and kills >= 1, (kills, [s["state"] for s in states])
    engine = fx.SyntheticEngine(ctx)
    for condition, family, replica in plan.runs():
        _, arrays = ctx.measure.load_run(ctx.root, condition, replica)
        _, expected = engine.execute(condition, family, replica)
        assert all(arrays[k].tobytes() == expected[k].tobytes() for k in expected)
    events = pg.journal_events(ctx)
    commits = [(e["condition"], e["replica"]) for e in events if e["event"] == "commit"]
    assert len(commits) == len(set(commits)) == plan.total()               # every run committed exactly once
    assert sum(e["event"] == "restart" for e in events) >= 1
    assert pg.require_complete(ctx)["complete"] == plan.total()
    # Every restart ran the same identity again, and the independent check reads the same history.
    starts = [e for e in events if e["event"] == "start"]
    assert len(starts) > plan.total() and max(e["attempt"] for e in starts) <= pg.MAX_ATTEMPTS
    for entry in starts:
        seeds = RC.run_seeds(fx.NAMESPACE, entry["condition"], entry["replica"])
        assert entry["seeds"] == seeds and entry["rng_input_sha256"] == RC.rng_input_sha256(seeds)
    check = VERIFY.check_attempt_log(ctx.root, fx.NAMESPACE, fx.audit_plan(plan))
    assert check["agrees"], check["discrepancies"][:5]
    assert check["interruptions"] == len(starts) - plan.total() and check["unresolved"] == []
