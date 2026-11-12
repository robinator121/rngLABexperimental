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


@pytest.mark.parametrize("policy", [None, "record_failure", "regenerate_same_id"])
def test_a_run_that_was_started_and_left_no_record_needs_the_authors_policy(tmp_path, policy):
    ctx = light(tmp_path)
    pg.journal(ctx, "start", condition="healthy", replica=1)       # the process died before any record
    if policy is None:
        with pytest.raises(pc.PilotError, match="interrupted-run policy"):
            pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
        assert not list((tmp_path / "results" / "runs").rglob("*.json")) if (tmp_path / "results" / "runs").exists() else True
        return
    pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0, test_policy=policy)
    phase = pc.read_json(ctx.record_path("healthy", 1))["failure"]["phase"]
    events = [e["event"] for e in pg.journal_events(ctx) if e.get("condition") == "healthy" and e.get("replica") == 1]
    if policy == "record_failure":
        assert phase == "interrupted" and events == ["start", "interrupted"]     # never executed again
    else:
        assert phase == "generation" and events == ["start", "interrupted", "restart", "start", "commit"]
        entry = next(e for e in pc.read_json(ctx.result("run_manifest.json"))["runs"]
                     if (e["condition"], e["replica"]) == ("healthy", 1))
        assert entry["attempts"] == 2
    assert pc.read_json(ctx.record_path("healthy", 0))["failure"]["phase"] == "generation"   # others are unaffected


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
        pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0, test_policy="regenerate_same_id")
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
    with pytest.raises(pc.PilotError, match="interrupted-run policy"):
        pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0)
    assert orphan.read_bytes() == b"half a trace"
    pg.generate(ctx, engine_spec=fx.FAILING_SPEC, workers=0, test_policy="record_failure")
    kept = list((tmp_path / "interrupted").rglob("00000.npz.*"))
    assert len(kept) == 1 and kept[0].read_bytes() == b"half a trace" and not orphan.exists()
    record = pc.read_json(ctx.record_path("fault/DUP", 0))
    assert record["failure"]["phase"] == "interrupted" and record["failure"]["trace_set_aside"].startswith("interrupted/traces/fault/DUP/")
    assert (tmp_path / ".staging").exists() is False


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
    with pytest.raises(pc.PilotError, match="no default"):
        pg.build_pre_generation_manifest(ctx, interrupted_run_policy=None, tdk_deadline=None)
    with pytest.raises(pc.PilotError, match="no default"):
        pg.build_pre_generation_manifest(ctx, interrupted_run_policy="retry", tdk_deadline=None)
    manifest = pg.build_pre_generation_manifest(ctx, interrupted_run_policy="record_failure", tdk_deadline=None, workers=6)
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
        pg.build_pre_generation_manifest(light(_Path("/nonexistent")), interrupted_run_policy="record_failure", tdk_deadline=None)


def test_generation_is_refused_without_the_authors_authorization(tmp_path, monkeypatch):
    for name in pc.THREAD_VARIABLES:
        monkeypatch.setenv(name, "1")
    ctx = pc.Context(pc.SMOKE_PLAN, tmp_path, None, au.dlib(), None)
    with pytest.raises(pc.PilotError, match="manifest.json is missing"):
        pg.require_generation_allowed(ctx)
    manifest = pg.build_pre_generation_manifest(ctx, interrupted_run_policy="record_failure", tdk_deadline="2027-04-30")
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
    assert pg.require_generation_allowed(ctx)["interrupted_run_policy"] == "record_failure"
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
    reauthorised({**manifest, "interrupted_run_policy": None})
    with pytest.raises(pc.PilotError, match="no interrupted-run policy"):
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


def test_command_line_stops_without_the_measurement_layer(capsys):
    if au.measure_path() is not None and au.measure_path().parent == au.REPO / "research/v019-d-pilot":
        pytest.skip("the measurement layer is in the tree")
    run_pilot = au.load("d_pilot_run_pilot_test", au.REPO / "research/v019-d-pilot/run_pilot.py")
    assert run_pilot.main(["status"]) == 2
    assert "measure.py is not in this tree" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        run_pilot.main(["calibrate", "--smoke"])     # no analysis of a smoke


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
fx.pg.generate(ctx, engine_spec=fx.ENGINE_SPEC, workers=0, test_policy="regenerate_same_id")
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
