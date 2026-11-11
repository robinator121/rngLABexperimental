"""Campaign integrity (research/v019-audit/integrity.py) and preflight (preflight.py).

The interruption test kills a writer process with SIGKILL at arbitrary moments
and checks that no half-written run is ever classified as complete, and that
the finished store equals an uninterrupted one array for array.
"""
import importlib.util as _ilu
import json
import os
import signal
import subprocess
import sys as _sys
import time
from pathlib import Path as _Path

import numpy as np
import pytest

_spec = _ilu.spec_from_file_location("_audit_loader", _Path(__file__).with_name("_audit_loader.py"))
au = _sys.modules.get("_audit_loader") or _ilu.module_from_spec(_spec)
if "_audit_loader" not in _sys.modules:
    _sys.modules["_audit_loader"] = au
    _spec.loader.exec_module(au)

ig = au.load("v019_audit_integrity", au.REPO / "research/v019-audit/integrity.py")
pf = au.load("v019_audit_preflight", au.REPO / "research/v019-audit/preflight.py")
PLAN = {"healthy": 3, "fault/DUP": 2}


def toy(condition, replica, size=64):
    rng = np.random.default_rng([replica, len(condition)])
    arrays = {"CV": rng.normal(size=size), "marks": np.cumsum(rng.integers(16, 33, size=size) * 8).astype(np.int64)}
    return {"condition": condition, "replica": replica, "status": "SUCCESS", "value": float(arrays["CV"].sum())}, arrays


def fill(root, plan=PLAN):
    for condition, count in plan.items():
        for replica in range(count):
            ig.commit_run(root, *toy(condition, replica))


def test_commit_scan_and_manifest_round_trip(tmp_path):
    assert ig.scan_store(tmp_path, PLAN)["counts"] == {"COMPLETE": 0, "MISSING": 5, "PARTIAL": 0, "CORRUPT": 0}
    fill(tmp_path)
    scan = ig.scan_store(tmp_path, PLAN)
    assert scan["all_complete"] and scan["counts"]["COMPLETE"] == 5 and ig.pending(tmp_path, PLAN) == []
    record_path, trace_path = ig.run_paths(tmp_path, "fault/DUP", 1)
    record = json.loads(record_path.read_text())
    assert record["trace_sha256"] == ig.sha256_file(trace_path) and record["trace_path"] == "traces/fault/DUP/00001.npz"
    assert record_path.read_bytes() == ig.canonical_json(record)
    with np.load(trace_path, allow_pickle=False) as archive:
        expected = toy("fault/DUP", 1)[1]
        assert all(archive[k].tobytes() == expected[k].tobytes() and archive[k].dtype == expected[k].dtype for k in expected)
    with pytest.raises(FileExistsError):
        ig.commit_run(tmp_path, *toy("fault/DUP", 1))  # a completed draw is not replaced
    manifest = ig.build_manifest(tmp_path, PLAN, {"protocol_sha256": "x"})
    assert ig.verify_manifest(tmp_path, manifest) == {"ok": True, "problems": [], "files_checked": 10}
    assert all(e["state"] == "COMPLETE" and e["trace"]["sha256"] for e in manifest["runs"]) and len(manifest["runs"]) == 5


@pytest.mark.parametrize("damage,state", [
    ("drop_record", "PARTIAL"), ("leftover_tmp", "PARTIAL"), ("truncate_record", "CORRUPT"),
    ("empty_record", "CORRUPT"), ("wrong_identity", "CORRUPT"), ("drop_trace", "CORRUPT"),
    ("flip_trace_byte", "CORRUPT"), ("empty_trace", "CORRUPT"), ("record_is_list", "CORRUPT")])
def test_every_damaged_state_is_classified_and_never_complete(tmp_path, damage, state):
    fill(tmp_path)
    record_path, trace_path = ig.run_paths(tmp_path, "healthy", 1)
    if damage == "drop_record":
        record_path.unlink()
    elif damage == "leftover_tmp":
        record_path.unlink()
        trace_path.unlink()
        (trace_path.parent / "00001.npz.4242.tmp").write_bytes(b"half")
    elif damage == "truncate_record":
        record_path.write_bytes(record_path.read_bytes()[:-20])
    elif damage == "empty_record":
        record_path.write_bytes(b"")
    elif damage == "wrong_identity":
        record_path.write_bytes(ig.canonical_json({**json.loads(record_path.read_text()), "replica": 2}))
    elif damage == "drop_trace":
        trace_path.unlink()
    elif damage == "flip_trace_byte":
        data = bytearray(trace_path.read_bytes())
        data[len(data) // 2] ^= 1
        trace_path.write_bytes(bytes(data))
    elif damage == "empty_trace":
        trace_path.write_bytes(b"")
    elif damage == "record_is_list":
        record_path.write_text("[1, 2]\n")
    found = ig.classify_run(tmp_path, "healthy", 1)
    assert found["state"] == state and found["reason"]
    scan = ig.scan_store(tmp_path, PLAN)
    assert not scan["all_complete"] and scan["counts"][state] == 1 and scan["counts"]["COMPLETE"] == 4
    assert ig.pending(tmp_path, PLAN) == [("healthy", 1, state)] or damage == "flip_trace_byte"  # pending skips hashing


def test_manifest_detects_change_loss_addition_and_self_alteration(tmp_path):
    fill(tmp_path)
    manifest = ig.build_manifest(tmp_path, PLAN)
    record_path, trace_path = ig.run_paths(tmp_path, "healthy", 0)
    original = record_path.read_bytes()
    record_path.write_bytes(original.replace(b'"SUCCESS"', b'"SUCCES5"'))
    assert ig.verify_manifest(tmp_path, manifest)["problems"] == [{"kind": "CHANGED", "path": "results/runs/healthy/00000.json"}]
    record_path.write_bytes(original)
    trace_path.rename(trace_path.with_suffix(".bak"))
    kinds = {(p["kind"], p["path"]) for p in ig.verify_manifest(tmp_path, manifest)["problems"]}
    assert kinds == {("MISSING", "traces/healthy/00000.npz"), ("UNLISTED", "traces/healthy/00000.bak")}
    trace_path.with_suffix(".bak").rename(trace_path)
    assert ig.verify_manifest(tmp_path, manifest)["ok"]
    assert not ig.verify_manifest(tmp_path, {**manifest, "counts": {**manifest["counts"], "COMPLETE": 4}})["ok"]
    bigger = {"healthy": 4, "fault/DUP": 2}
    partial_manifest = ig.build_manifest(tmp_path, bigger)
    assert [e["state"] for e in partial_manifest["runs"]].count("MISSING") == 1  # missing runs are explicit entries
    ig.commit_run(tmp_path, *toy("healthy", 3))
    assert {p["kind"] for p in ig.verify_manifest(tmp_path, partial_manifest)["problems"]} == {"APPEARED"}
    extra = tmp_path / "traces/healthy/99999.npz"
    extra.write_bytes(b"x")
    assert ig.scan_store(tmp_path, bigger)["unplanned_files"] == ["traces/healthy/99999.npz"]


WRITER = r"""
import importlib.util, sys, time
from pathlib import Path
import numpy as np
spec = importlib.util.spec_from_file_location("ig", sys.argv[1]); ig = importlib.util.module_from_spec(spec); spec.loader.exec_module(ig)
root, count, size = Path(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4])
for condition, replica, state in ig.pending(root, {"healthy": count}):
    rng = np.random.default_rng([replica, 7])
    arrays = {"CV": rng.normal(size=size), "marks": np.cumsum(rng.integers(16, 33, size=size) * 8).astype(np.int64)}
    ig.append_journal(root, "start", condition=condition, replica=replica, found=state)
    for stale in list((root / "traces" / condition).glob(f"{replica:05d}.*")) if (root / "traces" / condition).is_dir() else []:
        stale.unlink()  # a restart clears what an interrupted attempt left behind; it never touches a complete run
    ig.commit_run(root, {"condition": condition, "replica": replica, "status": "SUCCESS"}, arrays)
    ig.append_journal(root, "commit", condition=condition, replica=replica)
print("done")
"""


def test_killed_writer_never_leaves_a_run_that_passes_as_complete(tmp_path):
    script = tmp_path / "writer.py"
    script.write_text(WRITER)
    count, size = 40, 60000
    command = [_sys.executable, str(script), str(au.REPO / "research/v019-audit/integrity.py"), str(tmp_path / "store"),
               str(count), str(size)]
    kills = 0
    for attempt in range(40):
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            process.wait(timeout=0.25 + 0.05 * (attempt % 5))
        except subprocess.TimeoutExpired:
            process.send_signal(signal.SIGKILL)
            process.wait()
            kills += 1
        scan = ig.scan_store(tmp_path / "store", {"healthy": count})
        assert scan["counts"]["CORRUPT"] == 0, scan["not_complete"]
        for run in scan["runs"]:
            if run["state"] == "COMPLETE":  # whatever is accepted is exactly the planned content
                with np.load(ig.run_paths(tmp_path / "store", "healthy", run["replica"])[1], allow_pickle=False) as archive:
                    expected = np.random.default_rng([run["replica"], 7]).normal(size=size)
                    assert archive["CV"].tobytes() == expected.tobytes()
        if scan["all_complete"]:
            break
    assert scan["all_complete"] and kills >= 1, (kills, scan["counts"])
    journal = [json.loads(line) for line in (tmp_path / "store/results/journal.jsonl").read_text().splitlines()]
    commits = [(e["condition"], e["replica"]) for e in journal if e["event"] == "commit"]
    assert len(set(commits)) == len(commits) <= count  # no run is committed twice
    assert sum(e["event"] == "start" for e in journal) >= len(commits)
    assert not list((tmp_path / "store").rglob("*.tmp"))


def test_namespace_collisions_are_reported():
    rc = au.recount()
    assert ig.namespace_collisions(rc.namespace_seed, {"a": {"x": 50}, "b": {"x": 50}}, ("local", "fault")) == []
    constant = lambda namespace, condition, role, replica: 7 if replica == 3 else rc.namespace_seed(namespace, condition, role, replica)
    clashes = ig.namespace_collisions(constant, {"a": {"x": 5}}, ("local", "fault"))
    assert clashes == [(("a", "x", "local", 3), ("a", "x", "fault", 3))]


# --------------------------------------------------------------------------- preflight
def test_preflight_checks_on_the_current_tree(tmp_path):
    assert pf.check_protocol()["status"] == "PASS"
    golden = pf.check_golden()
    assert golden["status"] == "PASS" and golden["independent_replay_matches"] and golden["pilot_code_matches"]
    assert pf.check_seeds()["status"] == "PASS" and pf.check_seeds()["seeds_checked"] == 30445
    assert pf.check_store(tmp_path)["status"] == "PASS"
    code = pf.check_code()
    assert set(pf.HASHED) <= set(code["sha256"]) and code["sha256"]["research/v019-d-pilot/dlib.py"]
    assert pf.check_manifest(tmp_path, code)["status"] == "FAIL"  # no pre-generation manifest yet
    environment = pf.check_environment()
    assert environment["versions"]["numpy"] == np.__version__ and "libc" in environment
    resources = pf.check_resources(tmp_path, 8, 5999)
    assert resources["disk_needed_bytes"] > 30 * 2**30 and resources["memory_needed_bytes"] == 8 * 300 * 2**20
    assert pf.check_resources(tmp_path, 10**6, 5999)["status"] == "FAIL"


def test_preflight_manifest_and_store_rules(tmp_path):
    code = {"sha256": {"a.py": "1" * 64, "b.py": "2" * 64}}
    path = tmp_path / "pre_generation_manifest.json"
    path.write_text(json.dumps({"code_sha256": code["sha256"], "tdk_deadline": "2027-04-30"}))
    assert pf.check_manifest(tmp_path, code)["status"] == "PASS"
    path.write_text(json.dumps({"code_sha256": code["sha256"]}))
    out = pf.check_manifest(tmp_path, code)
    assert out["status"] == "WARN" and out["criterion_5"] == "INSUFFICIENT_INFORMATION"
    path.write_text(json.dumps({"code_sha256": {**code["sha256"], "a.py": "3" * 64}, "tdk_deadline": "2027-04-30"}))
    assert pf.check_manifest(tmp_path, code)["code_hashes_differing"] == ["a.py"]
    path.write_text("{not json")
    assert pf.check_manifest(tmp_path, code)["status"] == "FAIL"
    ig.commit_run(tmp_path, {"condition": "healthy", "replica": 0, "status": "SUCCESS"}, {"CV": np.zeros(4)})
    assert pf.check_store(tmp_path)["status"] == "WARN"       # generation has started
    ig.run_paths(tmp_path, "healthy", 0)[0].unlink()
    assert pf.check_store(tmp_path)["status"] == "FAIL"       # a trace without its record
    with pytest.raises(SystemExit):
        pf.run(tmp_path, "rnglab-v019-d-pilot-smoke", 8)


def test_preflight_golden_detects_a_different_platform(tmp_path, monkeypatch):
    golden = json.loads(pf.GOLDEN_PATH.read_text())
    golden["digests"]["LAG1"]["CV"] = "0" * 64
    altered = tmp_path / "golden.json"
    altered.write_text(json.dumps(golden))
    monkeypatch.setattr(pf, "GOLDEN_PATH", altered)
    assert pf.check_golden()["status"] == "FAIL"
    monkeypatch.setattr(pf, "GOLDEN_PATH", tmp_path / "absent.json")
    assert pf.check_golden()["status"] == "FAIL"
