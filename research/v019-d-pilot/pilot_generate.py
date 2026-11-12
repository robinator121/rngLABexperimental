"""Stage D pilot: generation driver, atomic storage, resume, manifests, I-4 replay. protocol.json governs.

What this module does not decide. The protocol does not say what happens to a run that
was started and left no record because its process died. `generate` therefore refuses to
continue past such a run unless the pre-generation manifest names a policy
(`interrupted_run_policy`), which only the author can choose:

    regenerate_same_id   run the same identity again; the draw is the same function of its
                         seed labels, so the result is the one the first attempt would have given
    record_failure       write a failure record with phase 'interrupted' and never run it again

A run with a record is never run again, whatever the record says (protocol.json failed_runs).
A record or trace that does not verify is reported and stops the driver; nothing is repaired.

Generation in the pilot or smoke namespace needs the pre-generation manifest and the
author's authorization record that names that manifest's hash (`require_generation_allowed`).
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Any, Mapping

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))
import pilot_common as pc  # noqa: E402

STATES = ("COMPLETE", "FAILED_RECORD", "MISSING", "PARTIAL", "CORRUPT")
POLICIES = ("regenerate_same_id", "record_failure")
# Fields that differ between two executions of the same run for reasons of timing only (I-4 ignores them).
TIMING_FIELDS = ("measurement_seconds", "generation_seconds", "check_seconds")
# Fields the storage layer adds to a record; a regenerated record has them only after it is written.
STORAGE_FIELDS = ("trace_path", "trace_sha256", "array_hashes", "compression")


# ---------------------------------------------------------------------------
# One run
# ---------------------------------------------------------------------------
def step_check(dlib: Any, run: Mapping[str, Any]) -> dict:
    """I-3 step test (protocol.json implementation_checks.I-3): on the recorded base stream and marks,
    the pilot adapter in step mode, the independent transformer and the Stage C adapter with
    onset_bytes = B0 must emit the same bytes."""
    h1lib = dlib.h1lib
    family, profile, seeds, marks = run["family"], run["profile"], run["seeds"], run["marks"]
    local_bytes, onset = int(run["local_bytes"]), run["B0"]
    if onset is None:
        return {"check": "I-3-step", "passed": False, "mismatches": [{"field": "B0", "actual": None}]}
    local = dlib.FaultyWords(dlib.PCG64Source(seeds["local"]), family, profile, seeds["fault"], step=True)
    clock = dlib.SweepClock(dlib.PCG64Source(seeds["centroid"]), local)
    for boundary in marks.tolist():
        while local.consumed_bytes < boundary:
            local.read_bytes(8)
        clock.read_bytes(16)
    pilot = bytes(local.emitted)
    stage_c = h1lib.FaultyWords(dlib.PCG64Source(seeds["local"]), family, float(profile["s_max"]), seeds["fault"], 0,
                                onset_bytes=onset)
    reference = b"".join(stage_c.read_bytes(8) for _ in range(local_bytes // 8))
    independent = dlib.reapply(family, profile, seeds, marks, local_bytes, step=True, base_stream=run["raw_base"])
    mismatches = [{"field": name} for name, stream in (("pilot_adapter", pilot), ("independent", independent))
                  if stream != reference]
    return {"check": "I-3-step", "passed": not mismatches, "mismatches": mismatches, "bytes": local_bytes}


def _partial_state(error: BaseException) -> dict:
    """Sweeps and bytes completed when the sampler raised, read from the driver's frame."""
    tb = error.__traceback__
    while tb is not None:
        if tb.tb_frame.f_code.co_name == "generate_run":
            local, clock = tb.tb_frame.f_locals.get("local"), tb.tb_frame.f_locals.get("clock")
            if local is not None and clock is not None:
                done = len(clock.marks)
                return {"completed_sweeps": done, "last_completed_sweep": done - 1 if done else None,
                        "local_bytes": int(local.consumed_bytes)}
        tb = tb.tb_next
    return {"completed_sweeps": None, "last_completed_sweep": None, "local_bytes": None}


def failure_record(plan: pc.Plan, condition: str, family: str, replica: int, phase: str,
                   error: BaseException | None, seeds: Mapping[str, int] | None, extra: Mapping[str, Any] | None = None) -> dict:
    """A record for a planned run that produced no trace (protocol.json failed_runs.recording)."""
    failure = {"phase": phase, "exception": None if error is None else {"class": type(error).__name__, "message": str(error)},
               **(dict(extra) if extra else {})}
    return {"namespace": plan.namespace, "condition": condition, "family": family, "replica": replica,
            "split": plan.split(condition, replica), "seeds": None if seeds is None else dict(seeds),
            "status": "FAILED", "failure": failure, "failures": [failure],
            "trace_path": None, "trace_sha256": None, "array_hashes": None}


class PilotEngine:
    """The pilot's own generation and measurement, with implementation checks I-1 to I-3."""

    def __init__(self, context: pc.Context):
        self.context = context

    def execute(self, condition: str, family: str, replica: int) -> tuple[dict, dict | None]:
        """(record, arrays) of a run that reached the measurement layer, else (failure record, None)."""
        context, plan = self.context, self.context.plan
        dlib, measure = context.dlib, context.measure
        seeds = dlib.seeds_for(plan.namespace, condition, replica)
        started = time.perf_counter()
        try:
            run = dlib.generate_run(plan.namespace, condition, family, replica)
        except Exception as error:  # a sampler or adapter exception is a recorded run failure
            return failure_record(plan, condition, family, replica, "generation", error, seeds, _partial_state(error)), None
        generated = time.perf_counter() - started
        try:
            run["I-2"] = dlib.check_i2(run)
            if family != "healthy":
                run["I-3"] = dlib.check_i3(run)
                if replica == 0:
                    run["I-3-step"] = step_check(dlib, run)
            run["failures"] = [{"phase": name, "reason": run[name]["mismatches"]}
                               for name in ("I-2", "I-3", "I-3-step") if name in run and not run[name]["passed"]]
            checked = time.perf_counter() - started - generated
            record, arrays = measure.measure_run(run)
        except Exception as error:
            hashes = {key: run.get(key) for key in ("emitted_sha256", "base_sha256", "prefix_sha256", "start_sha256")}
            return failure_record(plan, condition, family, replica, "measurement", error, seeds,
                                  {"completed_sweeps": int(len(run["marks"])), "local_bytes": int(run["local_bytes"]),
                                   "available_hashes": hashes}), None
        record["generation_seconds"], record["check_seconds"] = generated, checked
        return record, arrays


def load_engine(spec: str, context: pc.Context) -> Any:
    """'production', or '<file.py>:<factory>' where factory(context) returns an engine (tests)."""
    if spec == "production":
        return PilotEngine(context)
    if context.plan.kind != "test":
        raise pc.PilotError("only the production engine may run a pilot or smoke plan")
    path, _, name = spec.rpartition(":")
    return getattr(pc.load_module(f"d_pilot_engine_{pc.sha256_bytes(path.encode())[:12]}", Path(path)), name)(context)


# ---------------------------------------------------------------------------
# Storage
# ---------------------------------------------------------------------------
def _fsync_file(path: Path) -> None:
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def commit_run(context: pc.Context, record: dict, arrays: dict | None) -> dict:
    """Store one run so that an interruption leaves either the whole run or no record.

    The measurement layer writes its two files into a staging directory on the same file
    system; they are flushed and renamed into place, the record last. A failure record has
    no trace and is written in one step. An existing record or trace is never replaced.
    """
    condition, replica = record["condition"], int(record["replica"])
    trace_rel, record_rel = pc.run_paths(condition, replica)
    final_trace, final_record = context.root / trace_rel, context.root / record_rel
    if final_record.exists() or final_trace.exists():
        raise pc.PilotError(f"run {condition}/{replica:05d} already has stored files; they are not replaced")
    if arrays is None:
        digest = pc.atomic_write(final_record, pc.canonical_json(record))
        return {"condition": condition, "replica": replica, "status": record.get("status"),
                "files": {record_rel.as_posix(): digest}}
    staging = context.root / ".staging" / f"{os.getpid()}-{condition.replace('/', '_')}-{replica:05d}"
    try:
        entry = context.measure.write_run(record, arrays, staging)
        _fsync_file(staging / trace_rel)
        _fsync_file(staging / record_rel)
        for relative, final in ((trace_rel, final_trace), (record_rel, final_record)):  # the record goes last
            final.parent.mkdir(parents=True, exist_ok=True)
            os.replace(staging / relative, final)
            pc.fsync_directory(final.parent)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
        try:
            staging.parent.rmdir()      # only if no other process is staging a run
        except OSError:
            pass
    return entry


def classify(context: pc.Context, condition: str, replica: int, *, verify_hash: bool = True) -> dict:
    """State of one planned run on disk. Nothing is modified."""
    record_path, trace_path = context.record_path(condition, replica), context.trace_path(condition, replica)
    out: dict[str, Any] = {"condition": condition, "replica": replica, "reason": None}
    if not record_path.exists():
        return {**out, "state": "PARTIAL" if trace_path.exists() else "MISSING",
                "reason": "trace without a record" if trace_path.exists() else None}
    try:
        raw = record_path.read_bytes()
        record = json.loads(raw.decode("utf-8"))
        if not isinstance(record, dict):
            raise ValueError("record is not a JSON object")
    except (ValueError, OSError) as error:
        return {**out, "state": "CORRUPT", "reason": f"unreadable record: {type(error).__name__}: {error}"}
    out.update(record_sha256=pc.sha256_bytes(raw), status=record.get("status"))
    if (record.get("condition") != condition or record.get("replica") != replica
            or record.get("namespace") != context.plan.namespace):
        return {**out, "state": "CORRUPT", "reason": "record identity differs from its path or namespace"}
    if record.get("trace_path") is None:
        if record.get("status") != "FAILED" or not record.get("failure"):
            return {**out, "state": "CORRUPT", "reason": "record without a trace that is not a failure record"}
        if trace_path.exists():
            return {**out, "state": "CORRUPT", "reason": "failure record beside a trace file"}
        return {**out, "state": "FAILED_RECORD", "failure_phase": record["failure"].get("phase")}
    if not trace_path.exists():
        return {**out, "state": "CORRUPT", "reason": "record without its trace"}
    if verify_hash:
        actual = pc.sha256_file(trace_path)
        if actual != record.get("trace_sha256"):
            return {**out, "state": "CORRUPT", "reason": "trace hash differs from the record"}
        out["trace_sha256"] = actual
    return {**out, "state": "COMPLETE"}


def scan(context: pc.Context, *, verify_hash: bool = True) -> list[dict]:
    return [classify(context, condition, replica, verify_hash=verify_hash) for condition, _, replica in context.plan.runs()]


def journal(context: pc.Context, event: str, **fields: Any) -> None:
    """Append one event (start, commit, interrupted, restart, worker_error) as a single write."""
    path = context.result("journal.jsonl")
    path.parent.mkdir(parents=True, exist_ok=True)
    line = json.dumps({"time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "pid": os.getpid(), "event": event, **fields},
                      sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
    descriptor = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)
    try:
        os.write(descriptor, line.encode("utf-8"))
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def journal_events(context: pc.Context) -> list[dict]:
    path = context.result("journal.jsonl")
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line.strip()] if path.exists() else []


# ---------------------------------------------------------------------------
# Authorization and the pre-generation manifest
# ---------------------------------------------------------------------------
def code_hashes() -> dict[str, str]:
    """Every file the pilot executes: this directory's code and frozen documents, the imported
    Stage C modules, and every module under src/ that is loaded when the pilot's modules are."""
    paths = {p for p in HERE.glob("*.py")} | {HERE / "protocol.json", HERE / "design.json"}
    paths |= {pc.REPO / "research/v018-h1-detection" / name for name in ("h1lib.py", "analyze.py", "run_campaign.py")}
    source = (pc.REPO / "src").resolve()
    for module in list(sys.modules.values()):
        file = getattr(module, "__file__", None)
        if file and Path(file).resolve().is_relative_to(source):
            paths.add(Path(file).resolve())
    return {p.resolve().relative_to(pc.REPO.resolve()).as_posix(): pc.sha256_file(p) for p in sorted(paths)}


def i4_selection(context: pc.Context) -> dict[str, list[int]]:
    """protocol.json implementation_checks.I-4.sample."""
    from rnglab.seeds import namespace_seed
    chosen = {}
    for condition, _, count in context.plan.conditions:
        rng = np.random.default_rng(namespace_seed(context.plan.namespace, condition, "analysis", "I4-replay-sample"))
        chosen[condition] = np.sort(rng.choice(np.arange(count, dtype=np.int64), size=min(5, count), replace=False)).tolist()
    return chosen


def build_pre_generation_manifest(context: pc.Context, *, interrupted_run_policy: str, tdk_deadline: str | None,
                                  workers: int = pc.PROTOCOL_WORKERS) -> dict:
    """The record that must exist, and be pushed by the author, before the first draw.

    interrupted_run_policy has no default: it is the author's decision (see the module text).
    tdk_deadline may be None; stop criterion 5 is then INSUFFICIENT_INFORMATION, as the protocol says.
    """
    if interrupted_run_policy not in POLICIES:
        raise pc.PilotError(f"interrupted_run_policy must be one of {POLICIES}; it has no default")
    if context.plan.kind == "test":
        raise pc.PilotError("a test plan has no pre-generation manifest")
    plan = context.plan
    lock = pc.REPO / "research/phase2-20261030/environment/requirements-lock.txt"

    def git(*args: str) -> str | None:
        try:
            return subprocess.run(["git", "-C", str(pc.REPO), *args], capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return None

    dirty = git("status", "--porcelain")
    return {"schema": "rnglab-v019-d-pilot-pre-generation-manifest-v1", "namespace": plan.namespace, "plan_kind": plan.kind,
            "protocol_sha256": pc.protocol_sha256(), "design_sha256": pc.sha256_file(HERE / "design.json"),
            "code_sha256": code_hashes(), "environment_lock": {"path": lock.relative_to(pc.REPO).as_posix(),
                                                               "sha256": pc.sha256_file(lock) if lock.is_file() else None},
            "software": pc.software_versions(), "git": {"branch": git("rev-parse", "--abbrev-ref", "HEAD"),
                                                        "commit": git("rev-parse", "HEAD"),
                                                        "working_tree_clean": None if dirty is None else dirty == ""},
            "planned_runs": [{"condition": condition, "family": family, "replica": replica,
                              "split": plan.split(condition, replica),
                              "seed_labels": {role: f"{plan.namespace}/{condition}/{role}/{replica:05d}"
                                              for role in ("local", "centroid", "start", "fault", "profile")},
                              "seeds": context.dlib.seeds_for(plan.namespace, condition, replica)}
                             for condition, family, replica in plan.runs()],
            "i4_replay_selection": i4_selection(context),
            "workers": workers, "protocol_workers": pc.PROTOCOL_WORKERS,
            "workers_deviate_from_protocol": workers != pc.PROTOCOL_WORKERS,
            "thread_variables": {name: "1" for name in pc.THREAD_VARIABLES},
            "interrupted_run_policy": interrupted_run_policy, "tdk_deadline": tdk_deadline,
            "approved_amendments": pc.approved_amendments(),
            "written_at": time.strftime("%Y-%m-%dT%H:%M:%S%z")}


def require_generation_allowed(context: pc.Context) -> dict | None:
    """Stop unless this plan may draw. Returns the manifest for pilot and smoke plans."""
    plan = context.plan
    if plan.kind == "test":
        if plan.namespace in pc.RESERVED_NAMESPACES:
            raise pc.PilotError("a test plan must not use the pilot or smoke namespace")
        return None
    manifest_path, authorization_path = context.root / "pre_generation_manifest.json", context.root / "authorization.json"
    if not manifest_path.is_file():
        raise pc.PilotError("pre_generation_manifest.json is missing: no draw before it exists")
    manifest = pc.read_json(manifest_path)
    if manifest.get("namespace") != plan.namespace or manifest.get("protocol_sha256") != pc.protocol_sha256():
        raise pc.PilotError("pre-generation manifest is for another namespace or protocol")
    current = code_hashes()
    changed = sorted(p for p, h in manifest["code_sha256"].items() if current.get(p) != h)
    added = sorted(p for p in current if p not in manifest["code_sha256"])
    if changed or added:
        raise pc.PilotError(f"code differs from the pre-generation manifest: changed {changed}, not listed {added}")
    if manifest.get("interrupted_run_policy") not in POLICIES:
        raise pc.PilotError("the pre-generation manifest names no interrupted-run policy")
    if not authorization_path.is_file():
        raise pc.PilotError("authorization.json is missing: the author has not authorised generation")
    authorization = pc.read_json(authorization_path)
    if (authorization.get("pre_generation_manifest_sha256") != pc.sha256_file(manifest_path)
            or authorization.get("namespace") != plan.namespace
            or not authorization.get("authorized_by") or not authorization.get("authorized_on")):
        raise pc.PilotError("authorization.json does not authorise this manifest and namespace")
    for name in pc.THREAD_VARIABLES:
        if os.environ.get(name) != "1":
            raise pc.PilotError(f"{name} must be 1 before NumPy is imported (protocol.json storage.parallelism)")
    return manifest


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def run_one(context: pc.Context, engine: Any, condition: str, family: str, replica: int) -> dict:
    """Execute and store one planned run; every outcome becomes a record."""
    journal(context, "start", condition=condition, replica=replica)
    record, arrays = engine.execute(condition, family, replica)
    try:
        entry = commit_run(context, record, arrays)
    except pc.PilotError:
        raise
    except Exception as error:  # the measurement layer reports serialization errors as RuntimeError
        moved = None
        if context.trace_path(condition, replica).exists() and not context.record_path(condition, replica).exists():
            moved = _set_aside(context, condition, replica)
        record = failure_record(context.plan, condition, family, replica, "serialization", error, record.get("seeds"),
                                {"trace_set_aside": moved})
        entry = commit_run(context, record, None)
    journal(context, "commit", condition=condition, replica=replica, status=record.get("status"))
    return entry


def _set_aside(context: pc.Context, condition: str, replica: int) -> str:
    """Move a trace that has no record out of the result tree. It is kept, not deleted."""
    source = context.trace_path(condition, replica)
    target = context.root / "interrupted" / pc.run_paths(condition, replica)[0]
    target = target.with_name(f"{target.name}.{int(time.time() * 1000)}")
    target.parent.mkdir(parents=True, exist_ok=True)
    os.replace(source, target)
    return target.relative_to(context.root).as_posix()


def prepare(context: pc.Context, policy: str | None) -> list[tuple[str, str, int]]:
    """Decide what remains to be run. Raises on anything that does not verify."""
    shutil.rmtree(context.root / ".staging", ignore_errors=True)  # staging holds no results
    states = {(s["condition"], s["replica"]): s for s in scan(context)}
    corrupt = [s for s in states.values() if s["state"] == "CORRUPT"]
    if corrupt:
        listing = "; ".join(f"{s['condition']}/{s['replica']:05d}: {s['reason']}" for s in corrupt[:10])
        raise pc.PilotError(f"{len(corrupt)} stored runs do not verify and are not repaired: {listing}")
    events = journal_events(context)
    started = {(e["condition"], e["replica"]) for e in events if e["event"] == "start"}
    pending = []
    for condition, family, replica in context.plan.runs():
        state = states[(condition, replica)]["state"]
        if state in ("COMPLETE", "FAILED_RECORD"):
            continue
        interrupted = state == "PARTIAL" or (condition, replica) in started
        if interrupted:
            if policy not in POLICIES:
                raise pc.PilotError(f"run {condition}/{replica:05d} was started and left no record; an "
                                    "interrupted-run policy chosen by the author is required to continue")
            moved = _set_aside(context, condition, replica) if state == "PARTIAL" else None
            journal(context, "interrupted", condition=condition, replica=replica, trace_set_aside=moved, policy=policy)
            if policy == "record_failure":
                seeds = context.dlib.seeds_for(context.plan.namespace, condition, replica)
                commit_run(context, failure_record(context.plan, condition, family, replica, "interrupted", None, seeds,
                                                   {"trace_set_aside": moved}), None)
                continue
            journal(context, "restart", condition=condition, replica=replica)
        pending.append((condition, family, replica))
    return pending


def work(context: pc.Context, engine_spec: str, jobs: list[tuple[str, str, int]]) -> int:
    """Run the given planned runs in this process, in order."""
    engine = load_engine(engine_spec, context)
    for condition, family, replica in jobs:
        run_one(context, engine, condition, family, replica)
    return len(jobs)


def generate(context: pc.Context, *, engine_spec: str = "production", workers: int | None = None,
             worker_command: list[str] | None = None, test_policy: str | None = None) -> dict:
    """Run every planned run that has no record yet, then write the run manifest and the failure list.

    workers = 0 runs in this process. Otherwise `workers` operating-system processes are started
    with `worker_command` plus a shard argument; shard k takes every workers-th pending run, so
    the shards have the same mix of conditions. Results do not depend on the number of workers.
    """
    manifest = require_generation_allowed(context)
    policy = manifest["interrupted_run_policy"] if manifest else test_policy  # test_policy: test plans only
    if workers is None:
        workers = manifest["workers"] if manifest else 0
    context.root.mkdir(parents=True, exist_ok=True)
    pending = prepare(context, policy)
    started = time.time()
    if workers == 0 or not pending:
        work(context, engine_spec, pending)
        codes = []
    else:
        if not worker_command:
            raise pc.PilotError("a worker command is required to start worker processes")
        pc.atomic_write(context.root / ".pending.json", pc.canonical_json([list(p) for p in pending]), overwrite=True)
        environment = {**os.environ, **{name: "1" for name in pc.THREAD_VARIABLES}}
        processes = [subprocess.Popen([*worker_command, "--shard", f"{k}/{workers}"], env=environment)
                     for k in range(workers)]
        codes = [process.wait() for process in processes]
    summary = finalize(context, workers=workers, seconds=time.time() - started, manifest=manifest)
    if any(codes):
        raise pc.PilotError(f"worker processes ended with errors {codes}; stored runs are intact, see the journal")
    return summary


def shard(pending: list[tuple[str, str, int]], specification: str) -> list[tuple[str, str, int]]:
    k, _, n = specification.partition("/")
    return pending[int(k)::int(n)]


def finalize(context: pc.Context, *, workers: int, seconds: float | None, manifest: Mapping | None) -> dict:
    """results/run_manifest.json and results/failures.jsonl from what is stored (protocol.json storage.manifest)."""
    entries, failures = [], []
    events = journal_events(context)
    for condition, family, replica in context.plan.runs():
        state = classify(context, condition, replica)
        trace_rel, record_rel = pc.run_paths(condition, replica)
        entry = {"condition": condition, "family": family, "replica": replica, "split": context.plan.split(condition, replica),
                 "seeds": context.dlib.seeds_for(context.plan.namespace, condition, replica), "state": state["state"],
                 "reason": state["reason"], "status": state.get("status"),
                 "attempts": sum(1 for e in events if e["event"] == "start" and e["condition"] == condition and e["replica"] == replica),
                 "record": {"path": record_rel.as_posix(), "sha256": state.get("record_sha256")},
                 "trace": {"path": trace_rel.as_posix(), "sha256": state.get("trace_sha256")}}
        if state["state"] in ("COMPLETE", "FAILED_RECORD"):
            record = pc.read_json(context.record_path(condition, replica))
            entry.update(array_hashes=record.get("array_hashes"), failure=record.get("failure"),
                         completed_sweeps=record.get("completed_sweeps"),
                         completed_battery_checkpoints=record.get("completed_battery_checkpoints"),
                         completed_application_checkpoints=record.get("completed_application_checkpoints"))
            for item in record.get("failures") or []:
                failures.append({"condition": condition, "replica": replica, "disposition": "recorded, not replaced", **item})
        else:
            failures.append({"condition": condition, "replica": replica, "phase": "storage", "state": state["state"],
                             "reason": state["reason"], "disposition": "no valid record"})
        entries.append(entry)
    counts = {name: sum(1 for e in entries if e["state"] == name) for name in STATES}
    body = {"schema": "rnglab-v019-d-pilot-run-manifest-v1", "namespace": context.plan.namespace, "plan_kind": context.plan.kind,
            "protocol_sha256": pc.protocol_sha256(), "design_sha256": pc.sha256_file(HERE / "design.json"),
            "pre_generation_manifest_sha256": (pc.sha256_file(context.root / "pre_generation_manifest.json")
                                               if manifest is not None else None),
            "software": pc.software_versions(), "workers": workers,
            "thread_variables": {name: os.environ.get(name) for name in pc.THREAD_VARIABLES},
            "generation_wall_seconds_this_invocation": seconds, "counts": counts, "planned": len(entries), "runs": entries}
    digest = pc.atomic_write(context.result("run_manifest.json"), pc.canonical_json(body), overwrite=True)
    pc.atomic_write(context.result("failures.jsonl"), pc.canonical_jsonl(failures), overwrite=True)
    return {"counts": counts, "planned": len(entries), "failures": len(failures), "run_manifest_sha256": digest}


def require_complete(context: pc.Context) -> dict:
    """Every planned run has a verified record (success or failure) that is the one the run manifest
    lists. Analysis stages call this first; a record changed after generation stops them."""
    states = scan(context)
    bad = [s for s in states if s["state"] not in ("COMPLETE", "FAILED_RECORD")]
    if bad:
        listing = "; ".join(f"{s['condition']}/{s['replica']:05d}: {s['state']}" for s in bad[:10])
        raise pc.PilotError(f"{len(bad)} planned runs have no valid record: {listing}")
    manifest_path = context.result("run_manifest.json")
    if not manifest_path.is_file():
        raise pc.PilotError("results/run_manifest.json does not exist: generation has not been finalised")
    listed = {(e["condition"], e["replica"]): e for e in pc.read_json(manifest_path)["runs"]}
    changed = [s for s in states if listed.get((s["condition"], s["replica"]), {}).get("record", {}).get("sha256")
               != s.get("record_sha256")]
    if changed or len(listed) != len(states):
        listing = "; ".join(f"{s['condition']}/{s['replica']:05d}" for s in changed[:10])
        raise pc.PilotError(f"{len(changed)} run records differ from the run manifest: {listing}")
    return {"complete": sum(1 for s in states if s["state"] == "COMPLETE"),
            "failed_records": sum(1 for s in states if s["state"] == "FAILED_RECORD")}


# ---------------------------------------------------------------------------
# I-4 replay
# ---------------------------------------------------------------------------
def _comparable(record: Mapping[str, Any]) -> dict:
    return json.loads(json.dumps({k: v for k, v in record.items() if k not in TIMING_FIELDS + STORAGE_FIELDS}))


def replay_check(context: pc.Context, engine: Any, selection: Mapping[str, list[int]] | None = None) -> dict:
    """protocol.json implementation_checks.I-4: regenerate the selected runs and compare everything
    except timing fields. A selected run that has only a failure record is reported, not replaced."""
    selection = i4_selection(context) if selection is None else selection
    rows = []
    for condition, replicas in selection.items():
        family = context.plan.family(condition)
        for replica in replicas:
            state = classify(context, condition, replica)
            row: dict[str, Any] = {"condition": condition, "replica": replica, "state": state["state"], "mismatches": []}
            if state["state"] != "COMPLETE":
                row.update(passed=False, reason="no complete stored run to compare with")
                rows.append(row)
                continue
            stored, stored_arrays = context.measure.load_run(context.root, condition, replica)
            record, arrays = engine.execute(condition, family, replica)
            if arrays is None:
                row["mismatches"].append({"field": "regeneration", "actual": record.get("failure")})
            else:
                hashes = context.measure.array_hashes(arrays)
                for key, value in hashes.items():
                    if value != stored["array_hashes"].get(key):
                        row["mismatches"].append({"field": f"array:{key}"})
                a, b = _comparable(stored), _comparable(record)
                for key in sorted(set(a) | set(b)):
                    if a.get(key) != b.get(key):
                        row["mismatches"].append({"field": key})
            row["passed"] = not row["mismatches"]
            rows.append(row)
    return {"check": "I-4", "selection": {k: list(v) for k, v in selection.items()}, "runs": rows,
            "passed": bool(rows) and all(r["passed"] for r in rows),
            "ignored_fields": list(TIMING_FIELDS), "storage_fields_not_compared": list(STORAGE_FIELDS)}
