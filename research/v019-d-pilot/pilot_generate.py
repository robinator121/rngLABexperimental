"""Stage D pilot: generation driver, atomic storage, resume, manifests, I-4 replay. protocol.json governs.

Interrupted runs. The protocol does not say what happens to a run that was started and left
no record because its process died. The author decided it on 2026-11-13 (`INTERRUPTED_RUN_RULES`,
copied into the pre-generation manifest):

    regenerate_same_id   the same identity is run again: same run ID, seed labels, seeds,
                         generator configuration and protocol parameters. Before every attempt
                         the first outputs of the run's five generators are hashed and compared
                         with the manifest (`rng_input_sha256`). The hash verifies identity and
                         is used for nothing else. Nothing is assumed about floating-point
                         output being identical in another environment.
    three retries        a run may be attempted four times. After the fourth interruption it
                         gets a failure record with phase 'infrastructure' and the status
                         UNRESOLVED_INFRASTRUCTURE, is never run again, and the stop criteria
                         it touches cannot be PASS (pilot_analysis.criteria). The status says
                         that the run has no result; it does not say that the infrastructure
                         was the cause. The record keeps every attempt and its error evidence.
    two kinds of failure an interruption of the infrastructure (the process ended without a
                         record, an error of the operating system or of memory) is retried.
                         An exception of the sampler, a check, the measurement or the encoding
                         is experimental evidence: it gets a failure record and is never retried.

Every attempt, interruption and commit is one line of results/journal.jsonl. The file is only
appended to; every line carries the SHA-256 of the line before it.

A run with a record is never run again, whatever the record says (protocol.json failed_runs).
A record or trace that does not verify is reported and stops the driver; nothing is repaired.

Generation in the pilot or smoke namespace needs the pre-generation manifest and the
author's authorization record that names that manifest's hash (`require_generation_allowed`).
"""
from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shutil
import signal
import subprocess
import sys
import threading
import time
import traceback
from pathlib import Path
from typing import Any, Mapping

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.append(str(HERE))
import pilot_common as pc  # noqa: E402

STATES = ("COMPLETE", "FAILED_RECORD", "MISSING", "PARTIAL", "CORRUPT")
POLICY = "regenerate_same_id"
MAX_INFRASTRUCTURE_RETRIES = 3
MAX_ATTEMPTS = 1 + MAX_INFRASTRUCTURE_RETRIES
INFRASTRUCTURE_PHASE = "infrastructure"
UNRESOLVED = "UNRESOLVED_INFRASTRUCTURE"
RNG_ROLES = ("local", "centroid", "start", "fault", "profile")
RNG_INPUT_WORDS = 128
# The author's decision of 2026-11-13 on a question the frozen protocol leaves open. It changes no
# seed, sample size, threshold or criterion. The pre-generation manifest carries a copy.
INTERRUPTED_RUN_RULES = {
    "policy": POLICY, "decided_by": "author", "decided_on": "2026-11-13", "confirmed_on": "2026-11-14",
    "max_infrastructure_retries": MAX_INFRASTRUCTURE_RETRIES,
    "retry_identity": "same run ID, seed labels, seeds, generator configuration and protocol parameters",
    "infrastructure_interruptions": ["process_ended_without_record", "operating_system_error", "memory_error"],
    "experimental_failures": "an exception of the sampler, a check, the measurement or the encoding gets a failure "
                             "record (protocol.json failed_runs) and is never retried",
    "exhaustion": f"failure record with phase '{INFRASTRUCTURE_PHASE}' and infrastructure_status {UNRESOLVED}; "
                  "the stop criteria that use the run cannot be PASS",
    "unresolved_status": f"{UNRESOLVED} is an execution status: the run has no result after "
                         f"{1 + MAX_INFRASTRUCTURE_RETRIES} attempts. It does not establish that the infrastructure was "
                         "the cause. The record and the audit log keep every attempt and the error evidence of each",
    "rng_input": f"per role, SHA-256 of the first {RNG_INPUT_WORDS} 64-bit outputs, big-endian: the repository's PCG64 "
                 "source for local, centroid and start; numpy.random.default_rng(seed).bit_generator for fault and "
                 "profile. rng_input_sha256 is the SHA-256 of the compact sorted JSON of the five hashes",
    "rng_input_use": "identity verification only, computed when the manifest is written and before every attempt, "
                     "authorised by the author on 2026-11-14. It is not experimental sampling. No generated value "
                     "may influence parameters, selection, measurements or scientific decisions",
    "floating_point": "no bit-identical floating-point output is assumed across execution environments; every "
                      "attempt records the hash of its software versions",
}
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


class InfrastructureInterruption(pc.PilotError):
    """An attempt ended for a reason of the infrastructure. It left no record and may be retried."""

    def __init__(self, condition: str, replica: int, attempt: int, reason: Mapping[str, Any]):
        super().__init__(f"run {condition}/{replica:05d}, attempt {attempt}: interrupted by the infrastructure "
                         f"({reason['kind']}: {reason['detail']}); no record was written")
        self.reason = dict(reason)


def infrastructure_cause(error: BaseException) -> BaseException | None:
    """The operating-system or memory error behind an exception, if there is one.

    OSError covers a full or failing storage medium, a quota, a lost mount. MemoryError is the
    interpreter running out of memory. Neither is a property of the draw. Everything else
    (ArithmeticError, ValueError, a failed check, an encoding error) is the run's own failure.
    """
    seen: set[int] = set()
    while error is not None and id(error) not in seen:
        if isinstance(error, (OSError, MemoryError)):
            return error
        seen.add(id(error))
        error = error.__cause__ or error.__context__
    return None


def recorded_memory_error(value: Any) -> dict | None:
    """The first exception entry {"class": "MemoryError", ...} anywhere in a run record, if there is one.

    The measurement layer turns an exception inside a battery test or an endpoint into a recorded
    failure of that checkpoint. For a numerical or statistical failure that is the evidence. Running
    out of memory there says nothing about the draw, so such a record is not stored and the attempt
    counts as interrupted.
    """
    if isinstance(value, Mapping):
        if value.get("class") == "MemoryError":
            return dict(value)
        value = list(value.values())
    if isinstance(value, (list, tuple)):
        for item in value:
            found = recorded_memory_error(item)
            if found is not None:
                return found
    return None


def _evidence_text(text: str, root: Path | None) -> str:
    """Error text as raised, with the output root and the repository root written as <root> and <repo>,
    so that no path of the machine enters a result file."""
    for path, name in ((root, "<root>"), (pc.REPO, "<repo>")):
        if path is not None:
            for form in {str(Path(path).resolve()), str(path)}:
                text = text.replace(form, name)
    return text


def _frame_file(filename: str) -> str:
    path = Path(filename)
    try:
        return path.resolve().relative_to(pc.REPO.resolve()).as_posix()
    except ValueError:
        pass
    if "site-packages" in path.parts:
        return "/".join(path.parts[path.parts.index("site-packages"):])
    return path.name


def interruption_reason(error: BaseException, root: Path | None = None) -> dict:
    """What is known about an interruption: its kind, and the error as raised.

    The evidence is kept whole: every exception of the chain with its message, and the frames in
    which the operating-system or memory error was raised. Nothing is summarised away, because the
    kind is a classification and not a finding about the cause.
    """
    cause = infrastructure_cause(error)
    chain, seen, link = [], set(), error
    while link is not None and id(link) not in seen:
        chain.append({"class": type(link).__name__, "message": _evidence_text(str(link), root)})
        seen.add(id(link))
        link = link.__cause__ or link.__context__
    frames = [f"{_frame_file(frame.filename)}:{frame.lineno} in {frame.name}"
              for frame in traceback.extract_tb(cause.__traceback__)]
    return {"kind": "memory_error" if isinstance(cause, MemoryError) else "operating_system_error",
            "detail": _evidence_text(f"{type(cause).__name__}: {cause}", root),
            "exception_chain": chain, "traceback": frames,
            "paths": "the output root and the repository root are written as <root> and <repo>"}


def rng_input(dlib: Any, seeds: Mapping[str, int]) -> dict:
    """INTERRUPTED_RUN_RULES['rng_input']: what the run's five generators emit first, hashed per role.

    Integer arithmetic only, so the value is the same in every environment. It identifies the
    random input of an attempt without reference to any floating-point result.
    """
    def digest(role: str) -> str:
        if role in ("local", "centroid", "start"):
            data = dlib.PCG64Source(seeds[role]).read_bytes(8 * RNG_INPUT_WORDS)
        else:
            data = np.random.default_rng(seeds[role]).bit_generator.random_raw(RNG_INPUT_WORDS).astype(">u8").tobytes()
        return hashlib.sha256(data).hexdigest()

    roles = {role: digest(role) for role in RNG_ROLES}
    compact = json.dumps(roles, sort_keys=True, separators=(",", ":")).encode("ascii")
    return {"roles": roles, "sha256": hashlib.sha256(compact).hexdigest()}


def environment_sha256() -> str:
    """Hash of the software versions and machine type of this process."""
    return pc.sha256_bytes(json.dumps(pc.software_versions(), sort_keys=True, separators=(",", ":")).encode("utf-8"))


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
            if infrastructure_cause(error) is not None:
                raise                   # not the run's failure: run_one logs the interruption
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
            if infrastructure_cause(error) is not None:
                raise
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


def _journal_tail(descriptor: int) -> tuple[bytes | None, bool]:
    """(last complete line, whether an incomplete line follows it) of the open journal."""
    size = os.fstat(descriptor).st_size
    if size == 0:
        return None, False
    length = 1 << 16
    while True:                              # read backwards until the last complete line is whole
        block = os.pread(descriptor, min(size, length), max(0, size - length))
        if length >= size or block.count(b"\n") >= 2:
            break
        length *= 2
    torn = not block.endswith(b"\n")
    lines = block.split(b"\n")[:-1]          # what follows the last newline is not a line
    return (lines[-1] if lines else None), torn


def journal(context: pc.Context, event: str, **fields: Any) -> dict:
    """Append one event to results/journal.jsonl: the audit log of attempts, interruptions and commits.

    Events: generate_begin, start, commit, interrupted, restart, unresolved, trace_set_aside,
    stop_requested, worker_exit, generate_end. The file is opened for appending and locked for the duration of
    one write. Each entry has a sequence number and the SHA-256 of the previous line, so a
    removed, reordered or edited line is detectable. Nothing in this module rewrites the file.
    """
    path = context.result("journal.jsonl")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_RDWR | os.O_APPEND | os.O_CREAT, 0o644)
    try:
        fcntl.flock(descriptor, fcntl.LOCK_EX)
        last, torn = _journal_tail(descriptor)
        try:
            previous = json.loads(last) if last is not None else None
            if previous is not None and not isinstance(previous.get("seq"), int):
                raise ValueError("no sequence number")
        except (ValueError, AttributeError) as error:
            raise pc.PilotError(f"the last line of results/journal.jsonl is not a journal entry ({error}); "
                                "nothing is appended to a journal that does not verify") from None
        entry = {"seq": previous["seq"] + 1 if previous else 1,
                 "previous_sha256": pc.sha256_bytes(last) if last is not None else None,
                 "time": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "pid": os.getpid(), "event": event, **fields}
        if torn:    # a write that a power loss cut short: left in place, closed with a newline, and said so
            entry["after_incomplete_line"] = True
        line = json.dumps(entry, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n"
        os.write(descriptor, (b"\n" if torn else b"") + line.encode("utf-8"))
        os.fsync(descriptor)
    finally:
        os.close(descriptor)                 # closing releases the lock
    return entry


def journal_events(context: pc.Context) -> list[dict]:
    """Every entry of the journal, after checking its chain. A journal that does not verify stops the caller."""
    path = context.result("journal.jsonl")
    if not path.exists():
        return []
    events: list[dict] = []
    previous: bytes | None = None
    skipped = False
    for number, line in enumerate(path.read_bytes().split(b"\n"), 1):
        if not line:
            continue
        try:
            entry = json.loads(line)
            if not isinstance(entry, dict) or "seq" not in entry:
                raise ValueError("not a journal entry")
        except ValueError:
            skipped = True                   # acceptable only if the next entry declares it
            continue
        expected = pc.sha256_bytes(previous) if previous is not None else None
        if (entry["seq"] != len(events) + 1 or entry.get("previous_sha256") != expected
                or bool(entry.get("after_incomplete_line")) != skipped):
            raise pc.PilotError(f"results/journal.jsonl does not verify at line {number}: an entry was removed, "
                                "reordered or changed. It is not repaired")
        events.append(entry)
        previous, skipped = line, False
    return events


def attempt_history(events: list[dict]) -> dict[tuple[str, int], dict]:
    """Per run: number of starts, interruptions by attempt number, commits and the unresolved mark."""
    history: dict[tuple[str, int], dict] = {}
    for entry in events:
        if "condition" not in entry or "replica" not in entry:
            continue
        run = history.setdefault((entry["condition"], entry["replica"]),
                                 {"starts": 0, "started": [], "interrupted": {}, "commits": 0, "unresolved": False,
                                  "environments": []})
        if entry["event"] == "start":
            run["starts"] += 1
            run["started"].append({"attempt": entry.get("attempt"), "time": entry.get("time"), "pid": entry.get("pid"),
                                   "environment_sha256": entry.get("environment_sha256"), "journal_seq": entry.get("seq")})
            if entry.get("environment_sha256") and entry["environment_sha256"] not in run["environments"]:
                run["environments"].append(entry["environment_sha256"])
        elif entry["event"] == "interrupted":
            run["interrupted"][entry.get("attempt")] = {
                "reason": entry.get("reason"), "detected_by": entry.get("detected_by"), "time": entry.get("time"),
                "trace_set_aside": entry.get("trace_set_aside"), "journal_seq": entry.get("seq")}
        elif entry["event"] == "commit":
            run["commits"] += 1
        elif entry["event"] == "unresolved":
            run["unresolved"] = True
    return history


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


def build_pre_generation_manifest(context: pc.Context, *, tdk_deadline: str | None,
                                  workers: int = pc.PROTOCOL_WORKERS, interrupted_run_policy: str = POLICY) -> dict:
    """The record that must exist, and be pushed by the author, before the first draw.

    The interrupted-run policy is the author's decision of 2026-11-13 (INTERRUPTED_RUN_RULES); no
    other policy is implemented. Every planned run carries the hash of its random input.
    tdk_deadline may be None; stop criterion 5 is then INSUFFICIENT_INFORMATION, as the protocol says.
    """
    if interrupted_run_policy != POLICY:
        raise pc.PilotError(f"the author decided the interrupted-run policy {POLICY!r} on "
                            f"{INTERRUPTED_RUN_RULES['decided_on']}; {interrupted_run_policy!r} is not implemented")
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
    planned = []
    for condition, family, replica in plan.runs():
        seeds = context.dlib.seeds_for(plan.namespace, condition, replica)
        planned.append({"condition": condition, "family": family, "replica": replica, "split": plan.split(condition, replica),
                        "seed_labels": {role: f"{plan.namespace}/{condition}/{role}/{replica:05d}" for role in RNG_ROLES},
                        "seeds": seeds, "rng_input_sha256": rng_input(context.dlib, seeds)["sha256"]})
    return {"schema": "rnglab-v019-d-pilot-pre-generation-manifest-v1", "namespace": plan.namespace, "plan_kind": plan.kind,
            "protocol_sha256": pc.protocol_sha256(), "design_sha256": pc.sha256_file(HERE / "design.json"),
            "code_sha256": code_hashes(), "environment_lock": {"path": lock.relative_to(pc.REPO).as_posix(),
                                                               "sha256": pc.sha256_file(lock) if lock.is_file() else None},
            "software": pc.software_versions(), "git": {"branch": git("rev-parse", "--abbrev-ref", "HEAD"),
                                                        "commit": git("rev-parse", "HEAD"),
                                                        "working_tree_clean": None if dirty is None else dirty == ""},
            "planned_runs": planned,
            "i4_replay_selection": i4_selection(context),
            "workers": workers, "protocol_workers": pc.PROTOCOL_WORKERS,
            "workers_deviate_from_protocol": workers != pc.PROTOCOL_WORKERS,
            "thread_variables": {name: "1" for name in pc.THREAD_VARIABLES},
            "interrupted_run_policy": interrupted_run_policy, "interrupted_run_rules": dict(INTERRUPTED_RUN_RULES),
            "environment_sha256": environment_sha256(), "tdk_deadline": tdk_deadline,
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
    if manifest.get("interrupted_run_policy") != POLICY or manifest.get("interrupted_run_rules") != INTERRUPTED_RUN_RULES:
        raise pc.PilotError("the pre-generation manifest does not carry the interrupted-run policy the author decided")
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
def planned_inputs(manifest: Mapping | None) -> dict[tuple[str, int], dict] | None:
    """Seeds and random-input hash of every planned run, from the pre-generation manifest."""
    if manifest is None:
        return None
    return {(run["condition"], run["replica"]): {"seeds": run["seeds"], "rng_input_sha256": run["rng_input_sha256"]}
            for run in manifest["planned_runs"]}


def run_one(context: pc.Context, engine: Any, condition: str, family: str, replica: int, attempt: int = 1,
            planned: Mapping[tuple[str, int], Mapping] | None = None) -> dict:
    """One attempt at one planned run.

    It ends in a stored record (a result, or the run's own failure), or in an
    InfrastructureInterruption that left no record. Before anything is drawn, the seeds and
    the first outputs of the run's generators are compared with the pre-generation manifest.
    """
    seeds = context.dlib.seeds_for(context.plan.namespace, condition, replica)
    fingerprint = rng_input(context.dlib, seeds)["sha256"]
    expected = planned.get((condition, replica)) if planned is not None else None
    if planned is not None and (expected is None or dict(expected["seeds"]) != dict(seeds)
                                or expected["rng_input_sha256"] != fingerprint):
        raise pc.PilotError(f"run {condition}/{replica:05d}: seeds or random input differ from the pre-generation "
                            "manifest; nothing was drawn")
    journal(context, "start", condition=condition, replica=replica, attempt=attempt, seeds=dict(seeds),
            rng_input_sha256=fingerprint, environment_sha256=environment_sha256())
    try:
        try:
            record, arrays = engine.execute(condition, family, replica)
        except pc.PilotError:
            raise
        except Exception as error:  # an exception the engine did not classify is still the run's own failure
            if infrastructure_cause(error) is not None:
                raise
            record, arrays = failure_record(context.plan, condition, family, replica, "execution", error, seeds), None
        absorbed = recorded_memory_error(record)
        if absorbed is not None:    # not a property of the draw: the record is not stored
            raise MemoryError(f"recorded inside the run: {absorbed.get('message')}")
        try:
            entry = commit_run(context, record, arrays)
        except pc.PilotError:
            raise
        except Exception as error:  # the measurement layer reports serialization errors as RuntimeError
            if infrastructure_cause(error) is not None:
                raise               # the storage medium, not the encoding: no failure record
            moved = None
            if context.trace_path(condition, replica).exists() and not context.record_path(condition, replica).exists():
                moved = _set_aside(context, condition, replica)
            record = failure_record(context.plan, condition, family, replica, "serialization", error, record.get("seeds"),
                                    {"trace_set_aside": moved})
            entry = commit_run(context, record, None)
    except pc.PilotError:
        raise
    except Exception as error:
        if infrastructure_cause(error) is None:
            raise
        reason = interruption_reason(error, context.root)
        try:
            journal(context, "interrupted", condition=condition, replica=replica, attempt=attempt, reason=reason,
                    detected_by="attempt")
        except OSError:
            pass    # the medium that failed may hold the journal too; the next `generate` records the interruption
        raise InfrastructureInterruption(condition, replica, attempt, reason) from error
    journal(context, "commit", condition=condition, replica=replica, attempt=attempt, status=record.get("status"),
            record_sha256=pc.sha256_file(context.record_path(condition, replica)), detected_by="attempt")
    return entry


def _set_aside(context: pc.Context, condition: str, replica: int) -> str:
    """Move a trace that has no record out of the result tree. It is kept, not deleted."""
    source = context.trace_path(condition, replica)
    target = context.root / "interrupted" / pc.run_paths(condition, replica)[0]
    target = target.with_name(f"{target.name}.{int(time.time() * 1000)}")
    target.parent.mkdir(parents=True, exist_ok=True)
    os.replace(source, target)
    return target.relative_to(context.root).as_posix()


def prepare(context: pc.Context) -> list[tuple[str, str, int, int]]:
    """Decide what remains to be run: (condition, family, replica, attempt number). Raises on anything
    that does not verify.

    A run that was started and has no record was interrupted. It is logged, and run again with the
    same identity, at most MAX_INFRASTRUCTURE_RETRIES times; after that it gets its failure record
    with the status UNRESOLVED_INFRASTRUCTURE and is not run again.
    """
    shutil.rmtree(context.root / ".staging", ignore_errors=True)  # staging holds no results
    states = {(s["condition"], s["replica"]): s for s in scan(context)}
    corrupt = [s for s in states.values() if s["state"] == "CORRUPT"]
    if corrupt:
        listing = "; ".join(f"{s['condition']}/{s['replica']:05d}: {s['reason']}" for s in corrupt[:10])
        raise pc.PilotError(f"{len(corrupt)} stored runs do not verify and are not repaired: {listing}")
    history = attempt_history(journal_events(context))
    blank = {"starts": 0, "interrupted": {}, "commits": 0, "unresolved": False}
    # Every attempt is logged, and the log is flushed to disk, before anything is drawn. A trace
    # without a record and without a logged attempt was therefore not left by this driver.
    unknown = [s for s in states.values() if s["state"] == "PARTIAL"
               and not history.get((s["condition"], s["replica"]), blank)["starts"]]
    if unknown:
        listing = "; ".join(f"{s['condition']}/{s['replica']:05d}" for s in unknown[:10])
        raise pc.PilotError(f"{len(unknown)} traces have no record and no logged attempt; their origin is unknown "
                            f"and they are left as found: {listing}")
    pending = []
    for condition, family, replica in context.plan.runs():
        state, past = states[(condition, replica)], history.get((condition, replica), blank)
        if state["state"] in ("COMPLETE", "FAILED_RECORD"):
            if not past["commits"] and not past["unresolved"]:
                # The record was stored and the process ended before the journal said so.
                if state.get("failure_phase") == INFRASTRUCTURE_PHASE:
                    journal(context, "unresolved", condition=condition, replica=replica, attempts=past["starts"],
                            infrastructure_status=UNRESOLVED, record_sha256=state["record_sha256"], detected_by="resume")
                else:
                    journal(context, "commit", condition=condition, replica=replica, attempt=past["starts"] or None,
                            status=state.get("status"), record_sha256=state["record_sha256"], detected_by="resume")
            continue
        if state["state"] == "PARTIAL" or past["starts"]:
            moved = _set_aside(context, condition, replica) if state["state"] == "PARTIAL" else None
            if past["starts"] not in past["interrupted"]:
                journal(context, "interrupted", condition=condition, replica=replica, attempt=past["starts"],
                        reason={"kind": "process_ended_without_record",
                                "detail": "a trace without a record" if moved else "no record and no trace"},
                        detected_by="resume", trace_set_aside=moved)
                past = attempt_history(journal_events(context))[(condition, replica)]
            elif moved:
                journal(context, "trace_set_aside", condition=condition, replica=replica, attempt=past["starts"],
                        trace_set_aside=moved)
            if past["starts"] >= MAX_ATTEMPTS:
                seeds = context.dlib.seeds_for(context.plan.namespace, condition, replica)
                record = failure_record(
                    context.plan, condition, family, replica, INFRASTRUCTURE_PHASE, None, seeds,
                    {"infrastructure_status": UNRESOLVED,
                     "infrastructure_status_meaning": INTERRUPTED_RUN_RULES["unresolved_status"],
                     "attempts": past["starts"], "max_infrastructure_retries": MAX_INFRASTRUCTURE_RETRIES,
                     "attempts_started": list(past["started"]),
                     "interruptions": [{"attempt": number, **past["interrupted"][number]}
                                       for number in sorted(past["interrupted"], key=lambda n: (n is None, n))]})
                commit_run(context, record, None)
                journal(context, "unresolved", condition=condition, replica=replica, attempts=past["starts"],
                        infrastructure_status=UNRESOLVED, detected_by="prepare",
                        record_sha256=pc.sha256_file(context.record_path(condition, replica)))
                continue
            journal(context, "restart", condition=condition, replica=replica, attempt=past["starts"] + 1)
        pending.append((condition, family, replica, past["starts"] + 1))
    return pending


class _StopRequest:
    """SIGINT or SIGTERM once: finish the run in progress, then stop. A second signal acts as usual.

    A deliberate stop then costs no attempt: a run is either stored or was never started.
    """

    def __init__(self) -> None:
        self.requested: int | None = None
        self._previous: dict[int, Any] = {}

    def __enter__(self) -> "_StopRequest":
        if threading.current_thread() is threading.main_thread():
            for number in (signal.SIGINT, signal.SIGTERM):
                self._previous[number] = signal.signal(number, self._handle)
        return self

    def _handle(self, number: int, _frame: Any) -> None:
        self.requested = number
        for which, handler in self._previous.items():
            signal.signal(which, handler)
        self._previous = {}

    def __exit__(self, *_: Any) -> None:
        for which, handler in self._previous.items():
            signal.signal(which, handler)
        self._previous = {}


def work(context: pc.Context, engine_spec: str, jobs: list[tuple], planned: Mapping | None = None) -> int:
    """Run the given planned runs in this process, in order. Returns the number of runs attempted.

    An infrastructure interruption ends the work of this process: the cause (a full disk, no
    memory) would meet the next run too, and every start counts as an attempt.
    """
    engine = load_engine(engine_spec, context)
    done = 0
    with _StopRequest() as stop:
        for job in jobs:
            if stop.requested is not None:
                journal(context, "stop_requested", signal=int(stop.requested), runs_not_started=len(jobs) - done)
                break
            condition, family, replica = job[:3]
            run_one(context, engine, condition, family, replica, job[3] if len(job) > 3 else 1, planned)
            done += 1
    return done


def generate(context: pc.Context, *, engine_spec: str = "production", workers: int | None = None,
             worker_command: list[str] | None = None) -> dict:
    """Run every planned run that has no record yet, then write the run manifest and the failure list.

    workers = 0 runs in this process. Otherwise `workers` operating-system processes are started
    with `worker_command` plus a shard argument; shard k takes every workers-th pending run, so
    the shards have the same mix of conditions. Results do not depend on the number of workers.
    Only one `generate` may work on an output tree at a time.
    """
    manifest = require_generation_allowed(context)
    if workers is None:
        workers = manifest["workers"] if manifest else 0
    context.root.mkdir(parents=True, exist_ok=True)
    lock = os.open(context.root / ".generate.lock", os.O_RDWR | os.O_CREAT, 0o644)
    try:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise pc.PilotError("another generate is working on this output tree") from None
        journal(context, "generate_begin", workers=workers, policy=POLICY, environment_sha256=environment_sha256(),
                software=pc.software_versions(),
                pre_generation_manifest_sha256=(pc.sha256_file(context.root / "pre_generation_manifest.json")
                                                if manifest is not None else None))
        pending = prepare(context)
        started = time.time()
        interruption: InfrastructureInterruption | None = None
        codes: list[int] = []
        stopped = False
        if workers == 0 or not pending:
            try:
                stopped = work(context, engine_spec, pending, planned_inputs(manifest)) < len(pending)
            except InfrastructureInterruption as error:
                interruption = error
        else:
            if not worker_command:
                raise pc.PilotError("a worker command is required to start worker processes")
            pc.atomic_write(context.root / ".pending.json", pc.canonical_json([list(p) for p in pending]), overwrite=True)
            environment = {**os.environ, **{name: "1" for name in pc.THREAD_VARIABLES}}
            with _StopRequest() as stop:    # the workers get the signal too and finish their run in progress
                # The workers inherit the lock: the tree stays locked until the last of them has ended.
                processes = [subprocess.Popen([*worker_command, "--shard", f"{k}/{workers}"], env=environment,
                                              pass_fds=(lock,)) for k in range(workers)]
                codes = [process.wait() for process in processes]
                stopped = stop.requested is not None
            for k, code in enumerate(codes):
                if code:
                    journal(context, "worker_exit", shard=f"{k}/{workers}", code=code)
        summary = finalize(context, workers=workers, seconds=time.time() - started, manifest=manifest)
        summary["stopped_on_request"] = stopped
        journal(context, "generate_end", counts=summary["counts"], unresolved=summary["unresolved_infrastructure"],
                stopped_on_request=stopped)
    finally:
        os.close(lock)
    if interruption is not None:
        raise interruption
    if any(codes):
        raise pc.PilotError(f"worker processes ended with errors {codes}; stored runs are intact, see the journal")
    return summary


def shard(pending: list[tuple[str, str, int]], specification: str) -> list[tuple[str, str, int]]:
    k, _, n = specification.partition("/")
    return pending[int(k)::int(n)]


def finalize(context: pc.Context, *, workers: int, seconds: float | None, manifest: Mapping | None) -> dict:
    """results/run_manifest.json and results/failures.jsonl from what is stored (protocol.json storage.manifest)."""
    entries, failures = [], []
    history = attempt_history(journal_events(context))
    blank = {"starts": 0, "interrupted": {}, "environments": []}
    for condition, family, replica in context.plan.runs():
        state = classify(context, condition, replica)
        past = history.get((condition, replica), blank)
        trace_rel, record_rel = pc.run_paths(condition, replica)
        entry = {"condition": condition, "family": family, "replica": replica, "split": context.plan.split(condition, replica),
                 "seeds": context.dlib.seeds_for(context.plan.namespace, condition, replica), "state": state["state"],
                 "reason": state["reason"], "status": state.get("status"),
                 "attempts": past["starts"],
                 "interruptions": [{"attempt": number, **past["interrupted"][number]}
                                   for number in sorted(past["interrupted"], key=lambda n: (n is None, n))],
                 "attempt_environments": list(past["environments"]), "infrastructure_status": None,
                 "record": {"path": record_rel.as_posix(), "sha256": state.get("record_sha256")},
                 "trace": {"path": trace_rel.as_posix(), "sha256": state.get("trace_sha256")}}
        if state["state"] in ("COMPLETE", "FAILED_RECORD"):
            record = pc.read_json(context.record_path(condition, replica))
            entry.update(array_hashes=record.get("array_hashes"), failure=record.get("failure"),
                         completed_sweeps=record.get("completed_sweeps"),
                         completed_battery_checkpoints=record.get("completed_battery_checkpoints"),
                         completed_application_checkpoints=record.get("completed_application_checkpoints"))
            unresolved = (record.get("failure") or {}).get("phase") == INFRASTRUCTURE_PHASE
            if unresolved:
                entry["infrastructure_status"] = UNRESOLVED
            for item in record.get("failures") or []:
                failures.append({"condition": condition, "replica": replica,
                                 "disposition": ("unresolved: no result after the allowed attempts; the cause is not "
                                                 "established" if unresolved else "recorded, not replaced"), **item})
        else:
            failures.append({"condition": condition, "replica": replica, "phase": "storage", "state": state["state"],
                             "reason": state["reason"], "disposition": "no valid record"})
        entries.append(entry)
    counts = {name: sum(1 for e in entries if e["state"] == name) for name in STATES}
    unresolved_runs = [f"{e['condition']}/{e['replica']:05d}" for e in entries if e["infrastructure_status"] == UNRESOLVED]
    infrastructure = {"policy": POLICY, "max_infrastructure_retries": MAX_INFRASTRUCTURE_RETRIES,
                      "runs_with_more_than_one_attempt": sum(1 for e in entries if e["attempts"] > 1),
                      "interruptions": sum(len(e["interruptions"]) for e in entries), "unresolved": unresolved_runs,
                      "environments": sorted({h for e in entries for h in e["attempt_environments"]})}
    body = {"schema": "rnglab-v019-d-pilot-run-manifest-v1", "namespace": context.plan.namespace, "plan_kind": context.plan.kind,
            "infrastructure": infrastructure,
            "protocol_sha256": pc.protocol_sha256(), "design_sha256": pc.sha256_file(HERE / "design.json"),
            "pre_generation_manifest_sha256": (pc.sha256_file(context.root / "pre_generation_manifest.json")
                                               if manifest is not None else None),
            "software": pc.software_versions(), "workers": workers,
            "thread_variables": {name: os.environ.get(name) for name in pc.THREAD_VARIABLES},
            "generation_wall_seconds_this_invocation": seconds, "counts": counts, "planned": len(entries), "runs": entries}
    digest = pc.atomic_write(context.result("run_manifest.json"), pc.canonical_json(body), overwrite=True)
    pc.atomic_write(context.result("failures.jsonl"), pc.canonical_jsonl(failures), overwrite=True)
    return {"counts": counts, "planned": len(entries), "failures": len(failures), "run_manifest_sha256": digest,
            "unresolved_infrastructure": unresolved_runs}


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
