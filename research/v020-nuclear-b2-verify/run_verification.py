"""B2 verification against analytic laws and the frozen Layer A engine: the frozen study in protocol.json.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/run_verification.py [--smoke]

The script refuses to draw if a frozen file differs from protocol.sha256. The
main run refuses unless the tree is committed (only files of earlier attempts in
results/attempts/ may be uncommitted, and no ignored file other than bytecode
caches, egg-info metadata and smoke/ may sit in src/ or in this folder), PYTHONPYCACHEPREFIX names
a directory outside the repository, results/verification.json is absent and the
committed record calibration matches (record_calibration). It runs as a recorded
attempt (begin_attempt: reserved atomically before any draw, finished as
completed, failed or terminated; manifests go to results/attempts/verification-<k>/). --smoke uses the smoke namespace and reduced sizes, writes to smoke/
(not in Git) and is excluded from every result; before the freeze it runs on
the files' current hashes and records that it was not frozen.

Every test is a function of observations only (battery, path_tests,
ensemble_tests, comparison_tests, big_tests). calibration.py evaluates the
same functions on null data drawn from the exact laws. A test that cannot be
computed is recorded with its error and listed as a failed item; it never
reads as a passing test.
"""
from __future__ import annotations

import argparse
import copy
import datetime
import fcntl
import hashlib
import json
import math
import os
import signal
import socket
import subprocess
import sys
import time
from decimal import Context, Decimal, localcontext
from pathlib import Path

import numpy as np
from scipy.stats import binomtest

from rnglab.nuclear import (
    DecayScheme, NumpyPCG64Source, WordStream, decay_probabilities, simulate_decays, state_probabilities,
)
from rnglab.nuclear.aggregate import aggregate_manifest, simulate_aggregate, simulate_aggregate_reference
from rnglab.nuclear.binomial import PART_MAX, binomial_inverse
from rnglab.nuclear.checks import binomial_check, decay_invariants, holm, multinomial_check
from rnglab.nuclear.comparison import (
    allocation_check, binomial_z_check, consecutive_cells, fixed_cells_chisquare, homogeneity_check,
    planned_path_cells, pooled_cells, sample_covariance, scaled_chisquare_variance_check, summed_chisquare,
    variance_ratio_f_check,
)
from rnglab.nuclear.decay import Mode
from rnglab.nuclear.manifest import code_identity, environment, file_sha256
from rnglab.nuclear.moments import state_count_moments, state_sum_moments, two_time_sum_moments
from rnglab.nuclear.paths import bounds_hold, interval_law, tree_probabilities
from rnglab.seeds import namespace_seed
from rnglab14.sources import Xoshiro256PlusPlusSource

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
FROZEN = ("protocol.json", "design.json", "design.py", "run_verification.py", "calibration.py")
PHOTON_MIN_EXPECTED = 10.0      # per-interval photon cells: min(yield, 1 - yield) * expected decays
COMPARISON_MIN_EXPECTED = 5.0   # cells of the Layer A / B2 homogeneity tests, per engine


def frozen_hashes(*, smoke: bool = False) -> dict:
    if not (HERE / "protocol.sha256").exists():
        if smoke:
            return {"frozen": False, "files": {name: file_sha256(HERE / name) for name in FROZEN}}
        raise SystemExit("protocol.sha256 is missing: the protocol is not frozen; refusing to draw")
    listed = {}
    for line in (HERE / "protocol.sha256").read_text().splitlines():
        if line.strip():
            digest, name = line.split()
            listed[name] = digest
    if set(listed) != set(FROZEN):
        raise SystemExit("protocol.sha256 must list exactly the frozen files")
    for name in FROZEN:
        if file_sha256(HERE / name) != listed[name]:
            raise SystemExit(f"{name} differs from protocol.sha256; refusing to draw")
    return {"frozen": True, "files": listed}


ATTEMPT_SCHEMA = "rnglab-b2-attempt-v2"
STUDY = "research/v020-nuclear-b2-verify/"
RESULTS = STUDY + "results/"
ATTEMPTS = RESULTS + "attempts/"
CODE_KEYS = ("git_commit", "module_files_sha256", "source_files_sha256")
C1_LIMIT = 0.02


class Terminated(BaseException):
    """Raised in the main process on SIGTERM (a SLURM time limit, scancel), so that the attempt is recorded."""


def _raise_terminated(signum, frame):
    raise Terminated(f"signal {signum}")


def _utc() -> str:
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def _synced(path: Path, text: str, mode: str) -> None:
    with open(path, mode) as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())


def _sync_dir(folder: Path) -> None:
    """fsync a directory, so that a link or rename in it survives a crash of the node."""
    fd = os.open(folder, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def write_new(path: Path, text: str) -> None:
    """Publish a file atomically and never over an existing one (synced temporary file, then a hard link)."""
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    _synced(tmp, text, "x")
    try:
        os.link(tmp, path)  # FileExistsError if path exists: complete or absent, never overwritten
    finally:
        tmp.unlink()
    _sync_dir(path.parent)


def _replace(path: Path, text: str) -> None:
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    _synced(tmp, text, "w")
    os.replace(tmp, path)
    _sync_dir(path.parent)


def _git(repo: Path, *args: str, check: bool = True) -> subprocess.CompletedProcess:
    return subprocess.run(["git", "-C", str(repo), *args], capture_output=True, check=check)


def source_files_sha256(repo: Path) -> str:
    """SHA-256 over every tracked file under src/ (its path and its content's SHA-256): all the repository code
    a run can import, the seed derivation and rnglab14's multiple-test correction included."""
    digest = hashlib.sha256()
    for name in sorted(n for n in _git(repo, "ls-files", "-z", "--", "src").stdout.decode().split("\0") if n):
        digest.update(f"{name}\0{file_sha256(repo / name)}\n".encode())
    return digest.hexdigest()


def record_code(repo: Path) -> dict:
    """code_identity plus source_files_sha256: the code a record run is bound to, taken when it starts."""
    return {**code_identity(repo), "source_files_sha256": source_files_sha256(repo)}


def require_same_code(start: dict, repo: Path, what: str) -> None:
    """Refuse to publish a record run whose commit or code changed while it ran (checked when it ends)."""
    end = record_code(repo)
    changed = [key for key in CODE_KEYS if end.get(key) != start.get(key)]
    if changed:
        raise RuntimeError(f"the code of the {what} changed while it ran ({', '.join(changed)}); its results would "
                           "describe no single commit")


def require_pycache_outside(repo: Path) -> None:
    """Record runs read no bytecode from the tree: PYTHONPYCACHEPREFIX must name a directory outside it."""
    prefix = sys.pycache_prefix
    if not prefix or Path(prefix).resolve().is_relative_to(repo.resolve()):
        raise SystemExit("record runs need PYTHONPYCACHEPREFIX set to a directory outside the repository (a new, "
                         "empty one is recommended; protocol.json record_procedure.threads), so that no bytecode "
                         "in the tree is used")


def record_tree(repo: Path) -> dict:
    """Every change in the tree, split into files of earlier attempts (untracked, inside ATTEMPTS) and the rest.

    The rest includes ignored files in src/ and in this folder, other than bytecode caches, the egg-info
    metadata of an editable install (src/<name>.egg-info/, which cannot be imported) and smoke/: an ignored
    module beside the code (a sitecustomize.py, a module shadowing one of the study's) could change what a
    record run computes without appearing anywhere in Git."""
    out = _git(repo, "status", "--porcelain=v1", "-z", "--untracked-files=all", "--ignored=matching").stdout.decode()
    attempt_files, other, allowed = [], [], []
    for entry in out.split("\0"):
        if not entry:
            continue
        if entry.startswith("!! "):
            path = entry[3:]
            if not path.startswith(("src/", STUDY)):
                continue
            if ("__pycache__/" in path or path.startswith(STUDY + "smoke/")
                    or (path.startswith("src/") and path.split("/")[1].endswith(".egg-info"))):
                allowed.append(entry)       # recorded, so that a result shows what was there
                continue
        (attempt_files if entry.startswith("?? " + ATTEMPTS) else other).append(entry)
    return {"attempt_files": attempt_files, "other_changes": other, "ignored_allowed": allowed}


def require_record_tree(repo: Path, script: str) -> dict:
    """A record run starts from a committed tree; only files that earlier attempts left in ATTEMPTS may be
    uncommitted (they are evidence, committed with the results), so an interrupted attempt can be rerun."""
    tree = record_tree(repo)
    if tree["other_changes"]:
        raise SystemExit(f"the record {script} run needs a committed tree; only files of earlier attempts in "
                         f"{ATTEMPTS} may be uncommitted, and no ignored file may sit in src/ or {STUDY}, not: "
                         + "; ".join(tree["other_changes"][:5]))
    return tree


def _same_file(handle, path: Path) -> bool:
    """Whether path still names the file that handle has open."""
    held = os.fstat(handle.fileno())
    try:
        named = os.stat(path)
    except FileNotFoundError:
        return False
    return (named.st_dev, named.st_ino) == (held.st_dev, held.st_ino)


def _lock(folder: Path, script: str):
    """An exclusive flock on attempts/<script>.lock, held on the file the path still names once it is locked
    (a lock taken on a file that was deleted or replaced in the meantime would exclude nobody)."""
    lock_path = folder / f"{script}.lock"
    for _ in range(3):
        lock = open(lock_path, "a+")
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            lock.close()
            raise SystemExit(f"another {script} attempt is running: it holds {lock_path}") from None
        except OSError as error:
            lock.close()
            raise SystemExit(f"cannot lock {lock_path} ({error}); record runs need a file system with flock") from None
        if _same_file(lock, lock_path):
            return lock
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()
    raise SystemExit(f"{lock_path} keeps being deleted or replaced; never touch the lock file of a record run")


def _process_start(pid: int) -> int | None:
    """A process's start time in clock ticks since boot (field 22 of /proc/<pid>/stat), or None."""
    try:
        stat = Path(f"/proc/{pid}/stat").read_bytes()          # a command name need not be valid UTF-8
        return int(stat[stat.rindex(b")") + 2:].split()[19])    # fields after the command name start at field 3
    except (OSError, ValueError, IndexError):
        return None


def _alive_here(record: dict) -> bool:
    """Whether a record names a process that still runs on this host: the same process id and, where /proc
    gives it, the same start time (a process id reused by another program is not the attempt)."""
    try:
        pid = int(record["pid"])
    except (KeyError, TypeError, ValueError):
        return False
    if record.get("host") != socket.gethostname() or pid == os.getpid():
        return False
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        pass
    start = _process_start(pid)
    if start is None or record.get("process_start") is None:
        return True                     # cannot tell a reused process id apart: the safe answer
    return start == record["process_start"]


def begin_attempt(out_dir: Path, script: str, code: dict, *, deviation: str | None = None) -> dict:
    """Reserve the next attempt of a once-only record run before any random number is drawn.

    The attempt holds an exclusive flock on attempts/<script>.lock for as long as it runs; the kernel
    releases it when the process ends, however it ends. While it is held no other attempt of the script
    can run, so an earlier record still marked "started" belongs to a process that died (a node failure,
    SIGKILL): it is listed as interrupted in the new record. Refused: an attempt that holds the lock; a
    completed attempt (the run is never repeated); an earlier attempt on another commit, unless
    ``deviation`` gives the reason (recorded). Records are published atomically (complete or absent).
    A rerun on the same commit draws the same numbers (fixed seeds): it replays, it is not a new sample.
    """
    folder = out_dir / "attempts"
    folder.mkdir(parents=True, exist_ok=True)
    lock = _lock(folder, script)
    try:
        previous = []
        for path in sorted(folder.glob(f"{script}-*.json")):
            number = path.stem[len(script) + 1:]
            if number.isdigit():
                try:
                    previous.append((int(number), path, json.loads(path.read_text())))
                except (OSError, json.JSONDecodeError) as error:
                    raise SystemExit(f"{path} cannot be read ({error}); records are published atomically, so "
                                     "it was changed by hand: inspect it before any rerun") from None
        previous.sort()
        for number, path, record in previous:
            if record["status"] == "completed":
                raise SystemExit(f"{path} records a completed attempt; the run is never repeated")
            if record["status"] == "started" and _alive_here(record):
                raise SystemExit(f"attempt {number} still runs on this host (process {record['pid']}), although "
                                 f"{script}.lock was free: the lock file was deleted or replaced; never touch it")
            if record["code"].get("git_commit") != code.get("git_commit") and not deviation:
                raise SystemExit(f"attempt {number} ran on commit {record['code'].get('git_commit')}; a rerun on "
                                 "another commit is a protocol deviation: give --deviation '<reason>'")
        number = 1 + max((n for n, _, _ in previous), default=0)
        record = {"schema": ATTEMPT_SCHEMA, "script": script, "attempt": number, "status": "started",
                  "started_utc": _utc(), "host": socket.gethostname(), "pid": os.getpid(),
                  "process_start": _process_start(os.getpid()), "argv": sys.argv,
                  "code": code, "previous_attempts": [{"attempt": n, "status": r["status"],
                                                       "git_commit": r["code"].get("git_commit")}
                                                      for n, _, r in previous],
                  "interrupted_attempts": [n for n, _, r in previous if r["status"] == "started"],
                  "leftover_temporary_files": sorted(p.name for p in folder.glob(".*.tmp")),
                  "deviation": deviation}
        path = folder / f"{script}-{number}.json"
        write_new(path, json.dumps(record, indent=1) + "\n")
        (folder / f"{script}-{number}").mkdir()
        lock.seek(0)
        lock.truncate()
        lock.write(json.dumps({"attempt": number, "host": record["host"], "pid": record["pid"]}) + "\n")
        lock.flush()
    except BaseException:
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()
        raise
    return {"path": path, "dir": folder / f"{script}-{number}", "record": record, "lock": lock}


def complete_published_attempt(out_dir: Path, script: str, result_file: Path) -> bool:
    """Finish the record of an attempt that published its result but ended before recording that it completed
    (killed, or failed in the last steps of publication). The result names its attempt, and only that attempt
    can have published it (under the lock, after its code check); that attempt's record must show it started,
    terminated or failed on the result's commit. Done under the attempt lock. Returns True when a record was
    finished."""
    if not result_file.exists():
        return False
    data = json.loads(result_file.read_text())
    number = (data.get("attempt") or {}).get("number")
    folder = out_dir / "attempts"
    path = folder / f"{script}-{number}.json"
    if isinstance(number, bool) or not isinstance(number, int) or not path.is_file():
        raise SystemExit(f"{result_file} names no attempt record of the {script}; inspect it")
    lock = _lock(folder, script)
    try:
        record = json.loads(path.read_text())
        if record.get("status") == "completed":
            return False
        if (record.get("status") not in ("started", "terminated", "failed")
                or record["code"].get("git_commit") != (data.get("code") or {}).get("git_commit")):
            raise SystemExit(f"{path} ({record.get('status')}) does not describe {result_file}; inspect both")
        _replace(path, json.dumps({**record, "status": "completed", "exit_status": None, "result": result_file.name,
                                   "result_sha256": file_sha256(result_file),
                                   "completed_after_interruption": {"status_before": record["status"],
                                                                    "noted_utc": _utc()}},
                                  indent=1, default=str) + "\n")
        return True
    finally:
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()


def end_attempt(attempt: dict, status: str, **fields) -> None:
    attempt["record"] = {**attempt["record"], "status": status, "ended_utc": _utc(), **fields}
    _replace(attempt["path"], json.dumps(attempt["record"], indent=1, default=str) + "\n")


def release_attempt(attempt: dict) -> None:
    lock = attempt.pop("lock", None)
    if lock is not None:
        fcntl.flock(lock.fileno(), fcntl.LOCK_UN)
        lock.close()


def run_attempt(attempt: dict, body, result_file: Path) -> int:
    """Run body() as the reserved attempt; record how it ended (SIGTERM included) and the result's SHA-256."""
    handler = signal.signal(signal.SIGTERM, _raise_terminated)
    try:
        try:
            status = body()
        except BaseException as error:
            end_attempt(attempt, "terminated" if isinstance(error, (Terminated, KeyboardInterrupt)) else "failed",
                        error=f"{type(error).__name__}: {error}")
            raise
        end_attempt(attempt, "completed", exit_status=status, result=result_file.name,
                    result_sha256=file_sha256(result_file))
        return status
    finally:
        signal.signal(signal.SIGTERM, handler)
        release_attempt(attempt)


def clopper_pearson(k: int, n: int) -> list[float]:
    """The exact 95% interval of k successes in n trials; an impossible count is refused (scipy would return NaN)."""
    if not 0 <= int(k) <= int(n) or int(n) <= 0:
        raise ValueError(f"{k} of {n} is not a binomial count")
    interval = binomtest(int(k), int(n)).proportion_ci(confidence_level=0.95, method="exact")
    return [float(interval.low), float(interval.high)]


def r0_reading(families: dict, batteries: int, r1_event: dict | None = None, limit: float = C1_LIMIT) -> dict:
    """R0 as reading_rules.R0 words it. Criterion C1 is per family; the intervals are recomputed from the counts
    (stored intervals are not read), and an upper bound above the limit is not evidence of an excess. R1 fails
    when either family fails, so R1's own false-FAIL rate (the union) is reported beside them, not as a criterion."""
    out = {}
    for family, f in families.items():
        low, high = clopper_pearson(f["failed_batteries"], batteries)
        if high <= limit:
            statement = f"C1 holds: the false-FAIL rate is below {limit} (95% upper bound {high:.4f})"
        elif low > limit:
            statement = f"C1 fails: the false-FAIL rate exceeds {limit} (95% lower bound {low:.4f})"
        else:
            statement = (f"C1 fails: the calibration did not establish that the false-FAIL rate is below {limit} "
                         f"(95% interval {low:.4f} to {high:.4f})")
        out[family] = {"false_fail_rate": f["failed_batteries"] / batteries, "false_fail_ci95": [low, high],
                       "C1_holds": high <= limit, "excess_established": low > limit,
                       "statement": f"{family}: {statement}"}
    if r1_event is not None:
        k = r1_event["failed_batteries"]
        low, high = clopper_pearson(k, batteries)
        out["R1_event"] = {"failed_batteries": k, "false_fail_rate": k / batteries, "false_fail_ci95": [low, high],
                           "statement": f"R1 (a FAIL in NC1 or NC2): a correct engine FAILs in {k} of {batteries} "
                                        f"null batteries (rate {k / batteries:.4f}, 95% interval {low:.4f} to "
                                        f"{high:.4f}); reported beside C1, which is per family"}
    return out


def _committed(repo: Path, path: Path) -> bool:
    """Tracked, and identical in the working tree, the index and HEAD."""
    relative = path.resolve().relative_to(repo.resolve()).as_posix()
    return (_git(repo, "ls-files", "--error-unmatch", "--", relative, check=False).returncode == 0
            and _git(repo, "diff", "--quiet", "HEAD", "--", relative, check=False).returncode == 0)


def _mapping(value) -> dict:
    return value if isinstance(value, dict) else {}


def _count(value) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def record_calibration(protocol: dict, frozen: dict, code: dict, results: Path = HERE / "results",
                       repo: Path = REPO, design: dict | None = None) -> dict:
    """The completed record calibration that R0 reads. The verification refuses to start without one that ran
    on the same frozen files, module and source files and numerical libraries; whose attempt record (named as
    the record procedure names it) shows it completed with this file; that is committed, with its attempt
    record, exactly as it is; whose commit is an ancestor of this one, with nothing but results/ changed since;
    whose replay check found no mismatch; with exactly the families NC1 and NC2, each with the planned number of
    tests (design.json) and an interval, rate and C1 flag that follow from its count; a C1 that follows from
    them; and an R1 event whose count, rate and interval fit the family counts and the batteries. A malformed
    file is refused with a message, never raised on."""
    path = results / "calibration.json"
    if not path.exists():
        raise SystemExit("results/calibration.json is missing: run the record calibration first and commit it "
                         "(protocol.json record_procedure)")
    try:
        data = json.loads(path.read_text())
    except ValueError as error:
        raise SystemExit(f"results/calibration.json cannot serve R0: it is not JSON ({error})") from None
    if not isinstance(data, dict):
        raise SystemExit("results/calibration.json cannot serve R0: it is not a JSON object")
    spec, now, problems = protocol["calibration"], environment(), []
    recorded = _mapping(data.get("code"))
    if data.get("development") is not False:
        problems.append("it is a development run")
    if data.get("namespace") != spec["namespace"] or data.get("batteries") != spec["batteries"]:
        problems.append("its namespace or number of batteries differs from protocol.json")
    if data.get("protocol_files") != frozen:
        problems.append("it ran with other frozen files")
    if recorded.get("module_files_sha256") != code.get("module_files_sha256"):
        problems.append("it ran with other module files")
    if recorded.get("source_files_sha256") != code.get("source_files_sha256"):
        problems.append("it ran with other source files")
    commit = recorded.get("git_commit")
    if (not isinstance(commit, str) or not commit
            or _git(repo, "merge-base", "--is-ancestor", commit, "HEAD", check=False).returncode != 0):
        problems.append(f"its commit {commit} is not an ancestor of this one")
    else:
        changed = [n for n in _git(repo, "diff", "--name-only", "--no-renames", "-z", commit,
                                   "HEAD").stdout.decode().split("\0")
                   if n and not n.startswith(RESULTS)]
        if changed:
            problems.append("files outside results/ changed since it ran: " + ", ".join(changed[:5]))
    if _mapping(data.get("replay_check")).get("mismatches"):
        problems.append("its replay check found mismatches")
    batteries = data.get("batteries")
    if not _count(batteries) or batteries <= 0:
        problems.append("its number of batteries is not a positive integer")
        batteries = 0
    families = data.get("families")
    if not isinstance(families, dict) or not all(isinstance(f, dict) for f in families.values()):
        problems.append("its families are not a mapping of family records")
        families = {}
    if set(families) != {"NC1", "NC2"}:
        problems.append(f"its families are {sorted(families)}, not NC1 and NC2")
    planned = {f: len(t) for f, t in _mapping(_mapping(design).get("planned_tests")).items()}
    counts = []
    for family, f in families.items():
        k = f.get("failed_batteries")
        if not batteries or not _count(k) or not 0 <= k <= batteries:
            problems.append(f"its count for {family} is missing or impossible")
            continue
        counts.append(k)
        interval = clopper_pearson(k, batteries)
        if f.get("false_fail_ci95") != interval:
            problems.append(f"its stored interval for {family} does not follow from its counts")
        if f.get("false_fail_rate") != k / batteries or f.get("C1_holds") != (interval[1] <= C1_LIMIT):
            problems.append(f"its stored rate or C1 for {family} does not follow from its counts")
        if planned and f.get("tests") != planned.get(family):
            problems.append(f"its number of tests for {family} ({f.get('tests')}) differs from design.json "
                            f"({planned.get(family)})")
    if data.get("C1") is not all(f.get("C1_holds") is True for f in families.values()):
        problems.append("its C1 does not follow from its families")
    r1 = _mapping(data.get("R1_event"))
    k = r1.get("failed_batteries")
    if (len(counts) != 2 or not _count(k) or r1.get("batteries") != batteries
            or not max(counts) <= k <= min(sum(counts), batteries)
            or r1.get("false_fail_ci95") != clopper_pearson(k, batteries) or r1.get("false_fail_rate") != k / batteries):
        problems.append("its R1 event is missing, or its count, rate or interval does not fit its family counts "
                        "and batteries")
    environment_recorded = _mapping(data.get("environment"))
    for key in ("python", "numpy", "scipy"):
        if environment_recorded.get(key) != now.get(key):
            problems.append(f"it ran with {key} {environment_recorded.get(key)}, this run has {now.get(key)}")
    attempt = _mapping(data.get("attempt"))
    number, record_path = attempt.get("number"), None
    if not _count(number) or attempt.get("record") != f"attempts/calibration-{number}.json":
        problems.append("its attempt record is not named as the record procedure names it")
    else:
        record_path = results / "attempts" / f"calibration-{number}.json"
        try:
            record = json.loads(record_path.read_text()) if record_path.is_file() else None
        except ValueError:
            record = {}
        if record is None:
            problems.append("its attempt record is missing")
        elif (not isinstance(record, dict) or record.get("status") != "completed"
              or record.get("result_sha256") != file_sha256(path)):
            problems.append("its attempt record does not show it completed with this file")
    for item in (path, record_path):
        if item is not None and item.is_file() and not _committed(repo, item):
            problems.append(f"{item.name} is not committed as it is")
    if problems:
        raise SystemExit("results/calibration.json cannot serve R0: " + "; ".join(problems))
    return {"file": "results/calibration.json", "sha256": file_sha256(path), "attempt": number,
            "commit": commit, "deviation": attempt.get("deviation"), "replay_check": data.get("replay_check"),
            "tree": data.get("tree"), "R0": r0_reading(families, batteries, r1)}


_MEMO: dict = {}


def memo(function, scheme, *args):
    """Deterministic references (matrix-exponential laws, exact moments), computed once per argument set.

    The battery runs once in the study and many times in calibration.py; the values are the same either
    way. Arrays come back read-only.
    """
    key = (function.__module__, function.__qualname__, id(scheme), repr(args))
    hit = _MEMO.get(key)
    if hit is not None and hit[0] is scheme:
        return hit[1]
    value = function(scheme, *args)
    if isinstance(value, np.ndarray):
        value.setflags(write=False)
    _MEMO[key] = (scheme, value)
    return value


def evaluate(function, *args, **kwargs) -> dict:
    """A test's result, or {'error': ...} when it cannot be computed (a failed item, never a pass)."""
    try:
        return function(*args, **kwargs)
    except (ValueError, ArithmeticError) as error:
        return {"error": f"{type(error).__name__}: {error}"}


def make_source(generator: str, seed: int, label: str):
    if generator == "numpy-pcg64":
        return NumpyPCG64Source(seed, label=label)
    if generator == "xoshiro256pp":
        return Xoshiro256PlusPlusSource(seed=seed)
    raise ValueError(generator)


def groups_of(scheme: DecayScheme, initial: dict) -> list[tuple[str, list[int]]]:
    out = []
    for start in initial:
        reach, frontier = {scheme.index[start]}, [scheme.index[start]]
        while frontier:
            i = frontier.pop()
            for channel in scheme.nuclides[i].channels:
                j = scheme.index[channel.daughter]
                if j not in reach:
                    reach.add(j)
                    frontier.append(j)
        out.append((start, sorted(reach)))
    return out


def bookkeeping(run) -> dict[str, bool]:
    scheme = run.plan.scheme
    pops, decays, photons = run.populations, run.decays, run.photons
    total = sum(run.initial.values())
    inflow = np.zeros_like(pops[1:])
    yields = np.zeros(decays.shape[1:])
    neutron = np.zeros(decays.shape[1:], dtype=bool)
    for j, nuclide in enumerate(scheme.nuclides):
        for k, channel in enumerate(nuclide.channels):
            inflow[:, scheme.index[channel.daughter]] += decays[:, j, k]
            yields[j, k] = channel.photon_yield
            neutron[j, k] = channel.mode is Mode.BETA_DELAYED_NEUTRON
    outflow = decays.sum(axis=2)
    return {"conserved_and_nonnegative": bool(np.all(pops >= 0) and set(pops.sum(axis=1).tolist()) == {total}),
            "flow_balance": bool(np.array_equal(np.diff(pops, axis=0), inflow - outflow)),
            "outflow_bounded": bool(np.all(outflow <= pops[:-1] + inflow)),
            "photons_bounded": bool(np.all((photons >= 0) & (photons <= decays))),
            "photon_yield_extremes": bool(np.all(photons[:, yields == 0.0] == 0)
                                          and np.array_equal(photons[:, yields == 1.0], decays[:, yields == 1.0])),
            "neutrons_equal_beta_n_decays": bool(np.array_equal(run.neutrons, decays[:, neutron].sum(axis=1)))}


def b2_consumption(run, stream) -> dict:
    words = run.plan.deltas.size * run.plan.words_per_interval
    return {"words_expected": words, "words_read": run.words_read,
            "ok": run.words_read == words == stream.words_read and stream.source.consumed_bytes == 8 * words}


def chi2(counts, probabilities) -> dict:
    """multinomial_check, recording the pre-declared 'not applicable' case (fewer than two cells after pooling)."""
    try:
        return multinomial_check(counts, probabilities)
    except ValueError as error:
        if "fewer than two cells" in str(error):
            return {"not_applicable": str(error)}
        raise


def channel_label(nuclide, channel) -> str:
    return f"{nuclide.name}/{channel.mode.value}->{channel.daughter}"


def photon_interval_cells(ref: DecayScheme, case: dict) -> list[tuple[int, int, list[list[int]]]]:
    """(nuclide, channel, cells) of the per-interval photon tests, from expected decays only.

    For every channel with photon yield y strictly between 0 and 1, consecutive grid intervals are
    pooled until min(y, 1 - y) times the expected decays of the cell reaches PHOTON_MIN_EXPECTED;
    channels left with one cell get no per-interval test (the test of the totals covers them).
    """
    times = case["times"]
    expected: dict[tuple[int, int], np.ndarray] = {}
    for start, n in case["initial"].items():
        q = memo(decay_probabilities, ref, start, times)
        for j, nuclide in enumerate(ref.nuclides):
            for k, (channel, branch) in enumerate(zip(nuclide.channels, nuclide.channel_probabilities())):
                if 0.0 < channel.photon_yield < 1.0:
                    expected[(j, k)] = expected.get((j, k), 0.0) + n * branch * np.diff(q[:, j])
    out = []
    for (j, k), decays in sorted(expected.items()):
        y = ref.nuclides[j].channels[k].photon_yield
        cells = consecutive_cells(min(y, 1.0 - y) * decays, PHOTON_MIN_EXPECTED)
        if len(cells) >= 2:
            out.append((j, k, cells))
    return out


def neutron_law(ref: DecayScheme, case: dict, start: str) -> tuple[int, int, int, np.ndarray]:
    """BT-NEUTRON of a start nuclide: its index, its delayed-neutron channel, its initial count, and the cells
    of the test (a neutron in each interval, or none), whose counts are multinomial under the null."""
    i, n = ref.index[start], case["initial"][start]
    q = memo(decay_probabilities, ref, start, case["times"])[:, i]
    channels = ref.nuclides[i].channels
    k = next(c for c, ch in enumerate(channels) if ch.mode is Mode.BETA_DELAYED_NEUTRON)
    pn = ref.nuclides[i].channel_probabilities()[k]
    return i, k, n, np.append(pn * np.diff(q), 1.0 - pn * q[-1])


def photon_timing_test(ref: DecayScheme, j: int, k: int, cells, decays: np.ndarray, photons: np.ndarray) -> dict:
    """BT-PHOTON-T of one channel: its photons in the pre-declared cells of intervals, given its decays."""
    return evaluate(allocation_check, photons[:, j, k], decays[:, j, k], ref.nuclides[j].channels[k].photon_yield,
                    cells)


def battery(key: str, run, ref: DecayScheme, case: dict, prefix: str) -> dict:
    out = {}
    times = case["times"]
    pops, decays, photons = run.populations, run.decays, run.photons
    for start, members in groups_of(ref, case["initial"]):
        n = case["initial"][start]
        p = memo(state_probabilities, ref, start, times)
        for it, t in enumerate(times[1:], start=1):
            out[f"{prefix} {key} BT-STATE {start} t={t!r}"] = evaluate(chi2, pops[it, members], p[it, members])
        q_end = memo(decay_probabilities, ref, start, times[-1])[0]
        for j in members:
            if not ref.nuclides[j].stable:
                out[f"{prefix} {key} BT-DECAY {ref.nuclides[j].name}"] = evaluate(
                    binomial_check, decays[:, j, :].sum(), n, q_end[j])
        i = ref.index[start]
        q = memo(decay_probabilities, ref, start, times)[:, i]
        per_interval = decays[:, i, :].sum(axis=1)
        out[f"{prefix} {key} BT-BINS {start}"] = evaluate(
            chi2, np.append(per_interval, n - per_interval.sum()), np.append(np.diff(q), 1.0 - q[-1]))
        if key == "S4":
            i, k, n, probabilities = neutron_law(ref, case, start)
            neutrons = decays[:, i, k]
            out[f"{prefix} {key} BT-NEUTRON {start}"] = evaluate(
                chi2, np.append(neutrons, n - neutrons.sum()), probabilities)
    for j, nuclide in enumerate(ref.nuclides):
        total_j = int(decays[:, j, :].sum())
        if len(nuclide.channels) >= 2:
            out[f"{prefix} {key} BT-CHANNEL {nuclide.name}"] = evaluate(
                binomial_check, decays[:, j, 0].sum(), total_j, nuclide.channel_probabilities()[0])
        for k, channel in enumerate(nuclide.channels):
            if 0.0 < channel.photon_yield < 1.0:
                out[f"{prefix} {key} BT-PHOTON {channel_label(nuclide, channel)}"] = evaluate(
                    binomial_check, photons[:, j, k].sum(), decays[:, j, k].sum(), channel.photon_yield)
    for j, k, cells in photon_interval_cells(ref, case):
        nuclide = ref.nuclides[j]
        out[f"{prefix} {key} BT-PHOTON-T {channel_label(nuclide, nuclide.channels[k])}"] = photon_timing_test(
            ref, j, k, cells, decays, photons)
    return out


def path_tests(run, prefix: str, case: dict, min_expected: float) -> tuple[dict, dict]:
    """One test per source: Pearson statistics of the planned intervals on pre-declared cells, summed.

    Given the start count, the category counts of an interval are Multinomial(start count, q(delta))
    whatever happened before; with cells fixed in advance each interval's statistic is asymptotically
    chi-square with constant degrees of freedom given the past, so the sum is asymptotically
    chi-square(summed degrees of freedom). A statistic that cannot be computed is recorded as an
    error, which makes check D10 fail; it never reads as a passing test.
    """
    plan = run.plan
    scheme = plan.scheme
    (start, n0), = case["initial"].items()
    expected = [memo(state_probabilities, scheme, start, t)[0] * n0 for t in case["times"][:-1]]
    cells = planned_path_cells(plan, expected, min_expected=min_expected)
    out = {}
    for i, per_interval in cells.items():
        first, width = int(plan.first_category[i]), len(plan.categories[i])
        name = f"{prefix} S2 BT-PATH {scheme.nuclides[i].name}"
        try:
            parts = [fixed_cells_chisquare(run.path_counts[k, first:first + width],
                                           plan.source_plan(plan.deltas[k], i).law.as_float(),
                                           interval_cells)
                     for k, interval_cells in per_interval.items()]
            out[name] = summed_chisquare(parts)
            out[name]["intervals"] = sorted(per_interval)
            out[name]["interval_statistics"] = [part["statistic"] for part in parts]
        except (ValueError, ArithmeticError) as error:
            out[name] = {"error": f"{type(error).__name__}: {error}"}
    names = {scheme.nuclides[i].name: {str(k): c for k, c in per.items()} for i, per in cells.items()}
    return out, names


def naive_aggregate(scheme: DecayScheme, initial: dict, times, stream: WordStream) -> np.ndarray:
    """The deliberately wrong one-step aggregation (positive control): at most one decay per interval."""
    pops = np.array([initial.get(n.name, 0) for n in scheme.nuclides], dtype=np.int64)
    out = [pops.copy()]
    for delta in np.diff(times):
        new = np.zeros_like(pops)
        for i, nuclide in enumerate(scheme.nuclides):
            if nuclide.stable:
                new[i] += pops[i]
                continue
            stay = math.exp(-math.log(2.0) / nuclide.half_life * delta)
            probabilities = [stay] + [(1.0 - stay) * b for b in nuclide.channel_probabilities()]
            targets = [i] + [scheme.index[c.daughter] for c in nuclide.channels]
            remaining, tail = int(pops[i]), 1.0
            for k in range(len(probabilities) - 1):
                u = stream.uniforms(1)
                p = probabilities[k] / tail if tail > 0 else 0.0
                p = min(max(p, 0.0), 1.0)
                if p > 0.5:
                    x = remaining - int(binomial_inverse([remaining], [1.0 - p], u)[0])
                else:
                    x = int(binomial_inverse([remaining], [p], u)[0])
                new[targets[k]] += x
                remaining -= x
                tail -= probabilities[k]
            new[targets[-1]] += remaining
        pops = new
        out.append(pops.copy())
    return np.array(out)


def layer_a_observations(run, times):
    times = np.asarray(times)
    scheme = run.scheme
    pops = run.populations(times)
    ev = run.events
    interval = np.searchsorted(times, ev["time"], side="left") - 1
    at_zero = int(np.sum(ev["time"] <= times[0]))
    shape = (times.size - 1, len(scheme), max(scheme.max_channels, 1))
    decays = np.zeros(shape, dtype=np.int64)
    photons = np.zeros(shape, dtype=np.int64)
    keep = interval >= 0
    np.add.at(decays, (interval[keep], ev["parent"][keep], ev["channel"][keep]), 1)
    flag = keep & ev["photon"]
    np.add.at(photons, (interval[flag], ev["parent"][flag], ev["channel"][flag]), 1)
    return pops, decays, photons, at_zero


def law_checks(schemes: dict, cases: dict) -> dict:
    """D8. Every reference must be finite: a NaN would compare false and could never raise a maximum,
    so a non-finite reference is recorded and fails the check instead of reading as agreement."""
    worst_tree = worst_layer_a = 0.0
    worst_sum = Decimal(0)
    brackets, digits, laws, error, nonfinite = True, 0, 0, None, []
    try:
        for key, scheme in schemes.items():
            for delta in sorted(set(np.diff(cases[key]["times"]).tolist())):
                for i, nuclide in enumerate(scheme.nuclides):
                    law = interval_law(scheme, i, delta)
                    laws += 1
                    digits = max(digits, law.digits_used)
                    q = law.as_float()
                    tree = tree_probabilities(scheme, i, delta)
                    layer_a = state_probabilities(scheme, nuclide.name, delta)[0]
                    if not (np.all(np.isfinite(q)) and np.all(np.isfinite(tree)) and np.all(np.isfinite(layer_a))):
                        nonfinite.append(f"{key} {nuclide.name} delta={delta!r}")
                        continue
                    worst_tree = max(worst_tree, float(np.max(np.abs(q - tree))))
                    ends = np.zeros(len(scheme))
                    for category, value in zip(law.categories, q):
                        ends[category.end] += value
                    worst_layer_a = max(worst_layer_a, float(np.max(np.abs(ends - layer_a))))
                    with localcontext(Context(prec=300)):
                        worst_sum = max(worst_sum, abs(sum(law.probabilities, Decimal(0)) - 1))
                    brackets &= bounds_hold(scheme, law)
    except (ValueError, ArithmeticError) as caught:
        error = f"{type(caught).__name__}: {caught}"
    ok = (error is None and not nonfinite and worst_tree <= 1e-12 and worst_layer_a <= 1e-12
          and worst_sum <= Decimal("1e-30") and brackets)
    return {"laws_checked": laws, "max_abs_decimal_minus_tree": worst_tree,
            "max_abs_end_state_minus_layer_a": worst_layer_a, "max_abs_sum_minus_one": float(worst_sum),
            "brackets_hold": brackets, "max_digits": digits, "nonfinite_references": nonfinite, "error": error,
            "ok": bool(ok)}


def ensemble_tests(s2: DecayScheme, design: dict, reps: int, n_ens: int, times2, pop_b: np.ndarray) -> dict:
    """ENS-B2 against exact moments: ET-MEAN, ET-VAR, ET-VARSUM, ET-VARSUM2."""
    index_of_time = {repr(t): k for k, t in enumerate(times2)}
    out = {}
    for t in ("2.0", "8.0"):
        it = index_of_time[t]
        p = memo(state_probabilities, s2, "Xa-400", float(t))[0]
        for name in design["ens"]["eligible"][t]:
            j = s2.index[name]
            out[f"ENS ET-MEAN {name} t={t}"] = evaluate(binomial_check, pop_b[:, it, j].sum(), reps * n_ens, p[j])
            moments = memo(state_count_moments, s2, "Xa-400", n_ens, float(t), j, j)
            out[f"ENS ET-VAR {name} t={t}"] = evaluate(
                scaled_chisquare_variance_check, sample_covariance(pop_b[:, it, j], pop_b[:, it, j]),
                moments.sigma_ab, moments.sample_covariance_variance(reps))
    i2, i8 = index_of_time["2.0"], index_of_time["8.0"]
    for a, b in design["ens"]["pairs_t2"]:
        moments = memo(state_sum_moments, s2, "Xa-400", n_ens, 2.0, s2.index[a], s2.index[b])
        total = pop_b[:, i2, s2.index[a]] + pop_b[:, i2, s2.index[b]]
        out[f"ENS ET-VARSUM {a}+{b} t=2.0"] = evaluate(
            scaled_chisquare_variance_check, sample_covariance(total, total), moments.sigma_ab,
            moments.sample_covariance_variance(reps))
    for a, b in design["ens"]["pairs_two_time"]:
        moments = memo(two_time_sum_moments, s2, "Xa-400", n_ens, 2.0, 8.0, s2.index[a], s2.index[b])
        total = pop_b[:, i2, s2.index[a]] + pop_b[:, i8, s2.index[b]]
        out[f"ENS ET-VARSUM2 {a}@2+{b}@8"] = evaluate(
            scaled_chisquare_variance_check, sample_covariance(total, total), moments.sigma_ab,
            moments.sample_covariance_variance(reps))
    return out


def comparison_cells(s2: DecayScheme, times2, reps: int, n_ens: int) -> dict[str, list[list[int]]]:
    """Cells of the Layer A / B2 homogeneity tests, fixed before any data from expected counts per engine.

    Both engines draw reps * n_ens nuclei from the same law under the null, so the expected counts of
    each row are known; pooling follows comparison.pooled_cells (minimum COMPARISON_MIN_EXPECTED). Tests
    left with fewer than two cells are not applicable and are not planned.
    """
    pooled = reps * n_ens
    out = {}
    p = memo(state_probabilities, s2, "Xa-400", times2)
    for it, t in enumerate(times2[1:], start=1):
        out[f"COMP CT-HOMO t={t!r}"] = pooled_cells(pooled * p[it], COMPARISON_MIN_EXPECTED)
    q = memo(decay_probabilities, s2, "Xa-400", times2)
    for j, nuclide in enumerate(s2.nuclides):
        if not nuclide.stable:
            expected = pooled * np.append(np.diff(q[:, j]), 1.0 - q[-1, j])
            out[f"COMP CT-DECAY {nuclide.name}"] = pooled_cells(expected, COMPARISON_MIN_EXPECTED)
        for k, (channel, branch) in enumerate(zip(nuclide.channels, nuclide.channel_probabilities())):
            if 0.0 < channel.photon_yield < 1.0:
                decays = pooled * branch * q[-1, j]
                out[f"COMP CT-PHOTON {channel_label(nuclide, channel)}"] = pooled_cells(
                    [channel.photon_yield * decays, (1.0 - channel.photon_yield) * decays], COMPARISON_MIN_EXPECTED)
    return {name: cells for name, cells in out.items() if len(cells) >= 2}


def comparison_tests(s2: DecayScheme, design: dict, reps: int, n_ens: int, times2, pop_a, pop_b, dec_a, dec_b,
                     pho_a, pho_b) -> dict:
    """ENS-A against ENS-B2: CT-HOMO, CT-VAR, CT-VARSUM, CT-DECAY, CT-PHOTON."""
    index_of_time = {repr(t): k for k, t in enumerate(times2)}
    cells = comparison_cells(s2, times2, reps, n_ens)
    out = {}
    for t in times2[1:]:
        it = index_of_time[repr(t)]
        name = f"COMP CT-HOMO t={t!r}"
        if name in cells:
            out[name] = evaluate(homogeneity_check, pop_a[:, it, :].sum(axis=0), pop_b[:, it, :].sum(axis=0),
                                 cells=cells[name])
    i2 = index_of_time["2.0"]
    for name in design["ens"]["eligible"]["2.0"]:
        j = s2.index[name]
        moments = memo(state_count_moments, s2, "Xa-400", n_ens, 2.0, j, j)
        out[f"COMP CT-VAR {name} t=2.0"] = evaluate(
            variance_ratio_f_check, sample_covariance(pop_a[:, i2, j], pop_a[:, i2, j]),
            sample_covariance(pop_b[:, i2, j], pop_b[:, i2, j]), moments.sigma_ab,
            moments.sample_covariance_variance(reps))
    for a, b in design["ens"]["pairs_t2"]:
        ja, jb = s2.index[a], s2.index[b]
        moments = memo(state_sum_moments, s2, "Xa-400", n_ens, 2.0, ja, jb)
        total_a, total_b = pop_a[:, i2, ja] + pop_a[:, i2, jb], pop_b[:, i2, ja] + pop_b[:, i2, jb]
        out[f"COMP CT-VARSUM {a}+{b} t=2.0"] = evaluate(
            variance_ratio_f_check, sample_covariance(total_a, total_a), sample_covariance(total_b, total_b),
            moments.sigma_ab, moments.sample_covariance_variance(reps))
    for name in design["comparison"]["unstable"]:
        test = f"COMP CT-DECAY {name}"
        if test in cells:
            j = s2.index[name]
            va = dec_a[:, :, j, :].sum(axis=(0, 2))
            vb = dec_b[:, :, j, :].sum(axis=(0, 2))
            out[test] = evaluate(homogeneity_check, np.append(va, reps * n_ens - va.sum()),
                                 np.append(vb, reps * n_ens - vb.sum()), cells=cells[test])
    for j, nuclide in enumerate(s2.nuclides):
        for k, channel in enumerate(nuclide.channels):
            test = f"COMP CT-PHOTON {channel_label(nuclide, channel)}"
            if 0.0 < channel.photon_yield < 1.0 and test in cells:
                pa, da = pho_a[:, :, j, k].sum(), dec_a[:, :, j, k].sum()
                pb, db = pho_b[:, :, j, k].sum(), dec_b[:, :, j, k].sum()
                out[test] = evaluate(homogeneity_check, [pa, da - pa], [pb, db - pb], cells=cells[test])
    return out


def big_tests(s2: DecayScheme, times2, size: int, populations: np.ndarray, decays: np.ndarray) -> dict:
    """BIG-STATE (chi-square) and BIG-DECAY (normal approximation) of one huge-population run."""
    out = {}
    p_grid = memo(state_probabilities, s2, "Xa-400", times2)
    q_end = memo(decay_probabilities, s2, "Xa-400", times2[-1])[0]
    for it, t in enumerate(times2[1:], start=1):
        out[f"BIG {size} BIG-STATE t={t!r}"] = evaluate(chi2, populations[it], p_grid[it])
    for j, nuclide in enumerate(s2.nuclides):
        if not nuclide.stable:
            total_j = int(decays[:, j, :].sum())
            out[f"BIG {size} BIG-DECAY {nuclide.name}"] = evaluate(
                binomial_z_check, total_j, size, float(q_end[j]), "binomial_normal_approximation")
    return out


def family_members(found: dict) -> dict:
    return {name: value for name, value in found.items() if "p_value" in value}


def errors_of(found: dict) -> list[str]:
    return [name for name, value in found.items() if "error" in value]


def scaled(initial: dict, scale: float) -> dict:
    return {k: max(1, round(v * scale)) for k, v in initial.items()}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--deviation", default=None, help="reason for a rerun on another commit (recorded)")
    args = parser.parse_args()
    frozen = frozen_hashes(smoke=args.smoke)
    protocol = json.loads((HERE / "protocol.json").read_text())
    design = json.loads((HERE / "design.json").read_text())
    out_dir = HERE / ("smoke" if args.smoke else "results")
    if args.smoke:
        code, manifests = record_code(REPO), out_dir / "manifests"
        result, pending = study(args, frozen, protocol, design, code, out_dir, manifests, None, None, None)
        publish_verification(result, pending, out_dir, manifests, None, code)
        print_summary(result)
        return r1_status(result)
    require_pycache_outside(REPO)
    if complete_published_attempt(out_dir, "verification", out_dir / "verification.json"):
        raise SystemExit("the last verification attempt published results/verification.json but ended before its "
                         "record said so; the record now shows it completed: commit results/ "
                         "(record_procedure.recovery)")
    tree = require_record_tree(REPO, "verification")
    code = record_code(REPO)
    if (out_dir / "verification.json").exists():
        raise SystemExit("results/verification.json exists; results are never overwritten")
    calibration = record_calibration(protocol, frozen, code, out_dir, REPO, design)
    attempt = begin_attempt(out_dir, "verification", code, deviation=args.deviation)

    def body() -> int:
        manifests = attempt["dir"] / "manifests"
        result, pending = study(args, frozen, protocol, design, code, out_dir, manifests, attempt, calibration, tree)
        publish_verification(result, pending, out_dir, manifests, attempt, code)       # its last act
        return r1_status(result)
    status = run_attempt(attempt, body, out_dir / "verification.json")
    print_summary(json.loads((out_dir / "verification.json").read_text()))          # after the completed record
    return status


def r1_status(result: dict) -> int:
    return 0 if result["verdicts"]["R1"] == "PASS" else 1


def quiet_print(*lines: str) -> None:
    """Print and flush. An output error (a closed pipe, a full disk) is ignored, and the stream is pointed at
    /dev/null so that later prints and the interpreter's flush at exit cannot fail either: no record run is
    stopped by its output, a completed attempt stays completed, and its exit status is the study's."""
    try:
        for line in lines:
            print(line, flush=True)
    except OSError:                     # a closed pipe, a full disk, an output over quota
        try:
            target = sys.stdout.fileno()
        except (OSError, ValueError):
            return
        devnull = os.open(os.devnull, os.O_WRONLY)
        try:
            os.dup2(devnull, target)
        finally:
            os.close(devnull)


def print_summary(result: dict) -> None:
    """The verification's summary, printed only after publication (for the record run, after its attempt is
    recorded as completed), so nothing the study computes is shown before it completes."""
    quiet_print(json.dumps(result.get("summary", {"verdicts": result.get("verdicts")}), indent=1))


def study(args, frozen: dict, protocol: dict, design: dict, code: dict, out_dir: Path, manifests: Path,
          attempt: dict | None, calibration: dict | None, tree: dict | None) -> tuple[dict, dict]:
    """The study itself. It writes and prints nothing: it returns the result and the manifests, which main()
    publishes together (publish_verification) and only then summarises."""
    manifests.mkdir(parents=True, exist_ok=True)
    pending: dict[str, str] = {}
    ns = protocol["smoke_namespace"] if args.smoke else protocol["namespace"]
    smoke = protocol["smoke"]
    scale = smoke["scale"] if args.smoke else 1.0
    runs = protocol["runs"]
    experiment = protocol["id"] + ("-smoke" if args.smoke else "")
    schemes = {k: DecayScheme.from_dict(v) for k, v in protocol["schemes"].items()}
    if {k: s.sha256() for k, s in schemes.items()} != design["scheme_sha256"]:
        raise SystemExit("schemes differ from design.json")
    cases = {k: {"initial": scaled(c["initial"], scale), "times": c["times"]} for k, c in protocol["cases"].items()}
    clock = time.perf_counter()
    result = {"protocol": {"id": protocol["id"], "frozen_sha256": frozen}, "namespace": ns, "smoke": args.smoke,
              "code": code, "environment": environment(),
              "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(), "runs": {},
              "timings_seconds": {}}

    def seed(*parts):
        return namespace_seed(ns, *parts), "/".join([ns, *map(str, parts)])

    def b2(scheme, initial, times, generator, label_parts, **kw):
        s, label = seed(*label_parts)
        stream = WordStream(make_source(generator, s, label))
        t0 = time.perf_counter()
        run = simulate_aggregate(scheme, initial, times, stream, **kw)
        return run, stream, time.perf_counter() - t0, s, label

    def manifest(run_id, run, stream, label):
        pending[f"{run_id}.json"] = (json.dumps(
            aggregate_manifest(experiment=experiment, run_id=run_id, run=run, code=code,
                               extra={"seed_label": label}), indent=1) + "\n")

    deterministic: dict[str, bool] = {}
    result["laws"] = law_checks(schemes, protocol["cases"])
    deterministic["D8 interval laws"] = result["laws"]["ok"]

    tests = {"NC1": {}, "NC2": {}}
    errors: list[str] = []
    d1, d2 = True, True
    for family in ("NC1", "NC2"):
        generator = runs[family]["generator"]
        result["runs"][family] = {}
        for key in runs[family]["schemes"]:
            run, stream, seconds, s, label = b2(schemes[key], cases[key]["initial"], cases[key]["times"], generator,
                                                (family, key, "physics"))
            entry = {"seed": s, "seed_label": label, "seconds": seconds, "summary": run.summary(),
                     "stream_sha256": stream.source.stream_digest, "consumption": b2_consumption(run, stream),
                     "bookkeeping": bookkeeping(run),
                     "populations": run.populations.tolist(), "decays": run.decays.tolist(),
                     "photons": run.photons.tolist(), "path_counts": run.path_counts.tolist()}
            found = battery(key, run, schemes[key], cases[key], family)
            if key == "S2":
                path_found, path_cells = path_tests(run, family, cases[key], design["path_tests"]["min_expected"])
                found.update(path_found)
                if not args.smoke and path_cells != design["path_tests"]["cells"]:
                    raise SystemExit("BT-PATH cells differ from design.json")
            entry["tests"] = found
            tests[family].update(family_members(found))
            errors += [f"{family}: {name}" for name in errors_of(found)]
            d1 &= entry["consumption"]["ok"]
            d2 &= all(entry["bookkeeping"].values())
            result["runs"][family][key] = entry
            result["timings_seconds"][f"{family} {key}"] = seconds
            manifest(f"{family}_{key}", run, stream, label)
            if family == "NC1" and key == "S2":
                nc1_s2 = (run, stream, s, label)

    # Ensembles of S2: B2, Layer A, and the naive aggregation (positive control).
    s2, times2 = schemes["S2"], protocol["cases"]["S2"]["times"]
    reps = smoke["ens_replicates"] if args.smoke else runs["ENS-B2"]["replicates"]
    n_ens = runs["ENS-B2"]["initial"]["Xa-400"]
    pop_b, dec_b, pho_b, pop_a, dec_a, pho_a, pop_n = [], [], [], [], [], [], []
    seeds, digests = {"ENS-B2": [], "ENS-A": [], "ENS-NAIVE": []}, {"ENS-B2": [], "ENS-A": [], "ENS-NAIVE": []}
    a_invariants, a_consumption, a_zero = True, True, 0
    t0 = time.perf_counter()
    for r in range(reps):
        run, stream, _, s, _ = b2(s2, {"Xa-400": n_ens}, times2, "numpy-pcg64", ("ENS-B2", "S2", "physics", f"{r:04d}"))
        pop_b.append(run.populations)
        dec_b.append(run.decays)
        pho_b.append(run.photons)
        seeds["ENS-B2"].append(s)
        digests["ENS-B2"].append(stream.source.stream_digest)
        d1 &= b2_consumption(run, stream)["ok"]
        d2 &= all(bookkeeping(run).values())
        s, label = seed("ENS-A", "S2", "physics", f"{r:04d}")
        stream = WordStream(make_source("numpy-pcg64", s, label))
        layer_a = simulate_decays(s2, {"Xa-400": n_ens}, times2[-1], stream)
        pops, decays, photons, at_zero = layer_a_observations(layer_a, times2)
        pop_a.append(pops)
        dec_a.append(decays)
        pho_a.append(photons)
        a_zero += at_zero
        a_invariants &= all(decay_invariants(layer_a).values())
        a_consumption &= stream.words_read == n_ens * s2.words_per_nucleus
        seeds["ENS-A"].append(s)
        digests["ENS-A"].append(stream.source.stream_digest)
        s, label = seed("ENS-NAIVE", "S2", "physics", f"{r:04d}")
        stream = WordStream(make_source("numpy-pcg64", s, label))
        pop_n.append(naive_aggregate(s2, {"Xa-400": n_ens}, times2, stream))
        seeds["ENS-NAIVE"].append(s)
        digests["ENS-NAIVE"].append(stream.source.stream_digest)
    result["timings_seconds"]["ensembles"] = time.perf_counter() - t0
    pop_b, dec_b, pho_b = np.array(pop_b), np.array(dec_b), np.array(pho_b)
    pop_a, dec_a, pho_a, pop_n = np.array(pop_a), np.array(dec_a), np.array(pho_a), np.array(pop_n)
    d1 &= a_consumption
    d2 &= a_invariants
    deterministic["D7 ensemble streams"] = (all(len(set(v)) == reps for v in seeds.values())
                                            and all(len(set(v)) == reps for v in digests.values()) and a_zero == 0)
    names = [n.name for n in s2.nuclides]
    ens_tests = ensemble_tests(s2, design, reps, n_ens, times2, pop_b)
    comp_tests = comparison_tests(s2, design, reps, n_ens, times2, pop_a, pop_b, dec_a, dec_b, pho_a, pho_b)
    errors += [f"NC1: {name}" for name in errors_of(ens_tests) + errors_of(comp_tests)]
    result["ensemble"] = {"replicates": reps, "per_replicate": n_ens, "nuclides": names, "times": times2,
                          "populations_b2": pop_b.tolist(), "populations_layer_a": pop_a.tolist(),
                          "populations_naive": pop_n.tolist(), "layer_a_decays_at_time_zero": a_zero,
                          "tests": ens_tests, "comparison_tests": comp_tests}
    pending["ENSEMBLES_S2.json"] = (json.dumps(
        {"schema": "rnglab-nuclear-ensemble-manifest-v1", "experiment": experiment, "code": code,
         "seed_recipe": protocol["seed_recipe"], "seeds": seeds, "stream_sha256": digests}, indent=1) + "\n")
    tests["NC1"].update(family_members(ens_tests))
    tests["NC1"].update(family_members(comp_tests))

    # Huge populations.
    all_big, result["big"] = {}, {}
    d9 = True
    for size in (smoke["big_sizes"] if args.smoke else runs["BIG"]["initial_sizes"]):
        run, stream, seconds, s, label = b2(s2, {"Xa-400": size}, times2, "numpy-pcg64", ("BIG", "S2", "physics", size))
        d1 &= b2_consumption(run, stream)["ok"]
        d2 &= all(bookkeeping(run).values())
        d9 &= run.plan.parts == -(-size // PART_MAX) and set(run.populations.sum(axis=1).tolist()) == {size}
        found = big_tests(s2, times2, size, run.populations, run.decays)
        all_big.update(found)
        errors += [f"NC1: {name}" for name in errors_of(found)]
        result["big"][str(size)] = {"parts": run.plan.parts, "seconds": seconds, "summary": run.summary(),
                                    "populations": [[int(x) for x in row] for row in run.populations]}
        result["timings_seconds"][f"BIG {size}"] = seconds
        manifest(f"BIG_{size}", run, stream, label)
    if not args.smoke:
        d9 &= [result["big"][str(s)]["parts"] for s in runs["BIG"]["initial_sizes"]] == runs["BIG"]["parts"]
    deterministic["D9 big populations"] = d9
    tests["NC1"].update(family_members(all_big))
    result["big_tests"] = all_big

    # Positive controls.
    pc = runs["PC-B1"]
    s1_perturbed = copy.deepcopy(protocol["schemes"]["S1"])
    s1_perturbed["nuclides"][0]["half_life"] = pc["perturbation"]["Sa-400"]["half_life"]
    s1_perturbed["name"] += "-perturbed"
    pc1_case = {"initial": scaled(pc["initial"], scale), "times": protocol["cases"]["S1"]["times"]}
    run, stream, _, _, label = b2(DecayScheme.from_dict(s1_perturbed), pc1_case["initial"], pc1_case["times"],
                                  pc["generator"], ("PC-B1", "S1", "physics"))
    pc1_tests = battery("S1", run, schemes["S1"], pc1_case, "PC-B1")
    manifest("PC-B1_S1", run, stream, label)
    s_pc2, label_pc2 = seed("PC-B2", "S2", "physics")
    naive = naive_aggregate(s2, scaled(runs["PC-B2"]["initial"], scale), times2,
                            WordStream(make_source("numpy-pcg64", s_pc2, label_pc2)))
    p12 = state_probabilities(s2, "Xa-400", times2)
    pc2_tests = {f"PC-B2 S2 BT-STATE Xa-400 t={t!r}": chi2(naive[it], p12[it])
                 for it, t in enumerate(times2[1:], start=1)}
    pc2_cells = comparison_cells(s2, times2, reps, n_ens)[f"COMP CT-HOMO t={times2[-1]!r}"]
    pc2_comp = evaluate(homogeneity_check, pop_a[:, -1, :].sum(axis=0), pop_n[:, -1, :].sum(axis=0), cells=pc2_cells)

    # Software checks.
    xref = []
    for case in runs["XREF"]["cases"]:
        key = case["scheme"]
        initial = scaled(case["initial"], smoke["xref_scale"] if args.smoke else 1.0)
        s, label = seed("XREF", key, "physics", case["part_max"])
        fast_stream = WordStream(make_source("numpy-pcg64", s, label))
        slow_stream = WordStream(make_source("numpy-pcg64", s, label))
        fast = simulate_aggregate(schemes[key], initial, protocol["cases"][key]["times"], fast_stream,
                                  part_max=case["part_max"])
        t0 = time.perf_counter()
        slow = simulate_aggregate_reference(schemes[key], initial, protocol["cases"][key]["times"], slow_stream,
                                            part_max=case["part_max"])
        same = all(np.array_equal(getattr(fast, f), getattr(slow, f))
                   for f in ("populations", "path_counts", "decays", "photons", "neutrons"))
        xref.append({"scheme": key, "initial": initial, "part_max": case["part_max"], "parts": fast.plan.parts,
                     "digests_equal": fast.digest == slow.digest, "arrays_equal": bool(same),
                     "stream_digests_equal": fast_stream.source.stream_digest == slow_stream.source.stream_digest,
                     "reference_seconds": time.perf_counter() - t0})
        d1 &= b2_consumption(fast, fast_stream)["ok"] and b2_consumption(slow, slow_stream)["ok"]
    result["xref"] = xref
    deterministic["D3 xref"] = all(x["digests_equal"] and x["arrays_equal"] and x["stream_digests_equal"] for x in xref)
    run0, stream0, s0, label0 = nc1_s2
    replay, replay_stream, _, _, _ = b2(s2, cases["S2"]["initial"], times2, "numpy-pcg64", ("NC1", "S2", "physics"))
    deterministic["D4 replay"] = (replay.digest == run0.digest
                                  and replay_stream.source.stream_digest == stream0.source.stream_digest)
    checkpoints = []
    full, _, _, _, _ = b2(s2, cases["S2"]["initial"], times2, "numpy-pcg64", ("NC1", "S2", "physics"),
                          checkpoint_every=1, on_checkpoint=checkpoints.append)
    cut = runs["CHECKPOINT"]["resume_after_interval"]
    checkpoint = json.loads(json.dumps(checkpoints[cut - 1]))
    resumed = simulate_aggregate(s2, cases["S2"]["initial"], times2,
                                 WordStream(NumpyPCG64Source(s0, label=label0, start_word=checkpoint["words_consumed"])),
                                 resume_from=checkpoint)
    deterministic["D5 checkpoint"] = (resumed.digest == full.digest == run0.digest
                                      and np.array_equal(resumed.populations, run0.populations[cut:]))
    records = []
    lean, _, _, _, _ = b2(s2, cases["S2"]["initial"], times2, "numpy-pcg64", ("NC1", "S2", "physics"),
                          keep_outputs=False, sink=lambda k, record: records.append(record["populations"]))
    deterministic["D6 sink"] = lean.digest == run0.digest and np.array_equal(np.array(records), run0.populations[1:])
    deterministic["D1 consumption"] = bool(d1)
    deterministic["D2 bookkeeping and Layer A invariants"] = bool(d2)

    # Families, Holm, verdicts.
    alpha = protocol["families"]["alpha"]
    result["families"] = {}
    for family, members in tests.items():
        names_f = list(members)
        correction = holm([members[n]["p_value"] for n in names_f], alpha)
        result["families"][family] = {"tests": names_f, "p_values": [members[n]["p_value"] for n in names_f],
                                      "holm": correction,
                                      "rejected": [n for n, r in zip(names_f, correction["reject"]) if r]}
    if not args.smoke:
        deterministic["D10 planned tests ran"] = all(
            sorted(tests[f]) == sorted(design["planned_tests"][f]) for f in ("NC1", "NC2"))
    level = alpha / len(tests["NC1"])
    t1 = pc1_tests["PC-B1 S1 BT-STATE Sa-400 t=2.0"]["p_value"]
    t2 = pc2_tests["PC-B2 S2 BT-STATE Xa-400 t=12.0"]["p_value"]
    t3 = pc2_comp.get("p_value", float("nan"))
    result["positive_controls"] = {
        "level": level,
        "PC-B1": {"target": "BT-STATE Sa-400 t=2.0", "p_value": t1, "rejected": t1 < level, "all_tests": pc1_tests},
        "PC-B2 single run": {"target": "BT-STATE Xa-400 t=12.0", "p_value": t2, "rejected": t2 < level,
                             "all_tests": pc2_tests},
        "PC-B2 comparison": {"target": "CT-HOMO t=12.0, ENS-A against ENS-NAIVE", "p_value": t3,
                             "rejected": bool(t3 < level), "test": pc2_comp}}
    failed = [k for k, ok in deterministic.items() if not ok]
    failed += [f"test could not be computed in {name}" for name in errors]
    failed += [f"Holm rejection in {f}: {n}" for f, v in result["families"].items() for n in v["rejected"]]
    result["verdicts"] = {"deterministic": deterministic, "R1": "PASS" if not failed else "FAIL",
                          "R1_failed_items": failed,
                          "R2": "MET" if all(result["positive_controls"][k]["rejected"]
                                             for k in ("PC-B1", "PC-B2 single run", "PC-B2 comparison")) else "NOT MET"}
    result["calibration"] = calibration
    result["tree"] = tree
    result["verdicts"]["R0"] = calibration["R0"] if calibration else "not read: smoke run"
    result["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    result["elapsed_seconds"] = time.perf_counter() - clock
    result["summary"] = {"namespace": ns, "R1": result["verdicts"]["R1"], "failed": failed,
                         "R2": result["verdicts"]["R2"], "PC p-values": [t1, t2, t3], "level": level,
                         "family sizes": {f: len(v) for f, v in tests.items()},
                         "smallest p": {f: min(v["p_values"]) for f, v in result["families"].items()},
                         "elapsed_seconds": round(result["elapsed_seconds"], 1)}
    return result, pending


def publish_verification(result: dict, pending: dict, out_dir: Path, manifests: Path, attempt: dict | None,
                         code: dict) -> None:
    """Write the manifests and publish verification.json: the last act of an attempt (nothing is printed before
    the attempt is recorded as completed). A record run first checks that its code is still the code it took
    when it started."""
    if attempt is not None:
        require_same_code(code, REPO, "verification")
    manifests.mkdir(parents=True, exist_ok=True)
    for name, text in pending.items():
        (manifests / name).write_text(text)
    if attempt is None:
        (out_dir / "verification.json").write_text(json.dumps(result, indent=1, default=str) + "\n")
        return
    result["attempt"] = {"number": attempt["record"]["attempt"], "record": f"attempts/{attempt['path'].name}",
                         "manifests": f"attempts/{attempt['dir'].name}/manifests",
                         "previous_attempts": attempt["record"]["previous_attempts"],
                         "deviation": attempt["record"]["deviation"]}
    write_new(out_dir / "verification.json", json.dumps(result, indent=1, default=str) + "\n")


if __name__ == "__main__":
    sys.exit(main())
