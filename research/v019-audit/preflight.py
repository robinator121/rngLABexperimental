"""Campaign preflight for the v0.19 Stage D pilot. Generates no pilot data.

    PYTHONPATH=src python research/v019-audit/preflight.py \
        [--root research/v019-d-pilot] [--namespace rnglab-v019-d-pilot] [--workers 8] [--out preflight.json]

Checks, each PASS, WARN or FAIL with the evidence behind it:

    protocol      protocol.json matches protocol.sha256
    code          the pilot modules exist; hashes of all code the pilot executes
    git           branch, commit, and whether every hashed file is committed and unmodified
    environment   interpreter and library versions against the lock file; thread variables; libc
    golden        short engineering runs reproduce the recorded digests, through the pilot code
                  and through the independent replay (detects a platform that rounds differently)
    seeds         pilot, smoke and audit namespaces share no seed
    resources     free disk and memory against measured per-run needs
    store         state of the output tree: nothing partial or corrupt, nothing unplanned
    manifest      the pre-generation manifest exists, lists the same code hashes, records the
                  interrupted-run policy (regenerate_same_id, three retries) and the TDK deadline
                  (without it stop criterion 5 is INSUFFICIENT_INFORMATION); every planned run's
                  seeds and random-input hash are recomputed by the independent code

--namespace rnglab-v019-d-pilot-smoke checks the smoke (ten runs, root smoke/), which the freeze
procedure puts before the pilot's own manifest.

The engineering runs use the audit-only namespace and are 96 sweeps long. The
command refuses to use a pilot or smoke namespace for them. Exit status 0 means
no FAIL; it is not an authorisation to generate.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

_HERE = Path(__file__).resolve().parent
REPO = _HERE.parents[1]
AUDIT_NAMESPACE = "rnglab-v019-audit-smoke"
GOLDEN_PATH = _HERE / "golden.json"
PLAN = {"calibration": 999, "healthy": 2000, "fault/SHIFT": 1000, "fault/LAG1": 1000, "fault/DUP": 1000}
THREAD_VARIABLES = ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "BLIS_NUM_THREADS",
                    "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
# Code the pilot executes: its own modules, the imported Stage C modules, and the package modules
# named in protocol.json (application, endpoint, operational_alarm, data.seeds).
HASHED = ("research/v019-d-pilot/protocol.json", "research/v019-d-pilot/design.json", "research/v019-d-pilot/design.py",
          "research/v019-d-pilot/dlib.py", "research/v019-d-pilot/measure.py",
          "research/v018-h1-detection/h1lib.py", "research/v018-h1-detection/analyze.py",
          "research/v018-h1-detection/run_campaign.py", "src/rnglab/seeds.py", "src/rnglab11/statistics.py",
          "src/rnglab12/uncertainty.py", "src/rnglab14/models/pimc_harmonic.py", "src/rnglab14/sources/__init__.py",
          "src/rnglab14/sources/base.py", "src/rnglab14/sources/reference.py")
# Measured on full-length audit runs (docs/v019/IMPLEMENTATION_AUDIT.md, benchmark section).
BYTES_PER_RUN = 4_400_000        # compressed trace plus the indented run record
RSS_BYTES_PER_WORKER = 300 * 2**20
LOCK_FILE = "research/phase2-20261030/environment/requirements-lock.txt"
INTERRUPTED_RUN_POLICY = "regenerate_same_id"   # the author's decision of 2026-11-13
MAX_INFRASTRUCTURE_RETRIES = 3
PROTOCOL_WORKERS = 8                # protocol.json storage.parallelism
# Namespace -> (run plan, output root). The smoke comes first in the freeze procedure (step 5).
NAMESPACES = {"rnglab-v019-d-pilot": (PLAN, "research/v019-d-pilot"),
              "rnglab-v019-d-pilot-smoke": ({condition: 2 for condition in PLAN}, "research/v019-d-pilot/smoke")}


def _load(name: str, path: Path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def result(check: str, status: str, **detail: Any) -> dict[str, Any]:
    return {"check": check, "status": status, **detail}


def check_protocol() -> dict[str, Any]:
    recorded = (REPO / "research/v019-d-pilot/protocol.sha256").read_text().split()[0]
    actual = sha256_file(REPO / "research/v019-d-pilot/protocol.json")
    return result("protocol", "PASS" if recorded == actual else "FAIL", recorded=recorded, actual=actual)


def check_code(extra_pilot_modules: bool = True) -> dict[str, Any]:
    paths = list(HASHED)
    if extra_pilot_modules:
        paths += sorted(p.relative_to(REPO).as_posix() for p in (REPO / "research/v019-d-pilot").glob("*.py")
                        if p.relative_to(REPO).as_posix() not in HASHED)
    hashes = {p: (sha256_file(REPO / p) if (REPO / p).is_file() else None) for p in paths}
    missing = [p for p, h in hashes.items() if h is None]
    return result("code", "FAIL" if missing else "PASS", missing=missing, sha256=hashes)


def _git(*args: str) -> str | None:
    try:
        return subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        return None


def check_git(paths: list[str]) -> dict[str, Any]:
    head = _git("rev-parse", "HEAD")
    if head is None:
        return result("git", "FAIL", reason="not a git repository or git unavailable")
    status = _git("status", "--porcelain", "--", *paths) or ""
    dirty = sorted(line[3:] for line in status.splitlines() if line.strip())
    branch = (_git("rev-parse", "--abbrev-ref", "HEAD") or "").strip()
    return result("git", "FAIL" if dirty else "PASS", commit=head.strip(), branch=branch,
                  uncommitted_or_untracked=dirty,
                  note="every hashed file must be committed before the pre-generation manifest is written")


def check_environment() -> dict[str, Any]:
    import scipy
    versions = {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__}
    lock, mismatched = {}, []
    lock_path = REPO / LOCK_FILE
    if lock_path.is_file():
        for line in lock_path.read_text().splitlines():
            name, _, version = line.strip().partition("==")
            if version:
                lock[name.lower()] = version.split()[0]
        mismatched = [name for name in ("numpy", "scipy") if lock.get(name) not in (None, versions[name])]
    threads = {name: os.environ.get(name) for name in THREAD_VARIABLES}
    unset = [name for name, value in threads.items() if value != "1"]
    status = "FAIL" if mismatched else "WARN" if unset or not lock else "PASS"
    return result("environment", status, versions=versions, lock_file=LOCK_FILE if lock else None,
                  lock_mismatch=mismatched, thread_variables=threads, thread_variables_not_1=unset,
                  libc=" ".join(platform.libc_ver()), machine=platform.machine(), cpu_count=os.cpu_count(),
                  note="the sampler's traces depend on libm (pow, exp, log1p, cos): bit-for-bit replay is a "
                       "same-platform claim, which the golden check tests")


def _digest(arrays: dict[str, np.ndarray], stream_sha256: str) -> dict[str, str]:
    out = {key: hashlib.sha256(np.ascontiguousarray(arrays[key]).tobytes()).hexdigest() for key in ("CV", "TD", "marks")}
    out["emitted"] = stream_sha256
    return out


def golden_runs(use_pilot_code: bool) -> dict[str, dict[str, str]]:
    """Digests of four 96-sweep engineering runs (audit namespace, hand-built profile)."""
    replay = _load("v019_audit_replay", _HERE / "replay.py")
    caps = replay.rc.protocol_caps()
    dlib = _load("v019_preflight_dlib", REPO / "research/v019-d-pilot/dlib.py") if use_pilot_code else None
    out = {}
    for family in ("healthy", "SHIFT", "LAG1", "DUP"):
        condition = "healthy" if family == "healthy" else "fault/" + family
        profile = None
        if family != "healthy":
            cap = np.float64(caps[family])
            profile = {"t0": 3, "T_cap": np.float64(7.5), "r": np.float64(cap / np.float64(7.5)), "s_max": cap}
        if use_pilot_code:
            run = dlib.generate_run(AUDIT_NAMESPACE, condition, family, 2, sweeps=96, profile=profile)
        else:
            run = replay.replay(AUDIT_NAMESPACE, condition, family, 2, sweeps=96, profile=profile)
        out[family] = _digest(run, run["emitted_sha256"])
    return out


def check_golden() -> dict[str, Any]:
    if not GOLDEN_PATH.is_file():
        return result("golden", "FAIL", reason=f"{GOLDEN_PATH.name} is missing; write it once with --write-golden")
    golden = json.loads(GOLDEN_PATH.read_text())
    independent = golden_runs(False)
    pilot = golden_runs(True)
    ok = independent == golden["digests"] and pilot == golden["digests"]
    return result("golden", "PASS" if ok else "FAIL", recorded_environment=golden["environment"],
                  independent_replay_matches=independent == golden["digests"],
                  pilot_code_matches=pilot == golden["digests"])


def check_seeds() -> dict[str, Any]:
    integrity = _load("v019_audit_integrity", _HERE / "integrity.py")
    recount = _load("v019_audit_recount", _HERE / "recount.py")
    namespaces = {"rnglab-v019-d-pilot": PLAN, "rnglab-v019-d-pilot-smoke": {k: 2 for k in PLAN},
                  AUDIT_NAMESPACE: {k: 16 for k in PLAN}}
    clashes = integrity.namespace_collisions(recount.namespace_seed, namespaces,
                                             ("local", "centroid", "start", "fault", "profile"))
    return result("seeds", "FAIL" if clashes else "PASS", seeds_checked=5 * (5999 + 10 + 80), collisions=clashes[:5])


def check_resources(root: Path, workers: int, runs: int) -> dict[str, Any]:
    target = root if root.exists() else root.parent
    free = shutil.disk_usage(target).free
    need_disk = int(1.25 * runs * BYTES_PER_RUN)
    available = None
    try:
        for line in Path("/proc/meminfo").read_text().splitlines():
            if line.startswith("MemAvailable:"):
                available = int(line.split()[1]) * 1024
    except OSError:
        pass
    need_memory = workers * RSS_BYTES_PER_WORKER
    low_memory = available is not None and available < need_memory
    status = "FAIL" if free < need_disk or low_memory else "WARN" if available is None else "PASS"
    return result("resources", status, disk_free_bytes=free, disk_needed_bytes=need_disk,
                  memory_available_bytes=available, memory_needed_bytes=need_memory, workers=workers,
                  cpu_count=os.cpu_count(), basis={"bytes_per_run": BYTES_PER_RUN, "rss_per_worker": RSS_BYTES_PER_WORKER})


def check_store(root: Path, plan: dict[str, int] = PLAN) -> dict[str, Any]:
    integrity = _load("v019_audit_integrity", _HERE / "integrity.py")
    scan = integrity.scan_store(root, plan, verify_hash=False)
    counts = scan["counts"]
    bad = counts["PARTIAL"] + counts["CORRUPT"] + len(scan["unplanned_files"])
    status = "FAIL" if bad else "WARN" if counts["COMPLETE"] else "PASS"
    return result("store", status, counts=counts, unplanned_files=scan["unplanned_files"][:20],
                  not_complete_examples=[(r["condition"], r["replica"], r["state"], r["reason"])
                                         for r in scan["not_complete"] if r["state"] != "MISSING"][:20],
                  note="COMPLETE runs before authorisation mean generation has already started")


def planned_run_problems(manifest: dict[str, Any], namespace: str, plan: dict[str, int]) -> list[str]:
    """The manifest's planned runs against the plan of the namespace, with every run's seeds and
    random-input hash recomputed by the independent code. At most ten problems are listed."""
    recount = _load("v019_audit_recount", _HERE / "recount.py")
    problems: list[str] = []
    if manifest.get("namespace") != namespace:
        problems.append(f"the manifest is for the namespace {manifest.get('namespace')!r}")
    if manifest.get("protocol_sha256") != sha256_file(REPO / "research/v019-d-pilot/protocol.json"):
        problems.append("the manifest names another protocol hash")
    listed = {(run.get("condition"), run.get("replica")): run for run in manifest.get("planned_runs") or []}
    expected = {(condition, replica) for condition, count in plan.items() for replica in range(count)}
    if set(listed) != expected:
        problems.append(f"planned runs differ from the plan: {len(set(listed) - expected)} not in the plan, "
                        f"{len(expected - set(listed))} of the plan absent")
    for condition, replica in sorted(expected & set(listed)):
        if len(problems) >= 10:
            break
        seeds = recount.run_seeds(namespace, condition, replica)
        entry = listed[(condition, replica)]
        if entry.get("seeds") != seeds:
            problems.append(f"{condition}/{replica:05d}: seeds differ from the independent derivation")
        elif entry.get("rng_input_sha256") != recount.rng_input_sha256(seeds):
            problems.append(f"{condition}/{replica:05d}: rng_input_sha256 differs from the independent generators")
    return problems


def check_manifest(root: Path, code: dict[str, Any], namespace: str | None = None,
                   plan: dict[str, int] | None = None) -> dict[str, Any]:
    """With a namespace and plan, the planned runs are checked too (`planned_run_problems`)."""
    path = root / "pre_generation_manifest.json"
    if not path.is_file():
        return result("manifest", "FAIL", reason="pre_generation_manifest.json does not exist; the protocol requires "
                      "it to be pushed before the first pilot draw", tdk_deadline=None,
                      criterion_5="INSUFFICIENT_INFORMATION")
    try:
        manifest = json.loads(path.read_text(encoding="utf-8"))
    except ValueError as error:
        return result("manifest", "FAIL", reason=f"unreadable: {error}")
    deadline = manifest.get("tdk_deadline")
    listed = manifest.get("code_sha256") or manifest.get("sha256") or {}
    differing = sorted(p for p, h in code["sha256"].items() if p in listed and listed[p] != h)
    unlisted = sorted(p for p in code["sha256"] if p not in listed)
    # What happens to a run that was interrupted is the author's decision of 2026-11-13: the same
    # identity is run again, at most three times. The pipeline refuses to generate from a manifest
    # that does not record it.
    policy = manifest.get("interrupted_run_policy")
    rules = manifest.get("interrupted_run_rules") or {}
    policy_recorded = (policy == INTERRUPTED_RUN_POLICY and rules.get("policy") == INTERRUPTED_RUN_POLICY
                       and rules.get("max_infrastructure_retries") == MAX_INFRASTRUCTURE_RETRIES)
    identity = planned_run_problems(manifest, namespace, plan) if namespace is not None and plan is not None else []
    # The manifest is committed after it is written, so its commit is an ancestor of HEAD, not HEAD;
    # the code hashes above are the comparison that counts. A manifest written from a tree with
    # uncommitted changes names a commit that does not contain what it hashed.
    recorded_git = manifest.get("git") or {}
    unclean = bool(recorded_git) and recorded_git.get("working_tree_clean") is not True
    authorised, authorisation = False, root / "authorization.json"
    if authorisation.is_file():
        try:
            entry = json.loads(authorisation.read_text(encoding="utf-8"))
            authorised = bool(entry.get("pre_generation_manifest_sha256") == sha256_file(path)
                              and entry.get("authorized_by") and entry.get("authorized_on"))
        except ValueError:
            authorised = False
    deviates = manifest.get("workers") not in (None, PROTOCOL_WORKERS)
    status = ("FAIL" if differing or not policy_recorded or identity
              else "WARN" if deadline is None or unlisted or deviates or not authorised or unclean else "PASS")
    return result("manifest", status, tdk_deadline=deadline,
                  criterion_5=None if deadline else "INSUFFICIENT_INFORMATION",
                  code_hashes_differing=differing, code_hashes_not_listed=unlisted,
                  interrupted_run_policy=policy if policy_recorded else None,
                  max_infrastructure_retries=rules.get("max_infrastructure_retries") if policy_recorded else None,
                  planned_run_problems=identity, planned_runs_checked=namespace is not None and plan is not None,
                  manifest_commit=recorded_git.get("commit"), manifest_written_from_unclean_tree=unclean,
                  workers=manifest.get("workers"), workers_deviate_from_protocol=deviates,
                  authorised_for_this_manifest=authorised)


def run(root: Path | None, namespace: str, workers: int) -> dict[str, Any]:
    if namespace not in NAMESPACES:
        raise SystemExit("preflight checks the pilot namespace rnglab-v019-d-pilot or its smoke namespace only")
    plan, default_root = NAMESPACES[namespace]
    root = REPO / default_root if root is None else root
    code = check_code()
    checks = [check_protocol(), code, check_git([p for p, h in code["sha256"].items()]), check_environment(),
              check_golden(), check_seeds(), check_resources(root, workers, sum(plan.values())), check_store(root, plan),
              check_manifest(root, code, namespace, plan)]
    failed = [c["check"] for c in checks if c["status"] == "FAIL"]
    warned = [c["check"] for c in checks if c["status"] == "WARN"]
    return {"schema": "rnglab-v019-audit-preflight-v1", "root": str(root), "namespace": namespace,
            "generates_pilot_data": False, "checks": checks, "failed": failed, "warned": warned,
            "ready": not failed,
            "statement": "A clean preflight is a precondition, not an authorisation: the first pilot draw "
                         "needs the author's explicit go-ahead."}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--root", type=Path, default=None, help="default: the output root of the namespace")
    parser.add_argument("--namespace", default="rnglab-v019-d-pilot", choices=sorted(NAMESPACES))
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--write-golden", action="store_true",
                        help="record the golden digests of the engineering runs on this platform and exit")
    args = parser.parse_args(argv)
    if args.write_golden:
        independent, pilot = golden_runs(False), golden_runs(True)
        if independent != pilot:
            raise SystemExit("pilot code and independent replay disagree; golden digests not written")
        import scipy
        GOLDEN_PATH.write_text(json.dumps({
            "namespace": AUDIT_NAMESPACE, "runs": "replica 2, 96 sweeps, profile t0=3, T_cap=7.5, s_max=cap",
            "environment": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__,
                            "libc": " ".join(platform.libc_ver()), "machine": platform.machine()},
            "digests": independent}, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        print(f"wrote {GOLDEN_PATH}")
        return 0
    report = run(args.root, args.namespace, args.workers)
    text = json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(text, encoding="utf-8")
    for check in report["checks"]:
        print(f"{check['status']:4}  {check['check']}")
    print("READY (no FAIL); generation still needs explicit authorisation" if report["ready"]
          else "NOT READY: " + ", ".join(report["failed"]))
    return 0 if report["ready"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
