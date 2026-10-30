#!/usr/bin/env python3
"""Prepare an empty public Phase II replay tree without running simulations.

Uses only Python's standard library and a local Git object for the pinned
upstream revision. Neither Phase I inputs nor omitted trace archives are needed.
Frozen generation scripts and the two execution protocols are copied bytewise.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile


PINNED_COMMIT = "dd96adba23b8100aba9147d71fffc63932fc6b0d"
PUBLIC_PHASE = Path(__file__).resolve().parents[1]
DEFAULT_REPOSITORY = Path(__file__).resolve().parents[3]
PROTOCOLS = (
    Path("experiment/confirmation-v1/protocol.json"),
    Path("experiment/fault-study-v1/protocol.json"),
)
TEST_FILES = (
    "experiment/test_study.py",
    "experiment/test_fault_study.py",
    "baseline/test_independent_checker.py",
    "rus/test_counterexample.py",
    "analytic-dependence/test_analytic_dependence.py",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(2**20), b""):
            h.update(block)
    return h.hexdigest()


def git(repository: Path, *args: str) -> bytes:
    result = subprocess.run(
        ["git", "-C", str(repository), *args], capture_output=True, check=False
    )
    if result.returncode:
        detail = result.stderr.decode(errors="replace").strip()
        raise RuntimeError(f"Git command failed: {detail}")
    return result.stdout


def prepare(repository: Path, destination: Path) -> dict:
    repository, destination = repository.resolve(), destination.resolve()
    if destination.exists():
        raise FileExistsError("Destination must be new; existing trees are preserved")
    if destination.is_relative_to(PUBLIC_PHASE):
        raise ValueError("Choose a scratch destination outside the public payload")

    revision = git(repository, "rev-parse", "--verify", PINNED_COMMIT+"^{commit}").decode().strip()
    if revision != PINNED_COMMIT:
        raise RuntimeError("Local Git object does not resolve to the pinned revision")

    scripts = sorted(
        p for p in PUBLIC_PHASE.rglob("*.py")
        if "reproduction" not in p.relative_to(PUBLIC_PHASE).parts
        and "__pycache__" not in p.relative_to(PUBLIC_PHASE).parts
    )
    required = [PUBLIC_PHASE/name for name in TEST_FILES]
    required += [PUBLIC_PHASE/name for name in (
        "experiment/study.py", "experiment/fault_study.py",
        "analytic-dependence/analytic_dependence.py", "rus/rus_identifiability.py",
        "analytic-dependence/analyze_and_plot.py",
        "baseline/verify_baseline.py",
        "theory/harmonic_reference_checks.py", "experiment/plot_results.py",
        "experiment/plot_fault_results.py",
    )]
    required += [PUBLIC_PHASE/path for path in PROTOCOLS]
    missing = [str(p.relative_to(PUBLIC_PHASE)) for p in required if not p.is_file()]
    if missing:
        raise FileNotFoundError("Incomplete public replay payload: "+", ".join(missing))

    confirmation = json.loads((PUBLIC_PHASE/PROTOCOLS[0]).read_text())
    fault = json.loads((PUBLIC_PHASE/PROTOCOLS[1]).read_text())
    main_hash = sha256(PUBLIC_PHASE/"experiment/study.py")
    fault_hash = sha256(PUBLIC_PHASE/"experiment/fault_study.py")
    if (confirmation["script_sha256"] != main_hash
            or fault["study_code_sha256"] != main_hash
            or fault["code_sha256"] != fault_hash):
        raise ValueError("Published frozen script/protocol hashes do not match")

    # Read the fixed Git tree, independently of working-tree modifications.
    # The archive is held briefly in memory and is not shipped or retained.
    archive = git(repository, "archive", "--format=tar", PINNED_COMMIT)
    destination.mkdir(parents=True)
    upstream = destination/"work/rnglab-main-snapshot"
    upstream.mkdir(parents=True)
    with tarfile.open(fileobj=io.BytesIO(archive), mode="r:") as bundle:
        bundle.extractall(upstream, filter="data")

    phase = destination/"outputs/phase2-20261030"
    copied = {}
    for source in scripts+[PUBLIC_PHASE/path for path in PROTOCOLS]:
        relative = source.relative_to(PUBLIC_PHASE)
        target = phase/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        source_hash, copied_hash = sha256(source), sha256(target)
        if source_hash != copied_hash:
            raise RuntimeError("Bytewise copy verification failed: "+str(relative))
        copied[str(relative)] = {"sha256":copied_hash, "bytes":target.stat().st_size}

    # Explicit local configuration prevents the repository's testpaths and
    # pythonpath settings from redirecting the standalone Phase II test suite.
    configuration = """[pytest]
addopts = --import-mode=importlib
pythonpath =
    outputs/phase2-20261030/experiment
    outputs/phase2-20261030/baseline
    outputs/phase2-20261030/rus
    outputs/phase2-20261030/analytic-dependence
testpaths =
"""+"".join("    outputs/phase2-20261030/"+name+"\n" for name in TEST_FILES)
    (destination/"pytest.ini").write_text(configuration)
    manifest = {
        "preparer":"RNGLab public Phase II replay preparation v1",
        "timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "python_version":sys.version,
        "pinned_upstream_commit":PINNED_COMMIT,
        "git_archive_sha256":hashlib.sha256(archive).hexdigest(),
        "archived_source_file_count":sum(p.is_file() for p in upstream.rglob("*")),
        "copied_files":copied,
        "pytest_config_sha256":sha256(destination/"pytest.ini"),
        "phase2_test_cases_expected":35,
        "simulations_executed":False,
        "private_phase1_inputs_required":False,
        "source_payload_modified":False,
        "notes":[
            "Only source files and frozen confirmation/fault protocols were copied; results must be regenerated.",
            "The upstream materialization is a pinned Git archive, not the current checkout or an editable source install.",
            "Historical validate_artifacts.py also needs omitted Phase I inputs and original large trace archives.",
        ],
    }
    (destination/"public-replay-preparation.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return {"scratch_root":str(destination), "phase2_root":str(phase),
            "pinned_upstream_commit":PINNED_COMMIT, "copied_files":len(copied),
            "simulations_executed":False,
            "test_cwd":str(destination),
            "test_command":"python -m pytest -c pytest.ini -v"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, required=True,
                        help="new scratch root; must not already exist")
    parser.add_argument("--repository", type=Path, default=DEFAULT_REPOSITORY,
                        help="local Git repository containing pinned dd96adba (default: this checkout)")
    args = parser.parse_args()
    try:
        result = prepare(args.repository, args.destination)
    except (OSError, RuntimeError, ValueError, KeyError, tarfile.TarError) as exc:
        parser.exit(1, "Preparation failed: "+str(exc)+"\n")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
