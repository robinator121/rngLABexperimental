"""RNGLab v0.20 nuclear module: experiment manifests.

A run manifest records what is needed to repeat or audit one simulation: the
code (git commit, dirty flag, SHA-256 of every module file), the environment,
the random stream (source identity and seed label, offset, words and bytes
consumed, digest of the consumed bytes), the configuration and the output
checksums. A digest identifies a stream; it is not an entropy or independence
certificate.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any, Mapping

import numpy as np
import scipy

from .decay import DecayRun
from .streams import WordStream

MANIFEST_SCHEMA = "rnglab-nuclear-run-manifest-v1"
PACKAGE_DIR = Path(__file__).resolve().parent


def json_sha256(obj: Any) -> str:
    """SHA-256 of canonical JSON (sorted keys, no spaces, no NaN)."""
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()


def file_sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def code_identity(repo_root: Path | None = None) -> dict[str, Any]:
    """Git commit and dirty flag of repo_root (None where git is unavailable) and module file hashes."""
    files = {p.name: file_sha256(p) for p in sorted(PACKAGE_DIR.glob("*.py"))}
    identity: dict[str, Any] = {"module_files_sha256": files, "git_commit": None, "git_dirty": None}
    root = Path(repo_root) if repo_root is not None else PACKAGE_DIR.parents[2]
    try:
        commit = subprocess.run(["git", "-C", str(root), "rev-parse", "HEAD"], capture_output=True,
                                text=True, check=True).stdout.strip()
        status = subprocess.run(["git", "-C", str(root), "status", "--porcelain"], capture_output=True,
                                text=True, check=True).stdout
        identity.update(git_commit=commit, git_dirty=bool(status.strip()))
    except (OSError, subprocess.CalledProcessError):
        pass
    return identity


def environment() -> dict[str, str]:
    return {"python": sys.version.split()[0], "numpy": np.__version__, "scipy": scipy.__version__,
            "platform": platform.platform(), "machine": platform.machine()}


def run_manifest(*, experiment: str, run_id: str, run: DecayRun, stream: WordStream,
                 initial: Mapping[str, int], code: Mapping[str, Any], extra: Mapping[str, Any] | None = None
                 ) -> dict[str, Any]:
    """Manifest of one Layer A run (configuration, stream provenance, output checksums)."""
    if stream.words_read - run.first_word < run.words_used:
        raise ValueError("the stream did not supply the words this run used")
    return {
        "schema": MANIFEST_SCHEMA,
        "experiment": experiment,
        "run_id": run_id,
        "code": dict(code),
        "environment": environment(),
        "rng": {"role": "physics", "layout": run.layout, "words_per_nucleus": run.words_per_nucleus,
                "first_word_of_run": run.first_word, **stream.provenance()},
        "configuration": {"scheme": run.scheme.to_dict(), "scheme_sha256": run.scheme.sha256(),
                          "initial": dict(initial), "horizon": repr(run.horizon), "engine": run.engine},
        "outputs": run.summary(),
        "extra": dict(extra or {}),
    }
