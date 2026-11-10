"""Shared loaders for the v0.19 audit tests.

The independent audit code lives in research/v019-audit/ and is loaded by path,
like the other research modules. The Stage D measurement layer is loaded only
if it exists: research/v019-d-pilot/measure.py, or the file named by the
environment variable RNGLAB_DPILOT_MEASURE (its directory must hold dlib.py).
Tests that need it skip otherwise.
"""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
AUDIT_NAMESPACE = "rnglab-v019-audit-smoke"  # audit engineering runs only; never a pilot or smoke namespace


def load(name: str, path: Path):
    if name in sys.modules:
        return sys.modules[name]
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def recount():
    return load("v019_audit_recount", REPO / "research/v019-audit/recount.py")


def replay():
    return load("v019_audit_replay", REPO / "research/v019-audit/replay.py")


def dlib():
    return load("v019_audit_dlib", REPO / "research/v019-d-pilot/dlib.py")


def h1lib():
    return load("v019_audit_h1lib", REPO / "research/v018-h1-detection/h1lib.py")


def measure_path() -> Path | None:
    override = os.environ.get("RNGLAB_DPILOT_MEASURE")
    path = Path(override) if override else REPO / "research/v019-d-pilot/measure.py"
    return path if path.is_file() and path.with_name("dlib.py").is_file() else None


def measure():
    path = measure_path()
    if path is None:
        pytest.skip("Stage D measurement layer (measure.py) is not present in this tree")
    return load("v019_audit_measure", path)
