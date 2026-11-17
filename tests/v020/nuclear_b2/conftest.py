"""Shared synthetic schemes for the B2 (aggregated decay) tests.

Tests draw only from the namespace below, never from a reserved study
namespace. tests/v020/nuclear/test_layer_a_namespaces.py (frozen) scans every
test file for the top-level namespaces of every protocol; this package adds
(test_b2_guards.py) a scan for every namespace a protocol names at any depth
(calibration record and development namespaces included) and, below, a guard
at run time: while a B2 test runs, namespace_seed refuses every reserved
namespace, however a test reaches it (a script's default, a protocol field).
"""
import json
import math
import sys
from pathlib import Path

import pytest

import rnglab.seeds

from rnglab.nuclear import Channel, DecayScheme, Mode, Nuclide
from rnglab.nuclear.precursors import PrecursorGroup, precursor_scheme

B2_TEST_NAMESPACE = "rnglab-v020-nuclear-b2-test"
REPO = Path(__file__).resolve().parents[3]
_ORIGINAL_NAMESPACE_SEED = rnglab.seeds.namespace_seed


def reserved_namespaces() -> set[str]:
    """Every namespace any research/v020-nuclear-*/protocol.json names, at any depth."""
    found: set[str] = set()

    def walk(value):
        if isinstance(value, dict):
            for key, item in value.items():
                if "namespace" in key and isinstance(item, str):
                    found.add(item)
                walk(item)
        elif isinstance(value, list):
            for item in value:
                walk(item)

    for protocol in REPO.glob("research/v020-nuclear-*/protocol.json"):
        walk(json.loads(protocol.read_text()))
    return found


class ReservedNamespaceError(AssertionError):
    pass


@pytest.fixture(autouse=True)
def _refuse_reserved_namespaces(monkeypatch):
    """While a B2 test runs, namespace_seed raises for a reserved namespace, wherever it was imported."""
    reserved = reserved_namespaces()

    def guarded(namespace, *parts):
        if isinstance(namespace, str) and any(namespace == r or namespace.startswith(r + "/") for r in reserved):
            raise ReservedNamespaceError(f"a test asked for a seed in the reserved namespace {namespace!r}")
        return _ORIGINAL_NAMESPACE_SEED(namespace, *parts)

    for module in list(sys.modules.values()):
        if getattr(module, "namespace_seed", None) is _ORIGINAL_NAMESPACE_SEED:
            monkeypatch.setattr(module, "namespace_seed", guarded)
    yield


def branching_scheme() -> DecayScheme:
    n, c = Nuclide, Channel
    return DecayScheme("b2-branching", [
        n("Xa-400", 150, 400, 2.0, (c("Xb-396", Mode.ALPHA, 0.64, 0.25), c("Xc-400m", Mode.BETA_MINUS, 0.36, 0.3))),
        n("Xb-396", 148, 396, 0.5, (c("Xd-396", Mode.ELECTRON_CAPTURE, 1.0, 0.8),)),
        n("Xc-400m", 151, 400, 0.8, (c("Xc-400", Mode.ISOMERIC_TRANSITION, 1.0, 1.0),), level=1),
        n("Xc-400", 151, 400, 3.0, (c("Xe-396", Mode.ALPHA, 1.0),)),
        n("Xd-396", 147, 396, 1.25, (c("Xf-392", Mode.ALPHA, 1.0),)),
        n("Xe-396", 149, 396, 1.6, (c("Xg-396", Mode.BETA_MINUS, 0.9, 0.5),
                                    c("Xh-395", Mode.BETA_DELAYED_NEUTRON, 0.1))),
        n("Xf-392", 145, 392, math.inf),
        n("Xg-396", 150, 396, math.inf),
        n("Xh-395", 150, 395, math.inf),
    ])


def competing_scheme() -> DecayScheme:
    n, c = Nuclide, Channel
    return DecayScheme("b2-competing", [
        n("Pa-300", 145, 300, 1.0, (c("Pb-300", Mode.BETA_PLUS, 0.3, 0.2), c("Pb-300", Mode.ELECTRON_CAPTURE, 0.7, 0.6))),
        n("Pb-300", 144, 300, 0.4, (c("Pc-296", Mode.ALPHA, 0.5), c("Pd-300", Mode.BETA_PLUS, 0.3),
                                    c("Pe-299", Mode.BETA_DELAYED_NEUTRON, 0.2, 1.0))),
        n("Pc-296", 142, 296, math.inf),
        n("Pd-300", 143, 300, 2.5, (c("Pf-296", Mode.ALPHA, 1.0, 1.0),)),
        n("Pe-299", 145, 299, math.inf),
        n("Pf-296", 141, 296, math.inf),
    ])


@pytest.fixture
def ns():
    return B2_TEST_NAMESPACE


@pytest.fixture
def branching():
    return branching_scheme()


@pytest.fixture
def competing():
    return competing_scheme()


@pytest.fixture
def precursors():
    return precursor_scheme([PrecursorGroup("G1", 0.6, 0.3), PrecursorGroup("G2", 2.5, 0.1),
                             PrecursorGroup("G3", 9.0, 1.0)])
