"""No test may draw from a frozen study namespace of the nuclear module.

Frozen namespaces are read from every research/v020-nuclear-*/protocol.json;
none of them may occur in a test file, and the shared test namespace must
differ from all of them.
"""
import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]


def _frozen_namespaces():
    names = set()
    for protocol in REPO.glob("research/v020-nuclear-*/protocol.json"):
        data = json.loads(protocol.read_text())
        names.update(data[key] for key in ("namespace", "smoke_namespace") if key in data)
    return names


def test_frozen_namespaces_exist_and_differ_from_the_test_namespace(ns):
    frozen = _frozen_namespaces()
    assert "rnglab-v020-nuclear-a-verify" in frozen
    assert all(not ns.startswith(name) and not name.startswith(ns) for name in frozen)


def test_no_test_file_mentions_a_frozen_namespace():
    frozen = _frozen_namespaces()
    for path in (REPO / "tests").rglob("*.py"):
        if path.resolve() == Path(__file__).resolve():
            continue
        text = path.read_text()
        assert not any(name in text for name in frozen), path
