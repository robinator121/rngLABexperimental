"""End-to-end plumbing test for the Stage C campaign scripts (research/v018-h1-detection/).

A tiny campaign under the smoke namespace in a temporary folder: 4 calibration,
3 healthy and 2 runs in each of two fault cells, full length (10,240 sweeps).
It checks the scripts, not the science. run_all is never called with
verify=True, which would start the real campaign once its manifest exists.
"""
import importlib.util
import json
import sys
from pathlib import Path

import pytest

DIR = Path(__file__).resolve().parents[2] / "research" / "v018-h1-detection"
SMOKE = "rnglab-v018-h1-detection-smoke"
pytestmark = pytest.mark.filterwarnings("ignore::RuntimeWarning")


def load(name):
    spec = importlib.util.spec_from_file_location(f"v018_h1_{name}", DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


rc, an, rcnt = load("run_campaign"), load("analyze"), load("recount")
PLAN = [{"cell": "calibration", "condition": "calibration", "family": "healthy", "severity": 0, "replicas": 4},
        {"cell": "healthy", "condition": "test", "family": "healthy", "severity": 0, "replicas": 3},
        {"cell": "SHIFT-0", "condition": "test", "family": "SHIFT", "severity": 0.05, "replicas": 2},
        {"cell": "GFSR-0", "condition": "test", "family": "GFSR", "severity": [31, 13], "replicas": 2}]


@pytest.fixture(scope="module")
def campaign(tmp_path_factory):
    out = tmp_path_factory.mktemp("h1-campaign")
    rc.run_all(workers=4, out=out, ns=SMOKE, plan=PLAN, block_size=2, verify=False)
    replay = rc.replay(out, SMOKE, PLAN)
    gates = an.gates(out, ns=SMOKE, order_statistic=4, resamples=200)
    (out / "results" / "gates.json").write_text(json.dumps(gates))
    return out, replay, gates


def test_frozen_protocol_and_grid_build_the_planned_cells():
    protocol = rc.load_protocol()
    plan = rc.campaign_plan(protocol, json.loads((DIR / "grid.json").read_text()))
    assert len(plan) == 26 and sum(c["replicas"] for c in plan) == 999 + 4096 + 24 * 256


def test_merge_writes_complete_cells_and_blocks_are_not_regenerated(campaign):
    out, _, _ = campaign
    for cell in PLAN:
        records = (out / "results" / "cells" / f"{cell['cell']}.jsonl").read_text().splitlines()
        assert [json.loads(r)["replica"] for r in records] == list(range(cell["replicas"]))
    assert rc.run_block(SMOKE, PLAN[2], 0, out, 2).startswith("skip")


def test_replay_and_independent_fault_maps_agree(campaign):
    _, replay, gates = campaign
    assert all(v["I-4"] for v in replay.values())
    assert replay["SHIFT-0"]["I-3"] and replay["GFSR-0"]["I-3"]
    assert {k: v["result"] for k, v in gates["implementation"].items()} == dict.fromkeys(("I-1", "I-2", "I-3", "I-4"), "PASS")


def test_gates_have_every_cell_and_a_class(campaign):
    _, _, gates = campaign
    assert set(gates["cells"]) == {"SHIFT-0", "GFSR-0"}
    assert all(c["class"] in {"H1", "undetected_harm", "conventional_first", "both_detect_harmed", "harmless_quiet", "indeterminate"}
               for c in gates["cells"].values())
    assert set(gates["verdicts"]["H1-C"]) == {"SHIFT", "GFSR"}


def test_independent_recount_agrees(campaign):
    out, _, _ = campaign
    result = rcnt.recount(out, ns=SMOKE, resamples=200)
    assert result["all_agree"], [c for c in result["comparisons"] if not c["agree"]]


def test_secondaries_run_on_the_traces(campaign):
    out, _, gates = campaign
    s = an.secondary(out, gates)
    assert set(s["cells"]) == {"healthy", "SHIFT-0", "GFSR-0"}
    assert all(set(c["post_onset_alarms"]) == {"A1", "A2", "shape", "cusum_D"} for c in s["cells"].values())
