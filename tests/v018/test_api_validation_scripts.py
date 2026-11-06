"""Plumbing tests for the Phase 2 campaign scripts (research/v018-api-validation/).

Uses only the smoke namespace and a temporary output folder, at n = 4096 (8
batches of 512) with 4 replicas per condition. These tests check the scripts,
not coverage. They never call run_all with verify=True, which would start the
real campaign once its manifest exists.
"""
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

from rnglab14.models.pimc_harmonic import exact_finite_bead_energy

# With 4 replicas a bootstrap resample can repeat one replica, giving zero variance.
pytestmark = pytest.mark.filterwarnings("ignore:divide by zero:RuntimeWarning")

DIR = Path(__file__).resolve().parents[2] / "research" / "v018-api-validation"
SMOKE = "rnglab-v018-api-validation-smoke"
N, R = 4096, 4


def load(name):
    spec = importlib.util.spec_from_file_location(f"v018_api_{name}", DIR / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module  # worker processes unpickle tasks by module name
    spec.loader.exec_module(module)
    return module


rc, an, rcnt = load("run_campaign"), load("analyze"), load("recount")


def test_frozen_protocol_matches_script_constants():
    protocol = rc.load_protocol()
    assert protocol["data"]["namespace"] == "rnglab-v018-api-validation" != SMOKE
    assert protocol["data"]["replicas_per_condition"] == 4096


def test_recount_reference_is_bitwise_the_package_value():
    assert rcnt.reference_energy() == exact_finite_bead_energy(16, 2.0)


def test_local_and_refresh_share_metropolis_seeds_only():
    local, refresh, iid = (rc.seeds_for(SMOKE, c, 7) for c in ("local", "refresh", "iid"))
    assert set(local) == {"start", "proposal", "acceptance"}
    assert refresh == {**local, "centroid": refresh["centroid"]}
    assert len({*local.values(), refresh["centroid"], iid["paths"]}) == 5


@pytest.fixture(scope="module")
def campaign(tmp_path_factory):
    out = tmp_path_factory.mktemp("api-validation")
    rc.run_all(workers=2, out=out, ns=SMOKE, replicas=R, n=N, block_size=2, verify=False)
    rc.replay(out, SMOKE, N)
    gates = an.gates(out=out, ns=SMOKE, n=N, resamples=1000)
    (out / "results" / "gates.json").write_text(json.dumps(gates, indent=1))
    return out, gates


def test_blocks_merge_into_complete_tables_and_are_not_regenerated(campaign):
    out, _ = campaign
    for condition in rc.CONDITIONS:
        rows = an.load_rows(out / "results", condition)
        assert [int(r["replica"]) for r in rows] == list(range(R))
        before = (out / "traces" / condition / "block-0000.npz").stat().st_mtime_ns
        assert rc.run_block(SMOKE, condition, 0, R, out, N, 2).startswith("skip")
        assert (out / "traces" / condition / "block-0000.npz").stat().st_mtime_ns == before


def test_implementation_checks_pass_on_a_healthy_run(campaign):
    _, gates = campaign
    assert {k: v["result"] for k, v in gates["implementation"].items()} == dict.fromkeys(("I-1", "I-2", "I-3", "I-4"), "PASS")
    assert gates["verdicts"]["implementation_checks_pass"] is True
    assert gates["primary"]["refresh"]["units"] == R


def test_independent_recount_agrees_with_the_analysis(campaign):
    out, _ = campaign
    result = rcnt.recount(out, ns=SMOKE, resamples=1000)
    assert result["all_agree"], [c for c in result["comparisons"] if not c["agree"]]


def test_secondary_bm512_reproduces_the_api_flags(campaign):
    out, _ = campaign
    secondary = an.secondary(out, n=N)
    for condition in rc.CONDITIONS:
        assert secondary[condition]["cv_BM512_longrun_matches_api_flags"] == R
        assert sum(secondary[condition]["diagnostic"]["status_counts"].values()) == R


def test_gate_logic_on_synthetic_primaries():
    calibrated = {"cp95": [0.93, 0.96], "variance_ratio": {"ci95": [0.9, 1.1]}, "bias": {"ci95": [-0.001, 0.001]}}
    assert an.required_gate(calibrated)["result"] == "PASS"
    for key, bad in (("cp95", [0.91, 0.95]), ("variance_ratio", {"ci95": [0.7, 1.0]}), ("bias", {"ci95": [-0.02, 0.0]})):
        assert an.required_gate({**calibrated, key: bad})["result"] == "FAIL"


def test_bootstrap_is_reproducible_and_chunk_independent_of_total():
    x = np.random.default_rng(1).normal(size=50)
    first = an.bootstrap(lambda a: a.mean(axis=-1), [x], 123, 1000)
    assert first == an.bootstrap(lambda a: a.mean(axis=-1), [x], 123, 1000)
    assert first == rcnt.percentile_bootstrap(lambda a: a.mean(axis=-1), [x], 123, 1000)
