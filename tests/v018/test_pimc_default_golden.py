"""Bit-identical pin of the pre-v0.18 PIMC default behaviour.

The golden file tests/v018/data/pimc_default_golden.json was captured from the
committed, unpatched src/rnglab14/models/pimc_harmonic.py (its commit and
SHA-256 are recorded in the file) by running this module as a script:

    PYTHONPATH=src python tests/v018/test_pimc_default_golden.py --capture

Each case calls run_simulation without any v0.18 option. The test requires every
pre-existing PIMCResult field, the final model path, and each role source's byte
counter and stream digest to match the capture exactly (floats via repr, which
round-trips every bit). Fields added later are ignored here. Written by Claude
Code for the v0.18 Phase 2 API extension (research/v018-api-validation/).
"""
from dataclasses import asdict, fields
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

import numpy as np

from rnglab14.models.pimc_harmonic import PathIntegralHarmonicOscillator, PIMCResult
from rnglab14.sources import PCG64Source, StochasticChannels
from rnglab14.sources.faults import DegradedFaultSource
from rnglab14.sources.reference import Xoshiro256PlusPlusSource

GOLDEN = Path(__file__).resolve().parent / "data" / "pimc_default_golden.json"
MODULE = Path("src/rnglab14/models/pimc_harmonic.py")

LEGACY_FIELDS = (
    "beads_M", "temperature_T", "beta", "total_sweeps", "burn_in", "acceptance_rate",
    "virial_energy_mean", "virial_energy_std", "thermodynamic_energy_mean",
    "exact_quantum_energy", "energy_bias", "radius_of_gyration", "provenance_hash",
    "exact_finite_bead_energy", "discretization_bias", "virial_uncertainty",
    "thermodynamic_uncertainty", "sampler_protocol", "centroid_update",
)

# Model settings, run lengths and source layouts. "phase2_*" cases use the
# Phase 2 model (M=16, beta=2, delta=0.4); "applications" mirrors the call in
# src/rnglab15/applications.py.
CASES = {
    "single_source_M4": dict(beads=4, temperature=0.5, delta=0.5, sweeps=64, burn_in=16,
                             layout="single", centroid=None),
    "applications_M4": dict(beads=4, temperature=0.5, delta=0.5, sweeps=32, burn_in=0,
                            layout="channels", centroid=None),
    "phase2_channels": dict(beads=16, temperature=0.5, delta=0.4, sweeps=256, burn_in=64,
                            layout="channels", centroid=None),
    "phase2_channels_centroid": dict(beads=16, temperature=0.5, delta=0.4, sweeps=256, burn_in=0,
                                     layout="channels", centroid="separate"),
    "aliased_centroid": dict(beads=8, temperature=1.0, delta=0.5, sweeps=96, burn_in=32,
                             layout="single", centroid="aliased"),
    "faulted_mixed_sources": dict(beads=16, temperature=0.5, delta=0.4, sweeps=128, burn_in=32,
                                  layout="faulted", centroid="separate"),
}


def _seed(case, role):
    return int.from_bytes(hashlib.sha256(f"rnglab-v018-tests/pimc-default-golden/{case}/{role}".encode()).digest()[:16], "big")


def _sources(case, layout, centroid):
    if layout == "single":
        main = PCG64Source(_seed(case, "single"))
        roles = {"proposal": main, "acceptance": main, "initialization": main, "resampling": None}
        argument = main
    else:
        if layout == "channels":
            proposal, acceptance = PCG64Source(_seed(case, "proposal")), PCG64Source(_seed(case, "acceptance"))
        else:
            proposal = DegradedFaultSource(PCG64Source(_seed(case, "proposal")), bias_delta=0.02)
            acceptance = DegradedFaultSource(Xoshiro256PlusPlusSource(_seed(case, "acceptance") >> 64),
                                             autocorrelation_rho=0.1)
        argument = StochasticChannels(proposal, acceptance, PCG64Source(_seed(case, "initialization")),
                                      PCG64Source(_seed(case, "resampling")))
        roles = {"proposal": argument.proposal, "acceptance": argument.acceptance,
                 "initialization": argument.initialization, "resampling": argument.resampling}
    if centroid == "separate":
        roles["centroid"] = PCG64Source(_seed(case, "centroid"))
    elif centroid == "aliased":
        roles["centroid"] = roles["proposal"]
    else:
        roles["centroid"] = None
    return argument, roles


def _canonical(value):
    """JSON-ready copy; drops construction-time created_at stamps."""
    if isinstance(value, dict):
        return {str(k): _canonical(v) for k, v in value.items() if k != "created_at"}
    if isinstance(value, (list, tuple)):
        return [_canonical(v) for v in value]
    if isinstance(value, np.ndarray):
        return [_canonical(v) for v in value.tolist()]
    if isinstance(value, np.generic):
        return value.item()
    return value


def run_case(case):
    spec = CASES[case]
    argument, roles = _sources(case, spec["layout"], spec["centroid"])
    model = PathIntegralHarmonicOscillator(beads=spec["beads"], temperature=spec["temperature"],
                                           displacement_delta=spec["delta"])
    kwargs = {"centroid_source": roles["centroid"]} if roles["centroid"] is not None else {}
    result = model.run_simulation(argument, sweeps=spec["sweeps"], burn_in=spec["burn_in"], **kwargs)
    record = asdict(result)
    return _canonical({
        "result": {name: record[name] for name in LEGACY_FIELDS},
        "final_path": model.path,
        "sources_after": {role: (None if source is None else
                                 {"consumed_bytes": source.consumed_bytes, "stream_digest": source.stream_digest})
                          for role, source in roles.items()},
    })


def _dumps(value):
    return json.dumps(value, sort_keys=True, allow_nan=True)


def test_legacy_fields_keep_their_names_and_order():
    assert tuple(f.name for f in fields(PIMCResult))[:len(LEGACY_FIELDS)] == LEGACY_FIELDS


def test_default_calls_are_bit_identical_to_the_pre_v018_capture():
    golden = json.loads(GOLDEN.read_text())
    assert set(golden["cases"]) == set(CASES)
    for case in CASES:
        assert _dumps(run_case(case)) == _dumps(golden["cases"][case]), case


def _capture(force):
    if GOLDEN.exists() and not force:
        sys.exit(f"{GOLDEN} exists; pass --force to overwrite")
    root = Path(__file__).resolve().parents[2]
    git = lambda *args: subprocess.run(["git", *args], cwd=root, capture_output=True, text=True)
    if git("diff", "--quiet", "HEAD", "--", str(MODULE)).returncode != 0:
        sys.exit(f"{MODULE} has uncommitted changes; capture only from committed code")
    GOLDEN.parent.mkdir(exist_ok=True)
    GOLDEN.write_text(json.dumps({
        "captured_by": "PYTHONPATH=src python tests/v018/test_pimc_default_golden.py --capture",
        "source_commit": git("rev-parse", "HEAD").stdout.strip(),
        "module": str(MODULE),
        "module_sha256": hashlib.sha256((root / MODULE).read_bytes()).hexdigest(),
        "python": platform.python_version(), "numpy": np.__version__,
        "legacy_fields": list(LEGACY_FIELDS),
        "cases": {case: run_case(case) for case in CASES},
    }, indent=1, sort_keys=True, allow_nan=True) + "\n")
    print(f"wrote {GOLDEN}")


if __name__ == "__main__":
    if "--capture" in sys.argv:
        _capture(force="--force" in sys.argv)
    else:
        sys.exit("usage: PYTHONPATH=src python tests/v018/test_pimc_default_golden.py --capture [--force]")
