"""v0.18 namespace: re-exports are the original objects; legacy package intact."""
import importlib
import math

import pytest

import rnglab
from rnglab import faults, models, monitors, references, sources, stats

ORIGINS = {
    sources: ["rnglab14.sources.base", "rnglab14.sources.reference"],
    faults: ["rnglab14.sources.faults"],
    models: ["rnglab14.models.pimc_harmonic", "rnglab14.models.vmc_harmonic",
             "rnglab14.models.ising_2d", "rnglab14.models.decision_support", "rnglab16.core"],
    references: ["rnglab14.models.pimc_harmonic", "rnglab16.core"],
    monitors: ["rnglab14.monitoring.calibration", "rnglab14.standards.health_90b", "rnglab17.endpoints"],
    stats: ["rnglab12.uncertainty", "rnglab14.models.pimc_harmonic",
            "rnglab14.statistics.multitest", "rnglab17.endpoints"],
}


@pytest.mark.parametrize("module", list(ORIGINS), ids=lambda m: m.__name__)
def test_every_export_is_an_original_object(module):
    originals = [vars(importlib.import_module(name)) for name in ORIGINS[module]]
    for name in module.__all__:
        exported = getattr(module, name)
        assert any(exported is value for namespace in originals for value in namespace.values()), name


def test_legacy_package_metadata_unchanged():
    assert rnglab.__version__ == "0.1.0"
    assert rnglab.__doc__ == "Experimental framework V0.1. All default entropy is simulated."


def test_phase2_finite_bead_target():
    # research/phase2-20261030/experiment/confirmation-v1/summary.csv, column "target".
    assert math.isclose(references.exact_finite_bead_energy(16, 2.0), 0.6554741522981742, rel_tol=1e-14)


def test_v017_exact_susceptibility_reference():
    # SUSCEPTIBILITY_PROBE.json (v0.17), row ising_L4_T0, field "reference".
    exact = references.ising_exact_enumeration(4, 1.5)
    assert math.isclose(exact["susceptibility"], 10.403226593514693, rel_tol=1e-12)
