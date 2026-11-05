"""Public PIMC initialization, retained trace and audited role-ledger checks."""
from dataclasses import asdict
import hashlib
import math

import numpy as np
import pytest

from rnglab12.uncertainty import summarize_series
from rnglab14.models.pimc_harmonic import (
    PathIntegralHarmonicOscillator, harmonic_energy_estimators,
    sample_exact_finite_bead_paths,
)
from rnglab14.sources import PCG64Source, StochasticChannels


def sources(seed=18001):
    return StochasticChannels(PCG64Source(seed), PCG64Source(seed+1),
                              PCG64Source(seed+2), PCG64Source(seed+3))


@pytest.mark.parametrize("option", [
    {"initialization":"bad"}, {"initialization":None},
    {"record_trace":1}, {"record_trace":None},
    {"uncertainty_batch_sizes":None}, {"uncertainty_batch_sizes":2},
    {"uncertainty_batch_sizes":()}, {"uncertainty_batch_sizes":(2,2)},
    {"uncertainty_batch_sizes":(1,)}, {"uncertainty_batch_sizes":(True,)},
    {"uncertainty_batch_sizes":(2.,)}, {"uncertainty_batch_sizes":([2],)},
    {"initial_path":[0.,1.]}, {"initial_path":[[0.,1.,2.,3.]]},
    {"initial_path":[0.,1.,2.,math.nan]}, {"initial_path":[0.,1.,2.,math.inf]},
])
def test_invalid_options_fail_before_any_role_consumption(option):
    channel = sources()
    centroid = PCG64Source(19000)
    model = PathIntegralHarmonicOscillator(beads=4)
    with pytest.raises((ValueError, TypeError)):
        model.run_simulation(channel, sweeps=64, burn_in=16, centroid_source=centroid, **option)
    assert [s.consumed_bytes for s in (channel.proposal, channel.acceptance,
        channel.initialization, channel.resampling, centroid)] == [0]*5
    assert np.array_equal(model.path, np.zeros(4))


def test_exact_initialization_uses_existing_two_path_mapping_and_records_discard():
    channel = sources()
    independent_same_source = PCG64Source(18003)
    expected = sample_exact_finite_bead_paths(independent_same_source, 4, 2., 2)[0]
    result = PathIntegralHarmonicOscillator(beads=4).run_simulation(
        channel, sweeps=64, burn_in=16, initialization="exact_gaussian", record_trace=True)
    assert np.array_equal(result.trace["initial_path"], expected)
    assert channel.initialization.consumed_bytes == independent_same_source.consumed_bytes == 64
    assert channel.initialization.stream_digest == independent_same_source.stream_digest
    assert result.initialization_protocol["mode"] == "exact_gaussian"
    assert result.initialization_protocol["reference_paths_drawn"] == 2
    assert result.initialization_protocol["reference_paths_discarded"] == 1
    assert result.initialization_protocol["source_bytes_consumed"] == 64
    assert result.initialization_protocol["independence_verified"] is False


def test_exact_initialization_covariance_against_independent_harmonic_precision():
    n, beads, beta = 1024, 4, 2.
    initial = []
    for i in range(n):
        seed = int.from_bytes(hashlib.sha256(f"pimc-exact-init-test/{i}".encode()).digest()[:16],"big")
        result = PathIntegralHarmonicOscillator(beads=beads).run_simulation(
            sources(seed), sweeps=32, burn_in=0,
            initialization="exact_gaussian", record_trace=True, uncertainty_batch_sizes=(2,))
        initial.append(result.trace["initial_path"])
    initial = np.asarray(initial)
    eye = np.eye(beads)
    difference = np.roll(eye, -1, axis=0)-eye
    covariance = np.linalg.inv(beta/beads*eye+beads/beta*(difference.T @ difference))
    # Gaussian sample-mean and Wishart covariance standard errors; deliberately
    # loose six-SE bounds validate the initializer, not interval coverage.
    mean_se = np.sqrt(np.diag(covariance)/n)
    covariance_se = np.sqrt((covariance**2+np.outer(np.diag(covariance),np.diag(covariance)))/(n-1))
    assert np.all(abs(initial.mean(0)) < 6*mean_se)
    assert np.all(abs(np.cov(initial,rowvar=False)-covariance) < 6*covariance_se)


def test_explicit_path_is_copied_and_overrides_initializer_without_source_reads():
    path = np.array([-.6,.1,.4,.9])
    before = path.copy()
    channel = sources()
    model = PathIntegralHarmonicOscillator(beads=4)
    result = model.run_simulation(channel, sweeps=64, burn_in=16,
        initialization="exact_gaussian", initial_path=path, record_trace=True)
    assert np.array_equal(path, before)
    assert np.array_equal(result.trace["initial_path"], before)
    assert not np.shares_memory(result.trace["initial_path"], path)
    assert not np.shares_memory(result.trace["paths"], model.path)
    assert channel.initialization.consumed_bytes == 0
    assert channel.initialization.stream_digest == hashlib.sha256(b"").hexdigest()
    assert result.initialization_protocol["mode"] == "explicit_path"
    assert result.initialization_protocol["reference_paths_drawn"] == 0
    result.trace["initial_path"][:] = 100
    assert np.array_equal(path, before)


def test_full_retained_trace_rederives_estimators_geometry_and_requested_summaries():
    result = PathIntegralHarmonicOscillator(beads=4).run_simulation(
        sources(), sweeps=320, burn_in=64, centroid_source=PCG64Source(23000),
        initialization="exact_gaussian", record_trace=True,
        uncertainty_batch_sizes=(8,32,128))
    trace = result.trace
    assert set(trace) == {"initial_path","paths","TD","CV","centroid","rg2"}
    assert trace["paths"].shape == (256,4)
    for key in ("TD","CV","centroid","rg2"):
        assert trace[key].shape == (256,)
    td, cv = harmonic_energy_estimators(trace["paths"],2.)
    np.testing.assert_allclose(trace["TD"],td,rtol=0,atol=1e-14)
    np.testing.assert_allclose(trace["CV"],cv,rtol=0,atol=1e-14)
    np.testing.assert_allclose(trace["centroid"],trace["paths"].mean(1),rtol=0,atol=1e-15)
    np.testing.assert_allclose(trace["rg2"],((trace["paths"]-trace["centroid"][:,None])**2).mean(1),rtol=0,atol=1e-15)
    for values, summary in [(trace["TD"],result.thermodynamic_uncertainty),
                             (trace["CV"],result.virial_uncertainty)]:
        assert summary == summarize_series(values, reference=result.exact_finite_bead_energy, batch_sizes=(8,32,128))
        assert set(summary["batch_means"]) == {"8","32","128"}
    assert result.virial_energy_mean == round(float(trace["CV"].mean()),6)
    assert result.thermodynamic_energy_mean == round(float(trace["TD"].mean()),6)
    assert result.radius_of_gyration == round(float(np.sqrt(trace["rg2"]).mean()),6)


def test_recording_trace_does_not_change_any_transition_or_consumption():
    a, b = sources(), sources()
    common = dict(sweeps=64,burn_in=16,initialization="exact_gaussian")
    left = PathIntegralHarmonicOscillator(beads=4).run_simulation(a,**common)
    right = PathIntegralHarmonicOscillator(beads=4).run_simulation(b,record_trace=True,**common)
    assert left.trace is None
    x, y = asdict(left), asdict(right)
    x.pop("trace"); y.pop("trace")
    for record in (x,y):
        for entry in record["source_provenance"]["roles"].values():
            if entry is not None:
                entry["manifest"].pop("created_at", None)
    assert x == y


def test_role_ledger_preserves_prior_counts_unused_roles_and_shared_object_aliases():
    channel = sources()
    centroid = PCG64Source(23000)
    all_sources = {"proposal":channel.proposal,"acceptance":channel.acceptance,
        "initialization":channel.initialization,"resampling":channel.resampling,"centroid":centroid}
    for source in all_sources.values():
        source.read_bytes(3)
    before = {role:(source.consumed_bytes,source.stream_digest) for role,source in all_sources.items()}
    result = PathIntegralHarmonicOscillator(beads=4).run_simulation(
        channel,sweeps=64,burn_in=16,centroid_source=centroid)
    for role,source in all_sources.items():
        entry = result.source_provenance["roles"][role]
        assert entry["consumed_bytes_before"] == before[role][0] == 3
        assert entry["stream_digest_before"] == before[role][1]
        assert entry["consumed_bytes_after"] == source.consumed_bytes
        assert entry["stream_digest_after"] == source.stream_digest
        assert entry["consumed_bytes_delta"] == source.consumed_bytes-3
        identity = source.manifest.to_dict()
        identity.pop("created_at")
        assert entry["manifest"] == identity
        assert "created_at" not in entry["manifest"]
        assert entry["aliases"] == []
    assert result.source_provenance["roles"]["proposal"]["consumed_bytes_delta"] == 64*4*8
    assert result.source_provenance["roles"]["initialization"]["consumed_bytes_delta"] == 4*8
    assert result.source_provenance["roles"]["centroid"]["consumed_bytes_delta"] == 64*16
    assert result.source_provenance["roles"]["resampling"]["consumed_bytes_delta"] == 0
    assert result.source_provenance["independence_verified"] is False
    shared = PCG64Source(22399)
    aliased = PathIntegralHarmonicOscillator(beads=4).run_simulation(
        shared,sweeps=32,burn_in=0,centroid_source=shared)
    roles = aliased.source_provenance["roles"]
    assert roles["resampling"] is None
    for role in ("proposal","acceptance","initialization","centroid"):
        assert set(roles[role]["aliases"]) == {"proposal","acceptance","initialization","centroid"}-{role}
        assert roles[role]["consumed_bytes_delta"] == shared.consumed_bytes


def test_uniform_default_initial_path_and_default_no_trace():
    result = PathIntegralHarmonicOscillator(beads=4).run_simulation(sources(),sweeps=32,burn_in=0)
    assert result.trace is None
    assert result.initialization_protocol["mode"] == "uniform"
    assert result.initialization_protocol["source_bytes_consumed"] == 32
    assert result.source_provenance["roles"]["centroid"] is None
    assert result.source_provenance["roles"]["resampling"]["consumed_bytes_delta"] == 0
