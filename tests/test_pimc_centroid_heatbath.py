"""Analytical invariance and source-consumption checks for optional centroid Gibbs updates."""
from dataclasses import asdict
import hashlib
import math
import struct

import numpy as np
import pytest

from rnglab14.models import PathIntegralHarmonicOscillator, harmonic_centroid_heatbath
from rnglab14.sources import PCG64Source, RandomSource, SourceBudgetExceeded, SourceManifest, StochasticChannels


class FixedBytes(RandomSource):
    def __init__(self, raw, max_bytes=None):
        super().__init__(SourceManifest("fixture", "test_fixture", True, False), max_bytes=max_bytes)
        self.raw = raw

    def _read_bytes(self, n):
        result, self.raw = self.raw[:n], self.raw[n:]
        return result


def _channels():
    return StochasticChannels(PCG64Source(701), PCG64Source(702), PCG64Source(703))


def test_exact_box_muller_mapping_internal_coordinates_and_action_factorization():
    path = np.array([-2.0, -1.0, 0.0, 7.0])
    original = path.copy()
    beta = 2.0
    # First uniform exactly 1/2, second exactly zero: c_new=sqrt(log(4)/beta).
    raw = struct.pack(">QQ", 1 << 63, 0)
    source = FixedBytes(raw)
    updated = harmonic_centroid_heatbath(path, beta, source)
    old_c, new_c = path.mean(), math.sqrt(math.log(4)/beta)
    assert updated.mean() == pytest.approx(new_c, abs=2e-15)
    assert updated-updated.mean() == pytest.approx(path-old_c, abs=2e-15)
    assert np.roll(updated, -1)-updated == pytest.approx(np.roll(path, -1)-path, abs=2e-15)
    assert np.array_equal(path, original)
    def action(x):
        return len(x)/(2*beta)*np.sum((np.roll(x, -1)-x)**2) + beta/(2*len(x))*np.sum(x*x)
    assert action(updated)-action(path) == pytest.approx(beta/2*(new_c*new_c-old_c*old_c), abs=2e-13)
    assert source.consumed_bytes == 16
    assert source.stream_digest == hashlib.sha256(raw).hexdigest()


def test_zero_uniform_is_finite_and_source_budget_is_atomic():
    source = FixedBytes(struct.pack(">QQ", 0, (1 << 64)-1))
    updated = harmonic_centroid_heatbath([4.0, 6.0], 2.0, source)
    assert updated == pytest.approx([-1.0, 1.0])
    limited = FixedBytes(b"\0"*16, max_bytes=15)
    with pytest.raises(SourceBudgetExceeded):
        harmonic_centroid_heatbath([1.0, 2.0], 2.0, limited)
    assert limited.consumed_bytes == 0
    assert limited.raw == b"\0"*16


@pytest.mark.parametrize("path,beta", [([], 2.0), (1.0, 2.0), ([[1.0]], 2.0),
    ([math.nan], 2.0), ([math.inf], 2.0), ([1.0], 0), ([1.0], -1),
    ([1.0], math.nan), ([1.0], math.inf)])
def test_invalid_domains_fail_before_consumption(path, beta):
    source = PCG64Source(22)
    with pytest.raises(ValueError):
        harmonic_centroid_heatbath(path, beta, source)
    assert source.consumed_bytes == 0


def test_explicit_source_required_and_deterministic_replay():
    with pytest.raises(TypeError, match="RandomSource"):
        harmonic_centroid_heatbath([1.0], 2.0, None)
    a, b = PCG64Source(1901), PCG64Source(1901)
    x, y = np.array([-3.0, 0.0, 2.0]), np.array([-3.0, 0.0, 2.0])
    for _ in range(15):
        x = harmonic_centroid_heatbath(x, 3.0, a)
        y = harmonic_centroid_heatbath(y, 3.0, b)
        assert np.array_equal(x, y)
    assert a.consumed_bytes == b.consumed_bytes == 240
    assert a.stream_digest == b.stream_digest


@pytest.mark.parametrize("beta", [0.2, 2.0, 8.0])
def test_centroid_moments_against_analytic_gaussian_standard_errors(beta):
    n = 12000
    source = PCG64Source(391103)
    path = np.array([-2.0, -1.0, 1.0, 2.0])
    centroids = np.array([harmonic_centroid_heatbath(path, beta, source).mean() for _ in range(n)])
    # Under the explicit IID-source assumption, Var(c)=1/beta,
    # Var(c^2)=2/beta^2 and Var(c^4)=96/beta^4. Six-SE bounds are
    # deliberately loose regression checks, not uncertainty calibration claims.
    assert abs(centroids.mean()) < 6/math.sqrt(n*beta)
    assert abs(np.mean(centroids**2)-1/beta) < 6*math.sqrt(2/n)/beta
    assert abs(np.mean(centroids**4)-3/beta**2) < 6*math.sqrt(96/n)/beta**2


def test_heatbath_preserves_independent_gaussian_target_covariance():
    beads, beta, n = 4, 2.0, 10000
    eye = np.eye(beads)
    difference = np.roll(eye, -1, axis=0)-eye
    covariance = np.linalg.inv(beta/beads*eye+beads/beta*(difference.T @ difference))
    # Independent numerical construction of the finite-action Gaussian target.
    initial = np.random.default_rng(183310).multivariate_normal(np.zeros(beads), covariance, size=n)
    source = PCG64Source(985381)
    updated = np.array([harmonic_centroid_heatbath(x, beta, source) for x in initial])
    # Exact Gaussian and Wishart marginal standard errors; no sampler output
    # from the implementation is used to define either the target or tolerance.
    mean_se = np.sqrt(np.diag(covariance)/n)
    covariance_se = np.sqrt((covariance**2+np.outer(np.diag(covariance), np.diag(covariance)))/(n-1))
    assert np.all(np.abs(updated.mean(axis=0)) < 6*mean_se)
    assert np.all(np.abs(np.cov(updated, rowvar=False)-covariance) < 6*covariance_se)


def test_default_legacy_stream_consumption_and_optional_argument_are_unchanged():
    channels, explicit_channels = _channels(), _channels()
    default = PathIntegralHarmonicOscillator(beads=4).run_simulation(channels, sweeps=128, burn_in=32)
    explicit = PathIntegralHarmonicOscillator(beads=4).run_simulation(
        explicit_channels, sweeps=128, burn_in=32, centroid_source=None)
    assert asdict(default) == asdict(explicit)
    # Frozen from the unchanged dd96adba implementation before adding this option.
    assert default.virial_energy_mean == 0.828529
    assert default.thermodynamic_energy_mean == 0.990092
    assert default.provenance_hash == "7bf770bcaf05922c8ce3055fce4d43f0d083ce1cf25971531e37e94d69b98c74"
    assert [s.consumed_bytes for s in (channels.proposal, channels.acceptance, channels.initialization)] == [4096, 2384, 32]
    assert default.sampler_protocol == "LOCAL_METROPOLIS_V1"
    assert default.centroid_update == {"enabled": False}


def test_optional_simulation_schedule_replay_and_unverified_independence_metadata():
    first_model, second_model = PathIntegralHarmonicOscillator(beads=4), PathIntegralHarmonicOscillator(beads=4)
    a, b = PCG64Source(17100), PCG64Source(17100)
    first = first_model.run_simulation(_channels(), sweeps=64, burn_in=16, centroid_source=a)
    second = second_model.run_simulation(_channels(), sweeps=64, burn_in=16, centroid_source=b)
    assert np.array_equal(first_model.path, second_model.path)
    assert first.virial_energy_mean == second.virial_energy_mean
    assert first.virial_uncertainty == second.virial_uncertainty
    assert first.sampler_protocol == "LOCAL_METROPOLIS_CENTROID_HEATBATH_V1"
    meta = first.centroid_update
    assert meta["total_updates"] == 64  # includes the 16 burn-in sweeps
    assert meta["centroid_draw_bytes"] == meta["source_byte_counter_delta"] == 1024
    assert meta["source_digest_after"] == a.stream_digest == b.stream_digest
    assert meta["source_shared_with_local_channels"] is False
    assert meta["independence_verified"] is False
    assert first.virial_uncertainty["stationarity_verified"] is False
    shared = PCG64Source(1519)
    result = PathIntegralHarmonicOscillator(beads=4).run_simulation(
        shared, sweeps=32, burn_in=0, centroid_source=shared)
    assert result.centroid_update["source_shared_with_local_channels"] is True
    assert result.centroid_update["independence_verified"] is False
    assert result.centroid_update["source_byte_counter_delta"] == shared.consumed_bytes
    assert result.centroid_update["source_byte_counter_delta"] > result.centroid_update["centroid_draw_bytes"]


def test_invalid_optional_source_fails_before_any_local_consumption():
    source = PCG64Source(22)
    with pytest.raises(TypeError, match="RandomSource"):
        PathIntegralHarmonicOscillator(beads=4).run_simulation(source, centroid_source=object())
    assert source.consumed_bytes == 0
