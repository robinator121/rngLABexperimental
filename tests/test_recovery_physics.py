"""Independent analytical/negative controls for the v0.14.1 physics repair."""
import math

import numpy as np
import pytest

from rnglab12.uncertainty import summarize_series
from rnglab14.models.pimc_harmonic import (
    PathIntegralHarmonicOscillator, exact_finite_bead_energy, harmonic_energy_estimators,
    sample_exact_finite_bead_paths,
)
from rnglab14.models.vmc_harmonic import HarmonicOscillatorVMC
from rnglab14.monitoring.calibration import FrozenRankCalibration
from rnglab14.monitoring.dual_monitor import DualIntegrityMonitor, evaluate_alarm_trace
from rnglab14.sources.reference import PCG64Source
from rnglab14.sources.base import StochasticChannels


def _channels(seed):
    return StochasticChannels(PCG64Source(seed=seed), PCG64Source(seed=seed+1), PCG64Source(seed=seed+2))


@pytest.mark.parametrize("beta", [.2, 2., 8.])
@pytest.mark.parametrize("m", [2, 4, 16, 32])
def test_pimc_finite_bead_formula_matches_independent_gaussian_reference(beta, m):
    identity = np.eye(m)
    difference = np.roll(identity, -1, axis=0) - identity
    laplacian = difference.T @ difference
    precision = (beta/m)*identity + (m/beta)*laplacian
    covariance = np.linalg.inv(precision)
    # Two independent quadratic forms and their exact Gaussian expectations.
    td_quadratic = identity/(2*m) - m*laplacian/(2*beta**2)
    cv_quadratic = identity/m - np.ones((m, m))/(2*m*m)
    td_reference = m/(2*beta) + np.trace(td_quadratic @ covariance)
    cv_reference = 1/(2*beta) + np.trace(cv_quadratic @ covariance)
    reference = exact_finite_bead_energy(m, beta)
    assert td_reference == pytest.approx(reference, rel=2e-10, abs=2e-10)
    assert cv_reference == pytest.approx(reference, rel=2e-10, abs=2e-10)
    # Basis-vector polarization reconstructs the implementation's quadratic form.
    td, cv = harmonic_energy_estimators(identity, beta)
    assert td - m/(2*beta) == pytest.approx(np.diag(td_quadratic))
    assert cv - 1/(2*beta) == pytest.approx(np.diag(cv_quadratic))
    mixed = identity[0] + identity[1]
    td_mix, cv_mix = harmonic_energy_estimators(mixed, beta)
    assert td_mix - m/(2*beta) == pytest.approx(mixed @ td_quadratic @ mixed)
    assert cv_mix - 1/(2*beta) == pytest.approx(mixed @ cv_quadratic @ mixed)


def test_original_thermodynamic_coefficient_fails_deterministic_path_derivative():
    x = np.array([0., 1., 2., 3.])
    beta = 2.; m = len(x)
    springs = np.sum((np.roll(x, -1)-x)**2)
    old = m/(2*beta) - springs/(2*beta) + np.mean(.5*x*x)
    def action(b):
        return m*springs/(2*b) + b*np.mean(.5*x*x)
    h = 1e-5
    independent = m/(2*beta) + (action(beta+h)-action(beta-h))/(2*h)
    corrected, _ = harmonic_energy_estimators(x, beta)
    assert corrected == pytest.approx(independent, abs=1e-8)
    assert abs(old-independent) > 2.9


def test_discretization_converges_to_continuum_reference():
    beta = 2.
    continuum = .5/math.tanh(beta/2)
    errors = [abs(exact_finite_bead_energy(m, beta)-continuum) for m in [4, 8, 16, 32]]
    assert all(a > b for a, b in zip(errors, errors[1:]))
    assert errors[-2]/errors[-1] == pytest.approx(4., rel=.02)


def test_pimc_exposes_serial_uncertainty_and_finite_bead_target():
    result = PathIntegralHarmonicOscillator(beads=4).run_simulation(_channels(9001), sweeps=700, burn_in=100)
    assert result.virial_uncertainty["reference"] == result.exact_finite_bead_energy
    assert result.thermodynamic_uncertainty["reference"] == result.exact_finite_bead_energy
    assert result.virial_uncertainty["stationarity_verified"] is False
    assert result.discretization_bias < 0
    assert result.virial_uncertainty["n"] == 600


@pytest.mark.parametrize("kwargs", [{"beads": 1}, {"beads": 2.5}, {"temperature": 0},
    {"temperature": float("nan")}, {"displacement_delta": -1}])
def test_pimc_invalid_or_unsupported_inputs_fail(kwargs):
    with pytest.raises(ValueError):
        PathIntegralHarmonicOscillator(**kwargs)


def test_no_uncertainty_certificate_from_constant_energy():
    result = HarmonicOscillatorVMC().run_simulation(_channels(111), steps=256, burn_in=64)
    assert result.energy_mean == .5
    assert result.energy_uncertainty["status"] == "ABSTAIN_CONSTANT_TRACE"
    assert result.energy_uncertainty["geyer"]["ess"] is None
    assert result.observable_masking_detected is False
    assert result.masking_p_value is None
    assert result.monitoring_status == "UNCALIBRATED_EXPLORATORY_ONLY"


def test_vmc_observable_masking_identity_and_exact_moments():
    x = np.array([-1000., -1., 0., 2., 1000.])
    alpha = .5
    local_energy = alpha + x*x*(.5-2*alpha**2)
    assert np.all(local_energy == .5)
    # Independent Gaussian integral: density ∝ exp(-2 alpha x²).
    assert HarmonicOscillatorVMC(alpha=alpha).exact_x2 == 1/(4*alpha)


def test_rank_calibration_frozen_protocol_ties_and_finite_sample_rank():
    protocol = {"score": "fixed", "steps": 100}
    calibration = FrozenRankCalibration.fit(range(19), protocol, alpha=.05)
    assert calibration.p_value(19, protocol) == .05
    assert calibration.alarm(19, protocol)
    assert not calibration.alarm(18, protocol)
    with pytest.raises(ValueError, match="protocol"):
        calibration.p_value(19, {"score": "changed"})
    tied = FrozenRankCalibration.fit([1.]*39, protocol)
    assert tied.p_value(1.) == 1.
    assert not tied.alarm(1.)
    # Under all 20 equally likely ranks, exactly 1/20 is rejected.
    decisions = [FrozenRankCalibration.fit([s for s in range(20) if s != i], protocol).alarm(i)
                 for i in range(20)]
    assert sum(decisions) == 1


def test_vmc_calibration_cannot_be_reused_for_different_sampler_settings():
    model = HarmonicOscillatorVMC()
    calibration = FrozenRankCalibration.fit([.1]*39, model.monitoring_protocol(256, 64))
    with pytest.raises(ValueError, match="protocol"):
        model.run_simulation(_channels(111), steps=300, burn_in=64, masking_calibration=calibration)
    res = model.run_simulation(_channels(111), steps=256, burn_in=64, masking_calibration=calibration)
    assert res.masking_calibration_hash == calibration.digest
    assert res.monitoring_status == "CALIBRATED_EXCHANGEABLE_CHAIN_NULL"


def test_serial_uncertainty_detects_dependence_but_does_not_certify_mixing():
    rng = np.random.default_rng(8711)
    innovations = rng.normal(size=8000)
    values = np.zeros(len(innovations))
    for i in range(1, len(values)):
        values[i] = .8*values[i-1] + innovations[i]
    result = summarize_series(values[1000:])
    assert 5 < result["geyer"]["tau"] < 15
    assert result["geyer"]["mcse"] > 1.8*result["naive_mcse"]
    assert result["stationarity_verified"] is False


def test_monitor_does_not_invent_lead_time_for_censored_application():
    monitor = DualIntegrityMonitor()
    monitor.record_step(100, .2, .01)
    monitor.record_step(1000, .2, .02)
    result = monitor.summarize()
    assert result.mean_early_warning_lead_time_steps is None
    assert result.mean_detection_lag_steps is None
    assert result.application_failure_observed is False
    assert result.combined_monitoring_superior is None
    assert result.comparison_status.startswith("DESCRIPTIVE")


def test_monitor_nonuniform_times_and_source_missing_are_not_complete_events():
    monitor = DualIntegrityMonitor()
    monitor.record_step(10, .01, .2)
    monitor.record_step(5000, .01, .2)
    assert monitor.summarize().mean_detection_lag_steps is None
    monitor.record_step(5010, .2, .2)
    assert monitor.summarize().mean_detection_lag_steps == 5000.
    with pytest.raises(ValueError):
        monitor.record_step(5010, .1, .2)


def test_monitor_requires_finite_data_and_never_claims_superiority_from_quadrants():
    monitor = DualIntegrityMonitor()
    with pytest.raises(ValueError):
        monitor.record_step(1, float("nan"), .2)
    monitor.record_step(1, .2, .01)
    monitor.record_step(2, .2, .2)
    assert monitor.summarize().mean_early_warning_lead_time_steps == 1
    assert monitor.summarize().combined_monitoring_superior is None


def test_exact_reference_sampler_obeys_gaussian_covariance_and_separate_protocol():
    # Independent normal-mode precision, not a comparison of two samplers.
    beads=4; beta=2.
    source=PCG64Source(seed=189711)
    paths=sample_exact_finite_bead_paths(source,beads,beta,12000)
    indices=np.arange(beads)
    eigenvalues=beta/beads+4*beads/beta*np.sin(math.pi*indices/beads)**2
    covariance=np.array([[np.sum(np.cos(2*math.pi*indices*(j-k)/beads)/eigenvalues)/beads
                          for k in indices] for j in indices])
    assert np.mean(paths,axis=0)==pytest.approx(np.zeros(beads),abs=.025)
    assert np.cov(paths,rowvar=False)==pytest.approx(covariance,abs=.025)
    result=PathIntegralHarmonicOscillator(beads).run_exact_reference(PCG64Source(seed=189712),512)
    assert result["protocol"]=="EXACT_FINITE_BEAD_GAUSSIAN_REFERENCE_V1"
    assert result["independence_verified"] is False
    assert result["thermodynamic"]["status"]=="IID_REFERENCE_ASYMPTOTIC_T"
    assert result["channels_used"]==["supplied source only"]


def test_calibrated_monitor_uses_bonferroni_and_enforces_frozen_horizon():
    calibration=FrozenRankCalibration.fit(range(39),{"horizon":2,"score":"whole trajectory max"})
    monitor=DualIntegrityMonitor(source_calibration=calibration,app_calibration=calibration,horizon=2)
    # p=1/40=.025 triggers each channel's alpha/2 threshold.
    s=monitor.record_step(0,39.,0.)
    assert s.source_alarm_active
    assert not s.app_failure_active
    monitor.record_step(1,0.,39.)
    assert monitor.summarize().calibration_status.startswith("FROZEN")
    with pytest.raises(ValueError,match="horizon"):
        monitor.record_step(2,0.,0.)
    with pytest.raises(ValueError,match="horizon"):
        DualIntegrityMonitor(source_calibration=calibration,app_calibration=calibration,horizon=3)


def test_pre_fault_false_alarm_is_not_a_true_detection():
    result=evaluate_alarm_trace([True,False,False,False,False],fault_onset=2,failure_onset=3)
    assert result["first_alarm"]==0
    assert result["first_detection_from_fault_onset"] is None
    assert result["pre_fault_or_healthy_false_alarm"]
    assert result["missed_failure"] is True
    assert result["active_failure_step_hazard_fraction"]==1.


def test_active_gap_is_not_cumulative_first_alarm_gap():
    result=evaluate_alarm_trace([False,True,False,True,False],fault_onset=1,failure_onset=3)
    assert result["detection_delay_conditional_on_detection"]==-2
    assert result["failure_steps_before_first_true_detection"]==0
    assert result["active_failure_step_hazard_fraction"]==.5
    assert result["pre_fault_or_healthy_false_alarm"] is False
