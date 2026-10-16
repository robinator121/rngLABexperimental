"""Unit tests for RNGLab v0.13+ Scientific Models and Remaining Useful Samples Forecasting."""
import math
import unittest

from rnglab13.forecast import (
    CUSUMChangePointPredictor,
    ForecastingEvaluationReport,
    LinearExtrapolationPredictor,
    NoChangeBaselinePredictor,
    TrajectoryFamily,
    TrajectoryGenerator,
    evaluate_forecast_on_trajectory,
    summarize_forecasting_performance,
)
from rnglab13.models import (
    HarmonicOscillatorVMC,
    Ising2D,
    IsingObservableResult,
    PiMCResult,
    VMCResult,
    run_pi_monte_carlo,
)
from rnglab13.sources import MersenneTwisterSource, SyntheticTRNGSource


class ModelsAndForecastingV13Tests(unittest.TestCase):
    """Test scientific models and RUS forecasting subsystem."""

    def test_pi_monte_carlo_healthy(self):
        src = MersenneTwisterSource(seed=42)
        res: PiMCResult = run_pi_monte_carlo(src, samples=5000)
        self.assertEqual(res.total_samples, 5000)
        self.assertLess(abs(res.z_score), 3.0)  # Within 3 standard deviations
        self.assertAlmostEqual(res.true_pi, math.pi, places=5)

    def test_ising_2d_simulation(self):
        src = MersenneTwisterSource(seed=101)
        ising = Ising2D(L=8, temperature=2.269, seed=101)
        res: IsingObservableResult = ising.run_simulation(src, sweeps=200, burn_in=50)
        self.assertEqual(res.L, 8)
        self.assertEqual(res.total_sweeps, 200)
        self.assertTrue(res.healthy)
        self.assertGreater(res.specific_heat, 0.0)
        self.assertGreater(res.susceptibility, 0.0)
        self.assertGreater(res.effective_sample_size_energy, 10.0)

    def test_vmc_harmonic_oscillator_zero_variance_and_masking(self):
        # 1. Healthy source: alpha=0.5 yields exact energy 0.5 with near-zero variance
        src_healthy = MersenneTwisterSource(seed=777)
        vmc = HarmonicOscillatorVMC(alpha=0.5, delta=1.0)
        res_healthy: VMCResult = vmc.run_simulation(src_healthy, steps=2000, burn_in=500)
        self.assertAlmostEqual(res_healthy.energy_mean, 0.5, places=4)
        self.assertAlmostEqual(res_healthy.x2_mean, 0.5, places=1)
        self.assertFalse(res_healthy.observable_masking_detected)

        # 2. Defective source with thermal drift (synthetic TRNG):
        src_drift = SyntheticTRNGSource(seed=999, thermal_drift_rate=50.0)
        res_drift: VMCResult = vmc.run_simulation(src_drift, steps=2000, burn_in=500)
        # Energy remains ~ 0.5 because E_L(x) = 0.5 identically for alpha=0.5!
        self.assertAlmostEqual(res_drift.energy_mean, 0.5, places=3)
        # But position variance x^2 is shifted!
        self.assertNotEqual(round(res_drift.x2_mean, 2), 0.50)

    def test_all_trajectory_families_generate_paths(self):
        families = [
            TrajectoryFamily.CONSTANT_BIAS,
            TrajectoryFamily.LINEAR_DRIFT,
            TrajectoryFamily.ACCELERATING_DRIFT,
            TrajectoryFamily.GRADUAL_AUTOCORRELATION,
            TrajectoryFamily.STEPWISE_FAILURE,
            TrajectoryFamily.INTERMITTENT_FAULT,
            TrajectoryFamily.FAULT_RECOVERY,
            TrajectoryFamily.STOCHASTIC_DRIFT,
        ]
        for fam in families:
            gen = TrajectoryGenerator(family=fam, seed=42)
            path = gen.generate_path(total_steps=100)
            self.assertEqual(len(path), 100)
            self.assertEqual(path[0].t, 1)
            self.assertEqual(path[-1].t, 100)

    def test_forecasting_linear_drift_evaluation(self):
        # Generate linear drift trajectory: starts at 0, drifts by 0.001 per step, threshold at 0.50 -> crossing ~ step 500
        gen = TrajectoryGenerator(TrajectoryFamily.LINEAR_DRIFT, seed=42, start_val=0.0, drift_rate=0.001)
        path = gen.generate_path(total_steps=600)
        values = [p.parameter for p in path]

        pred = LinearExtrapolationPredictor(threshold=0.50, window=100)
        eval_times = [100, 200, 300, 400]
        forecasts = evaluate_forecast_on_trajectory(pred, values, threshold=0.50, evaluation_times=eval_times)
        self.assertEqual(len(forecasts), 4)

        report: ForecastingEvaluationReport = summarize_forecasting_performance(forecasts)
        self.assertEqual(report.total_evaluations, 4)
        # At step 200: value=0.20, threshold=0.50, remaining = 300. Linear extrapolation should be within 20 samples!
        f200 = next(f for f in forecasts if f.current_t == 200)
        self.assertLess(abs(f200.forecast_error), 20)
        # Lower bound should cover actual remaining
        self.assertTrue(f200.covered_by_lower_bound)
        self.assertFalse(f200.is_false_reassurance)

    def test_cusum_change_point_predictor(self):
        gen = TrajectoryGenerator(TrajectoryFamily.STEPWISE_FAILURE, seed=42, step_time=200, baseline=0.0, failure_level=0.45)
        path = gen.generate_path(total_steps=400)
        values = [p.parameter for p in path]

        pred = CUSUMChangePointPredictor(threshold=0.40, baseline_mean=0.0, k=0.01, h=0.05)
        point_h, lower_h = pred.predict(values[:250], current_t=250)
        # Value at 250 is 0.45 >= threshold 0.40 -> immediate crossing
        self.assertEqual(point_h, 0)
        self.assertEqual(lower_h, 0)


if __name__ == "__main__":
    unittest.main()
