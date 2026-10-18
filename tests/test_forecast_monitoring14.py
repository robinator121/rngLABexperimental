"""Unit tests for RNGLab v0.14 RUS Forecasting, Dual Monitoring, and Multi-Testing Statistics."""
import unittest

from rnglab14.forecast import (
    CUSUMChangePointPredictor,
    LinearExtrapolationPredictor,
    LipschitzStructuralBoundPredictor,
    NoChangeBaselinePredictor,
    TrajectoryFamily,
    TrajectoryGenerator,
    evaluate_forecast_on_trajectory,
    summarize_forecasting_performance,
)
from rnglab14.monitoring import (
    DualIntegrityMonitor,
    DualMonitoringReport,
    DualMonitorSnapshot,
    MonitoringQuadrant,
)
from rnglab14.statistics import (
    benjamini_hochberg_fdr,
    bonferroni_correction,
    holm_bonferroni_correction,
)


class ForecastAndMonitoringV14Tests(unittest.TestCase):
    def test_rus_forecasting_separates_censored_and_uncensored_events(self):
        # 1. Uncensored failing trajectory (linear drift crosses 0.40 around step 400)
        gen = TrajectoryGenerator(TrajectoryFamily.LINEAR_DRIFT, seed=42, start_val=0.0, drift_rate=0.001)
        path = gen.generate_path(total_steps=500)
        vals = [p.parameter for p in path]

        pred = LinearExtrapolationPredictor(threshold=0.40, window=50)
        fcasts = evaluate_forecast_on_trajectory(pred, vals, threshold=0.40, evaluation_times=[100, 200, 300])
        rep = summarize_forecasting_performance(fcasts)

        self.assertEqual(rep.uncensored_events_count, 3)
        self.assertEqual(rep.censored_events_count, 0)
        self.assertLess(rep.event_mae, 25.0)  # Accurate on linear drift
        self.assertFalse(rep.false_reassurance_count > 0)

        # 2. Censored non-failing trajectory (constant bias 0.10 never crosses 0.40)
        gen_censored = TrajectoryGenerator(TrajectoryFamily.CONSTANT_BIAS, seed=42, baseline=0.10)
        path_cens = gen_censored.generate_path(total_steps=500)
        vals_cens = [p.parameter for p in path_cens]

        fcasts_cens = evaluate_forecast_on_trajectory(pred, vals_cens, threshold=0.40, evaluation_times=[100, 200, 300])
        rep_cens = summarize_forecasting_performance(fcasts_cens)

        self.assertEqual(rep_cens.censored_events_count, 3)
        self.assertEqual(rep_cens.uncensored_events_count, 0)
        # No scored event predictions: error and calibration are undefined, not perfect.
        self.assertIsNone(rep_cens.event_mae)
        self.assertIsNone(rep_cens.is_well_calibrated)
        # Zero false reassurance on censored runs
        self.assertEqual(rep_cens.false_reassurance_count, 0)

    def test_uncertainty_abstention_on_flat_signals(self):
        # A flat zero signal has slope ~ 0. Model should abstain!
        vals = [0.0] * 100
        pred = LinearExtrapolationPredictor(threshold=0.40, window=50)
        point_h, lower_h, abstained = pred.predict(vals, current_t=100)
        self.assertTrue(abstained)
        self.assertIsNone(point_h)
        self.assertIsNone(lower_h)

    def test_lipschitz_bound_guarantees_zero_false_reassurance(self):
        gen = TrajectoryGenerator(TrajectoryFamily.LINEAR_DRIFT, seed=10, start_val=0.0, drift_rate=0.0008)
        path = gen.generate_path(total_steps=600)
        vals = [p.parameter for p in path]

        pred = LipschitzStructuralBoundPredictor(threshold=0.40, lipschitz_L=0.001)
        fcasts = evaluate_forecast_on_trajectory(pred, vals, threshold=0.40, evaluation_times=[100, 200, 300, 400])
        rep = summarize_forecasting_performance(fcasts)
        self.assertEqual(rep.false_reassurance_count, 0)
        self.assertEqual(rep.lower_bound_coverage_95, 1.0)

    def test_dual_integrity_monitor_quadrants(self):
        monitor = DualIntegrityMonitor(source_alarm_threshold=0.05, app_failure_threshold=0.10)

        # 1. Step 10: Normal
        s1 = monitor.record_step(10, source_defect_magnitude=0.01, app_observable_error=0.02)
        self.assertEqual(s1.quadrant, MonitoringQuadrant.Q1_NOMINAL)

        # 2. Step 20: Early warning (Source defect elevated, application not yet broken)
        s2 = monitor.record_step(20, source_defect_magnitude=0.08, app_observable_error=0.03)
        self.assertEqual(s2.quadrant, MonitoringQuadrant.Q2_EARLY_WARNING)

        # 3. Step 30: Synchronized failure
        s3 = monitor.record_step(30, source_defect_magnitude=0.12, app_observable_error=0.15)
        self.assertEqual(s3.quadrant, MonitoringQuadrant.Q3_SYNCHRONIZED_FAILURE)

        rep = monitor.summarize()
        self.assertEqual(rep.total_snapshots, 3)
        self.assertEqual(rep.q2_early_warning_count, 1)
        self.assertGreater(rep.mean_early_warning_lead_time_steps, 0.0)

    def test_multi_testing_error_rate_control(self):
        # 10 p-values, some small, some large
        p_vals = [0.001, 0.004, 0.012, 0.045, 0.060, 0.120, 0.250, 0.400, 0.700, 0.900]

        # Bonferroni: alpha / 10 = 0.005 -> only 0.001 and 0.004 rejected
        bonf = bonferroni_correction(p_vals, alpha=0.05)
        self.assertEqual(bonf.rejected_count, 2)

        # Holm-Bonferroni: step-down
        holm = holm_bonferroni_correction(p_vals, alpha=0.05)
        self.assertGreaterEqual(holm.rejected_count, 2)

        # Benjamini-Hochberg FDR at q=0.05
        fdr = benjamini_hochberg_fdr(p_vals, q_fdr=0.05)
        self.assertGreaterEqual(fdr.rejected_count, 2)


if __name__ == "__main__":
    unittest.main()
