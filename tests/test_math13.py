"""Unit tests for RNGLab v0.13+ Mathematical Research Work Packages (MATH-01 to MATH-06)."""
import math
import unittest
import numpy as np

from rnglab13.math import (
    MarkovEquivalenceExperiment,
    binary_entropy,
    compute_lipschitz_drift_bound,
    effective_sample_size,
    exact_variance_closed_form,
    exact_variance_summation,
    fundamental_equation_of_information_defect,
    markov_entropy_rate,
    markov_fundamental_matrix,
    predict_recovery_horizon,
    run_entropy_equivalence_experiment,
    run_indistinguishable_histories_counterexample,
    shannon_information_function,
    simulate_recovery_trajectory,
    test_hyers_ulam_stability_of_fei as check_hyers_ulam_fei,
    total_variation_distance,
    xor_composition_cauchy_identity,
)


class MathWorkPackagesV13Tests(unittest.TestCase):
    """Test mathematical identities, bounds, and counterexamples."""

    # --- MATH-01 Tests ---
    def test_finite_sample_variance_summation_matches_closed_form(self):
        test_cases = [
            (10, 0.2),
            (50, 0.5),
            (100, -0.3),
            (500, -0.1),
            (1000, 0.05),
        ]
        for n, rho in test_cases:
            var_sum = exact_variance_summation(n, rho)
            var_closed = exact_variance_closed_form(n, rho)
            self.assertAlmostEqual(
                var_sum,
                var_closed,
                places=9,
                msg=f"Summation and closed-form must match exactly for n={n}, rho={rho}",
            )

    def test_entropy_rate_is_strictly_symmetric_in_rho(self):
        """H_rate(+rho) == H_rate(-rho), proving entropy rate cannot detect correlation sign."""
        for rho in (0.05, 0.1, 0.3, 0.5, 0.7):
            h_pos = markov_entropy_rate(rho)
            h_neg = markov_entropy_rate(-rho)
            self.assertAlmostEqual(h_pos, h_neg, places=12)

    def test_variance_and_ess_diverge_radically_under_identical_entropy_rate(self):
        """Demonstrate that equal entropy rates yield radically different application variances."""
        n = 1000
        rho = 0.5
        h_pos = markov_entropy_rate(rho)
        h_neg = markov_entropy_rate(-rho)
        self.assertAlmostEqual(h_pos, h_neg, places=10)

        var_pos = exact_variance_closed_form(n, rho)
        var_neg = exact_variance_closed_form(n, -rho)
        # For rho=0.5, asymptotic variance ratio is ((1+0.5)/(1-0.5))^2 = 3^2 = 9
        ratio = var_pos / var_neg
        self.assertGreater(ratio, 8.5)
        self.assertLess(ratio, 9.5)

        ess_pos = effective_sample_size(n, rho)
        ess_neg = effective_sample_size(n, -rho)
        self.assertAlmostEqual(ess_pos, n / 3.0, delta=1.0)
        self.assertAlmostEqual(ess_neg, n * 3.0, delta=1.0)

    def test_monte_carlo_markov_experiment_matches_theory(self):
        res: MarkovEquivalenceExperiment = run_entropy_equivalence_experiment(
            n=200,
            rho=0.2,
            trials=500,
            base_seed=123,
        )
        self.assertLess(abs(res.mean_bias), 0.03)
        self.assertLess(res.relative_variance_error, 0.15)
        self.assertGreaterEqual(res.coverage_95, 0.90)

    # --- MATH-02 Tests ---
    def test_markov_fundamental_matrix_and_tv_bounds(self):
        P0 = np.array([[0.8, 0.2], [0.3, 0.7]])
        pi0, Z, gamma = markov_fundamental_matrix(P0)
        # Stationary check: pi * P = pi
        np.testing.assert_allclose(np.dot(pi0, P0), pi0, atol=1e-10)
        self.assertGreater(gamma, 0.0)

        # Perturbation
        Delta_P = np.array([[-0.1, 0.1], [0.1, -0.1]])
        P_perturbed = P0 + 0.1 * Delta_P
        pi_perturbed, _, _ = markov_fundamental_matrix(P_perturbed)
        tv = total_variation_distance(pi_perturbed, pi0)
        self.assertGreater(tv, 0.0)
        self.assertLess(tv, 0.1)

    # --- MATH-04 Tests ---
    def test_impossibility_counterexample(self):
        res = run_indistinguishable_histories_counterexample(t_history=500, t_future=200)
        self.assertEqual(res.mean_history_A, res.mean_history_B)
        self.assertIsNone(res.actual_crossing_A)
        self.assertIsNotNone(res.actual_crossing_B)
        self.assertIn("Proves Theorem 1", res.mathematical_conclusion)

    def test_lipschitz_drift_bound_calculation(self):
        # Current est 0.52, error margin 0.02, threshold 0.70, max drift 0.001/step
        # Safe distance = 0.70 - 0.54 = 0.16 -> 0.16 / 0.001 = 160 steps
        n_safe = compute_lipschitz_drift_bound(
            current_estimate=0.52,
            estimation_error_margin=0.02,
            drift_rate_bound_L=0.001,
            threshold=0.70,
        )
        self.assertEqual(n_safe, 160)

    # --- MATH-05 Tests ---
    def test_fei_exact_solution_has_zero_defect(self):
        defect = fundamental_equation_of_information_defect(shannon_information_function, 0.2, 0.3)
        self.assertAlmostEqual(defect, 0.0, places=12)

    def test_hyers_ulam_stability_of_fei(self):
        reports = check_hyers_ulam_fei(epsilon_values=[0.0, 0.01, 0.05])
        self.assertEqual(len(reports), 3)
        for r in reports:
            self.assertTrue(r.is_hyers_ulam_stable)

    def test_xor_cauchy_identity(self):
        psi_comb, psi_sum, defect = xor_composition_cauchy_identity(0.1, 0.05)
        self.assertAlmostEqual(defect, 0.0, places=12)

    # --- MATH-06 Tests ---
    def test_recovery_trajectory_memory_decay(self):
        traj = simulate_recovery_trajectory(t_healthy_pre=100, t_defect=100, t_healthy_post=300, defect_bias=0.25)
        self.assertEqual(len(traj), 500)

        # Before defect: healthy
        self.assertTrue(traj[50].is_source_healthy)
        self.assertEqual(traj[50].source_bias, 0.0)

        # During defect: bias active
        self.assertFalse(traj[150].is_source_healthy)
        self.assertEqual(traj[150].source_bias, 0.25)

        # Post defect: source restored healthy, but cumulative estimator retains memory
        self.assertTrue(traj[250].is_source_healthy)
        self.assertEqual(traj[250].source_bias, 0.0)
        self.assertGreater(traj[250].cumulative_estimator_bias, 0.05)
        # Reset estimator should have near zero bias
        self.assertIsNotNone(traj[250].reset_estimator_bias)
        self.assertLess(abs(traj[250].reset_estimator_bias), 0.10)

    def test_predict_recovery_horizon(self):
        # 100 defective samples with bias 0.20 -> total defective excess = 20
        # To get bias <= 0.01: need total samples >= 20 / 0.01 = 2000 -> 1900 healthy samples
        h = predict_recovery_horizon(tau_defect=100, defect_bias=0.20, acceptable_tolerance=0.01)
        self.assertEqual(h, 1900)


if __name__ == "__main__":
    unittest.main()
