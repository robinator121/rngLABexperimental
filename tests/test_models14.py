"""Unit tests for RNGLab v0.14 Computational Science Models."""
import math
import unittest

from rnglab14.models import (
    HarmonicOscillatorVMC,
    Ising2D,
    MonteCarloDecisionSupport,
    PathIntegralHarmonicOscillator,
    run_pi_monte_carlo,
)
from rnglab14.sources import (
    DegradedFaultSource,
    MersenneTwisterSource,
    PCG64Source,
    StochasticChannels,
    Xoshiro256PlusPlusSource,
)


class ModelsV14Tests(unittest.TestCase):
    def test_pi_monte_carlo_with_stochastic_channels(self):
        channels = StochasticChannels(
            proposal=PCG64Source(seed=101),
            acceptance=PCG64Source(seed=102),
            initialization=PCG64Source(seed=103),
        )
        res = run_pi_monte_carlo(channels, samples=5000)
        self.assertEqual(res.total_samples, 5000)
        self.assertAlmostEqual(res.estimated_pi, math.pi, places=1)
        self.assertLess(abs(res.z_score), 3.0)

    def test_ising_2d_with_decoupled_channels(self):
        channels = StochasticChannels(
            proposal=PCG64Source(seed=201),
            acceptance=Xoshiro256PlusPlusSource(seed=202),
            initialization=MersenneTwisterSource(seed=203),
        )
        ising = Ising2D(L=6, temperature=2.269)
        res = ising.run_simulation(channels, sweeps=200, burn_in=50)
        self.assertEqual(res.L, 6)
        self.assertTrue(res.healthy)
        self.assertGreater(res.specific_heat, 0.0)
        self.assertGreater(res.susceptibility, 0.0)
        self.assertGreater(res.effective_sample_size_energy, 5.0)

    def test_vmc_harmonic_oscillator_exact_zero_variance(self):
        # Ground state alpha = 0.5 has exact energy 0.5 with zero local variance
        channels = StochasticChannels(
            proposal=PCG64Source(seed=301),
            acceptance=PCG64Source(seed=302),
            initialization=PCG64Source(seed=303),
        )
        vmc = HarmonicOscillatorVMC(alpha=0.5, delta=1.0)
        res = vmc.run_simulation(channels, steps=2000, burn_in=500)
        self.assertAlmostEqual(res.energy_mean, 0.5, places=4)
        self.assertEqual(res.energy_std, 0.0)
        self.assertFalse(res.observable_masking_detected)

    def test_vmc_zero_variance_masking_contrast_with_imperfect_wavefunction(self):
        """Zero-variance property at alpha=0.5 causes an energetic blind spot to corrupted proposals,
        whereas an imperfect trial wavefunction (alpha=0.6) immediately detects the coordinate defect.
        """
        # Inject positive bias into proposal channel
        corrupted_prop_1 = DegradedFaultSource(PCG64Source(seed=42), bias_delta=0.25)
        channels_exact = StochasticChannels(
            proposal=corrupted_prop_1,
            acceptance=PCG64Source(seed=99),
            initialization=PCG64Source(seed=1),
        )
        vmc_exact = HarmonicOscillatorVMC(alpha=0.5, delta=1.0)
        res_exact = vmc_exact.run_simulation(channels_exact, steps=3000, burn_in=500)

        # At exact ground state (alpha=0.5), E_L(x)=0.5 identically, so variance is zero
        # despite severe coordinate sampling bias (<x^2> bias > 1.0)
        self.assertAlmostEqual(res_exact.energy_mean, 0.5, places=6)
        self.assertEqual(res_exact.energy_std, 0.0)
        self.assertGreater(abs(res_exact.x2_bias), 1.0)

        # With an imperfect trial wavefunction (alpha=0.6), E_L(x) depends on x^2
        corrupted_prop_2 = DegradedFaultSource(PCG64Source(seed=42), bias_delta=0.25)
        channels_imperfect = StochasticChannels(
            proposal=corrupted_prop_2,
            acceptance=PCG64Source(seed=99),
            initialization=PCG64Source(seed=1),
        )
        vmc_imperfect = HarmonicOscillatorVMC(alpha=0.6, delta=1.0)
        res_imperfect = vmc_imperfect.run_simulation(channels_imperfect, steps=3000, burn_in=500)

        # For alpha=0.6, energy variance is nonzero and energy bias shifts noticeably
        self.assertGreater(res_imperfect.energy_std, 0.2)
        self.assertGreater(abs(res_imperfect.energy_bias), 0.15)
        self.assertGreater(abs(res_imperfect.x2_bias), 0.8)

    def test_path_integral_monte_carlo_harmonic_oscillator(self):
        channels = StochasticChannels(
            proposal=PCG64Source(seed=401),
            acceptance=PCG64Source(seed=402),
            initialization=PCG64Source(seed=403),
        )
        # Beta = 2.0 (T = 0.5) -> E_exact = 0.5 * coth(1.0) ~ 0.656538
        pimc = PathIntegralHarmonicOscillator(beads=16, temperature=0.5, displacement_delta=0.4)
        res = pimc.run_simulation(channels, sweeps=1500, burn_in=400)
        self.assertEqual(res.beads_M, 16)
        self.assertAlmostEqual(res.exact_quantum_energy, 0.656538, places=4)
        # Virial energy estimate should be close to exact quantum energy
        self.assertAlmostEqual(res.virial_energy_mean, res.exact_quantum_energy, places=1)
        self.assertGreater(res.radius_of_gyration, 0.0)

    def test_monte_carlo_decision_support_system(self):
        src = PCG64Source(seed=501)
        dss = MonteCarloDecisionSupport()
        res = dss.evaluate_decisions(src, scenarios=3000)
        self.assertEqual(res.total_scenarios, 3000)
        # Under normal distributions, Growth_B (mean=120) should be optimal
        self.assertEqual(res.optimal_option, "Growth_B")
        self.assertTrue(res.ranking_correct)
        self.assertEqual(res.decision_regret, 0.0)
        self.assertGreater(res.confidence_spread, 5.0)


if __name__ == "__main__":
    unittest.main()
