"""Checks of independent quadrature against analytic and tabulated references."""
import math
import unittest
import verify_baseline as checker

class IndependentCheckerChecks(unittest.TestCase):
    def test_quadrature_matches_analytic_cauchy_integral(self):
        integral=checker.integrate(lambda x:1/(math.pi*(1+x*x)),0,1)
        self.assertAlmostEqual(integral,.25,places=11)

    def test_student_t_critical_matches_standard_df5_value(self):
        self.assertAlmostEqual(checker.tcrit(5),2.570581835636305,places=9)

    def test_wilson_boundary_has_exact_zero_lower_endpoint(self):
        lo,hi=checker.wilson(0,64)
        self.assertAlmostEqual(lo,0,places=14)
        self.assertGreater(hi,0)
        self.assertLess(hi,.1)

if __name__=="__main__": unittest.main()
