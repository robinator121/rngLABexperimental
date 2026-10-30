"""Deterministic checks of the synthetic RUS construction (stdlib unittest)."""
import math
import unittest
import rus_identifiability as model

class CounterexampleChecks(unittest.TestCase):
    def test_binomial_probability_at_one_half_matches_integer_count(self):
        exact_numerator=sum(math.comb(64,k) for k in range(65) if k<=16 or k>=48)
        self.assertAlmostEqual(model.failure_probability(.5), exact_numerator/(2**64), places=18)

    def test_complementary_bernoulli_laws_have_same_two_sided_failure(self):
        self.assertAlmostEqual(model.failure_probability(.1), model.failure_probability(.9), places=14)

    def test_failures_at_the_known_boundaries(self):
        self.assertEqual(model.first_failure([1]*2048),64)
        self.assertEqual(model.first_failure([0]*2048),64)
        self.assertIsNone(model.first_failure([0,1]*1024))
        # A count of47 is inside the acceptable interval;48 fails.
        self.assertIsNone(model.first_failure(([1]*47+[0]*17)*32))
        self.assertEqual(model.first_failure(([1]*48+[0]*16)*32),64)

    def test_finite_horizon_probabilities_separate_the_two_laws(self):
        q0=model.failure_probability(.5); q1=model.failure_probability(.9)
        self.assertLess(1-(1-q0)**32,.003)
        self.assertGreater(q1,.999)

    def test_seed_rule_is_unique_for_the_executed_replicas(self):
        seeds=[model.seed(f"{model.FUTURE_NAMESPACE}/{r:04d}")[0] for r in range(model.REPLICAS)]
        self.assertEqual(len(set(seeds)),model.REPLICAS)

if __name__=="__main__": unittest.main()
