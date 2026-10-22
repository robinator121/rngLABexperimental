# Mathematical mechanism
For s=a+b>0, d=delta_a+delta_b, and perturbed probabilities valid,
Delta pi1=(b delta_a-a delta_b)/(s(s+d)).
The linear term is (b delta_a-a delta_b)/s^2 and the exact remainder is -(b delta_a-a delta_b)d/[s^2(s+d)]. Quadratic order holds with s bounded away from zero and |d|<s. Scaling both rates leaves the stationary probability unchanged while changing mixing. At a=b=0 stationarity is nonunique; a=b=1 is periodic. Both lie outside the interior ergodic two-state simulation.

The nontrivial eigenvalue is lambda=1-a-b. A stationary binary observable has covariance pi(1-pi) lambda^k. The finite-N mean variance is pi(1-pi)/N [1+2 sum_(k=1)^(N-1)(1-k/N)lambda^k]. Negative lambda can reduce variance; slow mixing increases it. The simulations initialize in stationarity; nonequilibrium drift has no stationary variance certificate.

An 8-level driver with persistence rho requires augmented state (X,U). The driver marginal remains uniform. Exact augmented stationary probabilities of X=1 are:
[
  {
    "rho": 0,
    "pi_x1": 0.33333333333333337,
    "stationary_residual": 1.8041124150158794e-16,
    "power_difference": 2.220446049250313e-16,
    "naive_original_x_prediction": 0.3333333333333333
  },
  {
    "rho": 0.5,
    "pi_x1": 0.2857142857142857,
    "stationary_residual": 6.938893903907228e-18,
    "power_difference": 4.85722573273506e-17,
    "naive_original_x_prediction": 0.3333333333333333
  },
  {
    "rho": 0.9,
    "pi_x1": 0.25641025641025705,
    "stationary_residual": 4.163336342344337e-17,
    "power_difference": 3.594347042223944e-15,
    "naive_original_x_prediction": 0.3333333333333333
  }
]
The original-X IID transition model predicts 1/3 in every case and fails under persistence. These are established perturbation/dependence mechanisms, not a novelty claim.

The exact small-Ising IID acceptance-bias kernel has uphill acceptance max(exp(-DeltaE/T)-delta,0); downhill moves remain unconditional. At delta >= exp(-8/T), the two fully magnetized states become absorbing. Consequently a unique stationary state law need not exist. The forecast theory table is limited to delta<=.025 at Tc; it must not be extrapolated to strong faults or correlated driving.

References: [Mitrophanov 2005](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/26A9854BCB8D103B3A63A9B272616EC0/S0021900200001066a.pdf/sensitivity-and-convergence-of-uniformly-ergodic-markov-chains.pdf), [Murray & Elliott](https://arxiv.org/abs/1204.3187).


Fresh finite-N empirical/theoretical mean-variance ratios were1.053,1.112,.995 for the three stationary controls. The approximate bootstrap variance intervals all include their exact covariance references; see FINAL_TRACE_CHECKS.json. This is limited numerical agreement, not a theorem about arbitrary degraded sources.
