# Mechanism-to-consequence evidence

MECHANISM_EVIDENCE.json contains primary-cell detection rates with binomial intervals, paired detection differences and excess squared-loss bootstrap intervals. Raw observables remain separate in calibration_analysis.json.

- Marginal acceptance bias perturbs an IID transition kernel. Mild arms preserve positive transitions; structural arms deliberately eliminate uphill transitions and permit absorption.
- Proposal/acceptance coupling changes conditional draws even when marginals look plausible. In VMC it harms x²/x⁴ despite exact energy. No ordinary-X stationary perturbation formula is assumed for arbitrary correlated drivers.
- Temporal correlation uses an independent refresh selector and persistent previous uniform; driver marginal and algorithm dependence must be distinguished. Periodic faults are time-inhomogeneous and not assigned a single equilibrium reference law.
- Source anomalies in an UNUSED upstream channel leave paired Ising trajectories exactly unchanged. The NPZ stores diagnostic versus actually consumed source means separately. This is a causal negative control, not entropy certification. Quantum unused-source arms are exact healthy replays, not measured upstream-source anomalies.
- Healthy MCMC can have substantial uncertainty undercoverage without any defective RNG.

Nine mild-bias effects met significance/effect-size gates, but only seven also passed their own healthy false-alarm validation: ising_L2_T1, ising_L2_T2, ising_L3_T1, ising_L3_T2, ising_L4_T2, vmc, pimc4. H2 controlled-discrimination claims are limited to those seven. L4/Tc and PIMC16 effect significance cannot compensate for failed calibration. Near-critical susceptibility is not a universal RNG-vulnerability claim; finite size, channel,temperature and numerical severity matter.
