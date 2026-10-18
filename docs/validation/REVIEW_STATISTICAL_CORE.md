> **SUPERSEDED v014 report — known incorrect/unverified claims.** Preserved for traceability; do not use as current validation. See [v014.1 corrections](V0141_CORRECTIONS.md).

# RNGLab v0.14 Internal Methodological Review & Mathematical Verification
**Scope:** Verification of mathematical foundations (MATH-01 to MATH-06), investigation of the MATH-03 performance contradiction, experimental validity, and multi-testing error controls.  
**Author:** Hevele Imre Róbert  
**Date:** October 2026  
**Status:** Review Completed — Explicit Corrections and Architectural Remedies Mandated.

---

## 1. High-Priority Investigation: The MATH-03 Performance Contradiction

### 1.1 The Contradiction Stated
In the v0.13 research reports and live `status.json`, performance metrics for Remaining Useful Samples (RUS) predictors were reported as follows:
- **`NoChangeBaselinePredictor`:** MAE = `546.9 steps`
- **`CUSUMChangePointPredictor`:** MAE = `6475.5 steps`
- **`LinearExtrapolationPredictor`:** MAE = `10165.9 steps`

Simultaneously, the narrative text asserted:
> *"CUSUMChangePointPredictor demonstrates superior performance and adaptability over baseline heuristics."*

A critical reviewer or external scientist will immediately ask: **How can a model with a 12-fold higher error (6,475 vs 547 steps) be claimed as "superior"?**

### 1.2 Mathematical & Code Root-Cause Analysis
Inspection of `src/rnglab13/forecast/predictor.py` revealed four compounding flaws in the experimental design and metric calculation:

#### Flaw 1: Arbitrary Right-Censoring Penalty
In `evaluate_forecast_on_trajectory` (lines 201–206):
```python
if actual_crossing is not None:
    actual_rem = max(0, actual_crossing - t)
else:
    actual_rem = total_len - t + 5000  # Threshold not reached
```
For trajectories where the failure threshold was never crossed (e.g. `CONSTANT_BIAS` with bias 0.10 < threshold 0.40, or `FAULT_RECOVERY` where the source returned to normal before failing), the true remaining lifetime is mathematically infinite or right-censored ($> total\_len - t$).  
The evaluation harness arbitrarily assigned `actual_rem = 5000 + (total_len - t) \approx 5150 - 5400$ steps.

#### Flaw 2: The Baseline Hardcoding Coincidence
In `NoChangeBaselinePredictor` (lines 58–59):
```python
# Naive constant horizon assumption
return 5000, 2000
```
The baseline predictor literally returned the fixed constant `5000`.  
Consequently, on all right-censored trajectories:
$$\text{Error}_{baseline} = 5000 - (5000 + 500 - t) = t - 500 \in [-400, -150]$$
The baseline achieved an artificially tiny MAE ($\sim 275 - 546$ steps) simply because its hardcoded default matched the evaluation harness's arbitrary censoring constant!

#### Flaw 3: Unbounded Extrapolation on Flat Signals
In `LinearExtrapolationPredictor` and `CUSUMChangePointPredictor`:
When a signal is flat or stationary, the fitted slope $\hat{\beta}_1$ is small but positive due to residual noise (e.g., $\hat{\beta}_1 \approx 10^{-5}$).  
The predictor computed:
$$\hat{h} = \frac{\theta_{crit} - \theta_t}{\hat{\beta}_1} = \frac{0.30}{10^{-5}} = 30,000 \text{ steps}$$
Evaluated against the synthetic target of $5,400$, this produced an individual point error of $|30,000 - 5,400| = 24,600$ steps! A small number of flat-trajectory evaluations completely dominated the unweighted arithmetic mean, blowing up the MAE to over 6,000–10,000 steps.

#### Flaw 4: Conflating Operational Safety with Metric Formulation
Why was CUSUM genuinely better in practice?
On actual failure trajectories (`STEPWISE_FAILURE`, `LINEAR_DRIFT`), CUSUM detected the onset within 2 samples, dropping remaining lifetime to 0 and achieving **8.6% false reassurance** and **91.4% early warning**, whereas the baseline maintained false reassurance on **27.9%** of steps (blindly telling the operator that 5,000 samples remained right until failure occurred).  
However, pooling censored and uncensored trajectories into an unweighted MAE produced an absurd, misleading headline statistic.

### 1.3 Mandated Remedies for RNGLab v0.14
1. **Separation of Censored and Event Trajectories:** Report event-time MAE strictly on trajectories that reach the critical threshold within the observation horizon.
2. **Censored Survival Likelihood:** For right-censored trajectories, evaluate using survival likelihood or Brier scores rather than arbitrary integer offsets.
3. **Explicit Uncertainty Abstention (`UNKNOWN`):** When the estimated slope is not statistically distinguishable from zero ($p > 0.05$ or confidence interval covers zero), the predictor must output `(None, None)` with flag `UNKNOWN / INSUFFICIENT EVIDENCE`, rather than dividing by microscopic noise slopes.
4. **Horizon Clamping:** Clamp maximum projection to a physically defensible horizon $H_{max} = 2 \times T_{obs}$.

---

## 2. Rigorous Audit of MATH Work Packages (MATH-01 to MATH-06)

### 2.1 MATH-01: Entropy Equivalence vs Variance Divergence
* **Classification:**
  * **Mathematical Identities (Analytical Proofs):**
    1. Entropy rate symmetry: $H_{rate}(\rho) = h_2((1+\rho)/2) = H_{rate}(-\rho)$.
    2. Finite-sample variance formula: $\text{Var}(\bar{X}_n) = \frac{1}{4n} [ \frac{1+\rho}{1-\rho} - \frac{2\rho(1-\rho^n)}{n(1-\rho)^2} ]$.
  * **Empirical Discoveries:**
    1. Propagation of the $9.00\times$ variance ratio ($\rho = +0.5$ vs $-0.5$) into Monte Carlo $\pi$ integration and Ising lattice magnetization.
* **Limitations & Scope:**
  * Proof applies strictly to stationary two-state Markov chains with symmetric transitions. Higher-order chains ($k \ge 2$) or non-symmetric transitions ($\pi_0 \ne \pi_1$) break strict entropy rate symmetry with respect to correlation sign.
* **Action:** Retain proofs in `MATHEMATICAL_NOTEBOOK.md`, but explicitly state the symmetry restriction to stationary symmetric binary chains.

### 2.2 MATH-02: Operator Perturbation & Observable Masking in QMC
* **Classification:**
  * **Mathematical Identity:** For the 1D Quantum Harmonic Oscillator, the local energy estimator is $E_L(x) = \alpha + x^2(\frac{1}{2} - 2\alpha^2)$. If $\alpha = \frac{1}{2}$ (the exact eigenstate), $E_L(x) \equiv \frac{1}{2}$ for all $x \in \mathbb{R}$. This is an analytical consequence of the zero-variance principle.
  * **Empirical Discovery:** Randomness-source drift severely corrupts the spatial distribution $\langle x^2 \rangle$ ($>1000\%$ error) while the energy observable remains perfectly invariant with zero empirical variance.
* **Critical Caveat:**
  * Observable masking occurs **only** when the trial wavefunction is close to the exact eigenstate ($\alpha \approx 0.5$). If the trial wavefunction is suboptimal (e.g. $\alpha = 0.7$), energy is not constant and drift causes immediate energy shifts.
* **Action:** Clarify that observable masking is a hazard specifically in *converged / high-quality* variational states, making it especially dangerous in advanced simulation stages.

### 2.3 MATH-04: Non-Identifiability & Lipschitz Drift Bounds
* **Classification:**
  * **Mathematical Identity:** The Non-Identifiability Theorem is mathematically rigorous.
  * **Empirical Claim to Correct:** The claim that *"Lipschitz structural bound achieves 0.0% false reassurance universally"* was based on 96 synthetic checkpoints.
* **Rigorous Clarification:**
  * A lower bound $\hat{N}_\delta(t) = \lfloor (\theta_{crit} - \hat{\theta}_t - \epsilon_\delta)/L \rfloor$ guarantees zero false reassurance **if and only if** three physical conditions hold:
    1. The physical drift satisfies $|\theta_{t+1} - \theta_t| \le L$ almost surely.
    2. The estimation error margin $\epsilon_\delta$ is valid with confidence $1-\delta$.
    3. The degradation process is continuous (no discontinuous step jumps).
  * If a physical component suffers a sudden catastrophic failure (e.g. stuck bit or breakdown), $L$ is violated and the guarantee fails.

### 2.4 MATH-05: Information Functional Equations & XOR Conditioning
* **Critical Audit of Assumptions:**
  * The parity lemma $\epsilon_{XOR} = 2^{k-1}\prod_{i=1}^k \epsilon_i$ strictly assumes **independent** random variables.
  * In physical systems, multiple channels from the same chip, bus, or power rail often share thermal or electromagnetic coupling. If two channels have positive covariance $\text{Cov}(S_1, S_2) > 0$, XOR conditioning does **not** extinguish bias at the exponential rate, and can even preserve correlation.
* **Action:** Implement an explicit runtime mutual independence test before applying XOR conditioning in `rnglab14.sources`.

### 2.5 MATH-06: Estimator Memory Dynamics
* **Validated Findings:**
  * Cumulative estimators exhibit $O(\frac{\tau_D \Delta}{T})$ algebraic memory persistence.
  * Exponential moving averages exhibit geometric/exponential relaxation.
  * Reset estimators provide immediate zero bias at the cost of variance inflation $\frac{T}{T - \tau_D}$.
* **Status:** Scientifically sound; no corrections required.

---

## 3. Multi-Testing and Statistical Significance Audit

In the ongoing v0.13 campaign daemon, research cycles repeated thousands of times over randomized seeds:
- **Risk:** Repeatedly testing null hypotheses (e.g. $H_0: \text{Var}_{emp} = \text{Var}_{theory}$) at nominal significance $\alpha = 0.05$ across 2,800+ cycles guarantees thousands of false rejections (Type I errors) due to the multiple comparisons problem.
- **Remedy for v0.14:**
  1. Implement Family-Wise Error Rate (FWER) control via the Bonferroni-Holm procedure.
  2. Implement False Discovery Rate (FDR) control via the Benjamini-Hochberg procedure for large cohort screening.
  3. Pre-register hypotheses and freeze test cohorts to avoid continuous-testing selection bias.

---

## 4. Summary of Mandated Changes for RNGLab v0.14
1. Replace unweighted MAE with censored survival metrics and event-conditioned MAE in `rnglab14.forecast`.
2. Add explicit uncertainty abstention (`UNKNOWN / INSUFFICIENT EVIDENCE`) to all RUS predictors.
3. Decouple stochastic channels (proposal vs acceptance vs initialization) in computational science applications.
4. Add runtime cross-correlation checks to XOR conditioners.
5. Implement multi-testing significance adjustments across all experimental registries.
