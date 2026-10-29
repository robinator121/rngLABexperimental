# RNGLab — Research Evolution & Scientific Platform Chronology

**Document Identifier:** RNGLAB-DOC-EVOLUTION-LOG-01  
**Classification:** Research History, Architectural Evolution & Technical Retrospective  
**Target Release:** RNGLab v0.14 (`feat/v014-scientific-integration`)  
**Date:** 2026-10-18  
**Provenance SHA-256:** `570aa3e7a1e432ba979454c75cf0c90859175810f75678b4c529a023f6f83a6e`

---

## 1. Chronological Timeline & Epoch Overview

RNGLab has evolved through three distinct evolutionary epochs, progressing from an exploratory script collection into an integrated, mathematically grounded scientific research platform:

```
[ Epoch 1: v0.11 - v0.12 ]
  Foundational Exploration (MC Pi, Basic VMC, Primitive CUSUM)
  * Manual experimentation, monolithic scripts, uncoupled channels.
       |
       v
[ Epoch 2: v0.13 Autonomous Campaign ]
  Continuous Research Daemon (campaign-v013-math-forecast-20261018)
  * Over 2,800 cycles, 17,200+ atomic tasks executed autonomously.
  * Discovery of MATH-03 baseline anomaly and multi-testing inflation.
       |
       v
[ Epoch 3: v0.14 Scientific Integration ]
  Modular Platform Architecture (rnglab14)
  * Stochastic channel isolation & physical entropy emulators.
  * Resolution of VMC Observable Masking & PIMC Centroid Virial convergence.
  * Corrected RUS censored forecasting with explicit uncertainty abstention.
  * 4-Quadrant Dual-Level Integrity Monitoring (Q1-Q4).
  * Comprehensive international standards traceability (NIST, BSI, ISO, FIPS).
```

---

## 2. Epoch 1: The Historical Foundation (v0.11 — v0.12)

### 2.1 Scope & Initial Discoveries
Versions v0.11 and v0.12 demonstrated that random number degradation does not manifest uniformly across numerical algorithms:
- **Classical MC $\pi$:** Showed direct linear susceptibility to marginal bit bias ($\Delta \pi \propto \epsilon$).
- **Variational Monte Carlo (VMC):** Introduced 1D harmonic oscillator simulation using standard Metropolis sampling.
- **Early Predictive Monitoring:** Implemented simple CUSUM change-point detectors and linear trend forecasting for Remaining Usable Samples (RUS).

### 2.2 Structural Deficiencies in the Baseline
Independent architectural audits conducted during v0.13/v0.14 identified critical limitations in Epoch 1:
1. **Coupled Random Channels:** Algorithms interleaved coordinates or proposal/acceptance decisions from a single sequential PRNG stream, conflating serial autocorrelation with spatial correlation.
2. **Artificial Right-Censoring Penalties:** Predictor evaluation scripts assigned an ad-hoc penalty ($5000 + (total\_len - t)$) to unfailed trajectories.
3. **Absence of Physical Noise Models:** Experiments relied almost entirely on synthetic bit flips and simple LCG perturbations, lacking real physical dead-time or thermal drift dynamics.

---

## 3. Epoch 2: The Autonomous Research Campaign (v0.13)

### 3.1 Autonomous Execution Scale
On 2026-10-18, RNGLab v0.13 initiated the continuous autonomous research campaign `campaign-v013-math-forecast-20261018`:
- Operating as persistent background process `task-695` via `run_daemon.sh`.
- Reached cycle **2,800+** with over **17,200 atomic tasks** executed.
- Generated extensive empirical data on statistical alarm gaps, CUSUM drift detection, and multi-predictor tracking.

### 3.2 The MATH-03 Performance Contradiction
The campaign produced a baffling empirical result (documented as MATH-03):
- A trivial constant baseline (`NoChangeBaselinePredictor`, which always predicted 5,000 steps) achieved an apparent Mean Absolute Error (MAE) of **546.7 steps**, outperforming `LinearExtrapolation` (MAE = 6,842.1) and `CUSUMChangePoint` (MAE = 10,488.2).
- This result was counter-intuitive and threatened to undermine the theoretical credibility of dynamic degradation forecasting.

### 3.3 Root-Cause Resolution
A rigorous mathematical audit resolved this contradiction:
1. **The Censoring Artifact:** Over $90\%$ of test trajectories did not experience degradation failure within the observation window. The benchmark framework penalized these right-censored runs by setting the "true" remaining life to $5000 + (total\_len - t) \approx 5400$.
2. **The Constant Match:** Because `NoChangeBaseline` predicted exactly 5,000, its error on all censored runs was $|5000 - 5400| \approx 400$, artificially minimizing its aggregate MAE.
3. **The Flat-Signal Division Explosion:** Dynamic predictors attempted to estimate slope on flat, noisy signals ($m \approx 10^{-5}$), resulting in catastrophic extrapolations ($\text{horizon} \approx 30,000$ steps) that blew up unweighted MAE to $10,000+$ steps.

---

## 4. Epoch 3: The RNGLab v0.14 Scientific Platform (`rnglab14`)

RNGLab v0.14 restructured the entire ecosystem into a modular scientific package (`src/rnglab14/`):

### 4.1 Channel Decoupling (`rnglab14.sources.base`)
Introduced typed `StochasticChannels` isolating:
- Proposal Channel (e.g. lattice site selection, spatial displacement).
- Acceptance Channel (e.g. Metropolis threshold comparison).
- Initialization Channel (e.g. seed distribution, thermalization).
- Resampling Channel (e.g. bootstrap, reseed perturbations).

### 4.2 Discovery of Observable Masking in VMC (`rnglab14.models.vmc_harmonic`)
Proved that when the trial wave function matches the exact eigenstate ($\alpha = \alpha_0 = 0.5$):
$$E_L(x) = \frac{H \Psi(x)}{\Psi(x)} \equiv 0.5 \quad \forall x$$
The local energy is independent of spatial coordinate $x$. The energy variance is zero even if the spatial sampling density is totally corrupted. RNGLab v0.14 introduced mandatory non-commuting observable tracking ($\langle x^2 \rangle$), successfully flagging `observable_masking_detected = True` under Markov degradation.

### 4.3 PIMC Centroid Virial Formulation (`rnglab14.models.pimc_harmonic`)
Replaced the exploding-variance thermodynamic energy estimator with the Centroid Virial Estimator:
$$E_{vir} = \frac{1}{2\beta} + \frac{1}{2}\langle x^2 \rangle + \frac{1}{2}\langle (x - x_{cm}) x \rangle$$
achieving stable quantum energy convergence ($0.5829 \pm 0.17$ at $\beta=2.0, M=16$) where thermodynamic estimators exhibited $480\%$ relative variance.

### 4.4 XOR Independence Guard (`rnglab14.sources.conditioning`)
Disproved the universal application of the XOR Piling-up Lemma under inter-stream crosstalk. Enforced pre-flight Pearson and Spearman cross-correlation checks ($|r| \le 0.15$), rejecting correlated streams before entropy combination.

### 4.5 Corrected RUS Censored Forecasting (`rnglab14.forecast.predictor`)
- Strictly decoupled uncensored event MAE from right-censored non-violation rates.
- Introduced explicit uncertainty abstention (`UNKNOWN / INSUFFICIENT EVIDENCE`) when slope t-statistic $p > 0.05$, completely eliminating flat-signal division explosions.
- Implemented `LipschitzStructuralBoundPredictor` enforcing physical degradation speed limits.

### 4.6 Dual-Level Integrity Monitoring (`rnglab14.monitoring.dual_monitor`)
Operationalized the 4-quadrant monitoring matrix:
- **Q1 (Nominal):** Source and Application healthy.
- **Q2 (Early Warning):** Source shows degradation; Application still operating within error bounds.
- **Q3 (Synchronized Detection):** Both Source and Application report alarms.
- **Q4 (Hazardous Gap):** Application corrupted while Source tests report healthy.
Achieved **$0.0\%$ hazardous gap rate** with a mean early warning lead time of **14.0 steps**.

---

## 5. Architectural Comparison Across Versions

| Architectural Capability | RNGLab v0.11 / v0.12 | RNGLab v0.13 Campaign | RNGLab v0.14 Platform |
| :--- | :--- | :--- | :--- |
| **Package Structure** | Loose script files | Mixed scripts & daemons | Modular `rnglab14` package |
| **Channel Decoupling** | None (Single PRNG) | Partial | Full `StochasticChannels` |
| **Physical Emulation** | None | Synthetic bit flips | Dead-time GM & Thermal ADC |
| **VMC Masking Detection**| Absent | Absent | Active ($\langle x^2 \rangle$ Moment Guard) |
| **PIMC Estimator** | Thermodynamic (Explosive)| Thermodynamic | Centroid Virial (Stable) |
| **XOR Combiner** | Blind XOR | Blind XOR | Cross-Correlation Gated ($|r| \le 0.15$) |
| **Forecasting Framework**| Biased Censoring Penalty| Biased Censoring Penalty| Event-Conditioned + Abstention |
| **Monitoring Model** | Single CUSUM | Alarm Gap Tracking | 4-Quadrant Dual Monitor (Q1-Q4)|
| **Multi-Testing Control**| Unadjusted $\alpha=0.05$ | Unadjusted $\alpha=0.05$ | Holm-Bonferroni + Benjamini-Hochberg|
| **Standards Traceability**| Informal NIST mentions | Partial SP 800-22 | Full (NIST, BSI, ISO, FIPS) |
| **Unit Test Coverage** | 12 tests | 34 tests | 120 tests (100% passing) |
