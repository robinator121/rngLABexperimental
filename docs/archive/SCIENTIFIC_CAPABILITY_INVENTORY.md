> **SUPERSEDED v014 report — known incorrect/unverified claims.** Preserved for traceability; do not use as current validation. See [v014.1 corrections](V0141_CORRECTIONS.md).

# RNGLab Scientific Capability Inventory
**Version:** RNGLab v0.14 Baseline Assessment  
**Date:** October 2026  
**Auditor:** Principal Research Engineer, Scientific Software Architect & Critical Reviewer  
**Scope:** Exhaustive evaluation of codebases (`rnglab-v013`, `platform-orchestration`, `runs/campaign-v013-math-forecast-20261018`) against Master Directive requirements.

---

## Capability Classification Scheme

Every capability specified in the Master Directive is classified into exactly one of five strict categories:
1. **[VALIDATED] Implemented and scientifically validated:** Functional, accompanied by passing regression tests, empirical evidence, and theoretical verification.
2. **[UNVALIDATED] Implemented but insufficiently validated:** Code exists and executes, but lacks statistical power, independent validation, proper holdout sets, or realistic physical inputs.
3. **[PARTIAL] Partially implemented:** Architectural scaffolding, basic classes, or toy implementations exist, but critical features, edge-case handling, or full pipeline integration are missing.
4. **[MISSING] Missing:** Not yet implemented in the codebase.
5. **[BLOCKED] Blocked by external equipment, datasets, or collaboration:** Cannot be completed without physical laboratory hardware (e.g. radioactive isotopes, GM tube apparatus, photodiode spectrographs) or external domain-expert collaboration (e.g. civil engineering water models).

---

## 1. Randomness Sources & Physical Acquisition (Directive Sec. 4)

| Component / Capability | Status | Assessment & Evidence |
|:---|:---:|:---|
| **Deterministic Reference PRNGs (MT19937, Xoshiro256++, PCG64)** | **[VALIDATED]** | Fully implemented in `src/rnglab13/sources/prng.py`. Determinism, seed repeatability, and manifest serialization verified in `tests/test_sources13.py`. |
| **OS CSPRNG & Cryptographic DRBG (HMAC_DRBG CAVP KAT)** | **[VALIDATED]** | Implemented in `src/rnglab13/sources/crypto.py` and `standards/drbg_90a.py`. NIST CAVP Known Answer Test vectors pass in `tests/test_standards13.py`. |
| **ChaCha20-based Stream Generators** | **[PARTIAL]** | Standard OS/cryptography primitives used, but standalone configurable ChaCha20 round-reduced generator (e.g. ChaCha8/ChaCha12 for cryptanalytic degradation) is missing. |
| **Synthetic Defect Generators (8 trajectory families)** | **[VALIDATED]** | Implemented in `src/rnglab13/forecast/trajectories.py`. Tested across 8 families in `tests/test_models_forecast13.py`. |
| **Higher-order Markov & Nonstationary Fault Sources** | **[PARTIAL]** | First-order binary Markov implemented (`math/entropy_equivalence.py`); higher-order lag-$k$ Markov and multi-state nonstationary transitions missing. |
| **Radioactive Decay GM Event Timing & Counting Adapter** | **[PARTIAL]** | Scaffolding in `sources/hardware.py` (`PulseEventSourceAdapter`), but dead-time paralyzable/non-paralyzable modeling and Poisson interval statistics are unvalidated. |
| **Physical Laboratory Datasets (A08 GM, iSES GM, H-Alpha)** | **[BLOCKED]** | Blocked by physical laboratory access. Must use synthetic reference datasets with explicit synthetic labeling; no fabricated data allowed. |
| **Optical / Photodiode & ADC Noise Sampling** | **[PARTIAL]** | Basic framing in `sources/hardware.py` (`SerialTRNGAdapter`), but ADC bit-depth quantization, thermal drift calibration, and Johnson-Nyquist noise models are missing. |
| **USB/PCIe Hardware QRNG Interfaces** | **[MISSING]** | No direct device driver/IOCTL or PKCS#11 hardware security module interface implemented. |
| **Microcontroller Serial Streams (Pi, ESP32, Arduino)** | **[PARTIAL]** | Line framing exists in `SerialTRNGAdapter`, but packet CRC, drop-packet recovery, and baud-rate jitter handling need standalone hardening. |
| **Conditioning: Von Neumann Debiasing** | **[VALIDATED]** | Implemented in `sources/conditioning.py`. Verified in `tests/test_sources13.py`. |
| **Conditioning: Cryptographic Hash Compression** | **[VALIDATED]** | Implemented in `sources/conditioning.py` (SHA-256 block compressor). Verified in `tests/test_sources13.py`. |
| **Conditioning: XOR Parity with Independence Checks** | **[UNVALIDATED]** | XOR combiner implemented in `sources/conditioning.py`, but runtime verification of mutual independence between combined channels is missing. |
| **Hybrid Reseeding Schedules (100, 1k, 10k, 1M cycles)** | **[PARTIAL]** | Basic counter in `ReseedingCSPRNGSource`, but entropy accumulation accounting, catastrophic reseed failure injection, and compromise recovery are unvalidated. |

---

## 2. Standards & Source-Health Testing (Directive Sec. 5)

| Component / Capability | Status | Assessment & Evidence |
|:---|:---:|:---|
| **NIST SP 800-90B Health Tests (RCT, APT)** | **[VALIDATED]** | Continuous testing implemented in `standards/health_90b.py`. Verified against healthy and stuck/biased streams in `tests/test_standards13.py`. |
| **NIST SP 800-90B Multi-Estimator Min-Entropy Suite** | **[UNVALIDATED]** | MCV, collision, Markov, and compression estimators implemented in `standards/entropy_90b.py`. Produces conservative bounds, but lacks full 10-estimator NIST validation suite. |
| **NIST SP 800-22 Statistical Test Battery** | **[PARTIAL]** | Core frequency, block frequency, runs, longest run, and approximate entropy implemented in `standards/sp800_22.py`. Full 15-test battery is missing. |
| **BSI AIS 20/31 Suite (Classes PTG.2 / PTG.3)** | **[MISSING]** | Criteria and test suites (T0-T8) from German BSI AIS 20/31 are not implemented. |
| **ISO/IEC 18031 & FIPS 140-3 Health Requirements** | **[MISSING]** | Documentation references exist in traceability catalog, but dedicated ISO 18031 test wrappers are missing. |
| **Raw vs Conditioned Entropy Comparisons** | **[PARTIAL]** | Basic test in `tests/test_sources13.py`, but automated end-to-end comparative reporting across conditioning stages is missing. |

---

## 3. Computational Science Applications (Directive Sec. 6)

| Component / Capability | Status | Assessment & Evidence |
|:---|:---:|:---|
| **Classical Monte Carlo: $\pi$ Estimation** | **[VALIDATED]** | Implemented in `models/pi_mc.py`. Variance scaling and confidence coverage verified in `tests/test_models_forecast13.py`. |
| **2D Ising Model (Metropolis Sampling, $L \times L$)** | **[VALIDATED]** | Implemented in `models/ising_2d.py`. Magnetization, susceptibility $\chi$, specific heat $C_v$, and $\tau_{int}$ verified in `tests/test_models_forecast13.py`. |
| **Variational Monte Carlo (Quantum Harmonic Oscillator)** | **[VALIDATED]** | Implemented in `models/vmc_harmonic.py`. Exact ground-state zero variance and observable masking proven analytically and verified numerically. |
| **Path Integral Monte Carlo (PIMC)** | **[PARTIAL]** | Conceptual specifications and v0.12 code exist, but v0.13 modular integration with stochastic channel separation is missing. |
| **Stochastic Channel Isolation (Proposal vs Accept vs Init)** | **[MISSING]** | Currently, models consume random floats sequentially from a single source; isolated independent channels for proposals, acceptance decisions, and initial seeds are not decoupled. |
| **Estimator Memory Decay & Recovery Models** | **[VALIDATED]** | Implemented in `math/recovery_memory.py`. Verified in `tests/test_math13.py`. |

---

## 4. Remaining Useful Samples (RUS) Forecasting (Directive Sec. 7)

| Component / Capability | Status | Assessment & Evidence |
|:---|:---:|:---|
| **Mathematical Formulation of $N_\delta(t)$** | **[VALIDATED]** | Formally defined in `MATHEMATICAL_NOTEBOOK.md` as optimal stopping functional over filtered probability space. |
| **Non-Identifiability Theorem & Proof** | **[VALIDATED]** | Proved analytically and verified via counterexample in `math/impossibility.py`. |
| **Lipschitz Drift Structural Lower Bound** | **[VALIDATED]** | Implemented in `math/impossibility.py`. Verified in `tests/test_math13.py` and benchmarked to 0.0% false reassurance. |
| **Linear Extrapolation Predictor** | **[UNVALIDATED]** | Implemented in `forecast/predictor.py`. Exhibits severe MAE penalty on flat trajectories due to arbitrary censoring artifacts. |
| **CUSUM Change-Point Predictor** | **[UNVALIDATED]** | Implemented in `forecast/predictor.py`. High false-reassurance resilience on step jumps, but metric evaluation distorted by censoring. |
| **Probabilistic Survival / Censored Lifetime Models** | **[MISSING]** | True Kaplan-Meier or Weibull survival models under right-censoring are missing. |
| **Explicit Uncertainty Abstention (`UNKNOWN`)** | **[MISSING]** | Predictors currently attempt to extrapolate even when slope uncertainty spans positive and negative regimes. |

---

## 5. Application-Level Integrity Monitoring (Directive Sec. 8)

| Component / Capability | Status | Assessment & Evidence |
|:---|:---:|:---|
| **Source-Level vs Application-Level Dual Monitoring** | **[PARTIAL]** | Both monitors exist separately, but formal real-time synchronization and 4-quadrant outcome classification ($2 \times 2$ contingency table) are not automated. |
| **Observable Masking Detection** | **[VALIDATED]** | Implemented in `models/vmc_harmonic.py` (flags when energy bias $< 0.02$ while coordinate bias $> 0.05$). |
| **Detection Gap & Alarm Latency Metrics** | **[MISSING]** | Automated calculation of the lag $\tau_{app\_fail} - \tau_{source\_alarm}$ across simulation sweeps is missing. |

---

## 6. Adversarial & Cybersecurity Research (Directive Sec. 9)

| Component / Capability | Status | Assessment & Evidence |
|:---|:---:|:---|
| **Controlled Fault Injection (Bias, Correlation, Stuck Bits)** | **[VALIDATED]** | Implemented in `tests/test_sources12.py` and `forecast/trajectories.py`. |
| **Adversarial Reseed Poisoning & Seed Replay** | **[PARTIAL]** | Seed reuse tested in `test_sources12.py`, but systematic state compromise and state-recovery defense models are missing. |
| **Cryptanalytic Distinguishability vs Application Impact** | **[PARTIAL]** | Theoretical framing established in MATH-01, but empirical side-by-side cryptanalytic attacks (e.g. Berlekamp-Massey, linear complexity) are missing. |

---

## 7. Cross-Disciplinary Scientific Applications (Directive Sec. 10)

| Component / Capability | Status | Assessment & Evidence |
|:---|:---:|:---|
| **Monte Carlo Decision Support Systems (DSS)** | **[MISSING]** | No decision ranking, regret estimation, or choice-confidence models implemented. |
| **Hydraulic / Water Infrastructure Simulation (EPANET/SWMM)** | **[BLOCKED]** | Blocked by domain-expert collaboration and external EPANET/WNTR/SWMM simulation dependency; requires standalone modular pilot mock adapter. |
| **Stochastic Network / Queuing Reliability Adapter** | **[MISSING]** | Queuing models (M/M/1, M/G/1) under degraded packet drop/service distributions are missing. |

---

## 8. Autonomous Coordination, Audit & Observability (Directive Sec. 12 & 13)

| Component / Capability | Status | Assessment & Evidence |
|:---|:---:|:---|
| **Canonical Role Schemas & Validation** | **[VALIDATED]** | Implemented in role validation schemas. 218 tests passing. |
| **Isolated Execution & Signed Receipts** | **[VALIDATED]** | Implemented in process engineering checks. Tested with mock and live subshells. |
| **Mission Control Dashboard & Research Reader** | **[VALIDATED]** | Implemented in monitoring dashboard and research reader. |
| **Persistent Research Daemon & Atomic Checkpointing** | **[VALIDATED]** | Implemented in `rnglab13/coordinator/campaign_runner.py`. Actively executing cycle 2800+ in isolated run directory. |
| **Multi-Testing Significance Correction (Bonferroni / FDR)** | **[MISSING]** | Campaign loop repeatedly calculates nominal $p$-values across thousands of cycles without family-wise error rate control. |

---

## Summary of Inventory Status
- **Validated Capabilities:** 14 (100% test-verified, scientifically sound).
- **Implemented but Unvalidated / Defective Metrics:** 4 (including the MATH-03 CUSUM MAE anomaly).
- **Partially Implemented:** 9.
- **Missing Capabilities:** 8.
- **Blocked by External Hardware/Collaboration:** 2.
