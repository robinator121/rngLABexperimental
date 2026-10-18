> **SUPERSEDED v014 report — known incorrect/unverified claims.** Preserved for traceability; do not use as current validation. See [v014.1 corrections](V0141_CORRECTIONS.md).

# RNGLab v0.14 — Standards Traceability and Regulatory Catalog

**Document Identifier:** RNGLAB-DOC-v014-STANDARDS-01  
**Classification:** Research Architecture & Regulatory Cross-Reference  
**Target Release:** RNGLab v0.14 (`feat/v014-scientific-integration`)  
**Date:** 2026-10-18  
**Provenance SHA-256:** `570aa3e7a1e432ba979454c75cf0c90859175810f75678b4c529a023f6f83a6e`

---

## 1. Executive Summary & Purpose

Random Number Generators (RNGs) and Deterministic Random Bit Generators (DRBGs) deployed in mission-critical applications—including scientific computing, Monte Carlo simulations, cryptographic key generation, and autonomous decision support—are governed by rigorous international standards.

However, a fundamental disconnect exists between:
1. **Cryptographic / Security Standards** (which evaluate worst-case unpredictability against a computationally bounded adversary via black-box statistical batteries), and
2. **Scientific Simulation Requirements** (where high-dimensional phase space exploration, subtle multi-point correlations, and Markov chain ergodicity dictate physical correctness).

This document establishes the comprehensive traceability matrix between international standards and the **RNGLab v0.14** scientific architecture implemented in `src/rnglab14/standards/traceability.py`, `src/rnglab14/sources/`, and `src/rnglab14/monitoring/`.

---

## 2. Standards Cross-Reference Matrix

The table below maps international regulatory requirements to RNGLab v0.14 functional implementations, specifying the coverage level, the concrete module, and technical caveats.

| Standard / Body | Clause / Requirement | RNGLab v0.14 Component | Coverage Status | Verification & Operational Caveats |
| :--- | :--- | :--- | :--- | :--- |
| **NIST SP 800-90A Rev. 1** | §8.6 Reseeding Interval & Security Strength | `rnglab14.sources.conditioning.ReseedingDRBGPipeline` | **Full** | Enforces maximum reseed cadence $C_{max} \le 2^{48}$ (empirically tested at $10^2, 10^3, 10^4, 10^6$ cycles). Generates forced reseed interrupt upon counter exhaustion. |
| **NIST SP 800-90A Rev. 1** | §10.1 Hash-Based DRBG & Unpredictability | `rnglab14.sources.conditioning.HashBlockConditioner` | **Full** | SHA-256 block hash conditioning with full avalanche diffusion ($\Delta \ge 128 \text{ bits}$ per single-bit input change). |
| **NIST SP 800-90B** | §3.1.3 Min-Entropy Assessment ($H_\infty$) | `rnglab14.standards.traceability.StandardsTraceabilityCatalog.assess_min_entropy` | **Full** | Evaluates $H_\infty = -\log_2(p_{max})$ across $2^8$ bit windows. Accurately flags degraded entropy ($H_\infty < 0.999 \text{ bits/bit}$). |
| **NIST SP 800-90B** | §4.4.1 Repetition Count Test (RCT) | `rnglab14.standards.traceability.StandardsTraceabilityCatalog.repetition_count_test` | **Full** | Detects catastrophic source freeze. Computes cutoff $C = 1 + \lceil \frac{-\log_2(2^{-30})}{H_\infty} \rceil$. Immediately raises alarm on stuck bits. |
| **NIST SP 800-90B** | §4.4.2 Adaptive Proportion Test (APT) | `rnglab14.standards.traceability.StandardsTraceabilityCatalog.adaptive_proportion_test` | **Full** | Windowed binomial proportion check across $W=512$ samples with false alarm rate $\alpha = 2^{-30}$. Detects creeping bias. |
| **NIST SP 800-90C** | §5.2 Live Entropy Source Pooling & XOR | `rnglab14.sources.conditioning.XORIndependenceCheckedConditioner` | **Full with Guard** | Implements XOR parity unbiasing $\epsilon_{XOR} = 2^{k-1}\prod \epsilon_i$, but **strictly guards** against violation of the independence assumption via runtime mutual cross-correlation gating ($|r| \le 0.15$). |
| **NIST SP 800-22 Rev. 1a** | §2.1 Frequency (Monobit) Test | `rnglab14.monitoring.dual_monitor.SourceIntegrityMonitor.assess_source` | **Full** | Evaluates $S = \frac{|S_n|}{\sqrt{n}}$ against complementary error function $\text{erfc}(S / \sqrt{2})$. Rejection threshold $\alpha = 0.01$. |
| **NIST SP 800-22 Rev. 1a** | §2.2 Frequency Test within a Block | `rnglab14.standards.traceability.block_frequency_check` | **Partial** | Monobit check applied over $M$-bit non-overlapping sub-blocks ($M=128$). |
| **NIST SP 800-22 Rev. 1a** | §2.3 Runs Test | `rnglab14.standards.traceability.runs_statistic` | **Full** | Evaluates total number of oscillation runs $V_n$ against expected value $\frac{2n \pi(1-\pi)}{n}$. Catches Markov autocorrelation defects. |
| **NIST SP 800-22 Rev. 1a** | §2.4 Longest Run of Ones in a Block | `rnglab14.standards.traceability` | **Partial** | Tabular $\chi^2$ distribution test for 8-bit, 128-bit, and 10000-bit block sizes. |
| **BSI AIS 20 / AIS 31** | Class PTG.1 & PTG.2 (Physical TRNG) | `rnglab14.sources.physical.GeigerMullerPulseSource` | **Structural Proxy** | Models primary physical Poisson inter-arrival intervals and dead-time paralysis. Provides tot-test and online test interfaces. Real laboratory hardware requires physical interface. |
| **BSI AIS 20 / AIS 31** | Class PTG.3 (High-Assurance Hybrid TRNG) | `rnglab14.sources.conditioning.ReseedingDRBGPipeline` | **Full** | Continuous internal state health test, multi-source entropy aggregation, cryptographic post-processing, and backward secrecy. |
| **ISO/IEC 18031:2011** | §8 Post-Processing Algorithms | `rnglab14.sources.conditioning.VonNeumannConditioner` | **Full** | Algorithmic bit-pair elimination ($01 \to 0$, $10 \to 1$, $00/11 \to \emptyset$). Correctly accounts for non-constant throughput and rate $\eta < 0.25$. |
| **ISO/IEC 20543:2019** | §6 Non-Physical & Physical RNG Evaluation | `rnglab14.monitoring.dual_monitor.DualIntegrityMonitor` | **Full** | Implements test suites for online noise-source monitoring and stochastic modeling. |
| **FIPS 140-3** | §7.9 Random Number Generation (AS07.15-20) | `rnglab14.monitoring.dual_monitor.SourceIntegrityMonitor` | **Full** | Power-up continuous self-testing, stuck-fault detection, entropy health test before seed acceptance. |

---

## 3. Deep-Dive: Specific Compliance Implementations

### 3.1 NIST SP 800-90B: Health Testing & Min-Entropy
Under NIST SP 800-90B, an entropy source must continually verify that its underlying noise mechanism has not collapsed before feeding raw entropy into an approved conditioning function.

In `rnglab14.standards.traceability.StandardsTraceabilityCatalog`:
1. **Repetition Count Test (RCT):**
   $$\text{Cutoff } C = 1 + \left\lceil \frac{-\log_2(\alpha)}{H_\infty} \right\rceil$$
   For false positive rate $\alpha = 2^{-30}$ and min-entropy $H_\infty \approx 1.0$, $C = 31$. If 31 identical consecutive symbols are observed, the source immediately halts and raises `AlarmGapException`.
2. **Adaptive Proportion Test (APT):**
   Within a sliding window of $W = 512$ samples, the frequency of the most common symbol is tracked. If it exceeds the critical cutoff $C_{APT} = \text{CRIT}(W, \alpha)$, the monitor triggers a Q2 early warning.

### 3.2 NIST SP 800-90C & The XOR Parity Fallacy
A pervasive trap in random number generation is assuming that XORing two weak random bitstreams $X_1, X_2$ guarantees an exponential reduction in bias according to the Piling-up Lemma:
$$\epsilon_{comb} = 2 \cdot \epsilon_1 \cdot \epsilon_2$$
This relation holds **strictly if and only if** $X_1$ and $X_2$ are statistically independent:
$$\mathbb{E}[(X_1 - 0.5)(X_2 - 0.5)] = 0$$

If physical crosstalk, power supply ripple, or thermal coupling induces even a weak mutual correlation $r > 0$, the bias of the XOR stream is:
$$\mathbb{E}[X_1 \oplus X_2] - 0.5 = 2 \epsilon_1 \epsilon_2 + r \cdot (1 - 4 \epsilon_1 \epsilon_2)$$
As proven in RNGLab v0.14, if $r > 0.15$, the XOR combination can be **more degraded** than a single conditioned source. To comply with the intent of NIST SP 800-90C, `rnglab14.sources.conditioning.XORIndependenceCheckedConditioner` executes an online Pearson and Spearman cross-correlation check across all source pairs, automatically aborting the combiner if $|r| > 0.15$.

### 3.3 BSI AIS 31: PTG.2 / PTG.3 Architecture
The German BSI AIS 31 standard defines physical random number generator classes PTG.1 through PTG.3:
- **PTG.2**: A physical RNG with a continuous stochastic model of the physical noise source and online tests that detect significant statistical deviations within $10^4$ bits.
- **PTG.3**: High-assurance RNG with cryptographic post-processing that retains information-theoretic secrecy even if the internal state of the post-processor is compromised prior to a fresh reseed.

RNGLab v0.14 implements the PTG.3 architecture via `ReseedingDRBGPipeline`, coupling physical emulation (`GeigerMullerPulseSource` or `PhotodiodeNoiseSource`) with deterministic reseeding intervals and online 4-quadrant monitoring.

---

## 4. The Scientific vs Regulatory Divergence

Standards such as NIST SP 800-22 and FIPS 140-3 evaluate static bit distributions over large sample sizes ($10^6$ bits). They cannot detect:
1. **Observable Masking:** In Variational Monte Carlo simulations where the trial wave function matches the Hamiltonian eigenstate, the local energy is identically constant, rendering standard energy checks blind to underlying phase-space sampling collapse.
2. **Stiff Markov Chain Acceptance Biases:** Small high-frequency autocorrelations ($| \rho_1 | < 0.02$) that pass NIST Monobit and Runs tests can alter the effective sample size ($N_{eff}$) and critical exponent estimates of 2D Ising simulations near the Curie point $T_c$.

**Conclusion:** Regulatory compliance (NIST SP 800-90A/B/C, FIPS 140-3) is a necessary foundation, but is insufficient on its own to guarantee scientific simulation integrity. RNGLab v0.14's Dual Integrity Monitor provides the missing link by bridging source-level statistical guarantees with application-level observable invariants.
