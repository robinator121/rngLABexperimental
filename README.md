# RNGLab: Computational Randomness, Physical Entropy Acquisition & Stochastic Simulation Integrity

[![CI](https://github.com/robinator121/rngLABexperimental/actions/workflows/ci.yml/badge.svg)](https://github.com/robinator121/rngLABexperimental/actions/workflows/ci.yml)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](LICENSE)
[![Tests: 238 Passing](https://img.shields.io/badge/Tests-238%20Passing-brightgreen.svg)](tests/)

**Author:** Hevele Imre Róbert ([@robinator121](https://github.com/robinator121))  
**Contact:** [roberthevele@gmail.com](mailto:roberthevele@gmail.com)  
**Affiliation:** Faculty of Informatics, University of Debrecen  
**Project Nature:** Academic Research Prototype & BSc Thesis Codebase  
**Status:** Experimental Research Platform (v0.17 CALIBRATE)

---

## 1. Executive Summary

**RNGLab** is an experimental scientific research platform developed at the University of Debrecen to investigate how subtle, non-catastrophic degradations in pseudorandom and physical entropy streams propagate through stochastic simulations, and how such defects can be detected before corrupting physical or financial observables.

The project bridges three complementary domains:
1. **Physical Entropy Acquisition & Waveform Auditing:** Real-time hardware noise harvesting (avalanche breakdown in reverse-biased PN junctions, Geiger-Müller radiation pulse intervals, and cosmic ray / muon detector digitized waveforms) coupled with NIST SP 800-90B min-entropy evaluation.
2. **Computational Physics & Stochastic Modeling:** Markov Chain Monte Carlo (2D Ising Metropolis and Wolff cluster), Path Integral Monte Carlo (PIMC virial estimators), and Variational Monte Carlo (VMC harmonic oscillator) simulations evaluating observable masking and spectral gap collapse.
3. **Financial Risk & Model Risk Stress-Testing:** Quantitative financial engines (European options, Portfolio Value at Risk, and Expected Shortfall) designed to stress-test stochastic simulation stability under quantitative model validation principles (inspired by Federal Reserve SR 11-7, ECB TRIM, and Basel III guidelines).

> [!IMPORTANT]
> **Academic & Experimental Scope:** RNGLab is an exploratory academic prototype, not a commercial turnkey product. It implements methodological frameworks inspired by regulatory model risk guidance, but is **not officially certified by NIST, the Federal Reserve, the ECB, or other regulatory authorities**.

> [!WARNING]
> **Production Cryptography Disclaimer:** The conditioning routines implemented in this repository (von Neumann extraction, SHA-256 block hashing, ChaCha20 expansion) are research treatments designed for stochastic simulation benchmarking. They are **not certified or audited for production cryptographic systems or cryptographic key generation**. Passing statistical batteries (such as NIST SP 800-22) after cryptographic hashing demonstrates only the PRF diffusion of the hash function; true source unpredictability must be assessed directly on raw, unconditioned noise via NIST SP 800-90B min-entropy estimation.

---

## 2. Key Headline Findings & Experimental Results

1. **Hardware Digitizer Saturation & Clock Lattice Artifacts:**
   In the author's physical cosmic ray detector dataset (117,962 events, 120.7M samples), **91.96% of captured pulse peaks saturate the 16-bit digitizer ceiling** (code 65,535), and recorded timestamps lie on a discrete **80-unit clock lattice** ($\text{GCD} = 80$). Extracting raw timestamp LSBs yields identically zero bits, demonstrating why raw hardware inspection is mandatory before claiming Poisson entropy ([`docs/physical/MUON_DETECTOR_ASSESSMENT.md`](docs/physical/MUON_DETECTOR_ASSESSMENT.md), [`scripts/inspect_detector_export.py`](scripts/inspect_detector_export.py)).

2. **VMC Zero-Variance Blind Spot vs. Imperfect Wavefunctions:**
   Illustrates the zero-variance property and its blind spot: at the exact harmonic oscillator ground state ($\alpha = 0.5$), the local energy is constant ($E_L(x) \equiv 0.5$) with mathematically zero variance ($\sigma^2_{E} \equiv 0$). This completely masks severe proposal collapse ($\Delta \langle x^2 \rangle = 1.75$) behind a deceptively pristine energy bias $\Delta \langle E \rangle = 0.000000$. Contrasted against an imperfect trial wavefunction ($\alpha = 0.6$), local energy does depend on sampling position, immediately exposing the proposal corruption with an energy error $\Delta \langle E \rangle = -0.26$ and standard deviation $\sigma_E = 0.31$ ([`src/rnglab14/models/vmc_harmonic.py`](src/rnglab14/models/vmc_harmonic.py), [`tests/test_models14.py`](tests/test_models14.py)).

3. **Parallel Seed Collisions & False Confidence in Financial Risk:**
   When parallel Monte Carlo workers accidentally reuse seed states (parallel seed collision fault), multi-asset portfolio Value at Risk (95% VaR) standard error collapses artificially by a factor of up to $\sqrt{K}$ ($K$ parallel workers). This generates false statistical confidence intervals while masking severe risk underestimation under stress scenarios ([`src/rnglab_finance/faults/injection.py`](src/rnglab_finance/faults/injection.py), [`src/rnglab_finance/models/portfolio_var.py`](src/rnglab_finance/models/portfolio_var.py), [`tests/test_finance_faults.py`](tests/test_finance_faults.py)).

4. **Metropolis Acceptance Bias & Ising Ergodicity Breakdown:**
   In 2D Ising spin lattices, directional acceptance corruption $\delta \ge e^{-8/T}$ creates absorbing magnetized ground states, breaking ergodicity and causing perturbation approximations to fail near the critical temperature ($T_c \approx 2.269$) ([`docs/trident/REVIEW_TRIDENT_AUDIT.md`](docs/trident/REVIEW_TRIDENT_AUDIT.md), [`scripts/trident_report.py`](scripts/trident_report.py)).

---

## 3. Core Scientific Capabilities

### 3.1 Physical Entropy Harvesting & Cosmic Ray Acquisition (`rnglab15`, `hardware/`, `docs/physical/`)
- **Author's Physical Acquisition Dataset:** Continuous live acquisition conducted by Hevele Imre Róbert at the University of Debrecen (October 19, 2026, 14:05:19–14:27:23 CEST). Utilizes a particle detector bench with fast scintillation front-end, high-voltage bias, and a 16-bit digitizer exporting 1,023-sample waveforms over USB/serial. The 606.6 MB export contains 117,962 records and 120.7M waveform codes ([`docs/physical/MUON_DETECTOR_ASSESSMENT.md`](docs/physical/MUON_DETECTOR_ASSESSMENT.md)).
- **Microcontroller Hardware Interface:** Custom Arduino C++ firmware ([`hardware/arduino/RNGLabSensors.ino`](hardware/arduino/RNGLabSensors.ino)) interfacing with hardware noise diodes, optical PN junctions, and pulse-counting sensors via high-speed serial streaming.
- **NIST SP 800-90B Entropy & Health Testing:**
  - Online continuous health monitoring: **Repetition Count Test (RCT)** and **Adaptive Proportion Test (APT)** with exact binomial tail cutoffs.
  - Raw unconditioned min-entropy ($H_\infty$) estimation, collision tests, and Markov dependency bounds (`src/rnglab13/standards/entropy_90b.py`, `src/rnglab15/assessment.py`).

### 3.2 Computational Physics & Markov Chain Dynamics (`rnglab14`, `rnglab16` TRIDENT)
- **2D Ising Model:** Metropolis-Hastings local spin-flip dynamics alongside Wolff single-cluster algorithms on periodic $L \times L$ lattices ($L \in \{2, 3, 4, 8\}$). Validated against exact finite-size analytical enumeration states.
- **Variational Monte Carlo (VMC):** Quantum harmonic oscillator ground state simulation illustrating the zero-variance property blind spot and dual-level monitoring.
- **Path Integral Monte Carlo (PIMC):** Finite-temperature continuous-space quantum virial path integral simulations ($M$ time beads), evaluating observable masking and interval undercoverage.
- **Augmented Markov Chain Sensitivity Theory:** Analytical two-state and augmented-state transition kernels evaluating spectral gap degradation under directional acceptance bias.

### 3.3 Financial Risk & Model Risk Stress-Testing (`src/rnglab_finance/`, `docs/finance/`)
- **European Option Pricing:** Risk-neutral geometric Brownian motion simulation anchored against exact closed-form Black-Scholes analytical solutions. Enforces the risk-neutral martingale invariant ($\mathbb{E}^{\mathbb{Q}}[e^{-rT}S_T] = S_0$) and Put-Call Parity.
- **Portfolio Value at Risk (VaR) & Expected Shortfall (ES):** Multi-asset correlated portfolio risk engines computing parametric vs. Monte Carlo quantiles and coherent tail expectations.
- **Fault Injection & Model Risk Stress-Testing:** Controlled injection of non-catastrophic PRNG faults (parallel seed collisions, mean shift, bit truncation, Markovian autocorrelation, and periodic leakage) to evaluate risk capital bias under Fed SR 11-7 and ECB TRIM stress-testing principles.

### 3.4 Reliability Calibration & Statistical Identifiability Limits (`rnglab13`, `rnglab17` CALIBRATE)
- **Statistical Identifiability Limits:** Formal mathematical proof and counterexample showing the impossibility of causal distribution-free remaining-useful-sample (RUS) forecasting without structural drift constraints (e.g. Lipschitz bounds), implemented in [`src/rnglab13/math/impossibility.py`](src/rnglab13/math/impossibility.py).
- **Calibrated Operational Gating:** Holm-Bonferroni step-down and Benjamini-Hochberg FDR multi-testing controls over statistical alarm families.
- **Hardware & Precision Audits:** Detailed execution and precision benchmarks comparing multi-core CPU workers against NVIDIA L4 GPU acceleration across FP32, FP64, and BF16 execution modes ([`docs/benchmarks/RNGLAB_CPU_GPU_COMPARISON.md`](docs/benchmarks/RNGLAB_CPU_GPU_COMPARISON.md)).

---

## 4. Development Methodology & AI Assistance Disclosure

In accordance with modern academic research standards:

- **Author's Original Contributions:**
  - Formulated the research hypotheses: the zero-variance blind spot in VMC, PRNG parallel seed collisions in financial tail risk, and digitizer saturation artifacts in physical detectors.
  - Assembled and performed physical hardware measurements (detector bench acquisition, Arduino C++ firmware, and circuit wiring).
  - Derived the mathematical formulations, closed-form references, and verified all 238 automated unit, mathematical, and regression tests by hand.
- **AI Tooling Assistance:**
  - Large Language Model (LLM) coding assistants were used for boilerplate generation, test scaffolding, numerical verification script templates, and initial documentation drafting.
  - Review documents in `docs/validation/` and audit logs were co-drafted using LLM review prompts and subsequently checked, corrected, and verified against actual simulation and test suite outputs by the author.

---

## 5. Repository Architecture

```text
rngLABexperimental/
├── .github/
│   └── workflows/ci.yml   # GitHub Actions continuous integration workflow
├── src/
│   ├── rnglab11/          # Core statistical evaluation, DRBG pipelines, and CLI
│   ├── rnglab12/          # Experimental review and uncertainty estimation
│   ├── rnglab13/          # Identifiability bounds and multi-stream coordinator
│   ├── rnglab14/          # Abstract source hierarchy, VMC/PIMC models, RUS predictors
│   ├── rnglab15/          # Physical acquisition, conditioning, NIST tests, and pipeline
│   ├── rnglab16/          # TRIDENT: Vectorized MC engines, exact Ising, perturbation theory
│   ├── rnglab17/          # CALIBRATE: Out-of-distribution gates, reliability audit suites
│   └── rnglab_finance/    # Financial risk models (Black-Scholes, VaR, Expected Shortfall)
├── hardware/
│   ├── arduino/           # Microcontroller firmware (RNGLabSensors.ino)
│   ├── scripts/           # Real-time hardware serial acquisition daemons
│   └── docs/              # Circuit schematics and sensor inventory
├── docs/
│   ├── architecture/      # Standards traceability and regulatory mappings
│   ├── archive/           # Historical kickoff logs, capability inventories, corrections
│   ├── benchmarks/        # CPU vs. NVIDIA L4 GPU performance and precision audits
│   ├── calibrate/         # V017 calibration protocols and review documentation
│   ├── finance/           # Quantitative financial model specifications
│   ├── physical/          # Sensor guides, pipeline design, and muon detector assessment
│   ├── trident/           # V016 Ising studies, cross-model analyses, mathematical theory
│   ├── validation/        # Methodological reviews, critiques, and source comparisons
│   └── REPRODUCIBILITY.md # Verification protocol and provenance manifest
├── scripts/               # Benchmark runners, audit verifiers, report generators
├── tests/                 # 238 unit, mathematical, and regression tests (100% passing)
├── pyproject.toml         # Package build configuration
└── LICENSE                # GNU Affero General Public License v3 (AGPL-3.0)
```

---

## 6. Getting Started & Verification

### Prerequisites
- Linux x86_64 (recommended), macOS, or Windows WSL2.
- Python 3.11+ (tested on Python 3.11, 3.12, 3.14).

### Environment Setup

```bash
git clone https://github.com/robinator121/rngLABexperimental.git
cd rngLABexperimental

python3 -m venv .venv
source .venv/bin/activate
pip install numpy scipy cryptography matplotlib pytest
```

*(Optional dependencies for physical hardware capture: `pip install -r requirements-physical.lock`.)*

### Running the Complete Test Suite

Execute the full suite of **238 automated tests** (covering core statistical batteries, quantum/Ising models, physical acquisition, and financial risk engines):

```bash
PYTHONPATH=src pytest tests -q
```
*Expected result: 238 passed in ~9 seconds.*

### Executing Benchmark Suites

```bash
# Run TRIDENT cross-model and Ising simulation campaign
PYTHONPATH=src python -m rnglab16.campaign runs/trident_run/

# Run CALIBRATE forensic audit and failure predictors
PYTHONPATH=src python scripts/calibrate_audit.py runs/trident_run/

# Inspect digitized muon detector waveform export
python3 scripts/inspect_detector_export.py /path/to/meas_data.txt
```

---

## 7. Citation & Academic Contact

If you use this platform or its findings in academic research, please cite:

```bibtex
@misc{hevele2026rnglab,
  author       = {Hevele, Imre R{\'o}bert},
  title        = {{RNGLab: Computational Randomness, Physical Entropy Acquisition, and Stochastic Simulation Integrity}},
  year         = {2026},
  publisher    = {GitHub},
  howpublished = {\url{https://github.com/robinator121/rngLABexperimental}}
}
```

**Contact:** Hevele Imre Róbert ([roberthevele@gmail.com](mailto:roberthevele@gmail.com)), University of Debrecen.
