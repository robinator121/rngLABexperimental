# RNGLab: Computational Randomness, Physical Entropy Acquisition & Stochastic Simulation Integrity

**Author:** Hevele Imre Róbert ([@robinator121](https://github.com/robinator121))  
**Contact:** [roberthevele@gmail.com](mailto:roberthevele@gmail.com)  
**Affiliation:** Faculty of Informatics, University of Debrecen  
**Project Nature:** Academic Research Prototype & BSc Thesis Codebase  
**Status:** Experimental Research Platform (v0.17 CALIBRATE)

---

## 1. Executive Summary

**RNGLab** is an experimental scientific research platform developed at the University of Debrecen to investigate how subtle defects in pseudorandom and physical entropy streams propagate through stochastic simulations, and how such degradations can be detected before corrupting physical or financial observables.

The project bridges three complementary domains:
1. **Physical Entropy Acquisition:** Real-time hardware noise harvesting (avalanche breakdown in reverse-biased junctions, Geiger-Müller radiation pulse intervals, and cosmic ray / muon detector digitized waveforms) coupled with NIST SP 800-90B min-entropy evaluation.
2. **Computational Physics & Stochastic Modeling:** Rigorous Markov Chain Monte Carlo (2D Ising Metropolis and Wolff cluster), Path Integral Monte Carlo (PIMC virial estimators), and Variational Monte Carlo (VMC harmonic oscillator) simulations evaluating observable masking and spectral gap collapse.
3. **Financial Risk & Model Risk Stress-Testing:** Quantitative financial engines (European options, Portfolio Value at Risk, and Expected Shortfall) designed to stress-test stochastic simulation stability under quantitative model validation principles (inspired by Federal Reserve SR 11-7, ECB TRIM, and Basel III guidelines).

> [!IMPORTANT]
> **Academic & Experimental Scope:** RNGLab is an exploratory academic prototype, not a commercial turnkey product. It implements methodological frameworks inspired by regulatory model risk guidance, but is **not officially certified by NIST, the Federal Reserve, the ECB, or other regulatory authorities**.

> [!WARNING]
> **Production Cryptography Disclaimer:** The conditioning routines implemented in this repository (von Neumann extraction, SHA-256 block hashing, ChaCha20 expansion) are research treatments designed for stochastic simulation benchmarking. They are **not certified or audited for production cryptographic systems or cryptographic key generation**. Passing statistical batteries (such as NIST SP 800-22) after cryptographic hashing demonstrates only the PRF diffusion of the hash function; true source unpredictability must be assessed directly on raw, unconditioned noise via NIST SP 800-90B min-entropy estimation.

---

## 2. Core Scientific Capabilities

### 2.1 Physical Entropy Harvesting & Muon Detector Analysis (`rnglab15`, `hardware/`, `docs/physical/`)
- **Microcontroller Hardware Interface:** Custom Arduino C++ firmware ([`RNGLabSensors.ino`](hardware/arduino/RNGLabSensors.ino)) interfacing with hardware noise diodes, optical PN junctions, and pulse-counting sensors via high-speed serial streaming.
- **Cosmic Ray / Muon Detector Dataset Evaluation:** Full exploratory analysis of 606.6 MB of digitized detector waveforms (117,962 records, 120.7 million samples). Evaluates pulse shapes, ceiling clipping, timing gaps, baseline recovery, and raw physical entropy extraction feasibility ([`MUON_DETECTOR_ASSESSMENT.md`](docs/physical/MUON_DETECTOR_ASSESSMENT.md) and [`inspect_detector_export.py`](scripts/inspect_detector_export.py)).
- **NIST SP 800-90B Entropy & Health Testing:**
  - Online continuous health monitoring: **Repetition Count Test (RCT)** and **Adaptive Proportion Test (APT)** with exact binomial tail cutoffs.
  - Raw unconditioned min-entropy ($H_\infty$) estimation, collision tests, and Markov dependency bounds (`src/rnglab13/standards/entropy_90b.py`, `src/rnglab15/assessment.py`).

### 2.2 Computational Physics & Markov Chain Dynamics (`rnglab14`, `rnglab16` TRIDENT)
- **2D Ising Model:** Metropolis-Hastings local spin-flip dynamics alongside Wolff single-cluster algorithms on periodic $L \times L$ lattices ($L \in \{2, 3, 4, 8\}$). Validated against exact finite-size analytical enumeration states.
- **Variational Monte Carlo (VMC):** Quantum harmonic oscillator ground state simulation. Formally demonstrates **Observable Masking**—where global energetic observables remain invariant despite complete collapse in coordinate sampling.
- **Path Integral Monte Carlo (PIMC):** Finite-temperature continuous-space quantum virial path integral simulations ($M$ time beads), evaluating observable masking and interval undercoverage.
- **Augmented Markov Chain Sensitivity Theory:** Analytical two-state and augmented-state transition kernels evaluating spectral gap degradation under directional acceptance bias.

### 2.3 Financial Risk & Model Risk Stress-Testing (`src/rnglab_finance/`, `docs/finance/`)
- **European Option Pricing:** Risk-neutral geometric Brownian motion simulation anchored against exact closed-form Black-Scholes analytical solutions. Enforces the risk-neutral martingale invariant ($\mathbb{E}^{\mathbb{Q}}[e^{-rT}S_T] = S_0$) and Put-Call Parity.
- **Portfolio Value at Risk (VaR) & Expected Shortfall (ES):** Multi-asset correlated portfolio risk engines computing parametric vs. Monte Carlo quantiles and coherent tail expectations.
- **Fault Injection & Model Risk Stress-Testing:** Controlled injection of non-catastrophic PRNG faults (mean shift, bit truncation, Markovian autocorrelation, and periodic leakage) to evaluate risk capital bias under Fed SR 11-7 and ECB TRIM stress-testing principles.

### 2.4 Reliability Calibration & Benchmarks (`rnglab17` CALIBRATE, `docs/benchmarks/`)
- **Calibrated Operational Gating:** Holm-Bonferroni step-down and Benjamini-Hochberg FDR multi-testing controls over statistical alarm families.
- **Hardware & Precision Audits:** Detailed execution and precision benchmarks comparing multi-core CPU workers against NVIDIA L4 GPU acceleration across FP32, FP64, and BF16 execution modes ([`RNGLAB_CPU_GPU_COMPARISON.md`](docs/benchmarks/RNGLAB_CPU_GPU_COMPARISON.md)).

---

## 3. Repository Architecture

```text
rngLABexperimental/
├── src/
│   ├── rnglab11/          # Core statistical evaluation, DRBG pipelines, and CLI
│   ├── rnglab12/          # Experimental review and uncertainty estimation
│   ├── rnglab13/          # Impossibility theorems and multi-stream coordinator
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
│   ├── benchmarks/        # CPU vs. NVIDIA L4 GPU performance and precision audits
│   ├── calibrate/         # V017 calibration protocols and review documentation
│   ├── finance/           # Quantitative financial model specifications
│   ├── physical/          # Sensor guides, pipeline design, and muon detector assessment
│   ├── trident/           # V016 Ising studies, cross-model analyses, mathematical theory
│   ├── validation/        # Methodological reviews, critiques, and source comparisons
│   └── CODE_METRICS.md    # Syntactic line count audit (27,500+ lines of code)
├── scripts/               # Benchmark runners, audit verifiers, report generators
├── tests/                 # Comprehensive unit, mathematical, and regression test suite (237 tests)
├── pyproject.toml         # Package build configuration
└── LICENSE                # GNU Affero General Public License v3 (AGPL-3.0)
```

---

## 4. Codebase Metrics

A verified syntactic line count audit (excluding `.git`, `.venv`, and `__pycache__`) is maintained in [`docs/CODE_METRICS.md`](docs/CODE_METRICS.md):

| Language | Files | Code Lines | Comment Lines | Blank Lines | Total Lines |
| :--- | ---:| ---:| ---:| ---:| ---:|
| **Python 3** | 219 | 24,171 | 555 | 3,237 | 27,963 |
| **Arduino C++** | 1 | 166 | 9 | 54 | 229 |
| **Markdown Documentation** | 43 | 2,156 | 209 | 593 | 2,958 |
| **Configurations & Data** | 9 | 1,030 | 3 | 17 | 1,050 |
| **Total** | **272** | **27,523** | **776** | **3,901** | **32,200** |

---

## 5. Getting Started & Verification

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

Execute the full suite of **237 automated tests** (covering core statistical batteries, quantum/Ising models, physical acquisition, and financial risk engines):

```bash
PYTHONPATH=src pytest tests -q
```
*Expected result: 237 passed in ~9 seconds.*

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

## 6. Citation & Academic Contact

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
