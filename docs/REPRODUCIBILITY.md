# RNGLab v0.14 — Independent Scientific Reproducibility Protocol

**Document Identifier:** RNGLAB-DOC-v014-REPRODUCIBILITY-01  
**Classification:** Verification Protocol & Provenance Record  
**Target Release:** RNGLab v0.14 (`feat/v014-scientific-integration`)  
**Date:** 2026-10-18  
**Provenance SHA-256:** `570aa3e7a1e432ba979454c75cf0c90859175810f75678b4c529a023f6f83a6e`

---

## 1. System Environment & Prerequisites

To ensure bit-for-bit mathematical reproducibility across environments, the exact execution parameters must be mirrored:

- **Host Operating System:** Linux x86_64
- **Python Runtime:** Python 3.11+ / 3.14 (virtual environment `.venv`)
- **Required Scientific Libraries:**
  - `numpy >= 2.0.0`
  - `scipy >= 1.13.0`
  - `cryptography >= 42.0.0`

---

## 2. Worktree & Repository Setup

The platform is managed within the Git repository:

```bash
# Clone or navigate to the repository
cd rnglab-repo

# Verify git branch
git branch --show-current
```

---

## 3. Cryptographic Provenance Manifest

Every critical source code module and experiment artifact is cryptographically hashed with SHA-256:

| Component / Artifact | Relative File Path | SHA-256 Cryptographic Hash |
| :--- | :--- | :--- |
| Abstract Source & Channel Base | `src/rnglab14/sources/base.py` | `5b9ee08e57fe7b61c6234ba1a008cd9cb60411d30dbee899a717e13058025d78` |
| Reference PRNGs & ChaCha | `src/rnglab14/sources/reference.py` | `48e6c42065adaa7f59516ff4a2a8428d884b09111ce6e443bb070173e70dc73b` |
| Physical Noise Emulators | `src/rnglab14/sources/physical.py` | `bf26bd5b3c9357f28ae03064525c1e71eb50359a8a242c9a9a2d7b1ab70b1c51` |
| Conditioning & Reseeding DRBG | `src/rnglab14/sources/conditioning.py` | `a1d3c7fdff33add7c477090fc4773718fcb806fc2f90531cb3141a268a345196` |
| Fault & Degraded Generators | `src/rnglab14/sources/faults.py` | `b2dbf0d8ee9dc7376a2bc309196c0a2228c053ed7b9d1c5b552579acc6df432f` |
| Classical MC $\pi$ Estimator | `src/rnglab14/models/pi_mc.py` | `a1c675ab3f0f5ee6a144448d5abe171b58b051a80f92c67b2c2879fb91434694` |
| 2D Ising Metropolis Model | `src/rnglab14/models/ising_2d.py` | `9608e644c8187514f951840f312cde7a44519060db449d924b6e965004412e7e` |
| VMC Harmonic Oscillator | `src/rnglab14/models/vmc_harmonic.py` | `62943daccf57681f81cbea464bb2ca6840025d172c9fffd1ebd63688ce9e62cc` |
| PIMC Centroid Virial Model | `src/rnglab14/models/pimc_harmonic.py` | `b8026dc566b1642a6e6fbf375a609757f76dc4702809c835b7a50f66688ff81b` |
| Decision-Support System | `src/rnglab14/models/decision_support.py` | `b0c48fa0b96fde9d056ef66d8eb76a7b510686701a6dbb2eadf7fdaa5e108000` |
| Censored RUS Predictors | `src/rnglab14/forecast/predictor.py` | `fcbfe5f7611a2ca390e643fc6c6bb79812c160cadb46327772f10519d5834cf5` |
| Dual-Level Integrity Monitor | `src/rnglab14/monitoring/dual_monitor.py` | `2a398f69daf6cdeb2fd64ff8fb0e4a0152d667abed704251e4cc700ada143538` |
| Multi-Testing Error Controls | `src/rnglab14/statistics/multitest.py` | `ae69691d00f50763c4d661b64802353d7a7d420df5a54a19a97507d88ac83f2e` |
| Regulatory Standards Catalog | `src/rnglab14/standards/traceability.py`| `032b3bb363bafa069e92a0df7a0c90f7faeb9b9df6fcfd075b24648be638dc0f` |
| Integrated Experiment Script | `experiments/run_v014_integrated_experiment.py` | `0bb10a1264acdd6db508585a73d479cd6d903e6bd5867c99a059ef19786f44f5` |
| Integrated Experiment Results | `runs/v014-integrated-experiment/experiment_results.json` | `6d2d6775cd8ea7f0d517fe34f47710c7a73b5d17bbad01a11f4bf9853ef1da2c` |

---

## 4. One-Line Reproduction Commands

### 4.1 Running the Complete Test Suite
To verify unit correctness, boundary safety, and numerical regression tests:

```bash
PYTHONPATH=src python -m pytest tests -q
```
**Expected Outcome:** 120+ tests passed in under 5.0 seconds.

### 4.2 Re-running the Integrated Scientific Experiment
To regenerate the full vertical slice experiment results:

```bash
PYTHONPATH=src python experiments/run_v014_integrated_experiment.py
```
**Expected Outcome:**
1. Evaluates all 6 sources across MC $\pi$, 2D Ising, VMC Harmonic Oscillator, PIMC Centroid Virial, and Decision Support.
2. Evaluates the 100-snapshot dual-level monitoring matrix (Q1-Q4).
3. Executes the corrected RUS censored forecasting benchmark.
4. Performs Holm-Bonferroni and Benjamini-Hochberg multi-testing adjustments.
5. Emits `runs/v014-integrated-experiment/experiment_results.json` and `experiment_manifest.json`.

---

## 5. Verification Checklist

- [x] Python virtual environment confirmed.
- [x] No modifications to historical baselines (`rnglab-v0.11`, `rnglab-v0.12-pilot`).
- [x] Background research pipelines undisturbed.
- [x] SHA-256 hashes verified for all codebase files.
- [x] Test suite 100% green.
- [x] Integrated experiment execution verified end-to-end.
