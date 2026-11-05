# RNGLab: random-input faults and uncertainty failure in Monte Carlo simulation

[![CI](https://github.com/robinator121/rngLABexperimental/actions/workflows/ci.yml/badge.svg)](https://github.com/robinator121/rngLABexperimental/actions/workflows/ci.yml)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPL_v3-blue.svg)](LICENSE)

**Author:** Hevele Imre Róbert ([@robinator121](https://github.com/robinator121)), Faculty of Informatics, University of Debrecen
**Supervisor:** Dr. Herendi Tamás
**Status:** research prototype, v0.18 in progress on branch `feature/v018`. Not certified or validated for any production, cryptographic or regulatory use.

## Research question

How do defects in random inputs propagate through Monte Carlo simulations, and how can a sampler's own uncertainty failure be told apart from a fault in the random source?

The strongest evidence so far comes from one exact-reference workload: harmonic path-integral Monte Carlo (PIMC) with injected random-channel faults. Every number below names the file it comes from. [`docs/v018/AUDIT.md`](docs/v018/AUDIT.md) lists every quantitative claim in this repository as backed, partially backed or unbacked.

## Main result: Phase II, harmonic PIMC (October 2026)

Standalone harmonic PIMC at M = 16 beads, β = 2, proposal half-width δ = 0.4, exact finite-M Gaussian starts, no burn-in and 8,192 retained sweeps; nominal 95% interval on the centroid-virial (CV) energy from batch means with 16 batches of 512 (BM512) and a Student-t₁₅ critical value; 1,024 independent replicas per condition.

| Sampler (healthy random input) | Intervals covering the exact finite-M energy | Exact binomial 95% interval |
|---|---:|---:|
| Local Metropolis | 811 / 1,024 | 76.6–81.6% |
| Local Metropolis + exact centroid refresh | 964 / 1,024 | 92.5–95.5% |
| Independent Gaussian paths (control) | 969 / 1,024 | 93.1–95.9% |

Source: `research/phase2-20261030/experiment/confirmation-v1/summary.csv` (rows `CV,BM512`) and `primary_decision.json`.

- **Healthy sampler failure.** With healthy random numbers, local Metropolis under-covers. The failure persists from exact stationary starts, so it is not initialization bias. The slow centroid mode dominates the variance of the energy average and short batches miss it (estimated/empirical variance ratio 0.53 for local Metropolis, `primary_decision.json`). Refreshing the centroid from its exact conditional law restores coverage within the frozen gates for this workload.
- **Injected channel faults.** On the validated refreshed sampler, a paired 512-replica study covered the exact energy in 487 / 512 healthy runs, 0 / 512 under an acceptance-uniform shift and 13 / 512 under proposal/acceptance coupling (`research/phase2-20261030/experiment/fault-study-v1/summary.csv`, post-onset window). A marginal source-mean monitor missed the coupling, a joint-channel monitor detected it, and the application-shape monitor missed the acceptance shift.
- **Scope.** One harmonic workload, one interval method, two strong synthetic adapter faults. These are injected faults, not defects found in PCG64. No general monitor ranking follows. The coverage campaign used the standalone implementation preserved in `research/phase2-20261030`; validation through the package API is v0.18 Phase 2.
- **Known versus new.** That short batch means underestimate long-run variance is established (Flegal & Jones, [Ann. Statist. 38, 1034 (2010)](https://doi.org/10.1214/09-AOS735); Geyer, [Statist. Sci. 7, 473 (1992)](https://doi.org/10.1214/ss/1177011137)). Virial energy estimators for PIMC are established (Yamamoto, [J. Chem. Phys. 123, 104101 (2005)](https://arxiv.org/abs/physics/0505109)), and the centroid refresh samples the known exact conditional Gaussian law of the harmonic centroid; neither is new here. That generators can pass statistical tests yet bias an application is established (Ferrenberg, Landau & Wong, [Phys. Rev. Lett. 69, 3382 (1992)](https://doi.org/10.1103/PhysRevLett.69.3382); Coddington, [Int. J. Mod. Phys. C 7, 295 (1996)](https://doi.org/10.1142/S0129183196000235)). The Phase II contribution is a controlled separation, in one workload, of a healthy-sampler uncertainty failure from two random-channel contract violations, with paired effects and separately calibrated monitors.

Full reports, protocols and reproduction: [`research/phase2-20261030/README.md`](research/phase2-20261030/README.md).

## Negative and limited results

These carry the same weight as the positive results above.

- **Healthy PIMC undercoverage** (above): 811 / 1,024 for local Metropolis. Earlier short-block analyses were worse: 172 / 512 (BM32) and 311 / 512 (BM128) in v0.17 ([`docs/calibrate/V017_VMC_PIMC_VALIDATION.md`](docs/calibrate/V017_VMC_PIMC_VALIDATION.md); output off-repo, see the audit).
- **v0.17 CALIBRATE hypotheses:** H1 (healthy whole-path false-alarm ceiling) failed jointly, with 10 / 12 cells passing; the full H4 usefulness criterion and H5 (unseen-fault generalization) are unsupported ([`docs/calibrate/V017_NEGATIVE_FINDINGS.md`](docs/calibrate/V017_NEGATIVE_FINDINGS.md); outputs off-repo).
- **v0.16 TRIDENT forecasting:** theory features gave no Brier improvement; the forecasting endpoint confounded thermal fluctuation with reliability loss ([`docs/trident/V016_NEGATIVE_FINDINGS.md`](docs/trident/V016_NEGATIVE_FINDINGS.md)).
- **Crypto BLUE:** the frozen ≤ 10% healthy-alarm bound failed (8 / 128 healthy alarms, exact 95% interval 2.7–11.9%), and SHA256/ChaCha20 missed the preregistered 90% non-rejection criterion ([`research/crypto-blue/confirmation-reviewed/RESULTS.md`](research/crypto-blue/confirmation-reviewed/RESULTS.md)).
- **RUS forecasting:** a constructive counterexample shows two degradation processes with identical observed history and different futures. Unrestricted remaining-useful-sample forecasting is therefore not identifiable, and no successful general forecaster is claimed ([`research/phase2-20261030/rus/`](research/phase2-20261030/rus/), [`src/rnglab13/math/impossibility.py`](src/rnglab13/math/impossibility.py)).

## Supporting findings

**1. Acceptance-shift faults make the Ising chain absorbing (analytic).** In 2D Ising Metropolis with a clipped acceptance-uniform shift δ, an uphill move with energy cost ΔE is accepted with probability max(e^(−ΔE/T) − δ, 0). From a fully aligned state every single flip costs ΔE = 8, so for δ ≥ e^(−8/T) both magnetized states become absorbing and the chain loses ergodicity. This is structural damage, not a small stationary perturbation, and no stationary-response formula applies there. The v0.16 numerical illustration is in [`docs/trident/V016_ISING_STUDY.md`](docs/trident/V016_ISING_STUDY.md); its outputs are not in Git.

**2. A healthy chain gives a badly wrong susceptibility while its energy is right.** At L = 4, T = 1.5 the exact zero-field susceptibility is 10.403227. Across 512 healthy chains the within-chain plug-in estimate (variance of the signed magnetization) averaged 0.745828, because finite chains rarely visit both magnetized phases. A symmetry-aware estimator, N·mean(m²)/T, which uses E[m] = 0, averaged 10.403680 ([`docs/calibrate/V017_NEGATIVE_FINDINGS.md`](docs/calibrate/V017_NEGATIVE_FINDINGS.md)). This is an estimator and mixing failure with healthy random numbers, not an RNG defect. *Partially backed:* the script is `scripts/calibrate_susceptibility_probe.py`, and the numbers equal its output `SUSCEPTIBILITY_PROBE.json`, which is not yet in Git. v0.18 Phase 1 reproduces it as a regression test and validates both estimators on held-out data.

**3. Detector data must be inspected before any entropy claim.** The author's 606,632,259-byte detector export (University of Debrecen, 19 October 2026) contains 117,962 records of 1,023 waveform codes each. 108,473 records (91.96%) reach the 16-bit ceiling 65,535, and every timestamp is a multiple of 80, so its 4 lowest bits are always zero. Re-derived from the raw file: [`docs/v018/evidence/detector_inventory.json`](docs/v018/evidence/detector_inventory.json) (`scripts/inspect_detector_export.py`) and [`detector_ceiling.json`](docs/v018/evidence/detector_ceiling.json) (`scripts/v018/detector_ceiling_check.py`). Particle identity and entropy per sample are not established ([`docs/physical/MUON_DETECTOR_ASSESSMENT.md`](docs/physical/MUON_DETECTOR_ASSESSMENT.md)). The raw export is not in Git.

**4. Statistically plausible output can be fully predictable (crypto BLUE).** With a known 8-bit seed passed through a conditioner and an expander, the outputs usually pass three statistical diagnostics, and an attacker who knows the seed predicts all of them exactly ([`research/crypto-blue/confirmation-reviewed/RESULTS.md`](research/crypto-blue/confirmation-reviewed/RESULTS.md); 128 replicas, 256 operations each):

| Seed → conditioner → expander | Not rejected by monobit, block-frequency, runs (Holm, 1%) | Exactly predicted |
|---|---:|---:|
| SHA256 → HMAC-DRBG | 128 / 128 | 128 / 128 |
| SHA256 → ChaCha20 | 121 / 128 | 128 / 128 |
| HKDF → HMAC-DRBG | 128 / 128 | 128 / 128 |
| HKDF → ChaCha20 | 127 / 128 | 128 / 128 |

These are three diagnostics, not the full SP 800-22 battery, on synthetic data only. More results: [`docs/crypto/RESULTS.md`](docs/crypto/RESULTS.md).

## Demonstrations (known results, reproduced)

- **VMC zero-variance property.** At the exact harmonic trial function (α = 0.5) the local energy is constant, so a proposal-biased chain still reports energy 0.5 with zero variance while ⟨x²⟩ is wrong by 1.656. At α = 0.6 the same fault moves the energy by −0.248 (σ = 0.296). One deterministic run of the configuration in `tests/test_models14.py`, recorded by `scripts/v018/vmc_masking_values.py` in [`docs/v018/evidence/vmc_masking_values.json`](docs/v018/evidence/vmc_masking_values.json). This is the standard zero-variance property of variational Monte Carlo (Foulkes et al., [Rev. Mod. Phys. 73, 33 (2001)](https://doi.org/10.1103/RevModPhys.73.33); Assaraf & Caffarel, [Phys. Rev. Lett. 83, 4682 (1999)](https://doi.org/10.1103/PhysRevLett.83.4682)), shown here as a monitoring pitfall, not as a new hypothesis.
- **Reused seeds in parallel workers.** `rnglab_finance.faults.injection.ReusedSeedProvider` makes K workers produce identical paths (`tests/test_finance_faults.py`). The reported standard error is then too small by about √K, because only n/K paths are distinct. That magnitude is an analytic consequence and has not been measured in a saved experiment.

## Superseded, unbacked or not reproducible from this repository

Retained for traceability; do not cite. Details and evidence in [`docs/v018/AUDIT.md`](docs/v018/AUDIT.md).

- [`docs/validation/SOURCE_COMPARISON.md`](docs/validation/SOURCE_COMPARISON.md): the benchmark matrix has no producing script. Its throughput figures are 29–156 times what the repository code achieves for four of five sources (measured in [`docs/v018/evidence/source_throughput.json`](docs/v018/evidence/source_throughput.json)), and it reports "all 15" SP 800-22 tests although the code implements three.
- v0.14 "0.0% hazardous gap" and the other v0.14 registry results ([`docs/validation/EXPERIMENT_REGISTRY.md`](docs/validation/EXPERIMENT_REGISTRY.md)): one uncalibrated run, superseded by v0.14.1 and contradicted by Phase II's monitor blind spots. `APPLICATION_RELIABILITY.md`, `REVIEW_METHODOLOGY_CRITIQUE.md` and `docs/REPRODUCIBILITY.md` (stale file hashes) are also superseded although they carry no banner.
- CPU/GPU benchmarks ([`docs/benchmarks/`](docs/benchmarks/)): the author's measurements, but no benchmark script or raw timing is in Git.
- Earlier README statements that are withdrawn: regulatory framing (SR 11-7, ECB TRIM, Basel III), which no analysis used; "formal mathematical proof" for the RUS result, which is a constructive counterexample; physical SP 800-90B entropy evaluation, never run on a physical stream; exact-enumeration validation for L = 8, which exists only for L ≤ 4; and the previous VMC numbers, which do not reproduce.

## Repository layout

```text
src/rnglab/            v0.18 namespace (sources, faults, models, references, monitors, stats)
                       beside the legacy v0.1–v0.8 subpackages; see docs/v018/ARCHITECTURE.md
src/rnglab10 … 17/     historical snapshots, unchanged
src/rnglab_crypto/     crypto BLUE
src/rnglab_finance/    finance proof of concept
research/              phase1-20261030, phase2-20261030, crypto-blue (frozen study packages)
docs/v018/             v0.18 audit, architecture, evidence
docs/                  earlier version documentation (see the audit for status)
scripts/, scripts/v018/  campaign and audit scripts
hardware/              Arduino firmware and acquisition helpers
tests/                 automated test suite
```

## Running the tests

```bash
python3 -m venv .venv
.venv/bin/pip install -r research/phase2-20261030/environment/requirements-lock.txt
PYTHONPATH=src .venv/bin/pytest tests -q
```

Expected: **331 passed**, which is the 322 tests on `main` (257 core and 65 crypto BLUE) plus 9 v0.18 namespace tests ([`docs/v018/evidence/pytest_phase0.txt`](docs/v018/evidence/pytest_phase0.txt), Python 3.14.9). CI is configured for Python 3.11 and 3.12 but has never executed: GitHub has refused to start every job since the workflow was added on 29 October 2026 because of an account billing lock, so the CI badge shows failure without any test having run. The workflow's steps pass when run locally with its unpinned dependencies on Python 3.14.9 ([`docs/v018/evidence/ci_steps_local.txt`](docs/v018/evidence/ci_steps_local.txt)); Python 3.11 and 3.12 are untested. The 35 Phase II research checks run separately ([`research/phase2-20261030/REPRODUCE.md`](research/phase2-20261030/REPRODUCE.md)).

## How this work was produced (AI disclosure)

- **Author:** research question, design decisions and their approval, and the physical detector acquisition and Arduino firmware. *[Author to complete: which analyses, derivations and reviews you did yourself.]*
- **Codex (OpenAI):** Phase II computations, code, checks and report drafting (stated in [`research/phase2-20261030/README.md`](research/phase2-20261030/README.md)). The v0.16 TRIDENT, v0.17 CALIBRATE and crypto BLUE delivery folders on the author's machine are Codex project outputs. *[Author to confirm the extent, and which tools produced v0.10–v0.15.]* Review documents marked "drafted with AI-assisted code analysis tools" were AI-drafted.
- **Claude Code (Anthropic):** the v0.18 work on `feature/v018`, including the audit (`docs/v018/AUDIT.md`), the `src/rnglab` namespace modules, `scripts/v018/`, `docs/v018/` and this README.
- The test suite was not verified by hand. It is executed automatically, and its counts are recorded in the evidence files above. Independent review by a supervisor or another researcher is still required before thesis, TDK or publication use.

## Citation

```bibtex
@misc{hevele2026rnglab,
  author       = {Hevele, Imre R{\'o}bert},
  title        = {{RNGLab: random-input faults and uncertainty failure in Monte Carlo simulation}},
  year         = {2026},
  publisher    = {GitHub},
  howpublished = {\url{https://github.com/robinator121/rngLABexperimental}}
}
```

**Contact:** Hevele Imre Róbert ([roberthevele@gmail.com](mailto:roberthevele@gmail.com)), University of Debrecen.
