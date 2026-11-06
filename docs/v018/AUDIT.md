# v0.18 audit of quantitative claims in `README.md` and `docs/`

Audit date: 2026-11-05. Audited revision: `main` at `697c9436beb9f04c45365c4dae78eb54a90bbc12` (cloned fresh from GitHub). Audit produced by Claude Code (Anthropic, model Opus 5.5) at the request of the author; status judgements are the assistant's and are open to the author's review. No historical file was modified. Corrections belong in new documents (this one, the rewritten README), not in the audited files.

## Status definitions

| Status | Meaning |
|---|---|
| **backed** | A script in this repository produces an output in this repository that contains the number, and the output supports the sentence. |
| **partially backed** | One link is missing: the script exists but its output lives only outside Git (named below), or the output exists but the sentence says more than it shows, or only some numbers in a group trace. |
| **unbacked** | No script and no saved output produce the number, a re-run gives a different number, or the claim is of a kind no output can carry (e.g. "verified by hand"). |

**Superseded** is a separate annotation: the repository already marks the document as known-incorrect (`> SUPERSEDED v014 report` banner) or a later document corrects it.

## How the statuses were reached

1. Every document was read in full. Claims were listed one by one for `README.md` and grouped by table or paragraph for `docs/`.
2. `scripts/v018/audit_numbers.py` extracted every numeric token (3,008 tokens; 981 "discriminating", meaning precise enough that a chance match is unlikely) and searched for each in the text outputs of the campaign the document describes, at the quoted precision. Off-repo directories were searched separately. As a **decoy control**, every lookup was repeated with the value multiplied by 1.0371. When the decoy hit rate approaches the real hit rate, the matcher cannot discriminate and the claim was checked by hand instead. Output: `docs/v018/evidence/audit_numbers.csv` and `audit_numbers.summary.json`. A match only shows that the number occurs in an output of the right campaign; the status column is a judgement on top of that.
3. Five claims were re-derived from raw inputs or re-run with new scripts:

| Check | Script | Output |
|---|---|---|
| Test-suite baseline on unmodified `main` | `PYTHONPATH=src pytest tests -q` | `docs/v018/evidence/pytest_baseline.txt` (322 passed) |
| Detector inventory from the raw 606,632,259-byte export | `scripts/inspect_detector_export.py` (existing) | `docs/v018/evidence/detector_inventory.json` |
| Detector per-record ceiling and trigger counts | `scripts/v018/detector_ceiling_check.py` (new) | `docs/v018/evidence/detector_ceiling.json` |
| Source throughput of the repository generators | `scripts/v018/source_throughput.py` (new) | `docs/v018/evidence/source_throughput.json` |
| README VMC masking numbers | `scripts/v018/vmc_masking_values.py` (new) | `docs/v018/evidence/vmc_masking_values.json` |

The raw detector export (`~/rnglab-data/meas_data_2026-10-19 14.05.19.txt`) is on the author's machine, not in Git. Its SHA-256 (`docs/v018/evidence/detector_export.sha256`) equals the one recorded in `docs/physical/MUON_DETECTOR_ASSESSMENT.md`.

Off-repo evidence roots used (local machine only, not reproducible from GitHub):

- `v016`: `~/.codex/.chatgpt-projects/g-p-6ac60405d13c819197865c58f0a52c71/output/rnglab-v016-trident-20261022/` (147 MB)
- `v017`: `~/.codex/.chatgpt-projects/g-p-6ac60405d13c819197865c58f0a52c71/output/rnglab-v017-calibrate-20261022/` (933 MB)
- `phase2`: `~/Documents/Codex/2026-10-30/do-x20/outputs/phase2-20261030/` (1.9 GB, the full Phase II outputs including the 19 NPZ archives)

## 1. `README.md`, claim by claim

| # | Claim (README location) | Status | Evidence / problem |
|---|---|---|---|
| R1 | Badge and §6: "322 automated tests", "257 legacy plus 65 crypto" | **backed** | `pytest_baseline.txt`: 322 passed on `main`. Collected per file: the five `test_crypto_*.py` files hold 21+19+10+8+7 = 65. |
| R2 | Header: "BSc Thesis Codebase", "Status: v0.17 CALIBRATE" | n/a (not quantitative) | Outdated framing: v0.18 is research, not thesis work. |
| R3 | BLUE: "128 confirmation and 64 sensitivity replicas" | **backed** | `research/crypto-blue/confirmation-reviewed/protocol.json`, `sensitivity-reviewed/`. |
| R4 | §1: hardware harvesting from avalanche diodes, Geiger–Müller intervals and cosmic-ray waveforms "coupled with NIST SP 800-90B min-entropy evaluation" | **unbacked** as a result | Only the detector export has data. GM acquisition is a plan (`docs/physical/GM_NEXT_WEEK.md`). No 90B assessment of any physical stream was run; the detector assessment explicitly assigns no entropy. Estimators exist in code (`rnglab13/standards/entropy_90b.py`), which is a capability, not a result. |
| R5 | §1, §3.3: "inspired by Federal Reserve SR 11-7, ECB TRIM, and Basel III" | **unbacked** framing | No analysis in the repository applies these documents. Remove (brief rule 6). |
| R6 | Finding 1: 117,962 events, 120.7M samples, 91.96% at ceiling 65,535, GCD = 80, timestamp LSBs identically zero | **backed** (re-derived) | `detector_inventory.json`: 117,962 records, 120,675,126 values, field-1 GCD 80, 627,833 ceiling values. `detector_ceiling.json`: 108,473 records reach 65,535 (fraction 0.91956). GCD 80 = 16·5, so the 4 lowest timestamp bits are always zero. Wording fix: "records", not "events" (the assessment says particle identity is not established). |
| R7 | Finding 2: VMC zero-variance; "Δ⟨x²⟩ = 1.75", "ΔE = −0.26", "σ_E = 0.31" | **unbacked** numbers; framing overstated | No script prints these. The cited test (`tests/test_models14.py`) only asserts bounds. Re-running its exact configuration gives Δ⟨x²⟩ = 1.656373, ΔE = −0.247578, σ_E = 0.296491 (`vmc_masking_values.json`). The phenomenon is the standard zero-variance property of VMC (local energy is constant at an exact eigenfunction; see the VMC review of Foulkes et al., Rev. Mod. Phys. 73, 33 (2001), and the general principle of Assaraf & Caffarel, Phys. Rev. Lett. 83, 4682 (1999)). It must be presented as a demonstration, not as a hypothesis the author formulated (see R16). |
| R8 | Finding 3: parallel seed collision deflates VaR SE "by a factor of up to √K", "masking severe risk underestimation under stress scenarios" | **partially backed** | `tests/test_finance_faults.py::test_reused_seed_path_duplication` shows the K workers produce identical paths. The √K factor is an analytic consequence (K copies of n/K paths), stated in the `ReusedSeedProvider` docstring; no saved experiment measures it. "Severe risk underestimation under stress scenarios": no output. |
| R9 | Finding 4: acceptance shift δ ≥ e^(−8/T) makes the magnetised states absorbing | **backed** (analytic) / **partially backed** (numbers) | One-line derivation: from a fully aligned state every flip costs ΔE = 8, accepted with probability max(e^(−8/T) − δ, 0) = 0. `scripts/trident_audit.py` records the check (`strong_bias_creates_multiple_absorbing_states`), but its output is off-repo (`v016`). "Perturbation approximations fail near T_c": the theory table is restricted to δ ≤ 0.025 at T_c (`docs/trident/V016_MATHEMATICAL_THEORY.md`), so the failure is shown only as domain exclusion. |
| R10 | Phase II: 811/1,024, 964/1,024, 969/1,024 (CV BM512, M=16, β=2, 8,192 sweeps) | **backed** | `research/phase2-20261030/experiment/confirmation-v1/summary.csv` rows `CV,BM512`; gates in `primary_decision.json`. |
| R11 | Phase II faults: 487/512, 0/512, 13/512 | **backed** | `research/phase2-20261030/experiment/fault-study-v1/summary.csv` rows `CV,BM512,postonset`. |
| R12 | "35 separate research regression checks" | **backed** (historical) | `research/phase2-20261030/environment/pytest-phase2-consolidated.xml`: tests=35, failures=0. Not re-run in this audit (needs the prepared replay tree). |
| R13 | §3.1 dataset: 606.6 MB, 1,023-sample waveforms, 14:05:19–14:27:23 | **backed** (size, layout) / **partially backed** (time span) | Size and 1,023 values per row in `detector_inventory.json`. The start time is the file name; the end time is not derivable from the file (timestamp units unknown). |
| R14 | §3.1: RCT and APT "with exact binomial tail cutoffs" | **backed** (implementation) | `rnglab14/standards/health_90b.py`, tests in `tests/test_standards14.py`. Capability, no field result. |
| R15 | §3.2: Ising L ∈ {2, 3, 4, 8} "validated against exact finite-size analytical enumeration" | **partially backed** | Enumeration exists only for L ≤ 4 (`rnglab16/core.py::enumeration` raises otherwise). `docs/trident/V016_ISING_STUDY.md`: "L8 has no exact finite-size reference here". |
| R16 | §3.4: "Formal mathematical proof and counterexample showing the impossibility of … RUS forecasting" | **unbacked** as "formal proof" | `src/rnglab13/math/impossibility.py` states the theorem in a docstring and constructs a counterexample (`tests/test_math13.py::test_impossibility_counterexample`); Phase II adds `research/phase2-20261030/rus/`. A constructive counterexample does establish non-identifiability for that construction, but no written proof of the general theorem as stated exists. Call it a constructive non-identifiability counterexample. |
| R17 | §3.4: CPU vs NVIDIA L4 benchmarks across FP32/FP64/BF16 | **unbacked** in the repository | `docs/benchmarks/*` cite measurements, but no benchmark script, raw timing or telemetry is in Git (matcher: 0 of 67 numbers in any repo output). The author's own measurement; it may be correct, but cannot be reproduced from the repository. |
| R18 | §4: "Formulated the research hypotheses: the zero-variance blind spot in VMC, PRNG parallel seed collisions … digitizer saturation" | **unbacked** framing | The VMC property is textbook (R7). Seed-collision variance deflation is standard Monte Carlo practice knowledge. Digitizer saturation is an observation about one dataset, not a hypothesis. |
| R19 | §4: "verified all 238 automated … tests by hand" | **unbacked**; inconsistent | The suite has 322 tests (R1). 238 was the count at `dd96adb` (Phase II README). The Phase II README states its code and checks were produced with Codex assistance. The v0.16, v0.17 and BLUE delivery folders exist on this machine under a Codex project output directory (paths above). Manual verification of every test is not something any output can record. |
| R20 | §5 tree: "rnglab11 … rnglab17" and "tests/ 257 legacy plus 65 crypto" | **backed** (counts) | Omits `rnglab10`, `rnglab_crypto`, legacy `src/rnglab/`, `research/crypto-blue`, `research/phase1-20261030`. |
| R21 | §6: "tested on Python 3.11, 3.12, 3.14" | **partially backed** | 3.14 is backed: Phase II environment manifest (3.14.7) and this audit's runs (3.14.9). 3.11 and 3.12 are **unbacked**: the CI matrix lists them (`.github/workflows/ci.yml`), but GitHub has never started a CI job in this repository. Every run since the workflow was added (e.g. runs 40931415718 on 29 Oct, 41381077576 on 2 Nov, 42010260891 on 5 Nov) is annotated "The job was not started because your account is locked due to a billing issue." Added 2026-11-05 after PR #2's CI failed; the original Phase 0 entry said only "partially backed". The CI steps pass locally on 3.14.9 (`docs/v018/evidence/ci_steps_local.txt`). |

## 2. The known suspects

### 2.1 `docs/validation/SOURCE_COMPARISON.md`, benchmark matrix — **unbacked**

- The 10-row matrix (bias, lag-1, min-entropy, throughput, SP 800-22 status) matches no saved output. Of the document's 38 discriminating numbers, 28 are unmatched (matcher). The 10 matches are the ChaCha8/ChaCha20 application numbers in §5, which come from a single run in `runs/v014-integrated-experiment/experiment_results.json` (base seed 101). "Replicated", "10⁹ samples" and "zero detectable degradation" have no output.
- **Throughput is implausible for this code.** Measured through the same `read_bytes` boundary every model uses (`source_throughput.json`, Python 3.14.9, this machine), in million uint64/s: `ChaCha20Source(rounds=20)` 0.095 per-call / 0.099 bulk, against 14.8 claimed; `rounds=8` 0.216 / 0.240, against 28.6; `PCG64Source` 0.620 / 1.646, against 48.2; `Xoshiro256PlusPlusSource` 0.741 / 1.247, against 52.1; `MersenneTwisterSource` 2.72 / 42.8, against 31.4. For the four non-MT sources the claimed figures are 29–156 times what this code achieves (ratios computed from the same JSON). The ChaCha implementation is pure Python. Only the Mersenne Twister bulk read (backed by NumPy) exceeds its claimed value.
- **"Pass (All 15)" SP 800-22 is impossible here.** `rnglab14/standards/sp800_22.py::run_sp800_22_battery` runs three tests (monobit, block frequency, runs).
- The document is not marked superseded, although sibling v0.14 reports are. Its "Provenance SHA-256" is the same value printed on every v0.14 document and identifies nothing specific.

### 2.2 v0.14 "0.0% hazardous gap" (`docs/validation/EXPERIMENT_REGISTRY.md` EXP-05) — **partially backed, superseded**

The counts (Q1 49, Q2 14, Q3 37, Q4 0 of 100 snapshots; mean lead time 14.0 steps) are in `experiment_results.json → dual_monitoring_report`. They describe one deterministic sequence of 100 snapshots from one run with an uncalibrated monitor. Zero gaps in one sequence is not a 0.0% rate; with 0 of 100 the exact one-sided 95% upper bound is 1 − 0.05^(1/100) = 0.0295, and that bound assumes independent snapshots, which these are not. `docs/archive/V0141_CORRECTIONS.md` already supersedes it ("Held-out combined monitoring was not universally superior"). Phase II then showed opposite blind spots: the application-shape monitor missed acceptance bias and the marginal monitor missed coupling (`research/phase2-20261030/reports/Controlled_RNG_fault_report.md`). Other EXP-0x issues:

- EXP-01 says "replicated over 10 independent seeds"; the results file has one run per source.
- EXP-03 reports E_therm = 3.8098 and E_vir = 0.5829 as a healthy comparison. 0.58292 is the ChaCha8 PIMC row, and the thermodynamic estimator used the incorrect spring coefficient corrected in v0.14.1. The "exact" 0.6565 is the continuum value; the finite-bead target differs.
- EXP-07's throughput: see 2.1.

### 2.3 "verified all 238 tests by hand" versus 322 tests and Codex-produced Phase II code — **unbacked**

See R19. The correction belongs in an accurate AI disclosure in the rewritten README.

### 2.4 "Formal mathematical proof" for the RUS counterexample — **unbacked as worded**

See R16.

### 2.5 VMC zero-variance presented as an original hypothesis — **unbacked framing**

See R7 and R18. The demonstration itself (energy exact and zero-variance at α = 0.5 while ⟨x²⟩ is wrong) is correct and reproducible; only its framing and its three numbers need fixing.

## 3. `docs/` by folder

| Folder / document | Status | Notes |
|---|---|---|
| `docs/calibrate/*` (v0.17 CALIBRATE) | **partially backed** | Scripts are in the repo (`scripts/calibrate_*.py`, `src/rnglab17/`). The outputs they cite (`calibration_analysis.json`, `forecast_analysis.json`, `MECHANISM_EVIDENCE.json`, `SUSCEPTIBILITY_PROBE.json`, traces) are only in the off-repo `v017` root. Of 262 discriminating numbers, 65 occur in repo outputs (mostly the frozen protocol), 196 only off-repo, and 1 is unmatched. Spot check by hand: the susceptibility counterexample (exact 10.403227, plug-in 0.745828, bias interval [−9.859323, −9.455475], symmetry estimator 10.403680) equals `SUSCEPTIBILITY_PROBE.json` row `ising_L4_T0`. These documents are careful about scope and already report their negative results (H1, H4, H5 failed). |
| `docs/calibrate/V017_SCIENTIFIC_PROTOCOL.md`, `FROZEN_PROTOCOL.json` | **backed** | Protocol parameters; the JSON is in the repo. |
| `docs/trident/*` (v0.16 TRIDENT) | **partially backed** | Scripts in repo (`scripts/trident_*.py`, `src/rnglab16/`); all outputs (`analysis.json`, `CROSS_MODEL_TABLE.json`, `FORECAST_UNIT_TABLE.json`, `FINAL_TRACE_CHECKS.json`) off-repo in `v016`. Matcher: 180 of 181 numbers occur only off-repo, with a low decoy rate (the Ising table: 112 of 113, decoy 3), so the numbers do come from those outputs. |
| `docs/crypto/*` (BLUE) | **backed** | Outputs in `research/crypto-blue/` (reviewed confirmation and sensitivity, evidence logs). 122 of 129 discriminating numbers matched in repo; the decoy rate is high for `RESULTS.md` (large results file), so the headline counts (128/128, 121/128, 769 colliding paths, 8/128 healthy alarms and the failed ≤10% calibration bound) were checked by hand against the reviewed `RESULTS.md`/`results.json`. `PHYSICAL_INTEGRATION.md`: 3 numbers unmatched (fixture-level statements). |
| `docs/physical/MUON_DETECTOR_ASSESSMENT.md` | **partially backed** | Inventory and §2 figures re-derived (§1 above, R6). Later sections (baseline 3,341.25, noise proxies, recovery, gap statistics, plots) have no script in Git: 90 of 121 numbers unmatched, decoy 0. They were computed in a session whose analysis code was not committed. |
| `docs/physical/PIPELINE.md`, `REPRODUCE.md`, `GM_NEXT_WEEK.md` | n/a | Method descriptions and plans. `REPRODUCE.md` points at `examples/physical/local-20261019.json`, which is not in the repository. |
| `docs/benchmarks/*` | **unbacked** in the repository | See R17. `RNGLAB_CPU_GPU_COMPARISON.md` is consistent with the author's statement that the GPU is slower at small batches; its raw timings should be committed or archived before it is cited. |
| `docs/validation/EXPERIMENT_REGISTRY.md`, `REVIEW_STATISTICAL_CORE.md` | **superseded** (banner present) | See 2.2. 15 of 46 registry numbers unmatched (e.g. "33.7%", "14.53", "24.25", "14.2%"). |
| `docs/validation/APPLICATION_RELIABILITY.md` | **partially backed; not marked superseded** | 19 of 22 numbers trace to the v0.14 results file, but several are misattributed: 0.5829 "on healthy reference sources" is the ChaCha8 run; E_therm ≈ 3.81 used the incorrect pre-v0.14.1 coefficient; "+9800%" divides the degraded virial energy 65.47 by the continuum value. "1.615 → 0.758" susceptibility and "VaR underestimated by up to 18%" are unmatched. "A central theorem demonstrated" is not a theorem. |
| `docs/validation/REVIEW_METHODOLOGY_CRITIQUE.md` | **unbacked; not marked superseded** | 10 of 14 numbers unmatched. Claims "120/120 tests" and "all 15 NIST SP 800-22 tests", neither of which the code supports (see 2.1). |
| `docs/validation/SOURCE_COMPARISON.md` | **unbacked; not marked superseded** | See 2.1. |
| `docs/architecture/STANDARDS_TRACEABILITY.md`, `docs/archive/SCIENTIFIC_CAPABILITY_INVENTORY.md` | **superseded** (banner present) | 18 of 19 numbers unmatched. |
| `docs/archive/RESEARCH_EVOLUTION_LOG.md` | **partially backed** | Historical narrative; 5 of 9 numbers unmatched ("12 / 34 / 120 tests" per version). |
| `docs/archive/V0141_CORRECTIONS.md` | **partially backed** | "0/9 before → 9/9 after" refers to evidence under `../../evidence/physics`, which is not in the repository. |
| `docs/REPRODUCIBILITY.md` (v0.14) | **unbacked; stale** | Its file-hash table no longer matches the code: 10 of 14 listed `src/rnglab14` files have changed hashes (stale: `sources/base.py`, `reference.py`, `physical.py`, `conditioning.py`, `faults.py`, `models/vmc_harmonic.py`, `pimc_harmonic.py`, `forecast/predictor.py`, `monitoring/dual_monitor.py`, `standards/traceability.py`). Not marked superseded. |
| `docs/finance/FINANCIAL_MODELS.md` | n/a | Model specification (formulas). No empirical claims. |
| `docs/calibrate/V018_PRIORITIES.md` | n/a | Plan. |

## 4. Number-matcher results per document

Produced by `scripts/v018/audit_numbers.py`; counts copied from `docs/v018/evidence/audit_numbers.summary.json`. "Numbers" counts discriminating tokens only. "In repo", "Only off-repo" and "Unmatched" sum to "Numbers". "Decoy matches" counts tokens whose value × 1.0371 also matched; compare it to "In repo" + "Only off-repo" to judge how informative a match is.

| Document | Evidence groups searched | Numbers | In repo outputs | Only in off-repo outputs | Unmatched | Decoy matches |
|---|---|---:|---:|---:|---:|---:|
| `README.md` | v014, v016, v017, phase2, crypto, detector, benchmarks, tests | 36 | 34 | 2 | 0 | 31 |
| `docs/REPRODUCIBILITY.md` | v014 | 6 | 3 | 0 | 3 | 1 |
| `docs/architecture/STANDARDS_TRACEABILITY.md` | v014 | 13 | 0 | 0 | 13 | 0 |
| `docs/archive/RESEARCH_EVOLUTION_LOG.md` | v014 | 9 | 4 | 0 | 5 | 2 |
| `docs/archive/SCIENTIFIC_CAPABILITY_INVENTORY.md` | v014 | 6 | 1 | 0 | 5 | 1 |
| `docs/benchmarks/GPU_L4_AUDIT_SUMMARY.md` | benchmarks | 28 | 0 | 0 | 28 | 0 |
| `docs/benchmarks/RNGLAB_CPU_GPU_COMPARISON.md` | benchmarks | 39 | 0 | 0 | 39 | 0 |
| `docs/calibrate/ACADEMIC_SUMMARY.md` | v017, v016 | 17 | 6 | 11 | 0 | 14 |
| `docs/calibrate/REVIEW_CALIBRATE_AUDIT.md` | v017, v016 | 10 | 1 | 9 | 0 | 4 |
| `docs/calibrate/V017_FAILURE_DEFINITIONS.md` | v017, v016 | 1 | 0 | 1 | 0 | 1 |
| `docs/calibrate/V017_FORENSIC_BASELINE.md` | v017, v016 | 2 | 2 | 0 | 0 | 2 |
| `docs/calibrate/V017_HEALTHY_CALIBRATION.md` | v017, v016 | 23 | 3 | 20 | 0 | 19 |
| `docs/calibrate/V017_ISING_VALIDATION.md` | v017, v016 | 49 | 3 | 46 | 0 | 7 |
| `docs/calibrate/V017_NEGATIVE_FINDINGS.md` | v017, v016 | 9 | 1 | 8 | 0 | 3 |
| `docs/calibrate/V017_PREDICTIVE_RELIABILITY.md` | v017, v016 | 99 | 13 | 86 | 0 | 75 |
| `docs/calibrate/V017_REPRODUCIBILITY.md` | v017, v016 | 1 | 1 | 0 | 0 | 1 |
| `docs/calibrate/V017_SCIENTIFIC_PROTOCOL.md` | v017, v016 | 23 | 23 | 0 | 0 | 19 |
| `docs/calibrate/V017_SCIENTIFIC_REPORT.md` | v017, v016 | 11 | 2 | 9 | 0 | 5 |
| `docs/calibrate/V017_VMC_PIMC_VALIDATION.md` | v017, v016 | 17 | 10 | 6 | 1 | 13 |
| `docs/crypto/ARCHITECTURE_AUDIT.md` | crypto, tests | 8 | 8 | 0 | 0 | 2 |
| `docs/crypto/EXPERIMENT_PROTOCOL.md` | crypto, tests | 17 | 17 | 0 | 0 | 16 |
| `docs/crypto/FAULT_INTERFACE.md` | crypto, tests | 1 | 1 | 0 | 0 | 1 |
| `docs/crypto/PHYSICAL_INTEGRATION.md` | crypto, tests | 3 | 0 | 0 | 3 | 0 |
| `docs/crypto/REPRODUCE.md` | crypto, tests | 4 | 4 | 0 | 0 | 0 |
| `docs/crypto/RESULTS.md` | crypto, tests | 75 | 71 | 0 | 4 | 61 |
| `docs/crypto/SCIENTIFIC_REVIEW.md` | crypto, tests | 10 | 10 | 0 | 0 | 4 |
| `docs/crypto/STANDARDS_TRACEABILITY.md` | crypto, tests | 9 | 9 | 0 | 0 | 9 |
| `docs/crypto/THREAT_MODEL.md` | crypto, tests | 2 | 2 | 0 | 0 | 2 |
| `docs/finance/FINANCIAL_MODELS.md` | v014 | 7 | 0 | 0 | 7 | 0 |
| `docs/physical/GM_NEXT_WEEK.md` | detector | 1 | 0 | 0 | 1 | 0 |
| `docs/physical/MUON_DETECTOR_ASSESSMENT.md` | detector | 121 | 31 | 0 | 90 | 0 |
| `docs/physical/PIPELINE.md` | detector | 7 | 0 | 0 | 7 | 0 |
| `docs/physical/REPRODUCE.md` | detector | 5 | 0 | 0 | 5 | 1 |
| `docs/trident/EXECUTIVE_SUMMARY.md` | v016 | 12 | 0 | 12 | 0 | 1 |
| `docs/trident/V016_ISING_STUDY.md` | v016 | 113 | 0 | 112 | 1 | 3 |
| `docs/trident/V016_MATHEMATICAL_THEORY.md` | v016 | 2 | 0 | 2 | 0 | 2 |
| `docs/trident/V016_NEGATIVE_FINDINGS.md` | v016 | 12 | 0 | 12 | 0 | 1 |
| `docs/trident/V016_PREDICTIVE_RELIABILITY.md` | v016 | 18 | 0 | 18 | 0 | 7 |
| `docs/trident/V016_REPRODUCIBILITY.md` | v016 | 2 | 0 | 2 | 0 | 1 |
| `docs/trident/V016_SCIENTIFIC_PROTOCOL.md` | v016 | 4 | 0 | 4 | 0 | 3 |
| `docs/trident/V016_SCIENTIFIC_REPORT.md` | v016 | 18 | 0 | 18 | 0 | 7 |
| `docs/validation/APPLICATION_RELIABILITY.md` | v014 | 22 | 19 | 0 | 3 | 2 |
| `docs/validation/EXPERIMENT_REGISTRY.md` | v014 | 46 | 31 | 0 | 15 | 3 |
| `docs/validation/REVIEW_METHODOLOGY_CRITIQUE.md` | v014 | 14 | 4 | 0 | 10 | 1 |
| `docs/validation/REVIEW_STATISTICAL_CORE.md` | v014 | 11 | 1 | 0 | 10 | 0 |
| `docs/validation/SOURCE_COMPARISON.md` | v014 | 38 | 10 | 0 | 28 | 3 |
| **All** | | 981 | 325 | 378 | 278 | 328 |

Reading guide: `V016_ISING_STUDY.md` (112 off-repo matches, 3 decoy matches) is genuinely traceable to off-repo outputs. `V017_PREDICTIVE_RELIABILITY.md` (99 matches, 75 decoy matches) is not established by the matcher and relies on the hand checks above. `README.md` (36 matches, 31 decoy) was checked claim by claim in §1.

## 5. Other defects found during the audit

1. The "SUPERSEDED" banners in `docs/validation/` and `docs/architecture/` link to `V0141_CORRECTIONS.md` in the same folder; the file is in `docs/archive/`. Broken link.
2. `src/rnglab/` already exists as the legacy v0.1–v0.8 package (`__version__ = "0.1.0"`). The v0.18 namespace has to coexist with it; see `ARCHITECTURE.md`.
3. Three v0.14-era documents with superseded content carry no banner: `SOURCE_COMPARISON.md`, `APPLICATION_RELIABILITY.md`, `REVIEW_METHODOLOGY_CRITIQUE.md`. `docs/REPRODUCIBILITY.md` is stale. Fixing this means adding banners to old files, which the brief's rule 3 forbids, so the rewritten README lists them as superseded instead.
4. All v0.16/v0.17 outputs and the Phase II NPZ traces live only on this machine, in tool-specific directories (`~/.codex/...`, `~/Documents/Codex/...`). They are a single point of failure for every "partially backed" status above. The publication plan should archive them (e.g. Zenodo).

## 6. Summary

- **Backed and central:** Phase II coverage and fault counts; test counts; crypto BLUE; detector inventory and saturation; the Ising absorbing-state condition (analytically).
- **Partially backed:** all v0.16/v0.17 numbers (scripts in Git, outputs only on this machine); the healthy susceptibility counterexample (exact match to the off-repo output); the detector assessment beyond §2; the finance seed-collision magnitude.
- **Unbacked:** the SOURCE_COMPARISON matrix (throughput 29–156× above what the code achieves for four of five sources; 15-test SP 800-22 battery that does not exist); the README VMC numbers (a re-run gives 1.656 / −0.248 / 0.296); regulatory framing; "formal proof"; "verified all tests by hand"; the GPU benchmarks as far as the repository is concerned; physical 90B entropy results.
