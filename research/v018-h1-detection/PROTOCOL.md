# v0.18 Stage C protocol: application-aware detection at matched false alarms (H1)

**Status: FROZEN DESIGN, 2026-11-06.** Approved by the author as drafted (the eight choices in Section 10 take the first-listed option). The SHA-256 of `protocol.json` is in `protocol.sha256`. Grid values are not frozen here; they come from the pilot through the frozen rule (Section 3.4). No pilot or campaign data had been generated at freeze. The only runs before freezing were a runtime smoke under the separate namespace `rnglab-v018-h1-detection-smoke` (`smoke/smoke_timing.json`), which is excluded from every endpoint. Operational details fixed at freeze (byte-mark convention, onset sweep, LAG1 first word, DUP decision draws, GFSR state and warm-up, SP 800-22 parameters, health-test settings, storage, parallelism) are in `protocol.json` → `implementation_notes_added_at_freeze`. They change no gate or endpoint.

Machine-readable copy: [`protocol.json`](protocol.json). Design numbers come from [`design_calculations.py`](design_calculations.py) → [`design.json`](design.json). Claim scope follows the adopted novelty review (`docs/v018/RELATED_WORK.md`, Section 1) and the roadmap (`docs/v018/ROADMAP.md`, Stage C). Drafted by Claude Code (Opus 5.5); the author decides the design.

## 0. Freeze procedure (two steps, because the severity grids come from a pilot)

1. The author approves this draft. `protocol.json` is frozen with its SHA-256 recorded, as in Phases 1 and 2. This fixes the definitions, monitors, tests, gates, sample sizes, namespaces and **the pilot-to-grid rule** (Section 3.4), but not the grid values themselves.
2. The fault adapters, sweep clock, conventional tests and monitors are implemented. The added SP 800-22 tests are unit-tested against the worked examples in SP 800-22 Rev 1a, and the other parts against exact or synthetic inputs.
3. **Pilot**, under its own namespace and excluded from every endpoint. The frozen `make_grid.py` turns the pilot output into `grid.json` by the frozen rule.
4. A pre-generation manifest is committed and pushed before the first campaign draw. It records the SHA-256 of `grid.json`, every script, the code under test and the environment lock.
5. After that, any change needs a new, versioned protocol. A failed gate is a finding.

## 1. Question and hypotheses

**H1-C (per fault family).** At a whole-run false-alarm rate of 0.05 for both arms, there is a severity at which three things hold together:
1. the application's reported uncertainty interval is **wrong** in most runs;
2. a conventional test battery run online on **the same consumed byte stream** stays **quiet** in most runs;
3. reference-free application monitors **alarm more often** than the battery, paired by run.

Healthy sampler failure is controlled separately (gate V1).

Every severity cell is classified (Section 5.3). Three outcomes get the same prominence as an H1 cell:
- **conventional first:** the battery detects while the application is still unharmed;
- **undetected harm:** the application is wrong and neither arm alarms;
- **no harm in range.**

**Known versus new** (`RELATED_WORK.md`). That application tests catch generators which pass standard tests has been established since 1992. What is new here is the calibrated comparison:
- online, on the same stream;
- at matched whole-run false-alarm rates;
- with latencies in bytes and sweeps;
- with a validated healthy baseline.

The claim is never "application tests find what standard tests miss" in general.

## 2. Definitions

### 2.1 Application

- **Sampler.** The package API sampler validated in Phase 2: `PathIntegralHarmonicOscillator(beads=16, temperature=0.5, displacement_delta=0.4)` with the centroid refresh, an exact start (explicit `initial_path`), 10,240 sweeps and no burn-in.
- **Local source.** Proposal and acceptance words come from **one** local source, so "the consumed stream" is a single, unambiguous byte sequence. That source's initialization channel reads 0 bytes.
- **Centroid source.** A separate, healthy source. It is never faulted.

### 2.2 Sweep clock

The refresh reads the centroid source exactly once per sweep, after the local updates. A thin wrapper around the centroid source therefore records the local source's consumed-byte count at every sweep boundary. That gives an exact map from bytes to sweeps without changing the API.

### 2.3 Fault adapters

The local source is wrapped. Every emitted byte is recorded, and from the **onset byte offset B₀ = 2¹⁹ = 524,288** each 8-byte word is transformed. The application's `random_float` reads big-endian words and uses the top 53 bits, u = (w ≫ 11)/2⁵³. In the smoke, the onset fell at sweep 2,479–2,482.

| Family | Severity | Transformation from the onset word on |
|---|---|---|
| SHIFT | s | u′ = min(u + s, 1 − 2⁻⁵³); the low 11 bits are kept |
| LAG1 | ρ | Gaussian copula AR(1): z at the first transformed word = Φ⁻¹(u); afterwards z_t = ρ z_{t−1} + √(1−ρ²) Φ⁻¹(u_t + 2⁻⁵⁴), and u′ = min(Φ(z_t), 1 − 2⁻⁵³). Marginals stay uniform; consecutive words are correlated |
| DUP | c | With probability c, the word is replaced by the previously emitted word. The decisions come from a separate `fault` generator, never from the stream |
| GFSR | lag pair (p, q) | **Real weak generator (C3).** At onset the source switches to x_n = x_{n−p} ⊕ x_{n−q} on 32-bit words, two outputs per 64-bit word, with the state filled from the `gfsr` seed |
| Healthy | 0 | Identity; bytes are still recorded |

**GFSR lag pairs** are (31,13), (89,38), (127,63), (250,103) (R250, documented in [FLW92]), (521,168) and (1279,418), ordered from the shortest lag. A script must verify that each trinomial is irreducible before the pair is used. A pair that fails is dropped, and this is reported.

### 2.4 Application harm

Harm is assessed on the CV trace over a **fixed post-onset window**, sweeps 2,560–10,239 (7,680 sweeps):
- the interval is BM with 16 batches of 480 and t₁₅: `summarize_series(cv, reference=E_M, batch_sizes=(480,))`;
- a run is **wrong** if its interval excludes E_M = `exact_finite_bead_energy(16, 2.0)`;
- a status other than `ESTIMATED` counts as wrong and is reported separately;
- every run's onset sweep must be at most 2,560. A run violating this is flagged and reported, and its harm endpoint is excluded.

### 2.5 Conventional arm C, on the recorded local byte stream

**Windows.** Windows of 512 sweeps, mapped to byte ranges by the sweep clock, at checkpoints every 256 sweeps (window ends 512, 768, …, 10,240; 39 checkpoints). At the smoke's 211.4 bytes per sweep (`smoke/smoke_timing.json`), a 512-sweep window stays inside the 1,000,000-bit limit of the in-repository SP 800-22 implementation; a 1,024-sweep window would not.

**Tests in each window:**
- **SP 800-22 Rev 1a:** frequency, block frequency and runs (in-repository `rnglab11.evaluate_bits`), plus serial (m = 8), approximate entropy (m = 8), forward cumulative sums and longest run of ones in a block. The added four are implemented in research code and checked against the standard's worked examples.
- **Word-level** (Knuth-style), on the 64-bit words:
  - lag-1 serial correlation of the 53-bit uniforms, with the p-value from √n·r ~ N(0,1);
  - a count of adjacent identical words, with the p-value from Poisson((n−1)·2⁻⁶⁴).

**Arm statistic.**
- At each checkpoint, S_C = −log₁₀(min p over the 9 tests).
- The whole-run statistic is the maximum over all 39 checkpoints.
- The threshold τ_C is the 950th order statistic of 999 healthy calibration whole-run maxima.

**Health-test arm (descriptive, not calibrated).** The in-repository SP 800-90B RCT and APT on bytes, with the asserted min-entropy recorded and α = 2⁻²⁰.

### 2.6 Application arm A, reference-free

**Windows.** Windows of 1,024 sweeps at the same checkpoints (window ends 1,024 … 10,240; 37 checkpoints).

**Monitors:**
- **A1, estimator-pair invariant:** the window mean of D_t = TD_t − CV_t. Under the exact finite-bead target E[D] = 0, an identity that uses no reference energy (Phase II theory document). The score is z₁ = |mean D| / σ_D.
- **A2, dynamics:** the window fraction of uphill proposals, read from the byte ledger as (bytes in the sweep − 128)/(8·16). The score is z₂ = |f − μ_f| / σ_f.
- σ_D, μ_f and σ_f are estimated from the calibration runs.

**Arm statistic.** S_A = max(z₁, z₂) at each checkpoint. The whole-run statistic is the maximum over the checkpoints, and τ_A is the 950th order statistic of the healthy calibration maxima.

### 2.7 Detection and latency

- An arm **detects** in a run if any checkpoint whose window ends after the onset sweep exceeds its threshold.
- **Latency** is measured at the first such checkpoint:
  - in bytes: the local byte offset at the window end, minus B₀;
  - in sweeps: the window-end sweep minus the onset sweep.
- Checkpoints whose windows end at or before onset are false-alarm opportunities, exactly as in healthy runs.

## 3. Data generation

### 3.1 Namespaces

| Use | Namespace |
|---|---|
| Campaign | `rnglab-v018-h1-detection` (unused; searched 2026-11-06 in the repository and the local Codex folders) |
| Pilot | `rnglab-v018-h1-detection-pilot` |
| Smoke | `rnglab-v018-h1-detection-smoke` |

**Seed rule:** `namespace_seed(ns, condition, role, f"{i:05d}")`, giving a `PCG64Source` for the `local`, `centroid` and `start` roles. The `fault` and `gfsr` roles use NumPy `default_rng` with the same seed rule.

### 3.2 Conditions and pairing

- **Calibration:** condition key `calibration`, healthy, i = 0…998, with its own roles.
- **Test:** condition key `test`. Replica i uses the same local, centroid, start, fault and gfsr seeds in **every** test cell, so fault runs share their pre-onset prefix byte for byte with the healthy run of the same index.
  - Healthy test: i = 0…4,095.
  - Each fault cell: i = 0…255.

### 3.3 Cells

| Family | Severities | Replicas per cell |
|---|---|---|
| SHIFT | 6, from the pilot rule | 256 |
| LAG1 | 6, from the pilot rule | 256 |
| DUP | 6, from the pilot rule | 256 |
| GFSR | the 6 lag pairs (fewer if some fail the irreducibility check) | 256 |

### 3.4 Pilot and grid rule (frozen now; values computed later)

- **Pilot healthy runs:** 199, giving pilot-only thresholds at the 190th order statistic.
- **Pilot severities**, 32 replicas each:
  - SHIFT: 10 log-spaced values in [1e-4, 0.1];
  - LAG1: 10 log-spaced values in [0.003, 0.5];
  - DUP: 10 log-spaced values in [0.001, 0.5];
  - GFSR: all 6 lag pairs.
- **Rule for SHIFT, LAG1 and DUP.** Let w(s) be the pilot wrong-rate.
  - s_lo is the smallest pilot severity with w ≥ 0.10.
  - s_hi is the smallest pilot severity with w ≥ 0.90, or the largest pilot severity if none reaches it.
  - The campaign grid is 6 log-spaced severities from s_lo/2 to s_hi.
  - If no pilot severity reaches w ≥ 0.10, the grid is 6 log-spaced severities from s_max/32 to s_max, and the family is labelled "no harm observed in the pilot range".
- Pilot outcomes are reported descriptively only.

### 3.5 Saved data

- **Outside Git:** per run, the CV, TD and r_g² traces and the sweep byte marks.
- **Committed**, per run:
  - seeds;
  - local-stream SHA-256 and pre-onset-prefix SHA-256;
  - onset sweep;
  - harm endpoint;
  - per-checkpoint arm statistics and per-test p-values;
  - alarm flags and latencies.
- **Raw bytes are not saved.** They can be regenerated by deterministic replay.

## 4. Analysis rules

- **Unit.** One run; sweeps and windows are never units.
- **Rates.** Clopper–Pearson intervals: 95% for description, and the Bonferroni level 1 − 0.05/24 for the H1 cell criteria.
- **Paired comparisons** within a cell (application alarm minus conventional alarm, per run). Percentile bootstrap over runs, 10,000 resamples, at the same levels, with seeds `namespace_seed(ns, "analysis", "bootstrap", endpoint)`.
- **Latency.** Among runs where both arms detect, the median paired difference (A − C) in bytes and in sweeps, with a bootstrap 95% interval, plus the counts where only one arm detects or neither does.
- **No exclusions** beyond the flagged onset rule (Section 2.4).

## 5. Gates

### 5.1 Implementation checks (deterministic; a failure means stop, investigate and report)

- **I-1, sweep clock:** exactly 10,240 marks per run, monotone, with every sweep consuming 128 + 8k bytes.
- **I-2, prefix identity:** every fault run's pre-onset prefix SHA-256 equals that of the healthy test run with the same index.
- **I-3, fault realization:** for replica 0 of every cell, the faulty words regenerate exactly from the healthy words and the fault seed under the stated map.
- **I-4, replay:** replica 0 of every cell regenerates to identical trace and stream hashes.

### 5.2 Validity gates (if either fails, the H1 verdicts are reported as not interpretable)

- **V1, healthy application calibration:** the healthy test coverage on the post-onset window (4,096 runs) has its CP95 interval inside [0.92, 0.98].
- **V2, false alarms:** the healthy-test whole-run alarm rate of arm C and of arm A each has a CP95 upper bound ≤ 0.10. The paired difference between the two arms' healthy rates is reported.

### 5.3 Cell classification and the H1 verdict (Bonferroni level 1 − 0.05/24 for every bound)

| Class | Criterion |
|---|---|
| **H1 cell** | wrong-rate lower bound > 0.5 **and** C-alarm upper bound < 0.2 **and** paired (A − C) alarm-rate lower bound > 0 |
| Undetected harm | wrong-rate lower bound > 0.5 **and** C-alarm upper bound < 0.2 **and** A-alarm upper bound < 0.2 |
| Conventional first | C-alarm lower bound > 0.5 **and** wrong-rate upper bound < 0.2 |
| Both detect, harmed | wrong-rate lower bound > 0.5 **and** both alarm lower bounds > 0.5 |
| Harmless and quiet | wrong-rate upper bound < 0.2 **and** both alarm upper bounds < 0.2 |
| Indeterminate | anything else |

H1-C is **SUPPORTED for a family** if at least one of its cells is an H1 cell, and **NOT SUPPORTED** otherwise. No general H1 verdict is drawn across families.

### 5.4 Independent recount

A separate script recomputes every wrong flag from the saved CV traces, every alarm flag from the saved checkpoint statistics and thresholds, every cell rate, and the classification. All counts must agree.

## 6. Secondary endpoints (descriptive, no gates)

- Latency distributions, in bytes and in sweeps.
- Which conventional test fires first.
- A1 and A2 separately.
- The Phase II reference-based shape monitor, for continuity.
- A CUSUM of the standardized D [Pag54], with its threshold calibrated on the same healthy runs.
- The SP 800-90B health-test alarms.
- Mean CV bias and TD harm per cell.
- Pilot outcomes.

## 7. Sample sizes and runtime (`design.json`)

**Calibration.**
- 999 healthy runs, with the threshold at order statistic 950, give a marginal whole-run false-alarm rate of exactly 0.05.
- The rate conditional on the realized threshold spreads over 0.037–0.064 (2.5% and 97.5% quantiles). "Matched" therefore means matched in calibration, not identical in realization. V2 measures the realized rates.

**Precision at 256 runs per cell:**
- a 95% CP interval of [0.437, 0.563] at a rate of 0.5, and [0.067, 0.145] at 0.1;
- for the H1 rule (24 cells, Bonferroni level 0.99792), with true wrong-rate 0.7, C-alarm rate 0.1 and A − C gap 0.2, the three criteria hold with probabilities 0.9997, 0.888 and 0.976;
- at a wrong-rate of 0.6 the first criterion holds only with probability 0.507, so the rule is deliberately conservative.

**Healthy gate.** With 4,096 runs it passes with probability 0.9989 at a true coverage of 0.94.

**Runtime.**
- The smoke gives 2.62 s per run, including the battery.
- The campaign is 11,239 runs, about 8.2 single-core hours, or about 83 minutes on 8 workers (scaled by the 1.35 slowdown seen in Phase 2).
- The pilot is 1,351 runs, about 10 minutes.
- Peak memory is about 128 MiB per worker.
- Traces take about 2.6 GiB outside Git. CPU only.

## 8. What Stage C cannot show

- **Scope.** One application (harmonic PIMC at one parameter point), three synthetic word-level fault families and one real generator family. Results do not transfer to other applications or faults.
- **The battery.** Nine tests: a seven-test SP 800-22 subset plus two word-level tests, not the full SP 800-22 suite or TestU01. "Conventional tests pass" always names this battery, the window, the stream and the false-alarm rate (`ROADMAP.md` §8).
- **Matched false alarms are calibrated marginally,** not realized exactly.
- **Not H2 or H3.** No forecast or intervention is tested.
- **Synthetic faults.** Results for SHIFT, LAG1 and DUP say nothing about real generators. Only GFSR does.

## 9. Deliberately excluded

- **The Phase II channel-specific faults** (acceptance-only shift, proposal–acceptance coupling across separate sources). Single-stream batteries cannot see cross-channel coupling by construction, so comparing them on it would be unfair to the battery (Phase II report §3).
- TestU01 and PractRand: deferred, see choice 2.
- Physical sources, other models, time-varying severity profiles (Stage D), and interventions (Stage E).

## 10. Choices for the author (decided 2026-11-06: first-listed option in each)

1. **Fault families:** SHIFT, LAG1, DUP plus the GFSR lag family. Alternatives:
   - the Phase II channel-specific faults with separate channel sources;
   - a classic LCG (RANDU) as the real generator instead of, or as well as, GFSR.
2. **Conventional arm:** the nine-test in-Python battery, with TestU01 (Rabbit or Alphabit, version pinned) deferred to an external-battery replication in Stage F. Alternatives:
   - add TestU01 now, which needs your approval to download and build it;
   - the in-repository three tests only. Not recommended: the smoke showed they barely react to LAG1 or DUP, so the comparison would be a straw man.
3. **Application arm:** A1 and A2 combined. Alternatives: A1 alone; or add the Phase II shape monitor, which needs the known target, so the arm would no longer be reference-free.
4. **Harm:** the post-onset BM interval excludes E_M. Alternative: the Phase II rolling-window tolerance rule.
5. **Sizes:** 999 calibration runs, 4,096 healthy test runs, and 256 runs in each of 24 fault cells. Alternative: 512 runs per fault cell, which roughly doubles the fault cost.
6. **Grids:** the frozen pilot rule. Alternative: fix the grids now, without a pilot.
7. **H1 cell thresholds:** 0.5 and 0.2 at the Bonferroni level over 24 cells. Alternative: 95% per cell, with the multiplicity stated.
8. **Version label:** v0.18 Stage C. Alternative: v0.19 (roadmap decision 5).

## References (verified 2026-11-06; see `docs/v018/RELATED_WORK.md`)

- [FLW92] Ferrenberg, Landau & Wong, *Phys. Rev. Lett.* 69, 3382 (1992). doi:10.1103/PhysRevLett.69.3382
- [Pag54] Page, *Biometrika* 41, 100–115 (1954). doi:10.1093/biomet/41.1-2.100
- [SP800-22] L. Bassham, A. Rukhin, J. Soto, J. Nechvatal, M. Smid, S. Leigh, M. Levenson, M. Vangel, N. Heckert, D. Banks. *A Statistical Test Suite for Random and Pseudorandom Number Generators for Cryptographic Applications*. NIST SP 800-22 Rev 1a, April 2010 (checked 2026-11-06 against the NIST publication record and the legacy PDF at nvlpubs.nist.gov; no Rev 1a DOI confirmed). Its worked examples are the unit-test references for the added tests.
- [SP800-90B] Turan et al., NIST SP 800-90B (2018). doi:10.6028/NIST.SP.800-90B
- [LS07] L'Ecuyer & Simard, *ACM TOMS* 33(4), Art. 22 (2007). doi:10.1145/1268776.1268777
