# v0.18 Phase 1 protocol: estimator repair

**Status: FROZEN DESIGN, 2026-11-05.** Approved by the author as drafted (the six choices in Section 10 take the first-listed option). The SHA-256 of `protocol.json` is in `protocol.sha256`. No campaign data had been generated at freeze. The only simulations run before freezing were a runtime smoke under the separate namespace `rnglab-v018-estimators-smoke` (`smoke/smoke_timing.json`), which is excluded from every endpoint. Operational details fixed at freeze (trace storage, handling of an undefined Geyer τ̂, R̂ tie handling, parallelism) are listed in `protocol.json` → `implementation_notes_added_at_freeze`. They change no gate or endpoint.

Machine-readable copy: [`protocol.json`](protocol.json). Design numbers quoted below come from [`design_calculations.py`](design_calculations.py) → [`design.json`](design.json). Drafted by Claude Code (Opus 5.5) on 2026-11-05; the author decides the design.

## 0. Freeze procedure

1. The author approves this draft, with any changes.
2. `protocol.json` gets `"status": "FROZEN_DESIGN"`. It is committed, and its SHA-256 is written to `protocol.sha256` and quoted in the commit message. From then on definitions, gates, sample sizes and namespaces do not change.
3. The estimators and campaign scripts are implemented. Unit tests use only exact or synthetic inputs and the smoke namespace.
4. A pre-generation manifest (`manifest_pregeneration.json`) records the SHA-256 of every campaign script and of the `src/rnglab/` and `rnglab16` modules used, plus the environment lock. It is committed before the first campaign draw.
5. After that, any change is a new, separately versioned protocol file stating why. A failed gate is reported as failed. Nothing is re-run under new seeds to obtain a pass.

## 1. Questions and hypotheses

**H-S (susceptibility).** On healthy random-scan Metropolis Ising chains, the within-chain plug-in susceptibility fails at L=4, T=1.5 (reproducing v0.17). A symmetry-aware per-chain estimator and a generic multi-chain estimator both give calibrated intervals in all nine (L, T) cells with exact references.

**H-L (long-run variance).** On a stationary AR(1) process with known finite-n variance: batch means with short batches under-cover when the autocorrelation time is long (reproducing the Phase II mechanism); and BM, OBM, two spectral estimators and Geyer's estimator are calibrated when the run is long relative to τ.

**H-D (diagnostic).** A batch-length diagnostic based on Geyer's τ̂ either abstains or selects a batch size whose interval is calibrated, for every AR(1) setting. A multi-chain R̂ flags the Ising sign non-mixing that breaks the plug-in estimator.

Known versus new: all estimators are established (Section 2 cites them). What is new is their pre-registered validation, with exact references, on the two failures RNGLab already documented.

## 2. Definitions

### 2.1 Susceptibility

Lattice N = L² spins, periodic, J = k_B = 1, zero field; magnetization per spin m = M/N; temperature T.

- **Estimand.** χ = (N/T)·Var_π(m) = Var_π(M)/(N T), the zero-field susceptibility, which is exactly what `rnglab16.core.enumeration` returns.
- **P, per-chain plug-in (legacy; reproduce only).** χ̂_P = (N/T)(mean_t m_t² − (mean_t m_t)²) over the retained sweeps of one chain. Not a valid estimator when a chain does not visit both signs of m.
- **S, symmetry-aware, per chain.** χ̂_S = (N/T)·mean_t m_t². It uses E_π[m] = 0, which holds exactly at zero field on a finite periodic lattice by global spin-flip symmetry. Model-specific, not a general repair. Interval: batch means on the series (N/T)m_t² with 16 batches (b = 96 of the 1,536 retained sweeps) and t₁₅, mirroring the Phase II primary rule.
- **G, generic, ensemble of K = 16 independent chains.** With chain means A_i = mean_t m_t² and B_i = mean_t m_t (equal lengths), the pooled estimate is θ̂ = (N/T)(Ā − B̄²). The delete-one-chain jackknife gives θ̂_(−i) and the bias-corrected estimate θ̂_J = Kθ̂ − (K−1)·mean_i θ̂_(−i), with SE_J² = ((K−1)/K)·Σ_i(θ̂_(−i) − mean θ̂_(−·))². Interval: θ̂_J ± t₁₅·SE_J. The formula never uses E[m] = 0. Its validity needs the ensemble to be distributed as π. Here that holds exactly for the sign of m, because the initial law (independent uniform spins) and the kernel are both flip-symmetric, and holds approximately for |m| after burn-in. Outside symmetric designs there is no such guarantee, and the report must say so.
- **Secondary, different estimand.** χ_abs = (N/T)(E[m²] − E[|m|]²), per-chain plug-in, with its exact reference from enumeration. It is reported to show it is a different quantity, not as a repair.

### 2.2 Long-run variance of a scalar series x₁…x_n

The target is σ²_n = n·Var(x̄_n). For the AR(1) control it is exact; asymptotically it is σ² = Σ_k γ_k. Let γ̂_k be the biased (divide by n) sample autocovariances, computed by FFT.

| ID | Estimator | σ̂² | Interval for the mean |
|---|---|---|---|
| BM(b) | non-overlapping batch means, a = n/b batches | b/(a−1)·Σ_j(Ȳ_j − x̄)² | x̄ ± t_{a−1}·σ̂/√n |
| OBM(b) | overlapping batch means (Meketon & Schmeiser 1984) | n·b/((n−b)(n−b+1))·Σ_{j=1}^{n−b+1}(Ȳ_j(b) − x̄)² | t with ν = n/(b·∫w²) for the Bartlett window (∫w² = 2/3, so ν = 1.5n/b) |
| SV-B(b) | spectral variance, Bartlett window w(k) = 1 − abs(k)/b | Σ_{abs(k)<b} w(k)·γ̂_k | t with ν = n/(b·∫w²), ∫w² = 2/3 |
| SV-TH(b) | spectral variance, Tukey–Hanning w(k) = (1 + cos(πk/b))/2 | Σ_{abs(k)<b} w(k)·γ̂_k | t with ν = n/(b·∫w²), ∫w² = 3/4 |
| GEYER | initial positive then initial monotone sequence (Geyer 1992); existing `rnglab12.uncertainty._geyer` | −γ̂₀ + 2·Σ monotone pair sums | normal z (asymptotic) |

Batch sizes: b ∈ {32, 128, 512, 1024}. Reported with every estimate: τ̂ = σ̂²/γ̂₀ and ESS = n·γ̂₀/σ̂². The ν formula is the usual equivalent-degrees-of-freedom approximation for lag-window estimators. For OBM it gives 1.5·n/b, matching the approximation in Meketon & Schmeiser. Consistency of BM, OBM and spectral estimators under MCMC conditions: Flegal & Jones (2010).

### 2.3 Diagnostic

- **Batch length.** τ̂_G is Geyer's τ̂ on the same series. For a batch-type method (BM, OBM, SV) with batch or truncation length b, the status is **ADEQUATE** if b ≥ 10·τ̂_G and n/b ≥ 16; otherwise **WARN_BATCH_SHORT**.
- **Adaptive selection.** b* is the smallest power of two in [16, n/16] with b* ≥ 10·τ̂_G. If none exists, the status is **RUN_TOO_SHORT** and the method abstains: no interval is reported.
- **Why 10·τ̂.** For AR(1) the first-order relative bias of BM is −2φ/((1−φ)²·b·τ), about −5% at b = 10τ when τ is large (`design.json` → `bm_first_order_relative_bias`). Requiring 16 batches then needs n ≥ 160·τ. The constant is derived from AR(1); for other processes it is a heuristic, and testing it is part of H-D.
- **Convergence.** For each ensemble of 16 chains and each primitive series (m, m², e): rank-normalized split-R̂ (Vehtari et al. 2021). **NOT_CONVERGED** if R̂ > 1.01.

## 3. Data generation

Seed namespace: **`rnglab-v018-estimators`**, unused anywhere in the repository or the local Codex output folders (searched 2026-11-05). Seed rule: `int.from_bytes(sha256("rnglab-v018-estimators/<condition>/<role>/<unit>").digest()[:8], "big")`. The unit is a chunk; replicas within a chunk are `numpy.random.SeedSequence(seed).spawn(size)` children. The chunk-to-replica mapping is recorded in the manifest.

**A. Ising.** Use `rnglab16.core.ising` unchanged: random-scan Metropolis, fresh uniform random initial spins, healthy fault arm. Cells: L ∈ {2, 3, 4} × T ∈ {1.5, T_c = 2/ln(1+√2), 3.5}. Runs: 2048 sweeps, the first 512 discarded, 1,536 retained. **32,768 chains per cell**, in 128 chunks of 256 chains. Conditions are `ising_L{L}_T{1.5|Tc|3.5}`, role `chains`, unit `chunk-{0000..0127}`. Chain id = 256·chunk + child. Ensembles for G are consecutive groups of 16 chain ids (2,048 ensembles per cell). Exact references come from `rnglab16.core.enumeration`.

**B. AR(1).** x_t = φ·x_{t−1} + √(1−φ²)·ε_t with x₁ ~ N(0,1) and ε_t iid N(0,1) from NumPy `default_rng` (PCG64). Stationary from t = 1, mean 0, unit marginal variance. φ ∈ {0, 0.5, 0.9, 0.99, 0.997} (τ = 1, 3, 19, 199, 665.7); n ∈ {8192, 65536}. **R = 4,096 replicas** per (φ, n). Condition `ar1_phi{φ}_n{n}`, role `innovations`, unit `chunk-{c}` with 512 (n = 8192) or 64 (n = 65536) replicas per chunk. Exact Var(x̄_n) = [(1+φ)/(1−φ) − 2φ(1−φⁿ)/(n(1−φ)²)]/n.

**C. Regression replay (no new data).** Regenerate the v0.17 healthy confirmation cell `ising_L4_T0` from its own seed: `sha256("CALIBRATE-20261022-v1/ising_L4_T0/confirmation")`, 512 chains, `rnglab16.core.ising`. Require mean χ̂_P = 0.7458277684119011 and mean χ̂_S = 10.403679741753471, to relative 1e-9 (values from v0.17 `SUSCEPTIBILITY_PROBE.json`, row `ising_L4_T0`). This becomes a test in `tests/v018/`. It is a deterministic replay, not confirmation.

## 4. Analysis rules

- **Units.** One chain (S, P), one 16-chain ensemble (G), or one AR(1) replica. Sweeps within a run are never units.
- **Coverage.** Count of intervals that contain the exact reference, with a Clopper–Pearson 95% interval.
- **Relative mean bias.** For S and G: mean over units of (estimate/χ − 1), with a 95% t-interval.
- **Variance ratio.** Mean estimated variance of the estimator divided by the empirical variance of the estimates across units, with a percentile bootstrap over units (2,000 resamples, seed `rnglab-v018-estimators/analysis/bootstrap/<endpoint>`). For AR(1) it is also reported against the exact σ²_n.
- **Exclusions and retries.** No unit is excluded. Abstentions are counted and reported. Numerical failures (NaN or a non-positive σ̂²) count as non-coverage and are reported separately.
- **Pass rule.** A gate passes only if every criterion listed for it holds. No multiplicity adjustment: gates are acceptance requirements, not significance tests. At the planned counts, the chance that a method with true coverage 0.95 fails one coverage gate is 1.7e-4 at 2,048 units and below 1e-9 at 4,096 units (`design.json`).

## 5. Gates (frozen at approval)

### 5.1 Failure reproduction (expected to show the defect; non-reproduction is a finding)

| ID | Data | Criterion |
|---|---|---|
| S-F1 | A, L4 T1.5, estimator P | upper bound of the relative mean-bias 95% CI < −0.5 |
| L-F1 | B, φ = 0.99, n = 8192, BM(32) | coverage CP95 upper bound < 0.92 |
| L-F2 | B, φ = 0.997, n = 8192, BM(32) | coverage CP95 upper bound < 0.92 |
| D-F1 | A, L4 T1.5, split-R̂ of m over the 2,048 ensembles | fraction NOT_CONVERGED has CP95 lower bound ≥ 0.90 |

### 5.2 Control validity (a failure means a generator defect: stop and investigate, do not interpret estimators)

| ID | Data | Criterion |
|---|---|---|
| L-C | B, every (φ, n) | bootstrap 95% CI of empirical Var(x̄)/exact Var(x̄_n) contains 1; 95% t-CI of mean x̄ contains 0 |
| S-C | A, every cell | relative mean-bias 95% CI of mean_t e_t (energy per spin) against the exact energy lies within ±0.01 |

### 5.3 Validation gates (the hypotheses)

| ID | Data | Units | Criterion |
|---|---|---:|---|
| S-S[cell] ×9 | A, estimator S | 32,768 chains | coverage CP95 ⊂ [0.92, 0.98]; relative mean-bias CI ⊂ [−0.01, 0.01]; variance-ratio CI ⊂ [0.75, 1.25] |
| S-G[cell] ×9 | A, estimator G | 2,048 ensembles | same three criteria |
| L-A[φ, method] ×12 | B, n = 8192, b = 512, φ ∈ {0, 0.5, 0.9}, method ∈ {BM, OBM, SV-B, SV-TH} | 4,096 | coverage CP95 ⊂ [0.92, 0.98]; variance ratio against exact σ²_n: bootstrap CI ⊂ [0.75, 1.25] |
| L-G[φ] ×3 | B, n = 8192, GEYER, φ ∈ {0, 0.5, 0.9} | 4,096 | same two criteria |
| L-D[φ, n, family] ×30 | B, all φ and n; adaptive b* for family ∈ {BM, OBM, SV-TH} | non-abstaining replicas k | k ≥ 2,048: coverage CP95 ⊂ [0.92, 0.98] (pass/fail). 200 ≤ k < 2,048: **fail** if the CP95 upper bound < 0.92, otherwise inconclusive. k < 200: not evaluable. |

**Verdicts.** H-S is supported if S-F1 reproduces and every S-S and S-G gate passes. H-L is supported if L-F1 and L-F2 reproduce and every L-A and L-G gate passes. H-D is supported if D-F1 holds, no L-D gate fails, and at least one L-D gate per family is evaluable and passes. Every gate is reported individually whatever the verdict.

## 6. Secondary endpoints (descriptive, no gates)

- All methods and all b at φ ∈ {0.99, 0.997}: coverage, variance ratio, and median ESS-hat / true ESS, where true ESS = γ₀/Var(x̄_n).
- Distribution of τ̂_G/τ_n for each (φ, n); abstention and diagnostic-status rates per cell.
- Ising: S with OBM, SV-TH, GEYER and adaptive intervals; the diagnostic status of the m² series per cell; split-R̂ of m, m² and e in every cell (specificity is descriptive only); χ_abs against its exact value; the plug-in P bias in all nine cells.

## 7. Sample sizes and runtime

- **Replica counts.** CP95 half-width at true coverage 0.95: 0.0097 for 2,048 units, 0.0068 for 4,096, 0.0024 for 32,768. Probability of failing a gate when true coverage is 0.93: 0.63, 0.34 and 3.8e-7 respectively (`design.json`). Under-coverage of 2 percentage points is therefore usually detected at 2,048 to 4,096 units, and almost surely at 32,768.
- **Runtime.** From the smoke run, L4 costs 0.0127 s per chain for 2048 sweeps and L2 0.0032 s; L3 is scaled to 0.0071 s (not measured). Ising totals about 2,261 single-core seconds; AR(1) generation and FFT/cumsum work about 132 s, excluding the per-series Geyer loop. Planned execution: 3 worker processes, one BLAS thread each, 256-chain chunks (peak about 486 MiB per worker), comfortably inside the roughly 4 GiB of free RAM on this machine. Expected wall time is well under one hour. CPU only.

## 8. What this campaign cannot show

- Calibration for one model family (small Ising, AR(1)) at fixed run lengths. It does not transfer automatically to PIMC; Phase 2 carries these estimators as secondary endpoints.
- S's validity rests on an exact symmetry that most models lack. G's rests on a symmetric ensemble design.
- The 10·τ̂ constant is AR(1)-derived. Passing H-D on AR(1) does not prove it adequate for processes with several time scales, such as the PIMC centroid.
- Small lattices (L ≤ 4) and three temperatures: no thermodynamic-limit or critical-scaling statement.

## 9. Deliberately excluded

- Re-analysing the Phase II PIMC traces with the new estimators: those data have been seen, so it would be exploratory only. Phase 2 tests the estimators prospectively instead.
- Heat-capacity intervals: v0.17 errors were already small; this can be added later as its own protocol.
- Longer single-chain Ising runs as a third susceptibility route.

## 10. Choices for the author (decided 2026-11-05: first-listed option in each)

1. Generic estimator: the multi-chain jackknife G with K = 16. Alternatives: long single chains, or redefining the estimand via |m|. That would be χ_abs, kept secondary here because it is a different quantity.
2. Gate all nine Ising cells, or gate only L4 T1.5 and report the rest as secondary.
3. Diagnostic constants: 10·τ̂ and 16 batches, so n ≥ 160·τ̂.
4. S's primary interval: BM with 16 batches and t₁₅, consistent with Phase II. Alternative: the adaptive estimator.
5. Regression reference: the test can hold constants citing the off-repo file and its SHA-256, or `SUSCEPTIBILITY_PROBE.json` can be copied into `docs/v018/evidence/`.
6. Tolerances: relative bias ±1% (S, G) and variance ratio [0.75, 1.25] (as in Phase II).

## References (all checked 2026-11-05)

- Flegal, J. M. & Jones, G. L. (2010). Batch means and spectral variance estimators in Markov chain Monte Carlo. *Ann. Statist.* 38(2), 1034–1070. doi:10.1214/09-AOS735
- Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statist. Sci.* 7(4), 473–483. doi:10.1214/ss/1177011137
- Meketon, M. S. & Schmeiser, B. (1984). Overlapping batch means: something for nothing? *Proc. 1984 Winter Simulation Conference*, 227–230.
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. & Bürkner, P.-C. (2021). Rank-normalization, folding, and localization: an improved R̂ for assessing convergence of MCMC. *Bayesian Analysis* 16(2), 667–718. doi:10.1214/20-BA1221
