# v0.18 Phase 2 protocol: the Phase II core result through the package API

**Status: FROZEN DESIGN, 2026-11-05.** Approved by the author as drafted (the six choices in Section 10 take the first-listed option). The SHA-256 of `protocol.json` is in `protocol.sha256`. No campaign data had been generated at freeze. The only simulations run before freezing were a runtime smoke under the separate namespace `rnglab-v018-api-validation-smoke` (`smoke/smoke_timing.json`), which is excluded from every endpoint. Operational details fixed at freeze (trace blocks, path digest byte order, where check I-3 runs, IID secondary traces, an undefined Geyer τ̂, parallelism) are listed in `protocol.json` → `implementation_notes_added_at_freeze`. They change no gate or endpoint.

Machine-readable copy: [`protocol.json`](protocol.json). Design numbers come from [`design_calculations.py`](design_calculations.py) → [`design.json`](design.json). Phase II reference values come from `research/phase2-20261030/experiment/confirmation-v1/primary_decision.json`, via `design.json` → `phase2_reference`. Drafted by Claude Code (Opus 5.5); the author decides the design.

## 0. Freeze procedure

1. The author approves this draft, with any changes.
2. `protocol.json` gets `"status": "FROZEN_DESIGN"`. It is committed, and its SHA-256 is written to `protocol.sha256` and quoted in the commit message. From then on definitions, gates, sample sizes and namespaces do not change.
3. The campaign, analysis and recount scripts are written. Their tests use only the smoke namespace or synthetic inputs.
4. A pre-generation manifest (`manifest_pregeneration.json`) is committed before the first campaign draw. It records the SHA-256 of every campaign script and of the code under test, plus the environment lock. The code under test is `src/rnglab14/models/pimc_harmonic.py`, `src/rnglab14/sources/`, `src/rnglab12/uncertainty.py` and `src/rnglab/longrun.py`, at or after commit `4cf55a8`.
5. After that, any change is a new, separately versioned protocol file stating why. A failed gate is reported as failed. Nothing is re-run under new seeds to obtain a pass, and no other batch size is chosen after seeing the data.

## 1. Question and hypotheses

Phase II found its core result with a standalone NumPy/Numba implementation. With exact equilibrium starts at M = 16, β = 2, δ = 0.4 and n = 8192:
- local Metropolis under-covers: 811 of 1,024 CV BM512 intervals contain the reference;
- adding an exact centroid heat-bath refresh restores calibration: 964 of 1,024;
- an IID Gaussian-path control is calibrated: 969 of 1,024.

The remaining question is whether the repository's own package API reproduces this. The milestone that posed it is `research/phase2-20261030/reports/Research_advancement_and_next_milestone.md`, last section.

- **H-API (primary).** Through the package API, the refreshed sampler and the IID control both give calibrated CV BM512 intervals, meaning both pass the Phase II healthy gates (G-R and G-I).
- **H-REP (core-result replication).** In addition, plain local Metropolis under-covers (F-L), and the refresh improves paired coverage by more than 0.10 (P-RL). If H-API and H-REP both hold, the Phase II core result is reproduced through the API.

**Known versus new.** The sampler, the estimator and the gates are not new. What is new is an independent implementation path with fresh streams:
- the package API, run through the reviewed v0.18 run options;
- the package's pure-Python PCG64 sources;
- the package's acceptance-word consumption, read only for uphill moves where Phase II read one for every proposal;
- the package's 16-byte Box–Muller centroid draw.

The Markov kernel is the same in distribution as Phase II's: the same sequential scan, proposal, acceptance rule, centroid law and CV measurement point (checked against `study.py` `local_kernel`). Trajectories cannot match bit for bit, so the replication is distributional.

## 2. Definitions

- **Model.** `PathIntegralHarmonicOscillator(beads=16, temperature=0.5, displacement_delta=0.4)`: primitive action, m = ħ = ω = 1. Each sweep visits beads 0…15 in order. A proposal is x_j + δ(2u − 1), accepted by the Metropolis rule. The acceptance word is read only when ΔS > 0.
- **Refresh.** The model's `harmonic_centroid_heatbath`, applied once after every sweep. It draws the new centroid exactly from N(0, 1/β), using 16 bytes per draw.
- **Observable.** The centroid-virial energy CV = 1/(2β) + mean(x²) − c̄²/2, recorded after every sweep, including sweep 0. Its exact finite-bead mean is E_M = `exact_finite_bead_energy(16, 2.0)` = 0.6554741522981734 (`design.json` → `reference_energy`). The thermodynamic estimator TD is secondary.
- **Primary interval.** `summarize_series(cv, reference=E_M, batch_sizes=(512,))["batch_means"]["512"]`. This uses 16 non-overlapping batches of 512, with MCSE = sd(batch means)/4 and interval mean ± t₁₅,₀.₉₇₅ · MCSE. A replica is covered if E_M lies in the closed interval. This is the same rule as Phase II.
- **IID control.** Each replica draws 8,192 independent exact paths with `sample_exact_finite_bead_paths(source, 16, 2.0, 8192)` and computes CV per path with `harmonic_energy_estimators`. It is then analysed with the same BM512 interval.
- **Unit.** One whole-trajectory replica. Sweeps are never units.

## 3. Data generation

**Namespace.** `rnglab-v018-api-validation`. It is unused anywhere in the repository or the local Codex output folders (searched 2026-11-05). It is disjoint from Phase II (`rnglab-phase2-v1/…`, NumPy PCG64) and from Phase 1.

**Seed rule.** `namespace_seed("rnglab-v018-api-validation", condition, role, f"{i:05d}")`, which is the first 8 bytes of SHA-256, big-endian. Each seed goes into its own `PCG64Source(seed)`. Every role has its own source object, and the per-replica source ledger must show no aliasing.

| Condition key | Role | Used by | Bytes per replica |
|---|---|---|---:|
| `metropolis` | `start` | local and refresh (paired) | 256, before the run |
| `metropolis` | `proposal` | local and refresh (paired) | 1,048,576 |
| `metropolis` | `acceptance` | local and refresh (paired) | varies (uphill moves only) |
| `refresh` | `centroid` | refresh | 131,072 |
| `iid` | `paths` | IID control | 1,048,576 |

- **Start.** `sample_exact_finite_bead_paths(start_source, 16, 2.0, 2)[0]`, the same mapping as the API's `"exact_gaussian"` mode. It is passed as `initial_path`, with `initial_path_provenance` = {namespace, condition, role, replica, seed, sampler}. The start source is also passed as the initialization channel, which then reads 0 bytes.
- **Runs.**
  - Local: `run_simulation(channels, sweeps=8192, burn_in=0, initial_path=start, initial_path_provenance=…, uncertainty_batch_sizes=(512,), record_trace=True)`.
  - Refresh: the same call with `centroid_source`.
  - Pairing: local and refresh share each replica's start, proposal and acceptance seeds, as in Phase II. Replicas are independent across i under the healthy-PRNG assumption.
- **Replicas.** **R = 4,096 per condition** (Section 10, choice 1), i = 0…R − 1.
- **Saved.**
  - Traces, outside Git with their SHA-256 in the committed manifests: `traces/<condition>/block-XXXX.npz`, holding per replica the CV, TD, centroid and rg² series (float64, 8,192 each) and the start path.
  - Full paths are not saved. Each replica's path-array SHA-256 is recorded, and the paths are regenerable by deterministic replay.
  - Committed: `results/<condition>_replicas.csv`, with one row per replica giving seeds, start SHA-256, acceptance rate, CV mean, BM512 MCSE, interval, covered flag, the same for TD, and bytes and final digest per role.

## 4. Analysis rules

- **Coverage.** Covered count out of R, with a Clopper–Pearson 95% interval.
- **Variance ratio.** The mean over replicas of MCSE², divided by the between-replica variance of the CV means (ddof 1).
  - 95% percentile bootstrap over replicas, with 10,000 resamples.
  - Bootstrap seed: `namespace_seed(ns, "analysis", "bootstrap", endpoint)`.
- **Bias.** The mean over replicas of (CV mean − E_M), with a t_{R−1} 95% interval.
- **Paired difference.** d_i = covered_refresh,i − covered_local,i.
  - Its mean, with a 95% percentile bootstrap over replicas (10,000 resamples).
  - Bootstrap seed: `namespace_seed(ns, "analysis", "bootstrap", "paired")`.
- **Exclusions.** None. A replica whose BM512 status is not `ESTIMATED` counts as not covered and is reported separately.
- **Pass rule.** A gate passes only if every criterion listed for it holds. There is no multiplicity adjustment: these are acceptance requirements, not significance tests. Section 7 gives the chance-failure probabilities.

## 5. Gates (frozen at approval)

### 5.1 Implementation checks (deterministic; any failure is an implementation defect: stop, investigate, report)

- **I-1, byte ledger** (every replica): proposal 1,048,576; initialization 0; centroid 131,072 for refresh and absent for local; IID paths 1,048,576; resampling absent; no role aliases. Expected values are derived in `design.json` → `expected_bytes_per_replica`; the in-run values were observed in the smoke.
- **I-2, start pairing:** for each i, local and refresh record the same start SHA-256, equal to the digest of the start regenerated from its seed.
- **I-3, trace consistency** (every replica): CV recomputed from the recorded paths matches the CV trace to within 1e-12 absolute. The smoke maximum was 4.4e-16.
- **I-4, deterministic replay:** replica 0 of each condition, regenerated after the campaign, gives an identical CV-trace SHA-256.

### 5.2 Required validation gates (H-API; the milestone's gates)

| Gate | Condition | Coverage CP95 inside | Variance-ratio bootstrap CI95 inside | Bias CI95 inside |
|---|---|---|---|---|
| G-R | local + centroid refresh | [0.92, 0.98] | [0.75, 1.25] | [−0.01, 0.01] |
| G-I | IID Gaussian paths | [0.92, 0.98] | [0.75, 1.25] | [−0.01, 0.01] |

- H-API is **supported** if G-R and G-I both pass.
- Otherwise the API validation has **failed**, which is a finding. The milestone is then complete only when the discrepancy is localised and documented with a minimal reproducer.

### 5.3 Core-result replication (H-REP)

- **F-L, local under-coverage reproduces:** local Metropolis coverage CP95 upper bound < 0.92.
- **P-RL, refresh improvement:** the lower bootstrap 95% bound of the mean paired difference is > 0.10. This is Phase II's pre-specified margin.

H-REP is **supported** if F-L and P-RL both pass. The Phase II core result is reproduced through the API if H-API and H-REP are both supported.

### 5.4 Independent recount

A separate script recomputes every coverage event (3 × R) and every gate statistic from the saved CV traces. It uses its own NumPy/SciPy code, not `summarize_series`. All counts must agree. A disagreement blocks the report until it is explained.

## 6. Secondary endpoints (descriptive, no gates)

- **TD with BM512:** coverage, variance ratio and bias in all three conditions.
- **CV with the Phase 1 estimators** (`rnglab.longrun`): BM(32) and BM(128), which are the package's legacy defaults, plus OBM(512), SV-B(512), SV-TH(512) and Geyer. Report coverage and variance ratio per condition.
  - Geyer's pair theorem assumes reversibility, which a systematic sweep lacks (Phase II theory document §4), so Geyer is a sensitivity result only.
- **Batch-length diagnostic** on CV (10·τ̂ and at least 16 batches): status counts, and the coverage of the non-abstaining adaptive intervals.
  - Expectation recorded before generation, from the Phase 1 report: local Metropolis mostly abstains with RUN_TOO_SHORT, because Phase II's median τ̂ of 222.5 exceeds the 51.2 that n = 8192 allows. Refresh (τ̂ about 3.6) and IID (about 1.0) are mostly judged adequate (`design.json` → `phase2_reference`).
- **Phase II mechanism:** the between-replica variance of the CV means split into a centroid part (c̄²/2) and a shape part (rg²), together with twice their covariance.
- **Other distributions:** acceptance rates, Geyer τ̂ and ESS.
- **Comparison with Phase II, per condition:** the coverage difference (API minus Phase II) with a Wald 95% interval. Descriptive only, because the two runs used different implementations and streams.

## 7. Sample sizes and runtime

From `design.json`.

- **Coverage gate, exact pass probability.**

  | True coverage | R = 1,024 | R = 2,048 | R = 4,096 |
  |---:|---:|---:|---:|
  | 0.95 | 0.968 | 0.99983 | > 0.999999 |
  | 0.94 | 0.662 | 0.937 | 0.9989 |
  | 0.93 | 0.191 | 0.373 | 0.664 |

  Phase II's refresh point estimate was 0.941.
- **Variance-ratio gate.** Approximate pass probability at Phase II's refresh ratio of 0.944: 0.9989 at R = 1,024, and above 0.99999 at 2,048 and 4,096. This approximation reproduces Phase II's own bootstrap intervals: 0.863–1.032 against 0.865–1.036 for refresh.
- **Bias gate.** The CI95 half-width at Phase II's between-replica variance is about 0.0005 at R = 1,024 and 0.00025 at R = 4,096, far inside ±0.01.
- **Both required gates passing** (approximate):

  | Scenario | R = 1,024 | R = 2,048 | R = 4,096 |
  |---|---:|---:|---:|
  | Perfect calibration | 0.933 | 0.9997 | ≈ 1 |
  | Refresh at Phase II point estimates, IID perfect | 0.704 | 0.964 | 0.9997 |

- **P-RL power** at Phase II's difference of 0.149: 0.918, 0.997 and above 0.99999 for R = 1,024, 2,048 and 4,096.
- **Runtime.** Smoke means per replica: 0.92 s local, 1.10 s refresh, 0.083 s IID. At R = 4,096 that is about 8,600 single-core seconds, or about 36 minutes with 4 worker processes (one BLAS thread each, peak about 107 MiB each). CPU only.
- **Storage.** Traces need about 3.0 GiB at R = 4,096, outside Git.

## 8. What this campaign cannot show

- **One point only:** one model, one parameter point (M = 16, β = 2, δ = 0.4), one run length, exact starts and no burn-in. A pass does not transfer to other models, starts, run lengths or uncertainty methods.
- **Healthy sources are assumed, not certified.** The package PCG64 is a deterministic reference PRNG, and no entropy or independence claim is made.
- **Distributional, not bit-for-bit.** A pass means the package API reproduces the Phase II calibration in this domain. It does not make trajectories match Phase II bit for bit, and it does not prove the API correct in general. Regression tests are not coverage evidence.
- **No RNG-fault evidence.** This campaign tests none of H1, H2 or H3. It establishes the healthy API baseline on which the Phase 3 detection campaign (H1) builds.

## 9. Deliberately excluded

- Faulted sources: these belong to Phase 3, under its own protocol and namespace.
- Other bead counts, temperatures, step sizes, burn-in rules and non-stationary starts.
- Re-analysing the Phase II traces: those data have been seen.
- The API's `"exact_gaussian"` mode as a separate condition: it uses the same mapping as the explicit start and is covered by unit tests.

## 10. Choices for the author (decided 2026-11-05: first-listed option in each)

1. **Replicas per condition: 4,096**, 2,048, or 1,024 (the brief's value).
   - At 1,024 there is about a 30% chance that the required gates fail even if the API exactly reproduces Phase II, so a failure would say little about the API.
   - At 2,048 that chance is about 4%; at 4,096 it is about 0.03%.
   - 4,096 costs about 36 minutes and 3 GiB of traces outside Git.
2. **Gate H-REP (F-L and P-RL)**, or report them descriptively only. The milestone required only G-R and G-I.
3. **Start: an explicit `initial_path` from a dedicated start role**, which exercises the provenance API and makes pairing visible through the start digest. The alternative is the API's `"exact_gaussian"` mode.
4. **Pairing: local and refresh share start, proposal and acceptance seeds**, as in Phase II, or use independent seeds and an unpaired comparison.
5. **Bootstrap: 10,000 resamples for every interval**, or Phase II's 1,500 for the variance ratio. The criteria are the same either way; more resamples only reduce Monte Carlo noise in the interval endpoints.
6. **Storage: CV, TD, centroid and rg² traces plus starts, outside Git** (about 3 GiB), or CV only (about 0.75 GiB at R = 4,096). CV only would drop the TD and mechanism secondaries.

## References

- Flegal, J. M. & Jones, G. L. (2010). Batch means and spectral variance estimators in Markov chain Monte Carlo. *Ann. Statist.* 38(2), 1034–1070. doi:10.1214/09-AOS735 (checked 2026-11-05, Phase 1)
- Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statist. Sci.* 7(4), 473–483. doi:10.1214/ss/1177011137 (checked 2026-11-05, Phase 1)
- Meketon, M. S. & Schmeiser, B. (1984). Overlapping batch means: something for nothing? *Proc. 1984 Winter Simulation Conference*, 227–230. (checked 2026-11-05, Phase 1)
- Yamamoto, T. (2005). Path-integral virial estimator based on the scaling of fluctuation coordinates: application to quantum clusters with fourth-order propagators. *J. Chem. Phys.* 123, 104101. doi:10.1063/1.2013257 (checked 2026-11-05: arXiv physics/0505109 and the Kyoto University repository copy). Cited for the centroid-virial estimator.
