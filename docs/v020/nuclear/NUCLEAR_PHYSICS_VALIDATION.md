# RNGLab v0.20 nuclear module: physics validation report

**Status: Layer A only, 2026-11-16.** Results of the frozen study `research/v020-nuclear-a-verify`, run once at commit `c616a7b` on a clean tree (`results/verification.json`, `started_utc` 2026-11-16T12:17:07Z). Every number below is taken from `research/v020-nuclear-a-verify/results/verification.json` or from `results/tables.md`, which `report_tables.py` writes from it; design numbers come from `design.json`. Prepared for the author. Layers B and C have no validation yet.

## 1. What was verified, and what was not

**Question.** Does the Layer A event-level decay simulator, driven by healthy generators, reproduce exactly known decay laws for three synthetic schemes? And is the battery used for that check sensitive enough to notice small physics errors?

**Answer under the frozen reading rules: R1 PASS, R2 MET.**

| Law or property | Test (PROTOCOL.md) | Outcome |
|---|---|---|
| Exponential waiting times; survival probability 2⁻ᵗ | T1, T2 | not rejected (NC1, NC2) |
| Expected decay counts of each nuclide | T5 | not rejected |
| Daughter populations over time (Bateman, nine-state DAG and transient equilibrium) | T4, T9 | not rejected |
| Branching normalisation and channel probabilities | scheme validation, T6 | not rejected |
| Time-dependent activity in unit bins | T7 | not rejected |
| Photon yields (prompt γ, isomeric γ, X-ray after EC) | T3, T8 | not rejected (see the observation in Section 2.4) |
| Ensemble law and independence between replicate streams | T10a, T10b, D5 | not rejected; 400 distinct seeds and stream digests |
| Timestamp, lineage, (ΔZ, ΔA) and neutron bookkeeping; final states | D2 (`decay_invariants`) | holds in every run |
| Exact consumption: N·3·D words per run | D1 | holds in every run |
| Two engines agree bit for bit; exact replay | D3, D4 | hold |
| Two analytic routes agree | A1 | holds |

**Not verified here:** whether any test detects RNG defects, or how fast (Phases E and F); the aggregated (macroscopic) mode, which does not exist yet; Layers B and C; any agreement with nuclear data (the schemes are synthetic by design); behaviour at sample sizes or depths other than those declared.

| Status category | Layer A event engines |
|---|---|
| Implemented | yes |
| Tested | 115 unit tests in their own namespace (`tests/v020/nuclear`) |
| Independently verified | against analytic laws by the frozen study (R1 PASS, R2 MET); code and protocol reviewed twice by an independent reviewer before the freeze |
| Experimentally validated | no, and not applicable: the schemes have no experimental counterpart |

## 2. Results

### 2.1 Deterministic checks

All hold (D1–D5, A1). The largest |matrix exponential − Bateman path sum| is 5.33 × 10⁻¹⁵ (tolerance 10⁻¹¹), and the largest generator column sum is 2.08 × 10⁻¹⁷. XREF: 5,000 nuclei, 15,491 events in each engine, identical discrete columns, final states and stream digests, largest relative time difference 0. REPLAY with chunk 4,099 reproduces the NC1 S2 events digest and stream digest. NC1 S2 read exactly 1,200,000 words (9,600,000 bytes), as the layout requires.

### 2.2 Test families (Holm, family-wise α = 0.01)

| Family | Generator | Tests | Smallest p | Smallest Holm-adjusted p | Rejections |
|---|---|---:|---:|---:|---|
| NC1 | numpy PCG64 | 47 | 0.0017 (S1 T3) | 0.080 | none |
| NC2 | xoshiro256++ (pure Python) | 33 | 0.0074 (S1 T3) | 0.245 | none |

Ensemble (part of NC1): the pooled-count tests T10a have p from 0.31 to 0.94; the between-replicate variance ratios are 0.978–1.092 for the seven eligible states (smallest T10b p 0.20). Every test with its statistic, p-value and Holm-adjusted p-value is listed in `results/tables.md`.

### 2.3 Positive controls (level p < 0.01/47 = 2.13 × 10⁻⁴)

| Control | Perturbation | Target test | p | Rejected |
|---|---|---|---:|---|
| PC1 | half-life 1.02 instead of 1.0 | T2 at t = 2 | 9.48 × 10⁻⁹ | yes |
| PC2 | branching 0.62/0.38 instead of 0.64/0.36 | T6 Xa-400 | 2.54 × 10⁻³⁴ | yes |

The rest of the battery, outside any family, reacts where each perturbation acts and not elsewhere:

- PC1 is rejected by T1 and by T2 at t = 1, 2, 4 and 6. It is not rejected by T2 at t = 0.5 (p = 0.0018, where the survival difference is smallest) or by T3 (p = 0.68), since the photon yield was not changed.
- PC2 is rejected by T6 Xa-400, by T5 for every daughter, by T4 from t = 1 on, and by T7 for Xe-396. It is not rejected by tests that the branching of Xa-400 cannot affect: T5 and T7 for Xa-400 itself (p = 0.93 and 0.40), T6 Xe-396 (p = 0.88) and the five photon tests (p ≥ 0.19). T4 at t = 0.5, when few daughters exist yet, gives p = 0.003, which is not below the level.

### 2.4 Observation, after the run, not part of any rule

The smallest p-value of both families comes from the same test, S1 T3 (photons of the β⁻ channel, yield 0.6). The deviations have opposite signs: z = +3.14 with PCG64 (118,874 photons from 196,986 decays, expected 118,191.6) and z = −2.68 with xoshiro256++ (117,550 from 196,887, expected 118,132.2). A third, independent PCG64 stream with the same yield (PC1) gives z = −0.42, and the ten photon tests of S2 in NC1 and NC2 all have |z| ≤ 1.24. A systematic yield bias would push both families the same way, so this pattern is most consistent with chance. Neither value is a rejection after Holm. It is recorded rather than explained away; a fresh-namespace replication of T3 would settle it if it matters.

### 2.5 Descriptive (R3): transient equilibrium

From NC1 S3, Yb-410/Ya-410 decays in (20, 40]: 20,807/18,705 = 1.1124, approximate 95% interval [1.1059, 1.1189]. The standard error of the log ratio is 0.00296, and the per-nucleus window indicators correlate at 0.894. The exact expectation and the equilibrium limit λ_b/(λ_b − λ_a) are both 1.111111.

### 2.6 Single-machine timings (wall clock; not HPC measurements)

On the author's workstation (one Python process, Python 3.14.9, NumPy 2.5.3): the vectorised engine simulated 100,000 nuclei of S2 (309,469 events) in 0.049 s; the pure-Python xoshiro256++ source made the same run take 0.96 s; the 400 ensemble runs took 0.64 s; the whole study took 3.9 s. These figures characterise one machine and are not inputs to an HPC estimate until profiling (Phase G).

## 3. Numerical tolerances and known discrepancies

| Item | Value | Where it comes from | Effect |
|---|---|---|---|
| Uniform resolution | multiples of 2⁻⁵³ in [0, 1) | `streams.words_to_uniforms` (same as `RandomSource.random_float`) | the sampled waiting-time law differs from the exponential by at most 2⁻⁵³ in distribution function; waiting times are capped at 53·ln 2/λ (about 36.7 mean lives) |
| Branching normalisation | ratios must sum to 1 within 10⁻¹², never renormalised | `decay.BRANCHING_TOLERANCE` | sampler and analytic references use the same law, `Nuclide.channel_probabilities()`; channel probabilities are exact to within 2⁻⁵³ |
| Half-life range | finite half-lives in [10⁻³⁰⁰, 10³⁰⁰] | `decay.HALF_LIFE_RANGE` | rates and waiting times stay finite |
| Closed-form Bateman | refuses rates closer than 10⁻⁶ relative, subnormal denominators, and any evaluation whose rounding-error bound exceeds 10⁻¹⁰ | `bateman.DEGENERACY_RTOL`, `ROUNDING_LIMIT` | largest bound over every path of the three schemes 1.6 × 10⁻¹³ (`design.json`, `bateman_rounding_bound_max`); near-degenerate chains need the matrix exponential |
| Two analytic routes | must agree within 10⁻¹¹ (A1) | protocol | observed 5.33 × 10⁻¹⁵ |
| Engine agreement | discrete columns exact; times within 10⁻¹² relative (D3) | protocol | observed exact equality |
| χ² tests | Pearson χ², cells pooled below expected count 5 | `checks.multinomial_check` | p-values are asymptotic approximations |
| KS test | scipy `kstest`, automatic method | `checks.ks_check` | asymptotic p-value at these sample sizes |
| T10b | normal-theory χ²(R − 1) for a binomial variance | protocol | standard deviation of the statistic off by at most 0.87% (`design.json`) |
| R3 interval | delta method with the per-nucleus covariance | protocol | approximate; descriptive only |

## 4. Review before the freeze

Before anything was drawn from the study namespace, an independent reviewer read the code and the draft protocol and reported six defects. One was HIGH: the closed-form Bateman sum accepted rate sets for which it cancelled catastrophically. Five were MEDIUM or LOW: non-finite values could read as passing checks, the R3 interval assumed independent counts and was about 3.4 times too wide, a rate could overflow, a channel could be rounded away by the thresholds, and D1 did not cover XREF and REPLAY. Every numerical claim was reproduced before it was fixed. A confirmation review found two of the fixes incomplete (subnormal products, an invalid Holm level) and one new minor defect; those were fixed and covered by regression tests that reproduce the review's inputs. The list is in `research/v020-nuclear-a-verify/PROTOCOL.md`, section 7. None of the corrections changed what the sampler does.

## 5. Limitations and open questions

- The schemes are synthetic. Agreement with analytic laws says nothing about agreement with nuclear data.
- The check is algorithmically independent, but not organisationally: the samplers, the references and the tests come from one project. Two analytic routes, two engines, two generators and the outside review reduce the risk of a shared conceptual error without removing it. An independent recount from the stored counts is planned (Phase H).
- Three schemes of depth at most 4 were verified. Deeper or near-degenerate schemes rely on the matrix exponential alone.
- Only healthy generators were used. Whether this battery, or any application observable, detects RNG defects is the subject of Phases E and F.
- The S1 T3 observation (Section 2.4) is unresolved in the strict sense; a replication in a fresh namespace would close it.
- Timings come from one workstation.

## 6. Reproduction

```bash
PYTHONPATH=src .venv/bin/python research/v020-nuclear-a-verify/run_verification.py   # refuses: results exist
PYTHONPATH=src .venv/bin/python research/v020-nuclear-a-verify/report_tables.py      # rewrites tables.md from the results
```

A replication belongs in a new namespace under a new protocol version; the stored results are never overwritten.
