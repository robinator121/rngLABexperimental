# v0.20 nuclear module, Layer A: verification against analytic decay laws

**Status: FROZEN before generation, 2026-11-16.** `protocol.json` holds every parameter and reading rule and governs. `protocol.sha256` holds the SHA-256 of `protocol.json`, `design.json`, `design.py` and `run_verification.py`; the script refuses to draw if any of them differs. Nothing has been drawn from the study namespace `rnglab-v020-nuclear-a-verify`. A smoke run in its own namespace (`smoke/`, 5% sizes, not in Git) exercised the code and is excluded from every result.

## 1. What this study is, and what it is not

- **It is** a check that the Layer A event-level decay simulator (`src/rnglab/nuclear`) reproduces exactly known decay laws when it is driven by healthy generators, and that the battery used for that check is sensitive enough to notice small physics errors.
- **It is not** an RNG-sensitivity experiment, not evidence about branching processes, criticality or detectors, not a test of H1, H2 or H3, and not part of the Stage D pilot. A PASS shows agreement with the analytic laws for the three declared synthetic schemes and sample sizes; it does not make the module validated physics software.
- **Negative outcomes** are reported with the same prominence as positive ones. The study runs once.

## 2. Schemes and runs

All nuclides are synthetic, with Z ≥ 145, outside the chart of known nuclides. Time is in abstract units.

| Scheme | Content | Initial nuclei | Horizon |
|---|---|---|---|
| S1 | Sa-400 (T½ = 1, β⁻, photon yield 0.6) → stable | 200,000 | 6 |
| S2 | nine states: α/β⁻ split (0.64/0.36), electron capture, an isomeric transition, β⁻ and β-delayed-neutron split (0.9/0.1), three stable ends; longest path 4 decays | 100,000 | 12 |
| S3 | transient equilibrium: Ya-410 (T½ = 10) → Yb-410 (T½ = 1) → stable | 100,000 | 40 |

| Run | Generator | Purpose |
|---|---|---|
| NC1 | numpy PCG64 (`NumpyPCG64Source`) | negative control, primary healthy generator: S1, S2, S3 |
| NC2 | xoshiro256++ (`rnglab14`, pure Python) | negative control, a second, independently implemented healthy generator |
| ENS | numpy PCG64 | 400 replicates of S2 with 2,000 nuclei each, horizon 2: ensemble law and independence between replicate streams |
| PC1 | numpy PCG64 | positive control: S1 simulated with half-life 1.02 and tested against the reference 1.0 |
| PC2 | numpy PCG64 | positive control: S2 simulated with branching 0.62/0.38 and tested against the reference 0.64/0.36 |
| XREF | numpy PCG64 | software check: vectorised engine (chunk 997) against the pure-Python reference engine on the same stream, S2, 5,000 nuclei |
| REPLAY | numpy PCG64 | software check: NC1 S2 again with chunk 4,099 |

Seeds: `namespace_seed("rnglab-v020-nuclear-a-verify", run, scheme, "physics")`, and for ENS replicate r `namespace_seed(namespace, "ENS", "S2", "physics", f"{r:04d}")`. Layout `nuclear-decay-layout-v1`: three words per decay step, 3·D words per nucleus.

## 3. Tests and their exact laws

Nuclei are independent and the schemes are acyclic, so every law below is exact. The probabilities come from the matrix exponential of the decay generator (`rnglab.nuclear.bateman`), which draws no random numbers.

| Test | Scheme | Observable and law | Method |
|---|---|---|---|
| T1 | S1 | decay times against the exponential truncated at the horizon | Kolmogorov–Smirnov |
| T2 | S1 | Sa-400 count at t = 0.5, 1, 2, 4, 6 against Binomial(N, 2⁻ᵗ) | exact binomial |
| T3 | S1 | photons given the decays: Binomial(decays, 0.6) | exact binomial |
| T4 | S2 | counts in the nine states at t = 0.5, 1, 2, 4, 8, 12 against Multinomial(N, p(t)) | χ², pooling rule |
| T5 | S2 | decays of each of the six unstable nuclides within the horizon against Binomial(N, q_j) | exact binomial |
| T6 | S2 | first-channel count of Xa-400 and of Xe-396 given their decays | exact binomial |
| T7 | S2 | decay times of Xa-400 and of Xe-396 in twelve unit bins plus "no decay by 12" | χ², pooling rule |
| T8 | S2 | photons of the five channels with yield strictly between 0 and 1, given the channel's decays | exact binomial |
| T9 | S3 | counts in the three states at t = 1, 5, 10, 20, 40 | χ², pooling rule |
| T10a | ENS | count of each eligible state at t = 2 summed over 400 replicates: Binomial(800,000, p_j) | exact binomial |
| T10b | ENS | between-replicate variance of each eligible count against 2,000 p_j(1 − p_j) | χ²(399), two-sided; a normal-theory approximation |

Pooling rule (fixed in `rnglab.nuclear.checks.multinomial_check`): cells with expected count below 5 are pooled; a pooled cell still below 5 joins the smallest remaining cell. A count in a state of probability 0 gives p = 0. Eligible states for ENS are those with 2,000 p_j(1 − p_j) ≥ 20; `design.json` lists the seven, computed before generation. The normal-theory approximation in T10b ignores the binomial excess kurtosis (1 − 6pq)/(npq); `design.json` (`ens.T10b_normal_theory_error`) gives its largest value over the eligible states and the resulting change in the standard deviation of the statistic.

## 4. Deterministic checks (any failure is a FAIL)

- **D1:** each run, including both XREF engines, REPLAY and every ENS replicate, reads exactly N·3·D words, and its source's byte count is eight times that.
- **D2:** every NC1, NC2, PC1, PC2, XREF and ENS run satisfies every item of `decay_invariants`: times within the horizon, identifiers fixed by the layout, lineage links consistent and non-decreasing in time, first events leaving the initial state, channel defines daughter and mode, (ΔZ, ΔA) consistent with the mode, isomeric transitions lowering the level, neutrons only from the β-delayed-neutron channel, photon yields 0 and 1 respected, final states equal to the last daughters, populations at the horizon equal to the final states.
- **D3:** XREF engines have identical discrete columns, final states and stream digests, and a largest relative time difference of at most 10⁻¹².
- **D4:** REPLAY reproduces the NC1 S2 events digest and stream digest.
- **D5:** ENS uses 400 distinct seeds and produces 400 distinct stream digests.
- **A1 (no random numbers):** the matrix exponential and the closed-form Bateman path sum agree within 10⁻¹¹ for every scheme, start state and test time; generator column sums are within 10⁻¹⁵ of zero.

## 5. Families, controls and reading rules (fixed before generation)

- **Families.** NC1 = T1–T9 on NC1 plus T10a and T10b (47 tests); NC2 = T1–T9 on NC2 (33 tests). Holm step-down within each family at family-wise α = 0.01.
- **R1, agreement with the analytic laws.** PASS if and only if D1–D5 and A1 hold and Holm rejects no test in NC1 and none in NC2. Otherwise FAIL, with every failed item listed. Under correct code the chance of a false FAIL is at most about 2% (two families at 1%). A FAIL is reported as found and not explained away.
- **R2, sensitivity.** PC1 must reject its target (T2 at t = 2) and PC2 its target (T6 Xa-400), each at p < 0.01/47 = 2.13 × 10⁻⁴. `design.json` gives approximate power 0.9997 and 1.0 (normal approximation). If either control fails, R1 is reported together with the statement that the battery lacked the declared sensitivity. The other tests on PC1 and PC2 are reported descriptively, outside any family.
- **R3, descriptive, no verdict.** From NC1 S3: the ratio of Yb-410 to Ya-410 decays in (20, 40] with an approximate 95% interval (delta method on the log ratio; the two counts come from the same nuclei, so the variances and the covariance of the per-nucleus window indicators are estimated from the run), next to its exact expectation and to the transient-equilibrium limit λ_b/(λ_b − λ_a) = 10/9.
- **General.** The main run starts from a clean committed tree and refuses to overwrite results. No threshold, size or test changes after generation. A code change after the freeze is an amendment written before any affected number is regenerated.

## 6. Reproduction

```bash
PYTHONPATH=src .venv/bin/python research/v020-nuclear-a-verify/design.py              # analytic design numbers
PYTHONPATH=src .venv/bin/python research/v020-nuclear-a-verify/run_verification.py --smoke
PYTHONPATH=src .venv/bin/python research/v020-nuclear-a-verify/run_verification.py      # once, clean tree
```

Outputs: `results/verification.json` (every count, statistic, p-value, invariant, timing and verdict) and `results/manifests/` (one manifest per run, one for the ENS ensemble).

## 7. Corrections made before the freeze

A second, independent review of the code and of the draft of this protocol, made before anything was drawn from the study namespace, found six defects. Each was reproduced, corrected and covered by a regression test before the freeze:

1. The closed-form Bateman sum accepted close rates whose alternating sum cancels catastrophically (rates 1, 1.00001, …, 1.00004 gave −1280 instead of 0.0153). It now refuses any evaluation whose rounding-error bound exceeds 10⁻¹⁰, or whose result is not finite.
2. Non-finite values could read as passing checks (a NaN probability gave p = 1; Holm turned a NaN p-value into 1; A1 ignored NaN). Every check now raises on non-finite input or output, and A1 records a refused evaluation as a failure.
3. The R3 interval treated the parent and daughter counts as independent, although they come from the same nuclei. It was about 3.4 times too wide. It now uses the per-nucleus covariance.
4. A half-life of 10⁻³¹⁰ gave an infinite decay constant. Finite half-lives are now restricted to [10⁻³⁰⁰, 10³⁰⁰].
5. Branching ratios such as [1.0, 5 × 10⁻¹³] passed the normalisation tolerance, although the sampler could never choose the second channel while the analytic references kept it. One branching law (`Nuclide.channel_probabilities`) is now used everywhere, and a channel left without probability is refused.
6. D1 did not check the consumption of the XREF engines and of REPLAY. It now does.

A confirmation review of these corrections found two of them incomplete and one new minor defect: for tiny, close rates (about 10⁻⁷⁵, relative spacing 2 × 10⁻⁶) subnormal intermediate products made the rounding-error bound itself invalid (a result of −15.6 passed); Holm accepted a NaN family-wise level; and a stable nuclide reported a channel probability. The closed form is now evaluated on rates divided by the largest rate (the result is scale-invariant), refuses subnormal denominators, and is tested against the matrix exponential over rate scales from 10⁻²⁰⁰ to 10²⁰⁰ (each evaluation is refused or within 10⁻¹⁰); Holm refuses a level outside (0, 1); stable nuclides report no channel probabilities. These last corrections are covered by regression tests that reproduce the review's inputs; they were not reviewed a third time.

None of these changes what the sampler does: its thresholds were already the cumulative sums that `channel_probabilities` now exposes. The analytic references moved only in the last floating-point digit: `design.json` changed in the last digits of a few entries, and its eligibility list and family sizes are unchanged.
