# v0.20 nuclear module, milestone B2: verification of the aggregated decay engine

**Status: DRAFT, not frozen (2026-11-17), pending independent review of the hardened revision.** It becomes frozen when `protocol.sha256` is committed with the SHA-256 of `protocol.json`, `design.json`, `design.py` and `run_verification.py`; the script refuses to draw if any of them differs. That file is deliberately absent, so the study cannot be started. `protocol.json` governs. Nothing has been drawn from the study namespace `rnglab-v020-nuclear-b2-verify`. Smoke runs in their own namespace (`smoke/`, 5% sizes, not in Git) exercised the code and are excluded from every result. The Layer A study, its protocol, results and code are unchanged; this study only runs the frozen Layer A engine in its own namespace.

## 1. What this study is, and what it is not

- **It is** a check that the B2 aggregated decay engine (`src/rnglab/nuclear/aggregate.py`, with `binomial.py`, `paths.py`, `moments.py`, `comparison.py` and `precursors.py`) reproduces exactly known decay laws, including variances and covariances, and agrees statistically with the frozen Layer A event engine on every observable both can produce. It also checks that the battery notices a small physics error and a plausible aggregation error.
- **It is not** an RNG-sensitivity experiment, not evidence about branching processes, criticality or detectors, not a test of H1, H2 or H3, and not part of the Stage D pilot. Equivalence with Layer A is tested for healthy generators only. The two engines use random numbers in completely different ways, so how a defective generator distorts each of them is a separate question.
- **Negative outcomes** are reported with the same prominence as positive ones. The study runs once.

## 2. What B2 computes

Over an interval of length Δ every nucleus follows one decay-path category of its start nuclide. The categories' probabilities are closed-form Bateman terms evaluated in decimal arithmetic with adaptive precision (`paths.py`). The category counts of the n_i nuclei in nuclide i are Multinomial(n_i, q_i(Δ)), sampled as conditional binomials by exact inversion with one uniform per draw (`binomial.py`). Photons are Binomial(decays, yield); neutrons equal the β-delayed-neutron decays. Counts above 2⁵³ are split into parts, each with its own uniform. The layout `nuclear-aggregate-layout-v1` reads a fixed number of words per interval, independent of the population.

A uniform u = k·2⁻⁵³ stands for the cell [u, u + 2⁻⁵³); the binomial is inverted at the cell's midpoint u + 2⁻⁵⁴. Inverting at u itself would return 0 for u = 0, an outcome of probability 2⁻⁵³ that at n = 2⁵³, p = ½ would by itself double the variance (`research/v020-nuclear-b2-numerics`, `midpoint_convention`).

## 3. Schemes, cases and runs

S1, S2 and S3 are the frozen Layer A schemes (identical digests, checked by `design.py`). S4 holds three synthetic delayed-neutron precursor groups (half-lives 0.6, 2.5, 9; neutron probabilities 0.3, 0.1, 1). S5 is a stiff chain with a rate ratio of 10⁶ and two decay constants 10⁻³ apart. All nuclides are synthetic (Z ≥ 145).

| Run | Content |
|---|---|
| NC1, NC2 | S1–S5 with 10⁵ nuclei each (S4: 40,000 / 40,000 / 20,000) on non-uniform grids; numpy PCG64 and the pure-Python xoshiro256++ |
| ENS-B2, ENS-A | S2: 400 replicates of 2,000 nuclei in B2 and in the frozen Layer A engine |
| BIG | S2 with 10¹² nuclei (one part) and 2⁶⁰ + 12,345 nuclei (129 parts) |
| PC-B1 | S1 with half-life 1.02 instead of 1.0, 400,000 nuclei |
| PC-B2 | a deliberately wrong one-step aggregation (at most one decay per nucleus and interval) of S2: one run of 10⁵ nuclei and 400 replicates of 2,000 |
| XREF, REPLAY, CHECKPOINT, SINK | software checks: two engines on one stream; rerun; resume from a JSON checkpoint; streaming without kept outputs |

## 4. Tests

Every test, its null law and its method are in `protocol.json` (`tests`). `design.json` (`planned_tests`) names every test that runs: 215 in NC1 and 78 in NC2. Which tests run never depends on the data.

- **Against analytic laws (NC1, NC2).** State counts at every grid time (multinomial χ²); decays of each nuclide over the grid (exact binomial); decays per interval (multinomial χ²); channel splits and photon yields (exact binomial); for S4, neutrons per interval. For S2, one path test per source nuclide: the Pearson statistics of the planned intervals' path-category counts against the engine's decimal interval laws, on cells fixed in `design.json` (minimum expected count 10), summed and referred to χ² with the summed degrees of freedom. Counts of different intervals are not pooled, because later start counts depend on earlier outcomes.
- **Ensemble moments (ENS-B2).** Pooled means (exact binomial); sample variances of every eligible state at t = 2 and t = 8; sample variances of sums N_j + N_k of every eligible pair at t = 2 and of 15 two-time pairs (t = 2 with t = 8), which carry the covariances. The variance of each sample variance is exact (`moments.py`, from the per-nucleus categorical law), and the sample variance is referred to a scaled χ² matching its exact mean and variance, which keeps its right skew.
- **Layer A against B2 (ENS-A, ENS-B2).** Homogeneity χ² of pooled state counts at every grid time, of decays per interval and of photon yields; ratios of the two engines' sample variances, of every eligible state and of every eligible pair sum at t = 2, against F(ν, ν) with ν from the exact moments (the common scale cancels).
- **Huge populations (BIG).** State counts against the multinomial law, and decays against their expectations (normal approximation).

## 5. Deterministic checks (any failure is a FAIL)

| | Check |
|---|---|
| D1 | every B2 run reads intervals × words-per-interval words (8 bytes each); every Layer A run reads 3·D words per nucleus |
| D2 | every B2 run: populations conserved and nonnegative, changes equal inflow minus outflow, outflow bounded, photons between 0 and the channel's decays (equal for yield 1, zero for yield 0), neutrons equal the β-delayed-neutron decays; every Layer A run satisfies `decay_invariants` |
| D3 | the two B2 engines agree bit for bit on the same stream, including a forced split into 101 parts |
| D4 | replay reproduces the digest and the stream digest |
| D5 | a run resumed from a JSON checkpoint after interval 3 reproduces the digest and the remaining outputs |
| D6 | streaming with `keep_outputs=False` and a sink reproduces the digest |
| D7 | each ensemble has 400 distinct seeds and stream digests; no Layer A decay at time 0 |
| D8 | interval laws, no random numbers: the decimal law, the matrix exponential on the path tree and the Layer A matrix exponential agree within 10⁻¹², probabilities sum to 1 within 10⁻³⁰ and lie in their mean-value brackets |
| D9 | the BIG runs use 1 and 129 parts and conserve the population |
| D10 | exactly the tests named in `design.json` ran in each family |

## 6. Families, controls and reading rules

- **Families.** NC1: every test on the NC1 runs, ENS-B2, the Layer A / B2 comparison and BIG (215 tests). NC2: every test on the NC2 runs (78). Holm within each family at family-wise α = 0.01.
- **R1.** PASS iff D1–D10 hold and Holm rejects no test in NC1 and none in NC2. Otherwise FAIL, with every failed item listed.
- **How often a correct engine would FAIL.** Holm bounds the chance of a false rejection by 1% per family if every p-value is valid. The exact binomial tests are; the χ², scaled-χ² and normal tests are approximations. `calibration.py` estimates their null rejection rates down to the Holm level by sampling the exact laws with numpy, independently of both engines; its results are in `calibration.json` and summarised in section 8.
- **R2.** PC-B1 must reject its target (state counts at t = 2), and PC-B2 both of its targets (state counts at t = 12 of the single run; homogeneity at t = 12 between the naive ensemble and ENS-A), each at p < 0.01/215 = 4.65 × 10⁻⁵. `design.json` gives approximate power 1.0 for all three. If any is missed, R1 is reported with the statement that the battery lacked the declared sensitivity.
- **R3** (descriptive): run times, digits used by the decimal laws, the largest decimal-versus-matrix-exponential difference, and the parts and times of the BIG runs.
- The study runs once from a clean committed tree. No threshold, size or test changes after generation. A code change after the freeze is an amendment written before any affected number is regenerated.

## 7. Changes made before the freeze

From smoke runs and the numerical investigation (`research/v020-nuclear-b2-numerics`):

1. The binomial sampler evaluates its distribution function in two regimes (recurrence table for means up to 256, Boost's incomplete beta above), because Boost's incomplete beta has absolute errors up to 7.4 × 10⁻¹² (F) and 2.3 × 10⁻⁸ (its complement) for small means at 10⁶ ≤ n ≤ 10⁹.
2. The upper tail is consulted only where relative accuracy is needed (1 − u < 2⁻¹⁰) and only up to n = 10⁹, because Boost's upper tail loses accuracy near the centre at larger n.
3. The interval laws use the exact differences of the sampler's float thresholds as branching weights, so that they sum to exactly 1.
4. D8's tolerance is 10⁻¹², not 10⁻¹³: on S5 the float64 matrix exponential itself errs by 1.6 × 10⁻¹³, while the decimal law agrees with a 150-digit evaluation to 10⁻⁴⁶.

From a first independent review of the code and of the draft protocol (every reported defect was reproduced before it was corrected, and each correction has a regression test):

5. **Floating-point boundaries (HIGH).** The sampler inverted at u itself, so u = 0 returned 0. It now inverts at the cell midpoint (section 2). Adversarial tests compare it with exact integer arithmetic at the extreme uniforms, one ulp either side of distribution-function steps, both sides of the two regime switches, at 2⁵³ and for vanishing probabilities.
6. **Path tests.** Pooling path counts over intervals of equal length assumed a law that does not hold, because later start counts depend on earlier outcomes. The test is now a sum of per-interval statistics on pre-declared cells.
7. **Moment tests.** The normal approximation for sample variances understated their right tail at the Holm level. Variances are now referred to a scaled χ², and covariances are tested through variances of sums. For the Layer A / B2 comparison a z-test of the difference of two sample variances was tried and rejected: the calibration measured its null rejection rate at the Holm level at 1.66 times nominal (`calibration.json`, row "CT-VAR z-difference (not used)"). The comparison now uses variance-ratio F-tests.
8. **Neutron totals** could overflow 2⁶³ − 1 for schemes with several β-delayed-neutron steps. Such configurations are refused before any random number is read.
9. **State aliasing.** Sinks received live arrays and could alter a run without changing its digest. A sink now receives detached copies; the time grid is copied; plan arrays and returned arrays are read-only; checkpoints are built from fresh objects. Adversarial tests force writes on everything a caller can reach and require identical digests.
10. **Checkpoint integrity.** Checkpoints with fractional or mistyped fields were silently repaired. They are now refused; a checkpoint also carries the SHA-256 of its own contents and the identity of its stream (generator, seed, absolute word position), so an edited checkpoint and a stream at another position or with another seed are refused. The checksum is an integrity check, not an authentication: a forgery with a recomputed checksum resumes, and is found by comparing the final digest with an uninterrupted run (tested).
11. **Rules.** D10 was enforced by the script but missing from `protocol.json`, and this document had stale family sizes. All three now agree.
12. **Design.** The noncentrality of the PC-B2 comparison was halved by mistake (both powers round to 1 either way). Scalar inputs to the sampler raised an error.

## 8. Reviews and calibration before the freeze

- **Review status.** The corrections of section 7, items 5 to 12, and the adversarial tests were made after the first independent review. This hardened revision has not yet been reviewed independently. The study must not be run before such a review has passed; its outcome will be recorded here before the freeze.
- **Calibration (development run at this revision, to be repeated for the record before the freeze).** Null rejection rate divided by the nominal level, at the Holm level 4.65 × 10⁻⁵, with exact 95% intervals: ET-VAR 0.83 [0.64, 1.07]; ET-VARSUM 1.09 [0.89, 1.31]; ET-VARSUM2 1.00 [0.72, 1.36]; CT-VAR 1.38 [1.01, 1.85]; CT-VARSUM 0.96 [0.78, 1.18]; BT-PATH 1.13 [1.01, 1.26]; BT-STATE 1.02 [0.72, 1.40]; CT-HOMO 0.91 [0.53, 1.46]. The rejected z-difference alternative: 1.66 [1.25, 2.16]. Open point for the review: CT-VAR (seven tests) is above nominal; if all seven were at 1.38 the family-wise rate of NC1 would rise from 0.0100 to about 0.0101. Not calibrated: BT-BINS, BT-NEUTRON, CT-DECAY, CT-PHOTON (χ² tests with large cells) and BIG-DECAY (normal approximation).

## 9. Reproduction

```bash
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/design.py
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/calibration.py
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/run_verification.py --smoke
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/run_verification.py      # once, clean tree
```
