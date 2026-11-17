# v0.20 nuclear module, milestone B2: verification of the aggregated decay engine

**Status: DRAFT, not frozen (2026-11-17), pending independent review of the corrected revision.** It becomes frozen when `protocol.sha256` is committed with the SHA-256 of `protocol.json`, `design.json`, `design.py`, `run_verification.py` and `calibration.py`; both scripts refuse to draw for the record if any of them differs. That file is deliberately absent, so the study cannot be started. `protocol.json` governs. Nothing has been drawn from the study namespace `rnglab-v020-nuclear-b2-verify` or the record calibration namespace `rnglab-v020-nuclear-b2-calibration-record`. Smoke runs (`smoke/`, 5% sizes, not in Git) and development calibrations (namespace `rnglab-v020-nuclear-b2-calibration-dev`) exercised the code and are excluded from every result. The Layer A study, its protocol, results and code are unchanged; this study only runs the frozen Layer A engine in its own namespace.

## 1. What this study is, and what it is not

- **It is** a check of the B2 aggregated decay engine (`src/rnglab/nuclear/aggregate.py`, with `binomial.py`, `paths.py`, `moments.py`, `comparison.py` and `precursors.py`): whether its outputs can be told apart, at the stated sample sizes, from exactly known decay laws (counts, decays, channel splits, photons, neutrons, path categories, and the variances and covariances listed in section 4), and from the frozen Layer A event engine on the observables listed in section 4. It also checks that the battery notices a small physics error and a plausible aggregation error.
- **Not rejecting is not equivalence.** Section 4 states what is tested and how large a deviation the tests detect; deviations below that size, and features of the joint law that are not tested, can pass.
- **It is not** an RNG-sensitivity experiment, not evidence about branching processes, criticality or detectors, not a test of H1, H2 or H3, and not part of the Stage D pilot. The comparison with Layer A uses healthy generators only. The two engines use random numbers in completely different ways, so how a defective generator distorts each of them is a separate question.
- **Negative outcomes** are reported with the same prominence as positive ones. The study runs once.

## 2. What B2 computes

Over an interval of length Δ every nucleus follows one decay-path category of its start nuclide. The categories' probabilities are closed-form Bateman terms evaluated in decimal arithmetic with adaptive precision (`paths.py`). The category counts of the n_i nuclei in nuclide i are Multinomial(n_i, q_i(Δ)), sampled as conditional binomials by inversion with one uniform per draw (`binomial.py`). Photons are Binomial(decays, yield); neutrons equal the β-delayed-neutron decays. Counts above 2⁵⁰ are split into parts of at most 2⁵⁰ trials, each with its own uniform: Boost's incomplete beta returns NaN near the centre of the distribution when n is close to 2⁵³ (section 7, item 13). The layout `nuclear-aggregate-layout-v1` reads a fixed number of words per interval, independent of the population.

A uniform u = k·2⁻⁵³ stands for the cell [u, u + 2⁻⁵³); the binomial is inverted at the cell's midpoint u + 2⁻⁵⁴. Inverting at u itself would return 0 for u = 0, an outcome of probability 2⁻⁵³ that at n = 2⁵³, p = ½ would by itself double the variance (12.5% at the present part limit). What the draw realises is stated in `binomial.py`: with an exact distribution function it is within 2⁻⁵⁴ of Binomial(n, p) in Kolmogorov distance, outcomes much rarer than 2⁻⁵³ can be lost (item 17), and the distribution function is evaluated in floating point, so a draw can differ from exact inversion next to a step (item 18).

## 3. Schemes, cases and runs

S1, S2 and S3 are the frozen Layer A schemes (identical digests, checked by `design.py`). S4 holds three synthetic delayed-neutron precursor groups (half-lives 0.6, 2.5, 9; neutron probabilities 0.3, 0.1, 1). S5 is a stiff chain with a rate ratio of 10⁶ and two decay constants 10⁻³ apart. All nuclides are synthetic (Z ≥ 145).

| Run | Content |
|---|---|
| NC1, NC2 | S1–S5 with 10⁵ nuclei each (S4: 40,000 / 40,000 / 20,000) on non-uniform grids; numpy PCG64 and the pure-Python xoshiro256++ |
| ENS-B2, ENS-A | S2: 400 replicates of 2,000 nuclei in B2 and in the frozen Layer A engine |
| BIG | S2 with 10¹² nuclei (one part) and 2⁶⁰ + 12,345 nuclei (1,025 parts) |
| PC-B1 | S1 with half-life 1.02 instead of 1.0, 400,000 nuclei |
| PC-B2 | a deliberately wrong one-step aggregation (at most one decay per nucleus and interval) of S2: one run of 10⁵ nuclei and 400 replicates of 2,000 |
| XREF, REPLAY, CHECKPOINT, SINK | software checks: two engines on one stream; rerun; resume from a JSON checkpoint; streaming without kept outputs |

## 4. Tests

Every test, its null law and its method are in `protocol.json` (`tests`). `design.json` (`planned_tests`) names every test that runs: 222 in NC1 and 85 in NC2. Which tests run never depends on the data: pooling uses fixed totals and exact probabilities, and every other cell (path tests, per-interval photon tests, homogeneity tests) is fixed in `design.json` from expected counts. A planned test that cannot be computed is recorded with its error and is a failed item.

- **Against analytic laws (NC1, NC2).** State counts at every grid time (multinomial χ²); decays of each nuclide over the grid (exact binomial); decays per interval (multinomial χ²); channel splits and photon yields (exact binomial); photons per interval given the decays per interval, on cells of consecutive intervals (χ², BT-PHOTON-T), which tests how photons are spread over time; for S4, neutrons per interval. For S2, one path test per source nuclide: the Pearson statistics of the planned intervals' path-category counts against the engine's decimal interval laws, on cells fixed in `design.json` (minimum expected count 10), summed and referred to χ² with the summed degrees of freedom. Counts of different intervals are not pooled, because later start counts depend on earlier outcomes.
- **Ensemble moments (ENS-B2).** Pooled means (exact binomial); sample variances of every eligible state at t = 2 and t = 8; sample variances of sums N_j + N_k of every eligible pair at t = 2 and of 15 two-time pairs (t = 2 with t = 8), which carry the covariances. The variance of each sample variance is exact (`moments.py`, from the per-nucleus categorical law), and the sample variance is referred to a scaled χ² matching its exact mean and variance, which keeps its right skew.
- **Layer A against B2 (ENS-A, ENS-B2).** Homogeneity χ² on fixed cells of pooled state counts at every grid time, of decays per interval and of photon yields; ratios of the two engines' sample variances, of every eligible state and of every eligible pair sum at t = 2, against F(ν, ν) with ν from the exact moments (the common scale cancels).
- **Huge populations (BIG).** State counts against the multinomial law, and decays against their expectations (normal approximation).
- **Coverage of second moments.** Tested: the variance of every eligible state at t = 2 and t = 8, the covariances of eligible pairs at t = 2, and the two-time covariances of the 15 pairs in `design.json`; against Layer A, the variances of eligible states and pair sums at t = 2. Not tested: any other feature of the joint law across times. A deviation that keeps every marginal law and the tested sums (for example, moving probability between cells of the two-time law of t = 2 and t = 8 in a pattern that leaves those sums unchanged) passes. BT-PATH constrains each interval's law given its start counts, which covers dependence carried by the counts but not dependence carried by anything else, such as reused uniforms.
- **Sensitivity.** At the Holm level, a variance that is 10% too large is detected with probability below 0.005 by any one test (scaled χ²) and below 0.001 by the Layer A comparison (F test). Variance ratios detected with probability 0.9 are about 1.45 (ET tests) and 1.71 (CT tests) (`design.json`, `sensitivity`).

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
| D8 | interval laws, no random numbers: every reference is finite; the decimal law, the matrix exponential on the path tree and the Layer A matrix exponential agree within 10⁻¹², probabilities sum to 1 within 10⁻³⁰ and lie in their mean-value brackets |
| D9 | the BIG runs use 1 and 1,025 parts and conserve the population |
| D10 | exactly the tests named in `design.json` ran in each family |

## 6. Families, calibration, controls and reading rules

- **Families.** NC1: every test on the NC1 runs, ENS-B2, the Layer A / B2 comparison and BIG (222 tests). NC2: every test on the NC2 runs (85). Holm within each family at family-wise α = 0.01.
- **R1.** PASS iff D1–D10 hold and Holm rejects no test in NC1 and none in NC2. Otherwise FAIL, with every failed item listed.
- **How often a correct engine would FAIL.** The nominal family-wise error is 1% per family; it is a bound only if every p-value is valid. The exact binomial tests are valid; the χ², scaled-χ², F and normal tests are approximations. Their combined effect is estimated, not bounded, by `calibration.py`.
- **R0, calibration (run once after the freeze, before `run_verification.py`).** 50,000 independent null batteries in the record namespace: every run of NC1 and NC2 is drawn from the exact laws with numpy (no RNGLab sampler), and every planned test is evaluated with the functions `run_verification.py` uses. Criterion C1: for each family, the Clopper–Pearson 95% upper bound of the probability that a correct engine FAILs (Holm rejects a test, or a planned test cannot be computed) is at most 0.02. R0 is reported with R1 whichever way it falls; if C1 fails, an R1 FAIL caused only by approximate tests is reported together with the statement that the battery's false-FAIL rate exceeds its declared limit. No verdict is changed by R0.
- **R2.** PC-B1 must reject its target (state counts at t = 2), and PC-B2 both of its targets (state counts at t = 12 of the single run; homogeneity at t = 12 between the naive ensemble and ENS-A), each at p < 0.01/222 = 4.50 × 10⁻⁵. `design.json` gives approximate power 1.0 for all three. If any is missed, R1 is reported with the statement that the battery lacked the declared sensitivity.
- **R3** (descriptive): run times, digits used by the decimal laws, the largest decimal-versus-matrix-exponential difference, and the parts and times of the BIG runs.
- The study runs once from a clean committed tree. No threshold, size or test changes after generation. A code change after the freeze is an amendment written before any affected number is regenerated.

## 7. Changes made before the freeze

From smoke runs and the numerical investigation (`research/v020-nuclear-b2-numerics`):

1. The binomial sampler evaluates its distribution function in two regimes (recurrence table for means up to 256, Boost's incomplete beta above), because Boost's incomplete beta has absolute errors up to 7.4 × 10⁻¹² (F) and 2.3 × 10⁻⁸ (its complement) for small means at 10⁶ ≤ n ≤ 10⁹.
2. The upper tail is consulted only where relative accuracy is needed (1 − u < 2⁻¹⁰) and only up to n = 10⁹, because Boost's upper tail loses accuracy near the centre at larger n.
3. The interval laws use the exact differences of the sampler's float thresholds as branching weights, so that they sum to exactly 1.
4. D8's tolerance is 10⁻¹², not 10⁻¹³: on S5 the float64 matrix exponential itself errs by 1.6 × 10⁻¹³, while the decimal law agrees with a 150-digit evaluation to 10⁻⁴⁶.

From a first independent review of the code and of the draft protocol (every reported defect was reproduced before it was corrected, and each correction has a regression test):

5. **Floating-point boundaries (HIGH).** The sampler inverted at u itself, so u = 0 returned 0. It now inverts at the cell midpoint (section 2). Adversarial tests compare it with exact integer arithmetic at the extreme uniforms, one ulp either side of floating-point distribution-function steps, both sides of the two regime switches, at the part limit and for vanishing probabilities.
6. **Path tests.** Pooling path counts over intervals of equal length assumed a law that does not hold, because later start counts depend on earlier outcomes. The test is now a sum of per-interval statistics on pre-declared cells.
7. **Moment tests.** The normal approximation for sample variances understated their right tail at the Holm level. Variances are now referred to a scaled χ², and covariances are tested through variances of sums. For the Layer A / B2 comparison a z-test of the difference of two sample variances was tried and rejected (null rejection rate at the Holm level 1.66 times nominal in a development calibration); the comparison uses variance-ratio F-tests.
8. **Neutron totals** could overflow 2⁶³ − 1 for schemes with several β-delayed-neutron steps. Such configurations are refused before any random number is read.
9. **State aliasing.** Sinks received live arrays and could alter a run without changing its digest. A sink now receives detached copies; the time grid is copied; plan arrays and returned arrays are read-only; checkpoints are built from fresh objects. Adversarial tests force writes on everything a caller can reach and require identical digests.
10. **Checkpoint integrity.** Checkpoints with fractional or mistyped fields were silently repaired. They are now refused; a checkpoint also carries the SHA-256 of its own contents and the identity of its stream (generator, seed, absolute word position), so an edited checkpoint and a stream at another position or with another seed are refused. The checksum is an integrity check, not an authentication: a forgery with a recomputed checksum resumes, and is found by comparing the final digest with an uninterrupted run (tested).
11. **Rules.** D10 was enforced by the script but missing from `protocol.json`, and this document had stale family sizes. All three now agree.
12. **Design.** The noncentrality of the PC-B2 comparison was halved by mistake (both powers round to 1 either way). Scalar inputs to the sampler raised an error.

From a second independent review (numerical-correctness and statistical-methodology roles, read-only, on commit `50be70d`; every finding was reproduced before it was corrected, and the regression tests are in `tests/v020/nuclear_b2/test_b2_guards.py`):

13. **NaN from the incomplete beta (BLOCKING).** At n = 2⁵³ − 1, p = ½ Boost returned NaN next to the centre; a NaN compares false, so the bisection moved the draw by 3,842 (and by 8,649 with a wider bracket). The part limit is now 2⁵⁰ (BIG: 1,025 parts; XREF's large case: part limit 2⁵⁰), and any non-finite value of F or its complement raises an error instead of steering a draw. `numerics`, `cdf_finiteness`: no non-finite value at 14,001 points reachable by the bisection up to 2⁵⁰, 440 beyond it near 2⁵³ (development run, `research/v020-nuclear-b2-numerics/numerics_development.json`). Exact symmetry anchors up to 2⁵⁰ − 1 and the invariance of draws to the initial bracket are tested.
14. **Moments (BLOCKING).** The covariance was computed as E[fg] − E[f]E[g], which cancels for rare outcomes (−2,048 instead of 512 at N = 2⁶², p = 2⁻⁵³), and a NumPy integer count made N(N − 1) overflow (a negative variance at N = 10¹²). Moments now come from centred values, and counts are exact integers.
15. **D8 with a non-finite reference (BLOCKING).** A NaN reference never raised a maximum, so D8 read it as agreement (half-life 10⁻³⁰⁰, interval 10³⁰⁸). A non-finite decimal, path-tree or Layer A value now fails D8, and the engine's plan diagnostics count such comparisons as unavailable.
16. **Adaptive precision (MAJOR).** The precision growth could jump past the 1,000-digit maximum and refuse evaluations that 700 digits complete. The maximum is now tried before refusing.
17. **The midpoint variance bound (MAJOR).** The stated bound 2⁻⁵⁴(z_low² + z_high²) left out the tails: at n = 1, p = 2⁻⁵⁵ every draw is 0, a relative variance error of 1. It is replaced by a bound that includes the tails (`binomial.py`), evaluated exactly in `numerics` (`midpoint_convention`): it holds in all 12 cases, equals 1 at mean 2⁻⁵⁵ and is 7.9 × 10⁻¹⁵ for large means (development run, `numerics_development.json`).
18. **Floating-point distribution function (MINOR).** The documentation claimed exact inversion; at n = 11, p = ½, u = ½ the evaluated upper tail is ½ − 2⁻⁵³ and the draw is 5 where exact arithmetic gives 6. It now states that the sampler inverts the evaluated F; `numerics`, `table_boundaries`: 948 of 27,408 draws at uniforms adjacent to a step differ from exact inversion, none of 26,972 draws eight cells away (development run, `numerics_development.json`). The claimed truncation of the recurrence table (below 10⁻³⁰⁰) was wrong: the mass beyond it reaches 10⁻¹⁶²·⁴ (`table_truncation`, same file); the claim is now below 10⁻¹⁵⁰.
19. **Calibration (BLOCKING).** The calibration pooled dependent test instances into one binomial, so its "exact" intervals were not valid, several test types were not calibrated, and R1 claimed a false-FAIL bound that nothing established. `calibration.py` now simulates complete null batteries of both families (section 6, R0) with valid intervals for the family-wise rate and per test, R1 no longer claims a bound, and criterion C1 is declared before the freeze.
20. **Coverage and "equivalence" (MAJOR).** Photons were tested only as totals, so a wrong allocation over time passed every check. BT-PHOTON-T (7 tests per family) tests photons per interval. The tested second moments and the untested features of the joint law are stated (section 4), "equivalence" is replaced by non-rejection with a stated sensitivity, and `design.json` gives the power of the variance tests.
21. **Homogeneity cells (MINOR).** The homogeneity tests pooled cells by observed margins, so applicability depended on the data and an undefined comparison stopped the script. Their cells are now fixed in `design.json` from expected counts, a cell empty in both rows keeps its degree of freedom (conservative), and a test that cannot be computed is a failed item.

Also: `calibration.py` is now one of the frozen files; the smoke run no longer needs `protocol.sha256` and records that it ran unfrozen.

Third review (architecture and performance roles, on the revision with the development evidence), each finding reproduced before it was corrected; regression tests in `tests/v020/nuclear_b2/test_b2_integrity.py`:

22. **A plan could change under a run (BLOCKING).** A sink that emptied the plan's list of source nuclides made 4,999 nuclei vanish in the review's probe, and a plan array could be made writable again. Plans are now immutable: attributes cannot be rebound, containers are tuples and read-only mappings, arrays are views of immutable bytes whose write flag cannot be set again. A plan keeps its own snapshot of the nuclides and the scheme's SHA-256, and a run refuses a scheme that has changed since. Results are immutable in the same way. Conservation, signs and photons not exceeding decays are checked after every interval. Item 9's forced-write test now requires that the write flag cannot be set.
23. **Stream identity (BLOCKING).** A checkpoint named its generator and seed only, so a PCG64Source with another stream id resumed it. A checkpoint now records the exact generator class, every parameter that fixes its sequence (the seed, and for PCG64Source the stream id) and the word position. Checkpoint and resume accept NumpyPCG64Source, rnglab14's PCG64Source and Xoshiro256PlusPlusSource only, and refuse any other source before any draw: fault adapters, physical or OS entropy, conditioned sources, subclasses (checkpoint schema v2).
24. **Comparison arithmetic (BLOCKING as reported; not reachable with the declared cases).** `allocation_check` summed counts in int64: 2⁶² + 2⁶² wrapped negative, the cell was skipped and p was 1. A degree of freedom that underflowed to 0 gave NaN, which `min(1, NaN)` read as p = 1. Fractional counts were truncated silently. Counts are now exact integers and must be integers, degrees of freedom must be positive and finite (overflow included), and both tails must be finite before a two-sided p-value is formed. Otherwise the test raises, and in the study that is a failed item.
25. **Exclusive use of the source; manifests (MAJOR).** A sink that read the run's source shifted the physics stream unnoticed, and a manifest could be built with an unrelated stream. After every interval and callback the engine now checks that the source has emitted exactly the layout's words through the run's stream (RuntimeError otherwise). Manifests are built from the provenance the run records when it ends, and `aggregate_manifest` takes no stream (manifest schema v2).
26. **Once-only record runs (MAJOR).** Nothing durable was written before the end, an interruption left no trace, and two launches could race. The record runs of `calibration.py` and `run_verification.py` now reserve an attempt by exclusive creation before drawing and record how it ended (SIGTERM, which SLURM sends at a time limit, included). Their files are kept per attempt, and results are published without overwriting. Recovery is a rerun on the same commit, which repeats the same draws; a rerun on another commit is a recorded deviation (`protocol.json`, `record_procedure`). The calibration compares the digest of every finished chunk with earlier attempts.
27. **Failure leaves the stream mid-interval (MINOR).** The reference engine read its words slot by slot. It now reads each interval's words first, as the vectorised engine does, so both stop at an interval boundary.
28. **Namespace guard (MINOR; frozen test).** The frozen Layer A guard reads only top-level namespaces, so the calibration record namespace was not protected. A B2 guard reserves every namespace any nuclear protocol names, at any depth, scans every test file for them, and makes `namespace_seed` refuse them while a B2 test runs. The Layer A test is unchanged.
29. **Checkpoint cost (MAJOR, performance).** Every checkpoint serialised the whole grid, so checkpointing every interval cost O(I²). A checkpoint now holds the grid's SHA-256 and its number of intervals, O(nuclides) each.
30. **Memory (MAJOR, performance).** "Bounded memory" held for retained outputs only. The engine now states the plan's O(I) grid, its O(L·C) laws (L distinct interval lengths, C categories) and the tree check's dense K × K matrix. It counts categories exactly before enumerating them and refuses more than 1,024 for one start nuclide; the study's schemes have at most 9.
31. **Benchmark and calibration settings (MINOR, performance).** The benchmark took plan time from the first repeat. Every timing is now a median with its samples kept, repeats that disagree on a deterministic field raise an error, and parallel timings are labelled as including interpreter start-up. The calibration kept an inherited `OPENBLAS_NUM_THREADS=4`. It now sets one BLAS/OpenMP thread per worker whatever was inherited, records the inherited values, and by default uses the CPUs its affinity mask allows.

None of items 22–31 changes a draw or a statistic. The smoke run gives the same 13 run digests and the same p-values as the previous revision (whose smoke run was repeated from an exported copy), and `design.json` regenerates byte-identically.

## 8. Reviews and calibration before the freeze

- **Review status.** Items 5–12 answer the first review and items 13–21 the second, which covered the numerical and statistical roles only. Items 22–31 answer the third, in which the architecture and performance roles reviewed B2 for the first time; the numerical and statistical re-reviews of that round were refused by the reviewer's quota. The study must not be run before the corrected revision has passed an independent review in all four roles: the numerical and statistical re-reviews, and a check of items 22–31 by the architecture and performance roles. The outcome will be recorded here before the freeze.
- **Development calibration** (`calibration_development.json`: 4,000 batteries, development namespace; its `code` field names the commit it ran on; to be repeated for the record after the freeze with 50,000 batteries). A correct engine FAILed in 34 of 4,000 NC1 batteries (0.85%, 95% interval 0.59–1.19%) and in 39 of 4,000 NC2 batteries (0.97%, 0.69–1.33%); C1 holds for both. At the Holm level the per-type counts are small: the largest shares of the union bound were BT-NEUTRON in NC1 (5 batteries, 9.3 times its share, interval 3.0–21.6) and BT-PHOTON-T in NC2 (11 batteries, 3.3 times, 1.7–6.0); the same types in the other family had 1 and 0 batteries. The record run will report them with tighter intervals.

## 9. Reproduction

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1    # before Python starts (in a SLURM script too)
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/design.py
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/calibration.py --development --batteries 4000
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/run_verification.py --smoke
# The record, after the freeze commit, each step once, from a clean tree:
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/calibration.py --workers 8
#   commit results/calibration.json and results/attempts/ (the next step refuses a dirty tree)
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/run_verification.py
#   commit results/verification.json and results/attempts/
```

An attempt that did not complete is rerun with the same command on the same commit. If it left no end record (a node failure, SIGKILL), add `--previous-attempt-ended`. A rerun on another commit needs `--deviation "<reason>"` and is reported as a deviation.
