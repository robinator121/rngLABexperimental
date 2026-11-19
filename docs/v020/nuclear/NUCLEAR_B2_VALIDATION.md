# RNGLab v0.20 nuclear module, milestone B2: aggregated decay validation report

**Status: B2, 2026-11-19.** Results of the frozen study `research/v020-nuclear-b2-verify` (protocol frozen in commit `25ff551`), run once in the order the protocol prescribes, each from a committed tree: the record calibration at `25ff551` (`results/calibration.json`, started 2026-11-19T17:13:55Z), then the verification at `e0d3b09`, the commit that adds the calibration (`results/verification.json`, started 2026-11-19T17:43:04Z). Every number below is taken from those two files or from `results/tables.md`, which `report_tables.py` writes from them; design numbers come from `design.json`, numerical limits from `research/v020-nuclear-b2-numerics/results/numerics.json`, timings from `research/v020-nuclear-b2-bench/results/benchmark.json`, and development figures from the protocol (`PROTOCOL.md`, section 8). Prepared for the author. Layers B and C have no validation yet.

## 1. What was verified, and what was not

**Question.** Does the B2 aggregated (population-level) decay engine reproduce the exactly known laws of independent decay, and the frozen Layer A event engine, for five synthetic schemes, including populations far beyond what an event engine can simulate? And is the battery used for that check sensitive enough to notice a small physics error and a plausible aggregation error?

**Answer under the frozen reading rules: R1 PASS, R2 MET; R0: criterion C1 holds for both families.**

| Law or property | Test (PROTOCOL.md, section 4) | Outcome |
|---|---|---|
| State counts at every grid time (multinomial law of independent decay) | BT-STATE | not rejected (NC1, NC2) |
| Decays of each nuclide over the grid, and per interval | BT-DECAY, BT-BINS | not rejected |
| Channel splits; photon yields | BT-CHANNEL, BT-PHOTON | not rejected |
| Photons spread over time, given the decays per interval | BT-PHOTON-T | not rejected |
| β-delayed neutrons per interval (S4) | BT-NEUTRON | not rejected |
| Path-category counts of each interval against the decimal interval laws (S2) | BT-PATH | not rejected |
| Ensemble means, variances, covariances and two-time covariances (ENS-B2) | ET-MEAN, ET-VAR, ET-VARSUM, ET-VARSUM2 | not rejected (see Section 2.6) |
| B2 against the frozen Layer A engine: pooled counts, decays, photons, variance ratios | CT-HOMO, CT-DECAY, CT-PHOTON, CT-VAR, CT-VARSUM | not rejected (see Section 2.6) |
| Populations of 10¹² and 2⁶⁰ + 12,345 nuclei | BIG-STATE, BIG-DECAY, D9 | not rejected; 1 and 1,025 parts; population conserved |
| Interval laws: decimal Bateman terms against two matrix exponentials | D8 | agree within 1.57 × 10⁻¹³ (tolerance 10⁻¹²) |
| Fixed consumption; bookkeeping and Layer A invariants | D1, D2 | hold in every run |
| Two B2 engines bit for bit (also forced into 101 parts); replay; resume from a checkpoint; streaming | D3–D6 | hold |
| Distinct ensemble seeds and streams; exactly the planned tests ran | D7, D10 | hold |

**Not verified here:** whether any test detects RNG defects, or how fast (Phases E and F); Layers B and C; any agreement with nuclear data (the schemes are synthetic by design); features of the joint law across times other than the tested variances and covariances (protocol section 4, "Coverage of second moments"); behaviour at sample sizes, depths or grids other than those declared. Not rejecting is not equivalence: section 4 of the protocol states how large a deviation the tests detect.

| Status category | B2 aggregated decay engine |
|---|---|
| Implemented | yes: `aggregate.py`, with `binomial.py`, `paths.py`, `moments.py`, `comparison.py` and `precursors.py` |
| Tested | 202 unit tests in their own folder (`tests/v020/nuclear_b2`), beside Layer A's 115 |
| Independently verified | against analytic laws and the frozen Layer A engine by the frozen study (R1 PASS, R2 MET, C1 holds); code and protocol reviewed in eight rounds before the freeze |
| Experimentally validated | no, and not applicable: the schemes have no experimental counterpart |

## 2. Results

### 2.1 Calibration (R0): how often a correct engine would FAIL

50,000 independent null batteries in the record namespace `rnglab-v020-nuclear-b2-calibration-record` (attempt 1 at commit `25ff551`, the freeze commit; no deviation; no earlier attempt to compare; no battery with an error). Every run of NC1 and NC2 was drawn from the exact laws with numpy, without any RNGLab sampler, and every planned test was evaluated with the functions the verification uses.

| Family | Tests | Batteries that FAIL under the null | Rate | 95% interval (Clopper–Pearson) | C1 (upper bound ≤ 0.02) |
|---|---:|---:|---:|---|---|
| NC1 | 222 | 440 | 0.0088 | 0.0080–0.0097 | holds |
| NC2 | 85 | 399 | 0.0080 | 0.0072–0.0088 | holds |
| R1 event (NC1 or NC2 fails) | | 835 | 0.0167 | 0.0156–0.0179 | reported, not a criterion |

Under R0, each family's false-FAIL rate is below 2% (95% confidence). R1's own false-FAIL rate, the rate at which a correct engine would FAIL R1, is estimated at 1.67%: the two families add up, because a battery that fails both is rare (4 of 50,000). Both family rates are slightly below the nominal 1%.

Per test type, the calibration counts the batteries in which a type rejects at the family's Holm level (NC1: 0.01/222 = 4.50 × 10⁻⁵; NC2: 0.01/85 = 1.18 × 10⁻⁴), as a ratio to the union bound (the type's number of tests times the level). Most ratios lie between 0.7 and 1.4, and every interval but five contains 1 (`results/tables.md`). Below 1, as expected for conservative or positively dependent tests: ET-VARSUM2 in NC1 (0.50, 95% interval 0.29–0.81), BIG-DECAY in NC1 (0.55, 0.31–0.92) and the exact binomial BT-DECAY (NC1 0.67, 0.43–0.99; NC2 0.70, 0.54–0.89). Above 1: BT-NEUTRON in NC1, 17 batteries, 2.52 times the union bound (1.47–4.03). That is the largest of the 27 type ratios, and its interval is not adjusted for the selection. The development calibration had the same type as its largest ratio, and the development tail probe measured single BT-NEUTRON tests at 1.25, 1.07 and 1.07 times the NC1 Holm level; those values lie below this interval. In NC2, at its larger Holm level, BT-NEUTRON is at 0.96 (0.56–1.54). Section 5 lists this as an open question. C1 calibrates each family as a whole, this effect included, and no verdict depends on a type's ratio.

### 2.2 Deterministic checks

All hold (D1–D10). D8 compared 113 interval laws: the largest difference between the decimal law and the matrix exponential on the path tree is 1.57 × 10⁻¹³, and between the end states and the Layer A matrix exponential also 1.57 × 10⁻¹³ (tolerance 10⁻¹²). Probabilities sum to 1 within 5.91 × 10⁻⁴⁷ (tolerance 10⁻³⁰), every mean-value bracket holds, at most 50 digits were needed, and no reference was non-finite. XREF: S2 in one part and S4 forced into 101 parts give equal digests, arrays and stream digests in both B2 engines.

### 2.3 Test families (Holm, family-wise α = 0.01)

| Family | Generator | Tests | Smallest p | Smallest Holm-adjusted p | Rejections |
|---|---|---:|---:|---:|---|
| NC1 | numpy PCG64 | 222 | 0.0038 (ENS ET-VARSUM2, Xc-400m at t = 2 plus t = 8) | 0.84 | none |
| NC2 | xoshiro256++ (pure Python) | 85 | 0.030 (S4 BT-STATE G2 at t = 5) | 1 | none |

Every test with its statistic, p-value and Holm-adjusted p-value is listed in `results/tables.md`.

### 2.4 Positive controls (level p < 0.01/222 = 4.50 × 10⁻⁵)

| Control | Perturbation | Target test | p | Rejected |
|---|---|---|---:|---|
| PC-B1 | half-life 1.02 instead of 1.0 (S1, 400,000 nuclei) | BT-STATE Sa-400 at t = 2 | 2.19 × 10⁻²⁰ | yes |
| PC-B2, single run | one-step aggregation: at most one decay per nucleus and interval (S2, 10⁵ nuclei) | BT-STATE Xa-400 at t = 12 | 0 (below double precision) | yes |
| PC-B2, comparison | the same aggregation, 400 replicates of 2,000 nuclei, against ENS-A | CT-HOMO at t = 12 | 0 (below double precision) | yes |

Outside any family, PC-B1 is rejected by every state test (p from 1.1 × 10⁻¹¹ to 2.2 × 10⁻²⁰), by BT-DECAY (2.1 × 10⁻¹²) and by BT-BINS (1.5 × 10⁻²⁸), but not by the photon tests, since the photon yield was not changed (BT-PHOTON p = 0.69, BT-PHOTON-T p = 0.069). The aggregation error of PC-B2 is rejected at every grid time from t = 0.5 on.

### 2.5 Huge populations (BIG)

10¹² nuclei ran in one part (0.10 s), and 2⁶⁰ + 12,345 = 1,152,921,504,606,859,321 nuclei in 1,025 parts (21.3 s). Both conserve the population (D9). Their 24 tests, part of NC1, have p from 0.046 to 0.94.

### 2.6 Observation, after the run, not part of any rule

The three smallest p-values of NC1 all involve the ensemble variance of Xc-400m at t = 2: ET-VARSUM2 for Xc-400m at t = 2 plus t = 8 (p = 0.0038), CT-VARSUM for Xc-400m + Xe-396 at t = 2 (p = 0.0046) and CT-VAR for Xc-400m at t = 2 (p = 0.018). They use the same 400 replicates of each engine, so they are not independent evidence. B2's sample variances lie above the exact values (Xc-400m at t = 2: 160.7 against 143.1, ratio 1.12, ET-VAR p = 0.091; the two-time sum: ratio 1.22), and Layer A's below (126.8 at t = 2); the comparison tests see the sum of the two deviations. None is close to rejection after Holm (smallest adjusted p 0.84). By design these ensembles are weak for variances: a variance 10% too large is detected with probability below 0.005 by any one test, and the ratios detected with probability 0.9 are about 1.45 (ET tests) and 1.71 (CT tests) (`design.json`, `sensitivity`). So the ensembles neither show nor exclude a variance deviation of this size. A replication with more replicates in a fresh namespace would settle it if it matters.

### 2.7 Descriptive (R3) and single-machine timings (wall clock; not HPC measurements)

On the author's workstation (one Python process, Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1), each NC1 and NC2 run of 10⁵ nuclei took 0.006–0.041 s, the ensemble runs 17.2 s, the BIG runs 0.10 s and 21.3 s, and the whole verification 39.4 s. The decimal interval laws needed at most 50 digits (D8). Section 3 and `NUCLEAR_B2_BENCHMARK.md` give the numerical limits and the benchmarks.

## 3. Numerical limits

`NUCLEAR_B2_NUMERICAL_LIMITS.md` reports the numerical investigation run on the frozen code (`research/v020-nuclear-b2-numerics/results/numerics.json`, commit `d91cb45`). In brief:

| Item | Value | Effect |
|---|---|---|
| Binomial distribution function, as the sampler consults it | within 3.2 × 10⁻¹⁶ of exact sums for n up to 10⁹, 5.8 × 10⁻¹⁶ at small means up to n = 2⁵⁰, 9.3 × 10⁻¹⁴ at n from 10¹² + 1 to 2⁵⁰ − 1 | bounds how far a draw's law can be from the binomial law, beside the 2⁻⁵⁴ of the midpoint rule (protocol, section 2) |
| Part limit 2⁵⁰ | no non-finite value of Boost's F or S at 14,001 scanned points or along 1,008 traced draws below it; 440 non-finite values near 2⁵³ | populations above 2⁵⁰ are split into parts (1,025 for the BIG run) |
| Inversion | inverts the evaluated distribution function: next to a step, 948 of 27,408 draws differ from exact inversion; none elsewhere | a draw can differ from exact inversion only when its uniform lies next to a step of the distribution function |
| Midpoint convention | relative variance error at most 4.2 × 10⁻¹⁶ for means of 10 and more; outcomes much rarer than 2⁻⁵³ are lost | matters only for expected counts far below 1 |
| Interval laws | decimal Bateman terms within 1.1 × 10⁻¹⁶ (close rates) and 8.9 × 10⁻¹⁶ (extreme rates) of the matrix exponential; at most 100 digits | D8 observed 1.57 × 10⁻¹³ against a tolerance of 10⁻¹² |
| Populations | conserved from 10⁶ to 2⁶² nuclei; above 2⁶³ − 1 refused | the representable range is 64-bit signed integers |
| numpy's binomial (calibration only) | wrong for small expected counts at n ≥ 10¹⁵; no error visible at the BIG null's calls (smallest expected count 2.27 × 10¹²) | the calibration's BIG null is not exact, see Section 5 |

## 4. Review before the freeze

Before anything was drawn from the study or record calibration namespaces, the code and the draft protocol went through eight review rounds in four roles: numerical correctness, statistical methodology, scientific software architecture, and performance and memory. The reviews were read-only: reviewers inspected the code and ran their own checks but changed nothing. Every reported defect was reproduced before it was fixed, and the code fixes are covered by regression tests (among them the 66 of `test_b2_integrity.py` and the 24 of `test_b2_guards.py`). The 69 corrections (items 5–73) are listed in `research/v020-nuclear-b2-verify/PROTOCOL.md`, section 7, and the rounds in section 8. Twelve were blocking:

- **Numerics (second round).** Next to the centre of the distribution at n close to 2⁵³, Boost's incomplete beta returned NaN, which steered the bisection to a wrong draw (item 13; parts are now limited to 2⁵⁰ trials, and any non-finite value raises an error). Moments came from a cancelling formula and an overflowing integer product (item 14). D8 read a non-finite reference as agreement (item 15).
- **Statistics (second round).** The calibration pooled dependent tests into one binomial, so its intervals were not valid and R1 claimed a false-FAIL bound that nothing established; it now simulates complete null batteries (item 19).
- **Architecture and record machinery (third and fourth rounds).** A caller could change a plan or an array under a run (items 22 and 32; plans are now immutable and arrays are fresh read-only views). A checkpoint did not pin its generator's stream (item 23), and a source could be replaced during a run (item 33). The allocation statistic had int64 and cancellation errors (items 24 and 34; item 24 was not reachable with the declared cases; the statistic now uses exact integers and exact rationals). R0 inferred from an upper bound alone (item 36). The clean-tree rule blocked the documented rerun of an interrupted attempt (item 37).

The record machinery itself (attempt locks, atomic publication, code identity, the committed-tree rule, the calibration checks the verification makes before it starts) was reviewed in rounds four to eight. The last verdict of every role was "safe to freeze after the listed fixes"; those fixes are items 70–73. The author approved the freeze on 2026-11-19.

## 5. Limitations and open questions

- The schemes are synthetic. Agreement with analytic laws says nothing about agreement with nuclear data.
- The check is algorithmically independent, but not organisationally: the engines, the references, the tests and the reviews come from one project; they are not an external audit. Two analytic routes for the interval laws (decimal Bateman terms and the matrix exponential), two engines (B2 and Layer A), two generators and the calibration reduce the risk of a shared conceptual error without removing it.
- The calibration's null is drawn with numpy's samplers, which are not exact at large n (items 52 and 55). Every binomial call the BIG null makes above n = 10¹³ has an expected count of at least 2.3 × 10¹², where no error is visible; this is an empirical statement, not a proof.
- The decimal interval laws are evaluated in high-precision decimal arithmetic with adaptive precision, not in certified interval arithmetic.
- **Open: BT-NEUTRON in the calibration (Section 2.1).** In NC1, BT-NEUTRON rejected at the Holm level more often than the development tail probe measured for its single tests. C1 includes the effect. Two explanations remain: chance enlarged by selection, or a difference in the far tail between the calibration's null draws of S4 and the probe's direct multinomial draws. A probe of BT-NEUTRON through the calibration's null model at the NC1 Holm level would decide between them.
- **Open: the variance observation (Section 2.6)** is unresolved in the strict sense; a replication with more replicates in a fresh namespace would close it.
- Only healthy generators were used. Whether this battery, or any application observable, detects RNG defects is the subject of Phases E and F.
- Timings come from one workstation and are not HPC measurements.

## 6. Reproduction

```bash
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1    # before Python starts (in a SLURM script, before srun)
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/run_verification.py --smoke   # smoke sizes, not in Git
PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-verify/report_tables.py              # rewrites tables.md from the results
```

The record runs refuse to run again: their results exist and are never overwritten. A replication belongs in a new namespace under a new protocol version.
