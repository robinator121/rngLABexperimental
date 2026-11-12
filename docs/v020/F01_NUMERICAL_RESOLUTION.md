# F-01: the Stage D oracle probability returns NaN. Diagnosis and proposed resolution

**Status: for the author's decision, 2026-11-12.** Nothing in the frozen protocol or in the pilot's code was changed. Every number comes from `research/v019-audit/results/f01_oracle.json`, written by `research/v019-audit/f01_oracle.py`. Tests: `tests/v019/test_f01_oracle.py` (12).

## 1. Summary

| Question | Answer |
|---|---|
| Which operation produces NaN? | `scipy.special.nctdtr(15, δ, −c)`, the lower-tail term that `scipy.stats.nct.cdf(-c, df=15, nc=δ)` calls. The library reports "no result obtained" |
| Where? | Only for positive δ, at 1,554 of 12,001 points on [−60, 60] in steps of 0.01: between 8.53 and 43.61, at 44% of the points in that range. The upper term never fails |
| What kind of defect? | A numerical failure inside the library routine. The parameters are valid, and the formula is the intended one |
| How large is the failing term in truth? | At most Φ(−δ), which is 7.3 × 10⁻¹⁸ at 8.53. At 8.55 it is 3.6 × 10⁻²⁴, against q = 1.18 × 10⁻⁹ |
| Is the failure reproducible? | On one machine, for one double, yes. Across a rounding-level change of δ, no: 521 of 20,000 points change between NaN and a number when δ moves to the next double (section 3.1) |
| Is there a stable equivalent? | Yes, three. All agree with a 70-digit reference to better than 3 × 10⁻¹⁵ in absolute terms |
| Is an amendment needed? | **Yes, if the oracle is to be scored where the call fails.** The protocol fixes the call and says a non-finite value is unavailable |
| Is this the missing widening term (F-02)? | No. F-02 is a declared modelling choice. Fixing F-01 leaves it untouched |

## 2. The quantity

`protocol.json`, `oracle.probability`:

q(δ) = P(|T| ≤ c), with T noncentral t with 15 degrees of freedom and noncentrality δ, and c = `scipy.stats.t.ppf(0.975, df=15)` = 2.131449545559776.

It is evaluated as `nct.cdf(c, 15, δ) − nct.cdf(−c, 15, δ)`. `oracle.unavailable` says: "Nonfinite/out-of-range q is unavailable with reason, without clipping or substituting F0."

Two facts used below, both exact:
- **q is even in δ**, because −T under δ has the law of T under −δ.
- **For δ ≥ 0 the lower term is at most Φ(−δ)**: T ≤ −c < 0 requires Z + δ < 0.

## 3. Diagnosis

- **The call chain.** `nct.cdf(x, df, nc)` is `scipy.special.nctdtr(df, nc, x)`: the two return the same values and the same NaN at every grid point. With `scipy.special.errstate(all="raise")` the failing call raises `SpecialFunctionError: scipy.special/nctdtr: no result obtained`.
- **Only the lower term fails, and only for δ > 0.** On the grid the upper term is finite at all 12,001 points. For δ < 0 both terms are finite everywhere.
- **The failure is scattered.** At 9.5 the call fails, at 9.75 it returns 9.33 × 10⁻²⁸, at 10.0 it fails again.
- **It is not specific to this problem.** The same routine fails for 5 and 30 degrees of freedom and for arguments −0.5 and −5 (`failure_map.other_arguments`).

| Candidate cause | Verdict | Reason |
|---|---|---|
| Error in the pilot's implementation | no | the call is the one the protocol prescribes |
| Parameter outside the domain | no | df = 15, any real δ and x are valid; nearby points succeed |
| Wrong mathematical formulation | no | the formula is P(\|T\| ≤ c) for the stated law |
| Numerical failure in the library | **yes** | the routine reports that it obtained no result, for a tail probability below 10⁻²³ |

The library's source was not read. The diagnosis rests on its own error report and on the measured behaviour.

### 3.1 Whether a value exists depends on the last bits of δ

`rounding_sensitivity` in `results/f01_oracle.json`, on 20,000 equally spaced δ in [8.6, 43.5]:

| | Count |
|---|---|
| lower term not finite | 8,655 |
| changes between NaN and a number when δ moves to the next double | 521 |
| changes when δ is multiplied by 1 + 10⁻¹³ | 496 (2.5%) |

- **Consequence for the frozen rule.** "Non-finite q is unavailable" then makes the *set of runs with an oracle value* depend on rounding. Two correct programs that evaluate δ = β·S̄/se with S̄ summed in a different order agree on δ to 15 digits and can still disagree on availability.
- **Observed.** On a synthetic campaign of 207 runs, the pipeline and the independent recount agreed on 1,972 of 1,974 compared fields. The two differences were counts of unavailable oracle values: 3 against 2 (SHIFT, L = 20,480) and 7 against 8 (LAG1, L = 40,960). With option A of section 8 forced in memory on the same campaign, all 1,974 fields agreed and all nine oracle cells were scored (`docs/v020/evidence/synthetic_agreement.json`).
- **What this is not.** The synthetic campaign says nothing about how often the pilot's δ will fall in the failing range. It shows that when it does, exact agreement of two implementations cannot be demanded under the frozen rule.

**The frozen expression has a second, harmless weakness.** For negative δ it subtracts two numbers close to 1. At δ = −10 it returns the right value with a relative error of 2 × 10⁻⁴; at δ = −20 it returns 0 where the truth is 1.7 × 10⁻⁶⁰. The absolute error stays below 10⁻¹⁶.

## 4. Reference

`f01_oracle.py` computes q with Python's `decimal` module at 70 digits, without SciPy:

- F(t; δ) = ∫₀^∞ Φ(t·w − δ)·f(w) dw, with f the density of W = √(V/15). Both tails are integrated directly.
- The substitution w = eᵘ and the trapezoidal rule, with the step halved until two successive values agree to 10⁻⁴⁵. The largest last change was 4.8 × 10⁻⁴⁶.
- Φ from erfc: a positive-term series up to 2, a continued fraction above.
- c and δ enter as the exact values of their doubles; the reference is taken at |δ|.

Checks of the reference: π to 50 digits; erfc against `math.erfc` at 12 arguments; the two erfc branches agree at their switch point; q(0) equals the central-t value to 3 × 10⁻¹⁶.

## 5. Stable equivalents and their errors

| Name | Definition |
|---|---|
| bounded | the frozen expression at \|δ\|, with a non-finite lower term replaced by 0 |
| reflected | `nct.cdf(c, 15, δ) − nct.sf(c, 15, −δ)`, using F(−c; δ) = 1 − F(c; −δ) |
| mixture | E_W[Φ(c·W − \|δ\|) − Φ(−c·W − \|δ\|)] by Gauss–Legendre quadrature on fixed panels |

Against the reference at 33 values of δ from −20 to 100:

| | frozen | bounded | reflected | mixture |
|---|---|---|---|---|
| not finite | 8 | 0 | 0 | 0 |
| largest absolute error | 2.4 × 10⁻¹⁶ | 2.4 × 10⁻¹⁶ | 2.4 × 10⁻¹⁶ | 2.0 × 10⁻¹⁵ |
| largest relative error where q > 10⁻¹² | 1.5 × 10⁻¹⁴ | 7.4 × 10⁻¹⁵ | 1.5 × 10⁻¹⁴ | 6.1 × 10⁻¹⁵ |
| relative error at δ = 40 (q = 7.1 × 10⁻²⁵⁷) | 2.9 × 10⁻¹⁴ | 2.9 × 10⁻¹⁴ | 2.9 × 10⁻¹⁴ | 6.8 × 10⁻¹⁴ |
| at δ = 43.6 (q = 1.5 × 10⁻³⁰⁶) | NaN | absolute 8.5 × 10⁻¹⁵⁴ | absolute 8.5 × 10⁻¹⁵⁴ | relative 8.0 × 10⁻¹⁴ |

- **Absolute criterion.** Every equivalent is within 3 × 10⁻¹⁵ of the truth everywhere. A Brier score or a calibration bin cannot see that.
- **Relative criterion.** It holds down to about 10⁻²⁵⁰ for all three. Near the underflow limit the library's upper term is wrong in relative terms (8.5 × 10⁻¹⁵⁴ for 1.5 × 10⁻³⁰⁶); only the mixture stays accurate. No pilot score depends on that.
- **Beyond δ = 45** the truth is below the smallest positive double and every formulation returns 0.

## 6. Properties

On [0, 60] in steps of 0.01 (6,001 points):

| Property | bounded | reflected | mixture |
|---|---|---|---|
| inside [0, 1] | yes | yes | yes |
| value at 0 | 0.9500000000000002 | 0.9500000000000002 | 0.9500000000000007 |
| increases between neighbours (q must not increase in \|δ\|) | 52, largest 3.6 × 10⁻¹³⁵ | 52, largest 3.6 × 10⁻¹³⁵ | 0 |
| evenness, largest difference | 0 | 1.1 × 10⁻¹⁶ | 0 |
| at δ = 100, 10³, 10⁶ | 0 | 0 | 0 |

- **Unchanged where the frozen value exists.** At all 4,447 grid points where the frozen expression is finite, bounded and reflected return the same double, bit for bit.
- **Continuity.** The largest step between neighbours is 0.0037, at the steepest part of the curve; the filled values join the finite ones without a jump.
- **Convergence of the mixture.** With 25, 50, 100, 200 and 400 nodes per panel the values at δ = 0, 3, 8.55, 20 and 40 agree to 12 digits.

## 7. F-02 is a separate matter

The oracle formula assumes that the sixteen batch means share one mean. Under a ramp they do not, the interval widens, and the true probability is larger. That is F-02. The protocol declares it: "This ignores the widening of the interval by drift. It is a benchmark, not a bound."

| δ | frozen formula | with spread λ = 5 | with spread λ = 15 |
|---|---|---|---|
| 3 | 0.199 | 0.298 | 0.491 |
| 8.55 | 1.2 × 10⁻⁹ | 1.2 × 10⁻⁸ | 2.8 × 10⁻⁷ |

F-01 is about evaluating the frozen formula. F-02 is about whether that formula is the better model. Resolving F-01 does not touch F-02, and F-02 needs no amendment unless the author wants a different benchmark.

## 8. Proposed resolution (not applied)

**Proposed Amendment 01, option A (recommended): fill only what the frozen rule left unavailable.**

> q is computed by the frozen expression. If the lower term `nct.cdf(−c, 15, δ)` is not finite while δ > 0 and the upper term `nct.cdf(c, 15, δ)` is finite, the lower term is taken as 0, and the run is counted under `oracle_lower_term_replaced`. Any other non-finite or out-of-range value stays unavailable.

- It changes no value that the frozen rule can produce.
- The filled value differs from the exact probability by at most Φ(−δ) plus the error of the upper term: below 10⁻¹⁷ in total.
- It adds no new numerical method.

**Option B: evaluate q by the mixture quadrature everywhere.** More accurate in the far tail, exactly even and monotone, but it changes available values by up to 2.2 × 10⁻¹⁵ and puts a new algorithm into the protocol.

**Option C: no amendment.** Every evaluation run whose δ falls on a failing point has no oracle probability, and one such run makes the oracle's score for that family and length unavailable (`scores.brier`). Stop criteria 1 to 4 are unaffected. Two further costs follow from section 3.1:

- The count of unavailable oracle values is not a reproducible number: another SciPy build, another platform or another summation order can change it.
- The independent recount cannot be required to match that count exactly. Under option C the recount gate needs a declared exception for `oracle.unavailable` and `oracle.brier`, and that exception is itself a decision for the author.

How often option C would bite is **not known**: δ = β·S̄·√L/σ_eff depends on the bias slope and the noise scale, which the pilot estimates. With the scales used in the surrogate study (β = 0.06, σ_eff = 0.613), δ reaches about 9.3 at L = 40,960 for the fastest ramps, just inside the failing range. That figure is an illustration, not a prediction.

The pilot code on this branch implements the frozen rule. Option A is present behind a switch that is off unless an approved amendment record exists (`STAGE_D_IMPLEMENTATION_STATUS.md`).

## 9. Is a formal amendment required?

Yes for options A and B: they replace a rule that `protocol.json` fixes. The freeze procedure asks for "a versioned amendment, written before the affected analysis". The oracle is computed in the analysis, so the amendment can be written before or after generation, but before any oracle value is computed. Option C needs none.

## 10. Limitations

- One SciPy version (1.18.1) on one platform. Another version may fail elsewhere or not at all; the tests skip the failure-specific assertions in that case.
- The reference was computed at 33 values of δ. Between them the evidence is the agreement of the mixture with the other equivalents (largest difference 2.2 × 10⁻¹⁵ on 6,001 points).
- The library's upper term was validated, not proved. Option A relies on it; its absolute error at the reference points is below 2.4 × 10⁻¹⁶.
- A correction to `docs/v019/IMPLEMENTATION_AUDIT.md`: the quadrature used there as "stable" had large relative errors beyond δ ≈ 12 (its absolute error was below 10⁻¹⁹, so the statements there stand). It has been replaced by the mixture form in `research/v019-audit/recount.py`.

## 11. Reproduction

```bash
PYTHONPATH=src .venv/bin/python research/v019-audit/f01_oracle.py
```

```bash
PYTHONPATH=src .venv/bin/python -m pytest tests/v019/test_f01_oracle.py -q
```

The synthetic comparison of section 3.1 (needs the measurement layer; synthetic fixtures only):

```bash
PYTHONPATH=src RNGLAB_DPILOT_MEASURE=research/v019-d-pilot/measure.py .venv/bin/python tests/v019/synthetic_agreement.py
```

Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1, glibc 2.44.
