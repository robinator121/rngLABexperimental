# RNGLab v0.20 nuclear module, milestone B2: numerical limits

**Status: 2026-11-19.** Measured by `research/v020-nuclear-b2-numerics/numerics.py`, run once at commit `d91cb45` on a clean tree (`results/numerics.json`, started 2026-11-19T17:44:42Z, 246 s in total; Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1). It is a measurement, not a hypothesis test: the references are exact (rational arithmetic, or decimal arithmetic at 50 to 150 digits) or asymptotic with a stated error order, and the only random draws are B2 runs in the namespace `rnglab-v020-nuclear-b2-numerics`. Every number below comes from `results/numerics.json` or from `results/tables.md`, which `report_tables.py` writes from it. The development run of the same script (`numerics_development.json`, an earlier commit) informed items 13, 17, 18, 51, 52 and 55 of the B2 protocol; this run documents the frozen code.

## 1. Summary

| Component | Measured limit | Section |
|---|---|---|
| Binomial distribution function, as the sampler consults it (n from 10 to 10⁹) | largest absolute error 3.2 × 10⁻¹⁶ (Boost's F alone 3.2 × 10⁻¹⁴, its upper tail S alone 3.7 × 10⁻¹¹) | 2.1 |
| The same at large n with small means (n from 10⁶ to 2⁵⁰, means 1 to 1,000) | 5.8 × 10⁻¹⁶ (recurrence table 1.7 × 10⁻¹⁵; Boost's S alone 2.3 × 10⁻⁸, which is why means up to 256 use the table) | 2.1 |
| The same at n from 10¹² + 1 to 2⁵⁰ − 1, p = ½, against exact values | 9.3 × 10⁻¹⁴ | 2.1 |
| Finiteness of Boost's F and S up to the part limit 2⁵⁰ | no non-finite value at 14,001 scanned points or at the 7,995 evaluations of 1,008 traced draws; 440 non-finite values beyond the limit, near 2⁵³ | 2.2 |
| Inversion against exact arithmetic | no mismatch in 91,130 grid and near-step draws or 80 extreme-uniform draws; next to a step, 948 of 27,408 draws differ, none of 26,972 eight cells away | 2.3 |
| Midpoint convention | relative variance error at most 4.2 × 10⁻¹⁶ for means of 10 and more; outcomes much rarer than 2⁻⁵³ are lost (error 1 at mean 2⁻⁵⁵) | 2.4 |
| P(X ≥ 1) for p down to 10⁻³⁰⁰ | relative error at most 8.3 × 10⁻¹⁷ | 2.4 |
| Interval laws (decimal Bateman terms) | within 1.1 × 10⁻¹⁶ of the matrix exponential for close rates and 8.9 × 10⁻¹⁶ for extreme ones; 50 to 100 digits | 2.5 |
| Populations from 10⁶ to 2⁶² nuclei | conserved; 1 to 4,096 parts; above 2⁶³ − 1 refused | 2.6 |
| numpy's binomial (the calibration's null only) | wrong for small expected counts at n ≥ 10¹⁵; no error visible at the BIG null's calls | 3 |

## 2. The B2 sampler and its references

### 2.1 Binomial distribution function

Against exact sums in 60-digit decimal at nine points from −6 to +6 standard deviations (n from 10 to 10⁹, p from 10⁻⁶ to ½), the tail the sampler actually consults is within 3.2 × 10⁻¹⁶ of the exact value. Boost's F alone is within 3.2 × 10⁻¹⁴, and its upper tail S within 3.7 × 10⁻¹¹; the sampler consults the smaller of the two tails. Scipy's `bdtr` (cephes), which B2 does not use, is off by up to 0.34 at n = 10⁹.

At large n with small means, against a 50-digit probability recurrence (n from 10⁶ to 2⁵⁰, means 1 to 1,000), the sampler's tails are within 5.8 × 10⁻¹⁶ and the recurrence table within 1.7 × 10⁻¹⁵, while Boost's S alone is off by up to 2.3 × 10⁻⁸; B2 uses the table for means up to 256. The probability mass beyond the table is at most 10⁻¹⁶²·⁴.

From n = 10¹² + 1 to 2⁵⁰ − 1 at p = ½, against exact values (the centre value ½ for odd n plus exact probability sums), the tail the sampler consults is within 9.3 × 10⁻¹⁴; At the centre, Boost's F is exact for every n tested (10⁹ + 1 to 2⁵⁰ − 1), and its S from 10¹² + 1 on (off by 2.5 × 10⁻¹⁵ at 10⁹ + 1). One evaluation takes 0.04 ms at n = 10⁹ + 1 and about 0.5 ms from 10¹⁵ + 1 on.

### 2.2 Finiteness and the part limit

Boost's F and S were evaluated at 14,001 points up to the part limit 2⁵⁰ (the mean ± 9 standard deviations and the centre, nine probabilities): none is non-finite. Beyond the limit, near 2⁵³, 440 values are non-finite, which is why counts above 2⁵⁰ are split into parts. The 1,008 traced draws (n from 10⁶ to 2⁵⁰, extreme and central uniforms, narrow brackets that force the fallbacks) made 6,325 evaluations of F and 1,670 of S, reaching standardised points from −1.7 × 10⁷ to 1.6 × 10⁷: none non-finite, and no draw refused.

### 2.3 Inversion

Against exact rational inversion at the cell midpoint, the sampler has no mismatch in 90,045 draws on a grid of uniforms (n up to 55, five probabilities) or in 1,085 draws 10⁻¹² either side of every step. It has none in 80 draws at extreme uniforms (0, 2⁻⁵³, …, 1 − 2⁻⁵³; n from 50 to 10⁹, both regimes) against exact decimal quantiles. Four exact symmetry anchors (odd n from 10⁹ + 1 to 2⁵⁰ − 1, p = ½) hold, and the draws at 520 uniforms do not change with the initial bracket (half-widths 1, 4 and 64).

The sampler inverts the evaluated distribution function, not the exact one. At the uniforms adjacent to every step (n = 1 to 100 at four probabilities, and three cases in Boost's regime; 403 cases), 948 of 27,408 draws differ from exact inversion, because the evaluated F puts the step on the other side of the midpoint. None of 26,972 draws eight cells away differs.

### 2.4 The midpoint convention and very small probabilities

A uniform k·2⁻⁵³ is inverted at its cell's midpoint. Evaluated in high precision in 12 cases, the bound on the relative variance error, tails included, holds everywhere. For means of 10 and more the bound is about 7.9 × 10⁻¹⁵ and the exact error at most 4.2 × 10⁻¹⁶. For tiny means the draw loses outcomes much rarer than 2⁻⁵³: the relative variance error is 1 at mean 2⁻⁵⁵ (every draw is 0), 0.11 at mean 10⁻¹⁶, 8.0 × 10⁻⁴ at mean 10⁻¹⁴ and 1.6 × 10⁻⁸ at mean 10⁻⁸.

P(X ≥ 1) = 1 − (1 − p)ⁿ, as Boost evaluates it, is within a relative 8.3 × 10⁻¹⁷ of 80-digit decimal in all 28 cases (n = 1, 10⁶, 10¹² and 2⁵³; p from 10⁻³⁰⁰ to 10⁻⁵).

### 2.5 Interval laws

The decimal Bateman terms adapt their precision. For a four-member chain with relative rate spacing from 0.1 down to 10⁻¹⁴ they use 50 to 84 digits and stay within 1.1 × 10⁻¹⁶ of the matrix exponential, with every mean-value bracket holding; the float64 closed form of Layer A refuses every spacing from 10⁻³ down. For 25 cases with half-lives from 10⁻¹² to 10¹² and intervals from 10⁻¹⁵ to 10¹⁵, they stay within 8.9 × 10⁻¹⁶ of the matrix exponential, sum to 1 within 8.8 × 10⁻⁵⁰ and use at most 100 digits. For tiny intervals the deepest path of S2 (four decays) has probability 4.1 × 10⁻¹⁶ at Δ = 10⁻³ and 4.1 × 10⁻⁶⁴ at Δ = 10⁻¹⁵; the float64 matrix exponential gets it within a relative 5.7 × 10⁻¹³. On the stiff chain S5 the float64 matrix exponential is within 1.6 × 10⁻¹³ of 150-digit decimal laws, and the default precision within 3.4 × 10⁻⁴⁷.

### 2.6 Populations

B2 ran S2 over six intervals once per size, from 10⁶ to 2⁶² nuclei: one part up to 2⁵⁰, then 64 parts at 2⁵⁶, 1,024 at 2⁶⁰ and 4,096 at 2⁶², taking 0.04 s to 88 s. Every run conserves the population. Against the exact means, the largest |z| of the 54 values per size is between 1.9 and 2.7 (one replicate per size, descriptive). A total population above 2⁶³ − 1 is refused.

## 3. numpy's binomial (used by the calibration's null, not by B2)

With 4,000,000 draws per row (standard error of a variance ratio 0.0007), numpy's `binomial(n, p)` is wrong for small expected counts at large n:

| n | Expected count | Mean z | Variance ratio |
|---|---:|---:|---:|
| 10¹² | 115; 1,000; 10⁴ | −0.95; −0.37; 0.65 | 0.9986; 1.0005; 1.0006 |
| 10¹⁴ | 115; 1,000; 10⁴ | 2.86; −0.39; −1.39 | 1.0006; 1.0007; 1.0001 |
| 10¹⁵ | 115; 1,000; 10⁴ | 3.73; −7.93; 1.49 | 1.0084; 1.0049; 1.0007 |
| 2⁵³ − 1 | 115; 1,000; 10⁴ | 44.46; 26.70; 7.53 | 1.0483; 1.0167; 1.0034 |
| 2⁶⁰ + 12,345 | 115; 1,000; 10⁴; 10⁶; 10⁸; 10¹⁰ | −10.84; −93.93; 7.92; 0.59; −0.61; 0.15 | 1.2212; 0.9978; 1.0016; 1.0019; 0.9998; 0.9997 |

Every binomial call that the BIG null makes with n above 10¹³ has an expected count of at least 2.27 × 10¹² (`design.json`). At that smallest call's own n and p (n = 2.27 × 10¹³, p = 0.9), the mean z is −0.76 and the variance ratio 0.9994.

## 4. What this does not show

- The scans are not exhaustive. They cover the grids and traced draws listed, and they support the part limit 2⁵⁰ empirically; they do not prove that Boost's incomplete beta is finite everywhere below it.
- The decimal references are high-precision evaluations, not certified interval arithmetic.
- The numpy findings concern numpy 2.5.3's binomial, which only the calibration's null uses; B2 does not use it.
- Accuracy is measured for the probabilities, populations and interval lengths listed; other regimes rely on the same code without this evidence.
