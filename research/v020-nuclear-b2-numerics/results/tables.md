# B2 numerical limits: tables from `results/numerics.json`

Code `d91cb45a28c638f526655677c5bb8dc8466d60af` (dirty: False); started 2026-11-19T17:44:42.411066+00:00; 246 s. Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1.

## 1. Binomial distribution function against exact sums (60-digit decimal), n up to 1e9

Largest absolute error over nine points from -6 to +6 standard deviations.

| n | p | Boost F | Boost S | cephes `bdtr` (not used) | as consulted by the sampler |
|---:|---:|---:|---:|---:|---:|
| 1e+01 | 0.5 | 0 | 0 | 2.2e-16 | 0 |
| 1e+01 | 0.1 | 4.8e-17 | 4.7e-17 | 1.1e-16 | 4.3e-17 |
| 1e+01 | 0.001 | 2.5e-17 | 8.8e-19 | 0 | 8.5e-19 |
| 1e+01 | 1e-06 | 5.3e-17 | 8.8e-22 | 2.2e-16 | 8.2e-22 |
| 1e+02 | 0.5 | 3.0e-17 | 4.1e-16 | 3.3e-16 | 2.4e-17 |
| 1e+02 | 0.1 | 5.4e-17 | 6.2e-16 | 4.7e-16 | 4.1e-17 |
| 1e+02 | 0.001 | 5.5e-17 | 7.0e-18 | 1.1e-16 | 7.0e-18 |
| 1e+02 | 1e-06 | 2.4e-17 | 5.4e-21 | 2.9e-15 | 5.4e-21 |
| 1e+03 | 0.5 | 4.8e-17 | 1.5e-16 | 3.4e-13 | 4.8e-17 |
| 1e+03 | 0.1 | 4.2e-17 | 8.0e-16 | 5.3e-13 | 3.6e-17 |
| 1e+03 | 0.001 | 4.9e-17 | 2.3e-16 | 1.1e-13 | 4.9e-17 |
| 1e+03 | 1e-06 | 5.0e-17 | 1.9e-19 | 2.9e-14 | 2.4e-19 |
| 1e+04 | 0.5 | 5.4e-17 | 7.0e-16 | 3.8e-13 | 5.4e-17 |
| 1e+04 | 0.1 | 5.3e-17 | 8.9e-16 | 1.2e-11 | 5.3e-17 |
| 1e+04 | 0.001 | 5.3e-17 | 1.2e-15 | 3.3e-12 | 3.5e-17 |
| 1e+04 | 1e-06 | 3.2e-17 | 7.9e-19 | 2.8e-13 | 7.9e-19 |
| 1e+05 | 0.5 | 5.1e-17 | 3.0e-15 | 9.1e-11 | 5.1e-17 |
| 1e+05 | 0.1 | 5.0e-17 | 2.5e-15 | 2.2e-11 | 5.0e-17 |
| 1e+05 | 0.001 | 3.2e-17 | 3.7e-16 | 1.0e-10 | 1.4e-16 |
| 1e+05 | 1e-06 | 2.9e-17 | 6.9e-18 | 2.6e-12 | 2.1e-17 |
| 1e+06 | 0.5 | 4.6e-17 | 1.2e-14 | 1.2e-09 | 4.6e-17 |
| 1e+06 | 0.1 | 5.5e-17 | 2.1e-14 | 5.9e-10 | 2.3e-17 |
| 1e+06 | 0.001 | 5.5e-17 | 6.6e-16 | 2.0e-10 | 5.5e-17 |
| 1e+06 | 1e-06 | 5.3e-17 | 2.2e-16 | 5.6e-11 | 1.5e-17 |
| 1e+07 | 0.5 | 4.2e-17 | 3.4e-14 | 0.00136 | 3.4e-17 |
| 1e+07 | 0.1 | 3.5e-17 | 3.3e-14 | 8.6e-07 | 1.8e-17 |
| 1e+07 | 0.001 | 3.9e-17 | 9.5e-16 | 1.0e-08 | 3.9e-17 |
| 1e+07 | 1e-06 | 3.2e-14 | 3.7e-11 | 3.9e-09 | 5.7e-17 |
| 1e+08 | 0.5 | 6.8e-17 | 2.2e-14 | 0.11 | 6.8e-17 |
| 1e+08 | 0.1 | 5.0e-17 | 2.1e-13 | 0.012 | 2.1e-17 |
| 1e+09 | 0.5 | 2.2e-16 | 2.0e-13 | 0.341 | 3.2e-16 |
| 1e+09 | 0.1 | 5.5e-17 | 4.8e-13 | 0.199 | 1.8e-16 |

## 2. Large n with small means (reference: 50-digit pmf recurrence)

The sampler uses the recurrence table for means up to 256 and Boost above.

| n | mean | regime | Boost F | Boost S | recurrence table | as consulted by the sampler |
|---:|---:|---|---:|---:|---:|---:|
| 1e+06 | 1 | table | 5.3e-17 | 2.2e-16 | 2.4e-16 | 1.5e-17 |
| 1e+06 | 3 | table | 7.9e-15 | 1.0e-11 | 6.4e-16 | 5.0e-17 |
| 1e+06 | 10 | table | 1.2e-14 | 2.1e-11 | 3.1e-16 | 1.5e-16 |
| 1e+06 | 30 | table | 1.1e-14 | 1.2e-11 | 3.0e-16 | 8.5e-17 |
| 1e+06 | 100 | table | 5.5e-17 | 4.3e-16 | 5.8e-16 | 1.9e-16 |
| 1e+06 | 256 | table | 5.3e-17 | 1.5e-15 | 1.2e-15 | 1.5e-16 |
| 1e+06 | 300 | boost | 5.5e-17 | 5.2e-16 | - | 4.8e-17 |
| 1e+06 | 1000 | boost | 5.5e-17 | 1.1e-15 | - | 5.2e-17 |
| 1e+07 | 1 | table | 5.0e-17 | 3.5e-17 | 2.1e-16 | 3.5e-17 |
| 1e+07 | 3 | table | 3.1e-14 | 2.0e-10 | 5.6e-16 | 1.4e-16 |
| 1e+07 | 10 | table | 1.1e-13 | 1.3e-10 | 4.9e-16 | 5.7e-17 |
| 1e+07 | 30 | table | 8.9e-14 | 1.2e-10 | 6.4e-16 | 4.8e-17 |
| 1e+07 | 100 | table | 5.5e-17 | 4.4e-16 | 8.3e-16 | 4.3e-16 |
| 1e+07 | 256 | table | 5.5e-17 | 6.4e-16 | 5.4e-16 | 4.5e-16 |
| 1e+07 | 300 | boost | 5.6e-17 | 8.4e-16 | - | 5.3e-17 |
| 1e+07 | 1000 | boost | 5.5e-17 | 1.3e-15 | - | 5.5e-17 |
| 1e+08 | 1 | table | 5.5e-17 | 8.2e-17 | 3.7e-16 | 6.7e-17 |
| 1e+08 | 3 | table | 6.6e-13 | 1.7e-09 | 1.6e-16 | 9.7e-17 |
| 1e+08 | 10 | table | 1.1e-13 | 2.4e-09 | 4.9e-16 | 1.4e-16 |
| 1e+08 | 30 | table | 3.5e-13 | 2.2e-09 | 4.7e-16 | 1.9e-16 |
| 1e+08 | 100 | table | 5.5e-17 | 7.9e-16 | 7.8e-16 | 1.6e-16 |
| 1e+08 | 256 | table | 5.5e-17 | 4.7e-16 | 5.9e-16 | 4.2e-16 |
| 1e+08 | 300 | boost | 5.5e-17 | 1.4e-15 | - | 4.9e-17 |
| 1e+08 | 1000 | boost | 5.4e-17 | 1.9e-15 | - | 5.0e-17 |
| 1e+09 | 1 | table | 5.3e-17 | 9.4e-17 | 3.3e-16 | 5.2e-17 |
| 1e+09 | 3 | table | 3.0e-12 | 1.1e-08 | 1.9e-16 | 6.2e-17 |
| 1e+09 | 10 | table | 2.4e-12 | 2.3e-08 | 3.3e-16 | 1.8e-16 |
| 1e+09 | 30 | table | 7.4e-12 | 1.9e-08 | 5.4e-16 | 9.6e-17 |
| 1e+09 | 100 | table | 5.4e-17 | 6.8e-16 | 5.6e-16 | 4.4e-16 |
| 1e+09 | 256 | table | 4.6e-17 | 1.0e-15 | 1.7e-15 | 3.1e-16 |
| 1e+09 | 300 | boost | 5.5e-17 | 1.1e-15 | - | 5.3e-17 |
| 1e+09 | 1000 | boost | 5.5e-17 | 9.6e-16 | - | 4.2e-17 |
| 1e+10 | 1 | table | 5.0e-17 | 1.8e-16 | 3.2e-16 | 3.8e-17 |
| 1e+10 | 3 | table | 5.3e-17 | 2.8e-16 | 7.6e-16 | 1.1e-16 |
| 1e+10 | 10 | table | 5.5e-17 | 6.9e-16 | 4.9e-16 | 1.9e-16 |
| 1e+10 | 30 | table | 5.3e-17 | 1.6e-15 | 6.4e-16 | 9.7e-17 |
| 1e+10 | 100 | table | 5.6e-17 | 1.1e-15 | 7.8e-16 | 2.3e-16 |
| 1e+10 | 256 | table | 5.5e-17 | 1.2e-15 | 5.6e-16 | 3.1e-16 |
| 1e+10 | 300 | boost | 5.5e-17 | 1.6e-15 | - | 5.5e-17 |
| 1e+10 | 1000 | boost | 5.3e-17 | 2.2e-15 | - | 5.3e-17 |
| 1e+12 | 1 | table | 5.0e-17 | 7.3e-17 | 1.6e-16 | 3.8e-17 |
| 1e+12 | 3 | table | 4.6e-17 | 2.6e-16 | 2.2e-16 | 2.5e-17 |
| 1e+12 | 10 | table | 5.4e-17 | 1.2e-15 | 2.4e-16 | 6.8e-17 |
| 1e+12 | 30 | table | 5.4e-17 | 1.1e-15 | 4.7e-16 | 2.4e-16 |
| 1e+12 | 100 | table | 5.5e-17 | 2.8e-16 | 9.1e-16 | 4.1e-16 |
| 1e+12 | 256 | table | 5.5e-17 | 8.0e-16 | 1.7e-15 | 5.8e-16 |
| 1e+12 | 300 | boost | 5.5e-17 | 1.2e-15 | - | 5.5e-17 |
| 1e+12 | 1000 | boost | 5.5e-17 | 1.5e-15 | - | 5.5e-17 |
| 1e+15 | 1 | table | 5.3e-17 | 3.2e-17 | 4.7e-16 | 8.9e-18 |
| 1e+15 | 3 | table | 4.6e-17 | 2.0e-16 | 3.4e-16 | 7.3e-17 |
| 1e+15 | 10 | table | 5.3e-17 | 8.4e-16 | 3.3e-16 | 1.3e-16 |
| 1e+15 | 30 | table | 5.5e-17 | 3.2e-15 | 6.4e-16 | 4.5e-16 |
| 1e+15 | 100 | table | 5.5e-17 | 7.1e-16 | 7.8e-16 | 5.1e-16 |
| 1e+15 | 256 | table | 5.5e-17 | 5.1e-16 | 1.0e-15 | 6.6e-17 |
| 1e+15 | 300 | boost | 5.5e-17 | 3.5e-16 | - | 5.5e-17 |
| 1e+15 | 1000 | boost | 5.5e-17 | 1.7e-15 | - | 5.5e-17 |
| 1e+15 | 1 | table | 5.3e-17 | 1.9e-16 | 1.7e-16 | 6.5e-17 |
| 1e+15 | 3 | table | 5.5e-17 | 2.0e-16 | 5.6e-16 | 7.4e-17 |
| 1e+15 | 10 | table | 5.5e-17 | 1.1e-15 | 4.4e-16 | 3.1e-17 |
| 1e+15 | 30 | table | 5.5e-17 | 8.3e-16 | 5.9e-16 | 6.7e-17 |
| 1e+15 | 100 | table | 5.5e-17 | 6.7e-16 | 4.3e-16 | 6.2e-17 |
| 1e+15 | 256 | table | 5.5e-17 | 6.5e-16 | 6.1e-16 | 7.8e-17 |
| 1e+15 | 300 | boost | 5.5e-17 | 7.0e-16 | - | 5.5e-17 |
| 1e+15 | 1000 | boost | 5.5e-17 | 5.7e-16 | - | 5.5e-17 |

## 3. n from 1e12 to 2^53, p = 1/2: exact values (centre value 1/2 plus exact pmf sums)

| n | z | exact F | Boost F error | Boost S error | error of the tail the sampler consults | Edgeworth error |
|---:|---:|---:|---:|---:|---:|---:|
| 1e+12 | 1.0 | 0.841345 | -9.4e-16 | 7.8e-12 | -9.4e-16 | -5.2e-17 |
| 1e+12 | 3.0 | 0.998650 | -1.2e-16 | 2.2e-13 | -1.2e-16 | -8.9e-15 |
| 1e+15 | 0.1 | 0.539828 | 2.6e-14 | 2.3e-12 | 2.6e-14 | -1.5e-17 |
| 1e+15 | 0.25 | 0.598706 | 9.3e-14 | -3.1e-11 | 9.3e-14 | 4.4e-17 |
| 1.126e+15 | 0.05 | 0.519939 | 1.8e-14 | 3.9e-11 | 1.8e-14 | 4.4e-17 |

Centre anchor and Edgeworth comparison (error of the reference is of order 1/n):

| n | Boost F error at the centre | Boost S error at the centre | max abs(Boost F - Edgeworth), p = 0.5 | p = 0.1 | ms per evaluation |
|---:|---:|---:|---:|---:|---:|
| 1e+09 | 0 | 2.5e-15 | 2.7e-11 | 1.2e-10 | 0.04 |
| 1e+12 | 0 | 0 | 2.6e-14 | 1.2e-13 | 0.15 |
| 1e+15 | 0 | 0 | 3.4e-13 | 2.6e-13 | 0.53 |
| 1.126e+15 | 0 | 0 | 3.1e-13 | 7.2e-14 | 0.52 |

## 4. Inversion against exact arithmetic

- Grid of uniforms, n up to 55, five values of p: 0 mismatches in 90045 draws against exact rational inversion at the cell midpoint.
- Uniforms 1e-12 either side of every distribution-function step: 0 mismatches in 1085 draws.
- Extreme uniforms (0, 2^-53, 1e-15, ..., 1 - 2^-53), n up to 1e9, both regimes: 0 mismatches in 80 draws against exact decimal quantiles.

- Uniforms adjacent to every distribution-function step (n = 1..100 at four p, and three Boost-regime cases): 948 of 27408 draws differ from exact inversion (the evaluated F puts the step on the other side of the midpoint); 0 of 26972 draws eight cells away do.
- Exact symmetry anchors (odd n, p = 1/2, n up to PART_MAX - 1): all exact: True; differences between brackets of half-width 1, 4 and 64: 0.
- Boost's F and S at every point the bisection can reach up to PART_MAX: 0 non-finite of 14001 points; beyond PART_MAX (near 2^53): 440 non-finite (PART_MAX = 1125899906842624).
- Every evaluation made by 1008 traced draws (n from 1e6 to PART_MAX, extreme and central uniforms, narrow brackets that force the fallbacks): 6325 of F and 1670 of S, 0 non-finite; standardised points reached from -1.68e+07 to 1.64e+07; draws refused: 0.
- Mass beyond the recurrence table: at most 10^-162.4.

## 5. The discretisation of the midpoint rule (exact F; exact decimal evaluation)

| n | p | mean | z at u = 0 | z at u = 1 - 2^-53 | relative variance bound (with tails) | exact relative variance error | old bound without tails | variance added by inverting at u itself |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2.78e-17 | 2.78e-17 | -5.3e-09 | -5.3e-09 | 1 | 1 | 3.1e-33 | 3.1e-33 |
| 1 | 8.88e-16 | 8.88e-16 | -3.0e-08 | 3.4e+07 | 0.0625 | 0 | 0.0625 | 9.9e-32 |
| 10 | 1e-17 | 1e-16 | -1.0e-08 | 1.0e+08 | 0.555 | 0.11 | 0.555 | 1.1e-32 |
| 1e+06 | 1e-20 | 1e-14 | -1.0e-07 | 1.0e+07 | 0.00555 | 8.0e-04 | 0.00555 | 1.1e-30 |
| 1e+06 | 1e-14 | 1e-08 | -1.0e-04 | 1e+04 | 2.1e-08 | 1.6e-08 | 5.6e-09 | 1.1e-24 |
| 1e+06 | 1e-09 | 0.001 | -0.0316 | 126 | 2.7e-09 | 5.2e-14 | 8.9e-13 | 1.1e-19 |
| 1000 | 0.01 | 10 | -3.18 | 11.1 | 7.9e-15 | 1.0e-16 | 7.4e-15 | 1.1e-15 |
| 1000 | 0.5 | 500 | -8.22 | 8.22 | 7.8e-15 | 4.2e-16 | 7.5e-15 | 1.1e-13 |
| 1e+06 | 0.5 | 5e+05 | -8.29 | 8.29 | 7.9e-15 | 1.2e-16 | 7.6e-15 | 1.1e-10 |
| 1e+09 | 1e-06 | 1e+03 | -7.94 | 8.63 | 7.9e-15 | 3.3e-16 | 7.6e-15 | 1.1e-13 |
| 1e+08 | 0.5 | 5e+07 | -8.29 | 8.29 | 7.9e-15 | 1.5e-16 | 7.6e-15 | 1.1e-08 |
| 1e+09 | 0.5 | 5e+08 | -8.29 | 8.29 | 7.9e-15 | 1.6e-16 | 7.6e-15 | 1.1e-07 |

The bound holds in every row: True.


## 6. Very small probabilities

P(X >= 1) = betainc(1, n, p) against 80-digit decimal for n in 1, 1e6, 1e12, 2^53 and p from 1e-300 to 1e-5 (28 cases): largest relative error 8.3e-17.

## 7. Interval laws (decimal closed-form Bateman)

Close decay constants (four-member chain, half-lives 1 + i x spacing, interval 1):

| relative spacing | digits used | max abs(decimal - matrix exponential) | float64 closed form (Layer A `bateman_chain`) | brackets hold |
|---:|---:|---:|---|---|
| 1e-01 | 50 | 3.5e-18 | 0.00346313932128 (decimal: 0.0034631393212) | True |
| 1e-03 | 50 | 1.1e-16 | refused | True |
| 1e-06 | 60 | 5.6e-17 | refused | True |
| 1e-09 | 69 | 2.8e-17 | refused | True |
| 1e-12 | 78 | 5.6e-17 | refused | True |
| 1e-14 | 84 | 5.6e-17 | refused | True |

Extreme rate ratios (half-lives from 1e-12 to 1e12, intervals from 1e-15 to 1e15, 25 cases): largest abs(decimal - matrix exponential) 8.9e-16; largest abs(sum - 1) 8.8e-50; most digits used 100; brackets hold in all: True.

Tiny intervals (S2 from Xa-400, deepest path of four decays):

| interval | path probability (decimal) | relative error of the float64 matrix exponential |
|---:|---:|---:|
| 1e-03 | 4.0561e-16 | 5.7e-13 |
| 1e-06 | 4.0576e-28 | 0 |
| 1e-09 | 4.0576e-40 | 6.0e-16 |
| 1e-12 | 4.0576e-52 | 3.7e-16 |
| 1e-15 | 4.0576e-64 | 6.7e-16 |

Stiff, close-rate chain S5: largest abs(float64 matrix exponential - 150-digit decimal) 1.6e-13; largest abs(default-precision decimal - 150-digit decimal) 3.4e-47.

## 8. numpy's binomial at large n (the calibration's null; not used by the B2 sampler)

Each row: 4,000,000 draws of numpy's `binomial(n, p)`; mean z is (sample mean - np) / its standard error; the variance ratio is the sample variance over np(1 - p), with standard error 0.0007 under the binomial law.

| n | p | expected count | mean z | variance ratio | (ratio - 1) / se |
|---:|---:|---:|---:|---:|---:|
| 1,000,000,000,000 | 1.15e-10 | 115 | -0.95 | 0.9986 | -2.0 |
| 1,000,000,000,000 | 1e-09 | 1000 | -0.37 | 1.0005 | 0.7 |
| 1,000,000,000,000 | 1e-08 | 1e+04 | 0.65 | 1.0006 | 0.8 |
| 100,000,000,000,000 | 1.15e-12 | 115 | 2.86 | 1.0006 | 0.9 |
| 100,000,000,000,000 | 1e-11 | 1000 | -0.39 | 1.0007 | 1.0 |
| 100,000,000,000,000 | 1e-10 | 1e+04 | -1.39 | 1.0001 | 0.1 |
| 1,000,000,000,000,000 | 1.15e-13 | 115 | 3.73 | 1.0084 | 11.9 |
| 1,000,000,000,000,000 | 1e-12 | 1000 | -7.93 | 1.0049 | 6.9 |
| 1,000,000,000,000,000 | 1e-11 | 1e+04 | 1.49 | 1.0007 | 1.0 |
| 9,007,199,254,740,991 | 1.27676e-14 | 115 | 44.46 | 1.0483 | 68.3 |
| 9,007,199,254,740,991 | 1.11022e-13 | 1000 | 26.70 | 1.0167 | 23.7 |
| 9,007,199,254,740,991 | 1.11022e-12 | 1e+04 | 7.53 | 1.0034 | 4.7 |
| 1,152,921,504,606,859,321 | 9.97466e-17 | 115 | -10.84 | 1.2212 | 312.8 |
| 1,152,921,504,606,859,321 | 8.67362e-16 | 1000 | -93.93 | 0.9978 | -3.1 |
| 1,152,921,504,606,859,321 | 8.67362e-15 | 1e+04 | 7.92 | 1.0016 | 2.3 |
| 1,152,921,504,606,859,321 | 8.67362e-13 | 1e+06 | 0.59 | 1.0019 | 2.6 |
| 1,152,921,504,606,859,321 | 8.67362e-11 | 1e+08 | -0.61 | 0.9998 | -0.3 |
| 1,152,921,504,606,859,321 | 8.67362e-09 | 1e+10 | 0.15 | 0.9997 | -0.5 |
| 22,680,697,271,813 | 0.9 | 2.268e+12 | -0.76 | 0.9994 | -0.9 |

Smallest expected count among the BIG null's binomial calls with n above 1e13: 2.268e+12 (n = 2.26807e+13, p = 0.899999999998; interval 2, Xa-400, multinomial category 7, initial size 1152921504606859321). The last row above is that call's n and p.

## 9. Populations from 1e6 to 2^62 (S2, six intervals, one run each; descriptive)

| initial nuclei | parts | seconds | words read | conserved | max abs z against exact means (54 values) |
|---:|---:|---:|---:|---|---:|
| 1e+06 | 1 | 0.04 | 150 | True | 2.65 |
| 1e+09 | 1 | 0.05 | 150 | True | 2.28 |
| 1e+12 | 1 | 0.08 | 150 | True | 2.42 |
| 1.126e+15 | 1 | 0.13 | 150 | True | 1.93 |
| 7.206e+16 | 64 | 1.38 | 9600 | True | 2.69 |
| 1.153e+18 | 1024 | 21.92 | 153600 | True | 2.31 |
| 4.612e+18 | 4096 | 88.42 | 614400 | True | 2.54 |

A total population above 2^63 - 1 is refused: True.
