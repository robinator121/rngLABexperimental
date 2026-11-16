# Layer A verification: tables from `results/verification.json`

Namespace `rnglab-v020-nuclear-a-verify`; code `c616a7b2d3b2cc9b80d66c756f55315c1193faf2` (dirty: False); started 2026-11-16T12:17:07.447765+00:00; 3.9 s in total.

## Verdicts

- R1: **PASS**
- R2: **MET**

| Deterministic check | Holds |
|---|---|
| D1 consumption | True |
| D2 invariants | True |
| D3 xref | True |
| D4 replay | True |
| D5 ensemble streams | True |
| A1 analytic cross-check | True |

A1: largest |expm − Bateman path sum| = 5.33e-15; largest |generator column sum| = 2.08e-17.
XREF: 5000 nuclei, events 15491 and 15491, discrete columns equal: True, final states equal: True, stream digests equal: True, largest relative time difference 0.
REPLAY (chunk 4099): events digest equal: True, stream digest equal: True.

## Families (Holm, family-wise alpha 0.01)

| Family | Tests | Smallest p | Holm rejections |
|---|---:|---:|---|
| NC1 | 47 | 0.0017 | none |
| NC2 | 33 | 0.00742 | none |

### NC1: every test

| Test | Statistic (z for binomial tests) | p | Holm-adjusted p |
|---|---:|---:|---:|
| NC1 S1 T1 | 0.00151 | 0.757 | 1 |
| NC1 S1 T2 t=0.5 | 1.01 | 0.313 | 1 |
| NC1 S1 T2 t=1.0 | 0.604 | 0.548 | 1 |
| NC1 S1 T2 t=2.0 | -0.351 | 0.727 | 1 |
| NC1 S1 T2 t=4.0 | -2.23 | 0.026 | 1 |
| NC1 S1 T2 t=6.0 | -2 | 0.0454 | 1 |
| NC1 S1 T3 | 3.14 | 0.0017 | 0.0797 |
| NC1 S2 T4 t=0.5 | 5.46 | 0.486 | 1 |
| NC1 S2 T4 t=1.0 | 12.6 | 0.0838 | 1 |
| NC1 S2 T4 t=2.0 | 8.53 | 0.384 | 1 |
| NC1 S2 T4 t=4.0 | 6.16 | 0.63 | 1 |
| NC1 S2 T4 t=8.0 | 7.03 | 0.534 | 1 |
| NC1 S2 T4 t=12.0 | 4.73 | 0.786 | 1 |
| NC1 S2 T5 Xa-400 | 0.217 | 0.838 | 1 |
| NC1 S2 T5 Xb-396 | -1.12 | 0.262 | 1 |
| NC1 S2 T5 Xc-400m | 1.13 | 0.261 | 1 |
| NC1 S2 T5 Xc-400 | 0.831 | 0.407 | 1 |
| NC1 S2 T5 Xd-396 | -0.748 | 0.455 | 1 |
| NC1 S2 T5 Xe-396 | 0.854 | 0.393 | 1 |
| NC1 S2 T6 Xa-400 | -1.03 | 0.302 | 1 |
| NC1 S2 T6 Xe-396 | -0.414 | 0.674 | 1 |
| NC1 S2 T7 Xa-400 | 9.19 | 0.687 | 1 |
| NC1 S2 T7 Xe-396 | 6.08 | 0.912 | 1 |
| NC1 S2 T8 Xa-400/alpha | -1.24 | 0.217 | 1 |
| NC1 S2 T8 Xa-400/beta- | -0.275 | 0.786 | 1 |
| NC1 S2 T8 Xb-396/ec | -0.07 | 0.944 | 1 |
| NC1 S2 T8 Xc-400m/it | 0.165 | 0.88 | 1 |
| NC1 S2 T8 Xe-396/beta- | 0.0273 | 0.984 | 1 |
| NC1 S3 T9 t=1.0 | 0.835 | 0.659 | 1 |
| NC1 S3 T9 t=5.0 | 0.754 | 0.686 | 1 |
| NC1 S3 T9 t=10.0 | 1.28 | 0.527 | 1 |
| NC1 S3 T9 t=20.0 | 0.33 | 0.848 | 1 |
| NC1 S3 T9 t=40.0 | 2.06 | 0.356 | 1 |
| ENS T10a Xa-400 | 0.271 | 0.788 | 1 |
| ENS T10b Xa-400 | 392 | 0.812 | 1 |
| ENS T10a Xb-396 | -0.714 | 0.477 | 1 |
| ENS T10b Xb-396 | 403 | 0.882 | 1 |
| ENS T10a Xc-400m | 1.01 | 0.314 | 1 |
| ENS T10b Xc-400m | 390 | 0.773 | 1 |
| ENS T10a Xc-400 | -0.0815 | 0.938 | 1 |
| ENS T10b Xc-400 | 401 | 0.931 | 1 |
| ENS T10a Xd-396 | -0.149 | 0.883 | 1 |
| ENS T10b Xd-396 | 436 | 0.201 | 1 |
| ENS T10a Xe-396 | -0.366 | 0.718 | 1 |
| ENS T10b Xe-396 | 432 | 0.245 | 1 |
| ENS T10a Xf-392 | -0.272 | 0.787 | 1 |
| ENS T10b Xf-392 | 430 | 0.274 | 1 |

### NC2: every test

| Test | Statistic (z for binomial tests) | p | Holm-adjusted p |
|---|---:|---:|---:|
| NC2 S1 T1 | 0.00176 | 0.578 | 1 |
| NC2 S1 T2 t=0.5 | -0.881 | 0.378 | 1 |
| NC2 S1 T2 t=1.0 | -0.3 | 0.766 | 1 |
| NC2 S1 T2 t=2.0 | 0.501 | 0.616 | 1 |
| NC2 S1 T2 t=4.0 | 1.54 | 0.123 | 1 |
| NC2 S1 T2 t=6.0 | -0.216 | 0.836 | 1 |
| NC2 S1 T3 | -2.68 | 0.00742 | 0.245 |
| NC2 S2 T4 t=0.5 | 3.2 | 0.784 | 1 |
| NC2 S2 T4 t=1.0 | 10.1 | 0.183 | 1 |
| NC2 S2 T4 t=2.0 | 8.07 | 0.427 | 1 |
| NC2 S2 T4 t=4.0 | 10.8 | 0.214 | 1 |
| NC2 S2 T4 t=8.0 | 4.32 | 0.827 | 1 |
| NC2 S2 T4 t=12.0 | 5.65 | 0.686 | 1 |
| NC2 S2 T5 Xa-400 | -0.0382 | 0.969 | 1 |
| NC2 S2 T5 Xb-396 | -0.142 | 0.888 | 1 |
| NC2 S2 T5 Xc-400m | 0.118 | 0.908 | 1 |
| NC2 S2 T5 Xc-400 | 0.236 | 0.812 | 1 |
| NC2 S2 T5 Xd-396 | -0.198 | 0.843 | 1 |
| NC2 S2 T5 Xe-396 | 0.995 | 0.32 | 1 |
| NC2 S2 T6 Xa-400 | -0.126 | 0.9 | 1 |
| NC2 S2 T6 Xe-396 | -1.64 | 0.102 | 1 |
| NC2 S2 T7 Xa-400 | 17.3 | 0.139 | 1 |
| NC2 S2 T7 Xe-396 | 5.43 | 0.942 | 1 |
| NC2 S2 T8 Xa-400/alpha | -0.626 | 0.534 | 1 |
| NC2 S2 T8 Xa-400/beta- | -0.45 | 0.655 | 1 |
| NC2 S2 T8 Xb-396/ec | -1.05 | 0.294 | 1 |
| NC2 S2 T8 Xc-400m/it | -0.888 | 0.374 | 1 |
| NC2 S2 T8 Xe-396/beta- | -0.0136 | 0.995 | 1 |
| NC2 S3 T9 t=1.0 | 0.59 | 0.745 | 1 |
| NC2 S3 T9 t=5.0 | 0.564 | 0.754 | 1 |
| NC2 S3 T9 t=10.0 | 0.167 | 0.92 | 1 |
| NC2 S3 T9 t=20.0 | 1.93 | 0.381 | 1 |
| NC2 S3 T9 t=40.0 | 2.96 | 0.228 | 1 |

## Positive controls

Level: p < 2.13e-04.

| Control | Target | p | Rejected |
|---|---|---:|---|
| PC1 | T2 t=2.0 | 9.48e-09 | True |
| PC2 | T6 Xa-400 | 2.54e-34 | True |

PC1, every test (descriptive, outside the families):

| Test | p |
|---|---:|
| T1 | 7.81e-06 |
| T2 t=0.5 | 0.00184 |
| T2 t=1.0 | 1.07e-05 |
| T2 t=2.0 | 9.48e-09 |
| T2 t=4.0 | 1.80e-08 |
| T2 t=6.0 | 1.21e-07 |
| T3 | 0.675 |

PC2, every test (descriptive, outside the families):

| Test | p |
|---|---:|
| T4 t=0.5 | 0.00307 |
| T4 t=1.0 | 1.54e-06 |
| T4 t=2.0 | 1.16e-17 |
| T4 t=4.0 | 4.26e-24 |
| T4 t=8.0 | 7.03e-27 |
| T4 t=12.0 | 3.67e-30 |
| T5 Xa-400 | 0.929 |
| T5 Xb-396 | 1.05e-32 |
| T5 Xc-400m | 8.13e-32 |
| T5 Xc-400 | 3.60e-25 |
| T5 Xd-396 | 4.35e-29 |
| T5 Xe-396 | 7.17e-21 |
| T6 Xa-400 | 2.54e-34 |
| T6 Xe-396 | 0.875 |
| T7 Xa-400 | 0.4 |
| T7 Xe-396 | 9.05e-16 |
| T8 Xa-400/alpha | 0.614 |
| T8 Xa-400/beta- | 0.901 |
| T8 Xb-396/ec | 0.859 |
| T8 Xc-400m/it | 0.187 |
| T8 Xe-396/beta- | 0.894 |

## R3 (descriptive): transient equilibrium, NC1 S3

Decays in (20, 40]: Ya-410 18705, Yb-410 20807; ratio 1.1124, approximate 95% interval [1.1059, 1.1189] (standard error of the log ratio 0.00296; correlation of the per-nucleus window indicators 0.894); exact expectation 1.111111; equilibrium limit 1.111111.

## Single-machine timings (seconds, wall clock; not HPC measurements)

| Run | Seconds |
|---|---:|
| NC1 S1 | 0.039 |
| NC1 S2 | 0.049 |
| NC1 S3 | 0.022 |
| NC2 S1 | 0.475 |
| NC2 S2 | 0.961 |
| NC2 S3 | 0.503 |
| PC1 S1 | 0.020 |
| PC2 S2 | 0.049 |
| ENS | 0.636 |
| XREF reference engine | 0.042 |
