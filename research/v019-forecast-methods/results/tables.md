# Methods validation study: tables (frozen study)

Source: results.json, namespace `rnglab-v019-forecast-methods`, protocol SHA-256 `4b3b6a24d5e079dd5da539e7a7113ba4df14bdcd1d38d924bdc0426b15ff0dba`.

- alarm threshold 3.4038; estimated sigma_eff 0.6125, beta 0.05984, kappa 102.998
- healthy valid fraction per length: 10240: 0.9309, 20480: 0.9429, 40960: 0.9389

## A_declared

Runs 2000, alarms 2000, no alarm 0, pre-onset alarms 8.

**L = 10240**: eligible 1734, late alarm 266, valid fraction 0.8495

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.1279 | 0.4238 | -0.142 [-0.186, -0.099] | 0.78 | 0.007 | yes |
| alarm-time | 0.1120 | 0.3711 | 0.000 [0.000, 0.000] | 0.75 | 0.033 | yes |
| lite-logistic | 0.1068 | 0.3569 | 0.047 [0.026, 0.067] | 0.35 | 0.034 | yes |
| persistence | 0.1360 | 0.4763 | -0.214 [-0.266, -0.161] | 17.70 | 0.097 | no |
| bayes-ramp-no-widening | 0.1037 | 0.3473 | 0.074 [0.032, 0.115] | -0.28 | 0.012 | yes |
| bayes-ramp | 0.1034 | 0.3465 | 0.076 [0.044, 0.109] | 0.89 | 0.013 | yes |
| oracle | 0.0985 | 0.3320 | 0.120 [0.084, 0.158] | 0.61 | 0.011 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00026 [-0.00188, 0.00141]

bayes-ramp-minus-lite-logistic: -0.00334 [-0.00600, -0.00069]

lite-logistic-minus-alarm-time: -0.00522 [-0.00758, -0.00290]

**L = 20480**: eligible 2000, late alarm 0, valid fraction 0.5030

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.2500 | 0.6931 | -0.379 [-0.445, -0.320] | -0.07 | 0.001 | yes |
| alarm-time | 0.1812 | 0.5339 | 0.000 [0.000, 0.000] | -0.57 | 0.040 | yes |
| lite-logistic | 0.1531 | 0.4633 | 0.155 [0.123, 0.186] | -1.07 | 0.025 | yes |
| persistence | 0.4283 | 1.3354 | -1.363 [-1.519, -1.221] | 75.86 | 0.430 | no |
| bayes-ramp-no-widening | 0.1431 | 0.4423 | 0.210 [0.166, 0.254] | 3.55 | 0.045 | no |
| bayes-ramp | 0.1399 | 0.4316 | 0.228 [0.189, 0.268] | -0.20 | 0.022 | yes |
| oracle | 0.1153 | 0.3720 | 0.364 [0.317, 0.409] | -0.12 | 0.015 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00328 [-0.00561, -0.00092]

bayes-ramp-minus-lite-logistic: -0.01324 [-0.01805, -0.00841]

lite-logistic-minus-alarm-time: -0.02814 [-0.03386, -0.02227]

**L = 40960**: eligible 2000, late alarm 0, valid fraction 0.1830

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.1495 | 0.4760 | -0.574 [-0.690, -0.468] | -0.43 | 0.004 | yes |
| alarm-time | 0.0950 | 0.3145 | 0.000 [0.000, 0.000] | -1.06 | 0.023 | yes |
| lite-logistic | 0.0927 | 0.2961 | 0.025 [0.001, 0.049] | -0.06 | 0.012 | yes |
| persistence | 0.6508 | 1.8578 | -5.849 [-6.646, -5.173] | 105.64 | 0.712 | no |
| bayes-ramp-no-widening | 0.0813 | 0.2556 | 0.144 [0.107, 0.181] | 1.35 | 0.026 | yes |
| bayes-ramp | 0.0799 | 0.2520 | 0.160 [0.121, 0.196] | -1.27 | 0.018 | yes |
| oracle | 0.0577 | 0.1782 | 0.392 [0.337, 0.445] | -2.14 | 0.018 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00143 [-0.00260, -0.00027]

bayes-ramp-minus-lite-logistic: -0.01281 [-0.01667, -0.00900]

lite-logistic-minus-alarm-time: -0.00235 [-0.00466, -0.00007]

**Failure time**: 2000 runs with an alarm, 271 censored at the horizon, 241 already failed at issuance; Kaplan-Meier median failure sweep 17920.0, median remaining sweeps 11520.0.

| forecaster | RPS | 80% interval coverage [CP95] | upper end censored | median abs. error (sweeps) | both inside | forecast inside, outcome censored | forecast censored, outcome inside | both censored |
|---|---|---|---|---|---|---|---|---|
| empirical | 0.1385 | 0.845 [0.828, 0.861] | 801 | 5376.0 | 1729 | 271 | 0 | 0 |
| bayes-ramp | 0.1232 | 0.841 [0.824, 0.857] | 692 | 4352.0 | 1665 | 153 | 64 | 118 |

RPS skill of bayes-ramp over empirical: 0.111 [0.092, 0.129]

## B1_slow_rates

Runs 2000, alarms 2000, no alarm 0, pre-onset alarms 12.

**L = 10240**: eligible 99, late alarm 1901, valid fraction 0.9394

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.0639 | 0.2635 | 0.220 [-0.014, 0.411] | -2.36 | 0.083 | no |
| alarm-time | 0.0819 | 0.2901 | 0.000 [0.000, 0.000] | 0.40 | 0.010 | yes |
| lite-logistic | 0.0727 | 0.2719 | 0.112 [0.030, 0.240] | 0.49 | 0.026 | yes |
| persistence | 0.0570 | 0.2297 | 0.304 [0.067, 0.633] | 0.41 | 0.009 | yes |
| bayes-ramp-no-widening | 0.0574 | 0.2301 | 0.300 [0.118, 0.559] | -0.54 | 0.006 | yes |
| bayes-ramp | 0.0566 | 0.2273 | 0.309 [0.125, 0.585] | -0.41 | 0.008 | yes |
| oracle | 0.0570 | 0.2298 | 0.304 [0.059, 0.630] | 0.47 | 0.010 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00078 [-0.00274, 0.00130]

bayes-ramp-minus-lite-logistic: -0.01612 [-0.02828, -0.00606]

lite-logistic-minus-alarm-time: -0.00922 [-0.01635, -0.00279]

**L = 20480**: eligible 1368, late alarm 632, valid fraction 0.9393

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.2480 | 0.6892 | -2.600 [-3.347, -2.052] | -32.33 | 0.437 | no |
| alarm-time | 0.0689 | 0.3066 | 0.000 [0.000, 0.000] | 8.82 | 0.033 | no |
| lite-logistic | 0.0658 | 0.2934 | 0.045 [0.024, 0.069] | 8.79 | 0.032 | no |
| persistence | 0.0569 | 0.2284 | 0.173 [0.101, 0.250] | 1.33 | 0.008 | yes |
| bayes-ramp-no-widening | 0.0641 | 0.2442 | 0.070 [0.028, 0.111] | 0.69 | 0.001 | yes |
| bayes-ramp | 0.0625 | 0.2399 | 0.093 [0.050, 0.137] | 0.50 | 0.001 | yes |
| oracle | 0.0570 | 0.2293 | 0.172 [0.096, 0.248] | 1.32 | 0.008 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00160 [-0.00234, -0.00090]

bayes-ramp-minus-lite-logistic: -0.00327 [-0.00594, -0.00072]

lite-logistic-minus-alarm-time: -0.00312 [-0.00473, -0.00157]

**L = 40960**: eligible 2000, late alarm 0, valid fraction 0.9195

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.6109 | 1.5595 | -4.461 [-5.140, -3.902] | 84.09 | 0.733 | no |
| alarm-time | 0.1119 | 0.5316 | 0.000 [0.000, 0.000] | 19.36 | 0.071 | no |
| lite-logistic | 0.1099 | 0.5037 | 0.017 [0.003, 0.032] | 17.35 | 0.069 | no |
| persistence | 0.0740 | 0.2801 | 0.338 [0.280, 0.396] | 2.29 | 0.013 | yes |
| bayes-ramp-no-widening | 0.1030 | 0.3702 | 0.079 [0.053, 0.102] | -4.51 | 0.065 | no |
| bayes-ramp | 0.0995 | 0.3579 | 0.111 [0.087, 0.133] | -4.25 | 0.059 | no |
| oracle | 0.0731 | 0.2750 | 0.346 [0.292, 0.401] | 0.13 | 0.008 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00355 [-0.00412, -0.00300]

bayes-ramp-minus-lite-logistic: -0.01045 [-0.01378, -0.00718]

lite-logistic-minus-alarm-time: -0.00194 [-0.00354, -0.00031]

**Failure time**: 2000 runs with an alarm, 1378 censored at the horizon, 390 already failed at issuance; Kaplan-Meier median failure sweep None, median remaining sweeps None.

| forecaster | RPS | 80% interval coverage [CP95] | upper end censored | median abs. error (sweeps) | both inside | forecast inside, outcome censored | forecast censored, outcome inside | both censored |
|---|---|---|---|---|---|---|---|---|
| empirical | 0.1695 | 0.893 [0.879, 0.907] | 1959 | 24320.0 | 622 | 1378 | 0 | 0 |
| bayes-ramp | 0.1588 | 0.878 [0.863, 0.893] | 1973 | 13440.0 | 46 | 81 | 576 | 1297 |

RPS skill of bayes-ramp over empirical: 0.063 [0.044, 0.084]

## B2_late_onsets

Runs 2000, alarms 2000, no alarm 0, pre-onset alarms 22.

**L = 10240**: eligible 150, late alarm 1850, valid fraction 0.9467

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.0587 | 0.2505 | 0.069 [-0.293, 0.263] | -3.16 | 0.091 | no |
| alarm-time | 0.0630 | 0.2409 | 0.000 [0.000, 0.000] | 1.74 | 0.024 | yes |
| lite-logistic | 0.0578 | 0.2232 | 0.083 [0.024, 0.163] | 0.90 | 0.019 | yes |
| persistence | 0.0504 | 0.2077 | 0.200 [-0.066, 0.480] | 0.15 | 0.003 | yes |
| bayes-ramp-no-widening | 0.0541 | 0.2175 | 0.142 [-0.086, 0.371] | -0.21 | 0.000 | yes |
| bayes-ramp | 0.0530 | 0.2152 | 0.158 [-0.078, 0.403] | -0.19 | 0.001 | yes |
| oracle | 0.0505 | 0.2082 | 0.199 [-0.073, 0.464] | 0.14 | 0.003 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00104 [-0.00231, -0.00002]

bayes-ramp-minus-lite-logistic: -0.00473 [-0.01563, 0.00599]

lite-logistic-minus-alarm-time: -0.00526 [-0.00978, -0.00129]

**L = 20480**: eligible 1913, late alarm 87, valid fraction 0.8191

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.2486 | 0.6903 | -0.601 [-0.765, -0.465] | -27.72 | 0.317 | no |
| alarm-time | 0.1553 | 0.5725 | 0.000 [0.000, 0.000] | 30.49 | 0.131 | no |
| lite-logistic | 0.1415 | 0.4994 | 0.089 [0.047, 0.127] | 19.49 | 0.087 | no |
| persistence | 0.1629 | 0.5614 | -0.049 [-0.076, -0.023] | 24.14 | 0.126 | no |
| bayes-ramp-no-widening | 0.1524 | 0.5147 | 0.018 [-0.001, 0.038] | 19.30 | 0.106 | no |
| bayes-ramp | 0.1537 | 0.5195 | 0.010 [-0.007, 0.027] | 19.84 | 0.109 | no |
| oracle | 0.1205 | 0.3916 | 0.224 [0.182, 0.264] | 2.23 | 0.018 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: +0.00134 [0.00046, 0.00220]

bayes-ramp-minus-lite-logistic: +0.01222 [0.00539, 0.01919]

lite-logistic-minus-alarm-time: -0.01375 [-0.02064, -0.00701]

**L = 40960**: eligible 2000, late alarm 0, valid fraction 0.3110

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.2297 | 0.6643 | 0.491 [0.444, 0.534] | 14.26 | 0.124 | no |
| alarm-time | 0.4514 | 1.5253 | 0.000 [0.000, 0.000] | 92.70 | 0.525 | no |
| lite-logistic | 0.3501 | 1.1283 | 0.224 [0.208, 0.241] | 60.03 | 0.419 | no |
| persistence | 0.5752 | 1.7179 | -0.274 [-0.305, -0.246] | 98.66 | 0.606 | no |
| bayes-ramp-no-widening | 0.4101 | 1.1224 | 0.092 [0.082, 0.101] | 53.53 | 0.464 | no |
| bayes-ramp | 0.4242 | 1.1659 | 0.060 [0.050, 0.070] | 56.66 | 0.481 | no |
| oracle | 0.0850 | 0.2733 | 0.812 [0.792, 0.830] | -1.07 | 0.015 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: +0.01408 [0.01323, 0.01498]

bayes-ramp-minus-lite-logistic: +0.07405 [0.06519, 0.08311]

lite-logistic-minus-alarm-time: -0.10128 [-0.10959, -0.09309]

**Failure time**: 2000 runs with an alarm, 460 censored at the horizon, 342 already failed at issuance; Kaplan-Meier median failure sweep 27136.0, median remaining sweeps 13568.0.

| forecaster | RPS | 80% interval coverage [CP95] | upper end censored | median abs. error (sweeps) | both inside | forecast inside, outcome censored | forecast censored, outcome inside | both censored |
|---|---|---|---|---|---|---|---|---|
| empirical | 0.1891 | 0.905 [0.892, 0.918] | 1987 | 13312.0 | 1540 | 460 | 0 | 0 |
| bayes-ramp | 0.2104 | 0.890 [0.875, 0.903] | 1918 | 10752.0 | 297 | 8 | 1243 | 452 |

RPS skill of bayes-ramp over empirical: -0.113 [-0.126, -0.099]

## C_step_hold

Runs 2000, alarms 1804, no alarm 196, pre-onset alarms 12.

**L = 10240**: eligible 1514, late alarm 290, valid fraction 0.8692

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.1138 | 0.3886 | 0.042 [-0.011, 0.095] | -1.46 | 0.013 | yes |
| alarm-time | 0.1189 | 0.3938 | 0.000 [0.000, 0.000] | -4.02 | 0.068 | no |
| lite-logistic | 0.1755 | 0.5635 | -0.476 [-0.608, -0.360] | 5.95 | 0.210 | no |
| persistence | 0.1113 | 0.3818 | 0.064 [0.003, 0.128] | 9.25 | 0.058 | no |
| bayes-ramp-no-widening | 0.1858 | 0.5523 | -0.563 [-0.688, -0.451] | 6.96 | 0.218 | no |
| bayes-ramp | 0.1539 | 0.4704 | -0.294 [-0.380, -0.215] | 1.41 | 0.162 | no |
| oracle | 0.0928 | 0.3206 | 0.220 [0.169, 0.272] | 0.00 | 0.006 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.03196 [-0.03507, -0.02883]

bayes-ramp-minus-lite-logistic: -0.02159 [-0.03169, -0.01157]

lite-logistic-minus-alarm-time: +0.05658 [0.04416, 0.06917]

**L = 20480**: eligible 1682, late alarm 122, valid fraction 0.7140

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.2490 | 0.6912 | 0.206 [0.171, 0.238] | -17.37 | 0.212 | no |
| alarm-time | 0.3136 | 0.8795 | 0.000 [0.000, 0.000] | 32.57 | 0.322 | no |
| lite-logistic | 0.3772 | 1.4896 | -0.203 [-0.243, -0.164] | 70.82 | 0.424 | no |
| persistence | 0.1985 | 0.6031 | 0.367 [0.302, 0.431] | 20.12 | 0.162 | no |
| bayes-ramp-no-widening | 0.4329 | 1.9297 | -0.380 [-0.430, -0.333] | 97.38 | 0.480 | no |
| bayes-ramp | 0.3938 | 1.4858 | -0.256 [-0.299, -0.214] | 70.84 | 0.434 | no |
| oracle | 0.1066 | 0.3530 | 0.660 [0.623, 0.695] | -0.50 | 0.020 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.03904 [-0.04178, -0.03641]

bayes-ramp-minus-lite-logistic: +0.01662 [0.00951, 0.02369]

lite-logistic-minus-alarm-time: +0.06363 [0.05153, 0.07591]

**L = 40960**: eligible 1803, late alarm 1, valid fraction 0.6112

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.4178 | 1.1060 | -0.008 [-0.044, 0.024] | 46.25 | 0.424 | no |
| alarm-time | 0.4144 | 1.5681 | 0.000 [0.000, 0.000] | 88.12 | 0.428 | no |
| lite-logistic | 0.4208 | 1.9351 | -0.016 [-0.023, -0.008] | 122.61 | 0.452 | no |
| persistence | 0.2027 | 0.5997 | 0.511 [0.459, 0.558] | 15.73 | 0.181 | no |
| bayes-ramp-no-widening | 0.4295 | 3.1476 | -0.037 [-0.049, -0.024] | 122.48 | 0.452 | no |
| bayes-ramp | 0.4214 | 2.9788 | -0.017 [-0.030, -0.003] | 116.72 | 0.442 | no |
| oracle | 0.0972 | 0.3248 | 0.765 [0.738, 0.791] | 0.20 | 0.024 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00809 [-0.00900, -0.00721]

bayes-ramp-minus-lite-logistic: +0.00062 [-0.00410, 0.00536]

lite-logistic-minus-alarm-time: +0.00642 [0.00321, 0.00971]

**Failure time**: 1804 runs with an alarm, 812 censored at the horizon, 199 already failed at issuance; Kaplan-Meier median failure sweep 32000.0, median remaining sweeps 23808.0.

| forecaster | RPS | 80% interval coverage [CP95] | upper end censored | median abs. error (sweeps) | both inside | forecast inside, outcome censored | forecast censored, outcome inside | both censored |
|---|---|---|---|---|---|---|---|---|
| empirical | 0.2505 | 0.585 [0.562, 0.608] | 525 | 7168.0 | 992 | 812 | 0 | 0 |
| bayes-ramp | 0.2752 | 0.450 [0.427, 0.473] | 449 | 6656.0 | 917 | 595 | 75 | 217 |

RPS skill of bayes-ramp over empirical: -0.098 [-0.118, -0.079]

## D_random_walk

Runs 2000, alarms 2000, no alarm 0, pre-onset alarms None.

**L = 10240**: eligible 1962, late alarm 38, valid fraction 0.9332

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.0683 | 0.2745 | 0.140 [0.100, 0.178] | -9.74 | 0.077 | no |
| alarm-time | 0.0794 | 0.2963 | 0.000 [0.000, 0.000] | -8.20 | 0.088 | no |
| lite-logistic | 0.0791 | 0.2945 | 0.004 [-0.017, 0.024] | -8.06 | 0.085 | no |
| persistence | 0.0624 | 0.2456 | 0.214 [0.155, 0.277] | 2.51 | 0.013 | yes |
| bayes-ramp-no-widening | 0.0989 | 0.3355 | -0.246 [-0.305, -0.193] | -4.34 | 0.106 | no |
| bayes-ramp | 0.0858 | 0.3045 | -0.081 [-0.118, -0.046] | -5.25 | 0.080 | no |
| oracle | 0.0618 | 0.2404 | 0.222 [0.172, 0.277] | -0.86 | 0.005 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.01311 [-0.01461, -0.01161]

bayes-ramp-minus-lite-logistic: +0.00671 [0.00451, 0.00899]

lite-logistic-minus-alarm-time: -0.00029 [-0.00190, 0.00134]

**L = 20480**: eligible 2000, late alarm 0, valid fraction 0.6750

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.2492 | 0.6916 | 0.163 [0.133, 0.192] | -15.45 | 0.173 | no |
| alarm-time | 0.2978 | 0.8222 | 0.000 [0.000, 0.000] | 25.38 | 0.266 | no |
| lite-logistic | 0.3461 | 0.9991 | -0.162 [-0.193, -0.132] | 41.60 | 0.313 | no |
| persistence | 0.2831 | 0.9070 | 0.049 [-0.032, 0.124] | 44.48 | 0.254 | no |
| bayes-ramp-no-widening | 0.4087 | 1.3392 | -0.372 [-0.416, -0.330] | 65.89 | 0.390 | no |
| bayes-ramp | 0.3628 | 1.0929 | -0.218 [-0.252, -0.185] | 49.24 | 0.336 | no |
| oracle | 0.1719 | 0.5191 | 0.423 [0.385, 0.459] | 0.75 | 0.017 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.04587 [-0.04927, -0.04255]

bayes-ramp-minus-lite-logistic: +0.01676 [0.01014, 0.02353]

lite-logistic-minus-alarm-time: +0.04827 [0.03942, 0.05682]

**L = 40960**: eligible 2000, late alarm 0, valid fraction 0.1630

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.1370 | 0.4465 | 0.052 [0.023, 0.078] | -2.73 | 0.024 | yes |
| alarm-time | 0.1445 | 0.5221 | 0.000 [0.000, 0.000] | 15.52 | 0.088 | no |
| lite-logistic | 0.1467 | 0.5820 | -0.015 [-0.035, 0.005] | 22.45 | 0.101 | no |
| persistence | 0.6620 | 1.8744 | -3.580 [-4.113, -3.127] | 105.05 | 0.726 | no |
| bayes-ramp-no-widening | 0.1543 | 0.8682 | -0.067 [-0.093, -0.042] | 37.82 | 0.125 | no |
| bayes-ramp | 0.1545 | 0.8120 | -0.069 [-0.096, -0.042] | 33.74 | 0.119 | no |
| oracle | 0.0862 | 0.2726 | 0.403 [0.358, 0.447] | 0.47 | 0.010 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: +0.00022 [-0.00053, 0.00096]

bayes-ramp-minus-lite-logistic: +0.00784 [0.00474, 0.01107]

lite-logistic-minus-alarm-time: +0.00212 [-0.00063, 0.00507]

**Failure time**: 2000 runs with an alarm, 214 censored at the horizon, 244 already failed at issuance; Kaplan-Meier median failure sweep 22016.0, median remaining sweeps 15872.0.

| forecaster | RPS | 80% interval coverage [CP95] | upper end censored | median abs. error (sweeps) | both inside | forecast inside, outcome censored | forecast censored, outcome inside | both censored |
|---|---|---|---|---|---|---|---|---|
| empirical | 0.1867 | 0.690 [0.669, 0.710] | 561 | 8448.0 | 1786 | 214 | 0 | 0 |
| bayes-ramp | 0.2047 | 0.522 [0.500, 0.544] | 307 | 8704.0 | 1775 | 210 | 11 | 4 |

RPS skill of bayes-ramp over empirical: -0.096 [-0.115, -0.078]

## E_continue

Runs 2000, alarms 2000, no alarm 0, pre-onset alarms 9.

**L = 10240**: eligible 1710, late alarm 290, valid fraction 0.8409

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.1340 | 0.4390 | -0.144 [-0.186, -0.103] | 1.78 | 0.015 | yes |
| alarm-time | 0.1171 | 0.3848 | 0.000 [0.000, 0.000] | 1.74 | 0.026 | yes |
| lite-logistic | 0.1121 | 0.3731 | 0.043 [0.021, 0.065] | 1.24 | 0.027 | yes |
| persistence | 0.1435 | 0.4999 | -0.225 [-0.275, -0.175] | 19.09 | 0.105 | no |
| bayes-ramp-no-widening | 0.1073 | 0.3562 | 0.084 [0.042, 0.125] | -0.08 | 0.014 | yes |
| bayes-ramp | 0.1073 | 0.3564 | 0.084 [0.053, 0.114] | 1.10 | 0.015 | yes |
| oracle | 0.1032 | 0.3454 | 0.119 [0.084, 0.154] | 1.03 | 0.012 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: +0.00006 [-0.00172, 0.00191]

bayes-ramp-minus-lite-logistic: -0.00473 [-0.00736, -0.00214]

lite-logistic-minus-alarm-time: -0.00507 [-0.00779, -0.00247]

**L = 20480**: eligible 2000, late alarm 0, valid fraction 0.5205

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.2499 | 0.6930 | -0.370 [-0.437, -0.306] | -1.63 | 0.018 | yes |
| alarm-time | 0.1824 | 0.5475 | 0.000 [0.000, 0.000] | 0.55 | 0.056 | no |
| lite-logistic | 0.1563 | 0.4841 | 0.143 [0.113, 0.174] | 0.56 | 0.035 | yes |
| persistence | 0.4125 | 1.2864 | -1.261 [-1.412, -1.115] | 72.55 | 0.412 | no |
| bayes-ramp-no-widening | 0.1442 | 0.4509 | 0.210 [0.167, 0.252] | 4.83 | 0.059 | no |
| bayes-ramp | 0.1392 | 0.4343 | 0.237 [0.199, 0.275] | 0.86 | 0.026 | yes |
| oracle | 0.1208 | 0.3889 | 0.338 [0.291, 0.383] | 2.19 | 0.021 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00504 [-0.00741, -0.00274]

bayes-ramp-minus-lite-logistic: -0.01715 [-0.02150, -0.01276]

lite-logistic-minus-alarm-time: -0.02614 [-0.03182, -0.02035]

**L = 40960**: eligible 2000, late alarm 0, valid fraction 0.1885

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.1530 | 0.4841 | -0.566 [-0.684, -0.457] | 0.20 | 0.002 | yes |
| alarm-time | 0.0977 | 0.3229 | 0.000 [0.000, 0.000] | -1.04 | 0.022 | yes |
| lite-logistic | 0.0954 | 0.3045 | 0.023 [-0.001, 0.047] | -0.16 | 0.010 | yes |
| persistence | 0.6449 | 1.8424 | -5.602 [-6.340, -4.976] | 104.73 | 0.707 | no |
| bayes-ramp-no-widening | 0.0866 | 0.2734 | 0.114 [0.075, 0.151] | 1.30 | 0.023 | yes |
| bayes-ramp | 0.0863 | 0.2723 | 0.117 [0.076, 0.156] | -0.80 | 0.019 | yes |
| oracle | 0.0652 | 0.2036 | 0.333 [0.273, 0.387] | -0.72 | 0.021 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.00029 [-0.00145, 0.00085]

bayes-ramp-minus-lite-logistic: -0.00910 [-0.01316, -0.00508]

lite-logistic-minus-alarm-time: -0.00229 [-0.00467, 0.00010]

**Failure time**: 2000 runs with an alarm, 268 censored at the horizon, 239 already failed at issuance; Kaplan-Meier median failure sweep 18176.0, median remaining sweeps 11264.0.

| forecaster | RPS | 80% interval coverage [CP95] | upper end censored | median abs. error (sweeps) | both inside | forecast inside, outcome censored | forecast censored, outcome inside | both censored |
|---|---|---|---|---|---|---|---|---|
| empirical | 0.1373 | 0.855 [0.838, 0.870] | 841 | 5376.0 | 1732 | 268 | 0 | 0 |
| bayes-ramp | 0.1226 | 0.846 [0.829, 0.862] | 746 | 4352.0 | 1649 | 147 | 83 | 121 |

RPS skill of bayes-ramp over empirical: 0.107 [0.088, 0.125]

## E_freeze

Runs 2000, alarms 2000, no alarm 0, pre-onset alarms 9.

**L = 10240**: eligible 1710, late alarm 290, valid fraction 0.9386

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.0644 | 0.2649 | 0.165 [0.128, 0.202] | -9.72 | 0.083 | no |
| alarm-time | 0.0771 | 0.2923 | 0.000 [0.000, 0.000] | -6.30 | 0.088 | no |
| lite-logistic | 0.0812 | 0.3022 | -0.053 [-0.082, -0.024] | -5.07 | 0.083 | no |
| persistence | 0.0577 | 0.2319 | 0.251 [0.195, 0.314] | 1.44 | 0.008 | yes |
| bayes-ramp-no-widening | 0.1031 | 0.3449 | -0.337 [-0.414, -0.269] | -0.92 | 0.078 | no |
| bayes-ramp | 0.0887 | 0.3108 | -0.151 [-0.199, -0.107] | -2.64 | 0.058 | no |
| oracle | 0.0579 | 0.2327 | 0.249 [0.193, 0.310] | 0.62 | 0.005 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.01439 [-0.01616, -0.01268]

bayes-ramp-minus-lite-logistic: +0.00755 [0.00512, 0.00997]

lite-logistic-minus-alarm-time: +0.00406 [0.00191, 0.00626]

**L = 20480**: eligible 2000, late alarm 0, valid fraction 0.9190

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.2481 | 0.6894 | 0.195 [0.164, 0.224] | -37.28 | 0.417 | no |
| alarm-time | 0.3084 | 0.8593 | 0.000 [0.000, 0.000] | 31.04 | 0.420 | no |
| lite-logistic | 0.3358 | 1.0380 | -0.089 [-0.117, -0.062] | 43.70 | 0.424 | no |
| persistence | 0.0741 | 0.2798 | 0.760 [0.724, 0.793] | 2.34 | 0.013 | yes |
| bayes-ramp-no-widening | 0.3934 | 1.4304 | -0.276 [-0.315, -0.238] | 63.54 | 0.449 | no |
| bayes-ramp | 0.3454 | 1.1361 | -0.120 [-0.152, -0.088] | 48.28 | 0.407 | no |
| oracle | 0.0740 | 0.2790 | 0.760 [0.727, 0.791] | -0.29 | 0.015 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.04807 [-0.05072, -0.04547]

bayes-ramp-minus-lite-logistic: +0.00952 [0.00398, 0.01508]

lite-logistic-minus-alarm-time: +0.02748 [0.01956, 0.03573]

**L = 40960**: eligible 2000, late alarm 0, valid fraction 0.8620

| forecaster | Brier | log score | skill over alarm-time [95%] | Spiegelhalter Z | weighted gap | calibrated |
|---|---|---|---|---|---|---|
| base-rate | 0.5749 | 1.4749 | 0.052 [0.033, 0.071] | 77.49 | 0.675 | no |
| alarm-time | 0.6064 | 2.0630 | 0.000 [0.000, 0.000] | 120.11 | 0.710 | no |
| lite-logistic | 0.6090 | 2.2724 | -0.004 [-0.011, 0.003] | 131.05 | 0.700 | no |
| persistence | 0.1143 | 0.3839 | 0.812 [0.789, 0.833] | 4.12 | 0.032 | no |
| bayes-ramp-no-widening | 0.6379 | 3.4428 | -0.052 [-0.064, -0.041] | 163.85 | 0.705 | no |
| bayes-ramp | 0.6143 | 3.1798 | -0.013 [-0.025, -0.001] | 149.81 | 0.669 | no |
| oracle | 0.1080 | 0.3654 | 0.822 [0.803, 0.840] | -0.98 | 0.016 | yes |

bayes-ramp-minus-bayes-ramp-no-widening: -0.02356 [-0.02492, -0.02225]

bayes-ramp-minus-lite-logistic: +0.00532 [-0.00073, 0.01128]

lite-logistic-minus-alarm-time: +0.00256 [-0.00185, 0.00682]

**Failure time**: 2000 runs with an alarm, 1235 censored at the horizon, 239 already failed at issuance; Kaplan-Meier median failure sweep None, median remaining sweeps None.

| forecaster | RPS | 80% interval coverage [CP95] | upper end censored | median abs. error (sweeps) | both inside | forecast inside, outcome censored | forecast censored, outcome inside | both censored |
|---|---|---|---|---|---|---|---|---|
| empirical | 0.2959 | 0.512 [0.490, 0.535] | 841 | 11520.0 | 765 | 1235 | 0 | 0 |
| bayes-ramp | 0.3098 | 0.446 [0.424, 0.468] | 746 | 11520.0 | 701 | 1095 | 64 | 140 |

RPS skill of bayes-ramp over empirical: -0.047 [-0.063, -0.031]

## Identical histories, two continuations (E_continue, E_freeze)

Alarms identical: True.

| L | forecasts identical | mean p (continue) | mean p (freeze) | bound on worst-world excess | bound holds for every forecaster |
|---|---|---|---|---|---|
| 10240 | True | 0.8527 | 0.9415 | 0.00790 | True |
| 20480 | True | 0.5116 | 0.9152 | 0.07006 | True |
| 40960 | True | 0.1900 | 0.8513 | 0.12682 | True |

Excess Brier score over the path oracle (mean of (q - p)^2), per world:

| L | forecaster | continue | freeze | worse of the two |
|---|---|---|---|---|
| 10240 | alarm-time | 0.01144 | 0.01809 | 0.01809 |
| 10240 | base-rate | 0.02744 | 0.00748 | 0.02744 |
| 10240 | bayes-ramp | 0.00342 | 0.02785 | 0.02785 |
| 10240 | bayes-ramp-no-widening | 0.00486 | 0.04194 | 0.04194 |
| 10240 | lite-logistic | 0.00765 | 0.02090 | 0.02090 |
| 10240 | persistence | 0.03495 | 0.00011 | 0.03495 |
| 20480 | alarm-time | 0.06693 | 0.22476 | 0.22476 |
| 20480 | base-rate | 0.13905 | 0.17261 | 0.17261 |
| 20480 | bayes-ramp | 0.02450 | 0.25752 | 0.25752 |
| 20480 | bayes-ramp-no-widening | 0.02761 | 0.30497 | 0.30497 |
| 20480 | lite-logistic | 0.03983 | 0.24843 | 0.24843 |
| 20480 | persistence | 0.30924 | 0.00142 | 0.30924 |
| 40960 | alarm-time | 0.03418 | 0.48157 | 0.48157 |
| 40960 | base-rate | 0.08675 | 0.45520 | 0.45520 |
| 40960 | bayes-ramp | 0.01995 | 0.48599 | 0.48599 |
| 40960 | bayes-ramp-no-widening | 0.02051 | 0.50965 | 0.50965 |
| 40960 | lite-logistic | 0.03232 | 0.48343 | 0.48343 |
| 40960 | persistence | 0.57667 | 0.00968 | 0.57667 |
