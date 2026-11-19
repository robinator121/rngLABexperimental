# Layer A and B2 benchmarks: tables from `results/benchmark.json`

Code `d91cb45a28c638f526655677c5bb8dc8466d60af` (dirty: True); started 2026-11-19T17:48:50.131180+00:00; 581 s in total.

Machine: AMD Ryzen 7 7435HS, 16 logical CPUs, 15.3 GiB RAM; Linux-7.2.14-arch1-1-x86_64-with-glibc2.44; Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1. One workstation; none of this is an HPC measurement.

Every configuration ran in fresh processes: times are medians of three untraced processes; "Python peak" is tracemalloc's peak from a fourth process; "peak RSS" includes the interpreter and the imported libraries (baseline: 64 MiB, import only). "samples" gives the smallest and largest of the three wall times. Fields that differ by seed (Layer A's events and output size) are medians of the three processes.

## 1. Measured: Layer A (S2, horizon 12, events observed on the six-interval grid)

| nuclei | wall s | samples | CPU s | events | events/s | words/s | output MiB | Python peak MiB | peak RSS MiB |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1e+03 | 0.0044 | 0.002106-0.007988 | 0.0038 | 3087 | 7.05e+05 | 2.74e+06 | 0.1 | 0.5 | 66 |
| 1e+04 | 0.0091 | 0.008956-0.01037 | 0.0090 | 30927 | 3.4e+06 | 1.32e+07 | 1.3 | 4.8 | 67 |
| 1e+05 | 0.0927 | 0.08393-0.09446 | 0.0916 | 309195 | 3.34e+06 | 1.3e+07 | 12.7 | 40.3 | 104 |
| 1e+06 | 1.0582 | 0.8875-1.243 | 1.0459 | 3093913 | 2.92e+06 | 1.13e+07 | 126.9 | 365.6 | 434 |
| 2e+06 | 2.4710 | 2.242-2.601 | 2.4392 | 6187478 | 2.5e+06 | 9.71e+06 | 253.7 | 731.5 | 797 |

## 2. Measured: B2 against the population (S2, six intervals, plan construction included)

| nuclei | parts | words/interval | wall s | samples | CPU s | plan s | s per interval | Python peak MiB | peak RSS MiB |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1000 | 1 | 25 | 0.0410 | 0.03454-0.07024 | 0.3380 | 0.0286 | 0.00208 | 0.17 | 84 |
| 1e+05 | 1 | 25 | 0.0499 | 0.04883-0.05367 | 0.4136 | 0.0277 | 0.00370 | 0.18 | 84 |
| 1e+06 | 1 | 25 | 0.0512 | 0.04609-0.05395 | 0.4230 | 0.0256 | 0.00428 | 0.18 | 85 |
| 1e+09 | 1 | 25 | 0.0667 | 0.06478-0.07146 | 0.5647 | 0.0276 | 0.00651 | 0.15 | 84 |
| 1e+12 | 1 | 25 | 0.1058 | 0.09826-0.1213 | 0.9060 | 0.0278 | 0.01301 | 0.15 | 85 |
| 1.126e+15 | 1 | 25 | 0.3207 | 0.3027-0.4726 | 1.1522 | 0.0245 | 0.04936 | 0.15 | 85 |
| 9.007e+15 | 8 | 200 | 0.3932 | 0.2749-0.4181 | 1.2276 | 0.0249 | 0.06138 | 0.15 | 85 |
| 7.206e+16 | 64 | 1600 | 1.2266 | 1.223-1.366 | 2.0847 | 0.0301 | 0.19942 | 0.19 | 85 |
| 1.153e+18 | 1024 | 25600 | 22.8234 | 20.77-22.83 | 23.3784 | 0.0285 | 3.79915 | 0.92 | 84 |

## 3. Measured: B2 against the number of intervals (S2, 1e6 nuclei, uniform grid on [0, 12])

| intervals | outputs kept | wall s | samples | CPU s | plan s | ms per interval | output MiB | Python peak MiB | peak RSS MiB |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | True | 0.067 | 0.06487-0.07219 | 0.560 | 0.033 | 3.389 | 0.01 | 0.21 | 85 |
| 10 | False | 0.072 | 0.06786-0.08711 | 0.603 | 0.033 | 3.859 | 0.00 | 0.20 | 85 |
| 100 | True | 0.283 | 0.2791-0.2871 | 1.196 | 0.034 | 2.486 | 0.06 | 0.44 | 85 |
| 100 | False | 0.310 | 0.268-0.3312 | 1.239 | 0.039 | 2.711 | 0.00 | 0.24 | 85 |
| 1000 | True | 2.087 | 2.031-2.177 | 3.212 | 0.062 | 2.025 | 0.57 | 2.67 | 85 |
| 1000 | False | 1.940 | 1.928-2.413 | 3.069 | 0.060 | 1.880 | 0.00 | 0.37 | 84 |
| 10000 | True | 13.991 | 12.98-14.47 | 15.066 | 0.082 | 1.391 | 5.72 | 24.50 | 106 |
| 10000 | False | 13.231 | 12.85-13.23 | 14.270 | 0.071 | 1.316 | 0.00 | 1.58 | 77 |

## 4. Measured: checkpoints (B2, 1e6 nuclei, 2,000 uniform intervals, outputs not kept)

| checkpoint every | wall s | samples | CPU s | ms per interval | Python peak MiB | peak RSS MiB |
|---:|---:|---:|---:|---:|---:|---:|
| never | 3.485 | 3.481-3.619 | 4.623 | 1.707 | 0.47 | 76 |
| 100 | 3.195 | 3.107-3.247 | 4.283 | 1.569 | 0.48 | 68 |
| 1 | 3.225 | 3.155-3.243 | 4.369 | 1.580 | 0.48 | 70 |

## 5. Measured: plan construction against distinct interval lengths (B2, 1e6 nuclei, 200 intervals)

One set of interval laws serves every interval of a uniform grid; a grid whose 200 lengths all differ needs 200 sets.

| grid | plan s | plan samples | wall s | Python peak MiB | peak RSS MiB |
|---|---:|---:|---:|---:|---:|
| uniform | 0.0328 | 0.03133-0.03978 | 0.4279 | 0.29 | 68 |
| all lengths distinct | 0.9645 | 0.9146-1.13 | 1.3625 | 3.72 | 70 |

## 6. Measured: thread settings (B2, 1e6 nuclei, 1,000 intervals)

| setting | wall s | CPU s | CPU/wall |
|---|---:|---:|---:|
| default (unset) | 1.671 | 2.779 | 1.66 |
| 1 thread | 2.085 | 2.063 | 0.99 |

## 7. Measured: independent processes started together (single-thread settings)

Wall time runs from the start of the first process to the last result, so it includes interpreter start-up and imports; the median task time is each worker's own timing, which excludes them.

| workload | processes | wall s | runs/s | median task wall s |
|---|---:|---:|---:|---:|
| layer_a_n_1e5 | 1 | 1.404 | 0.71 | 0.105 |
| layer_a_n_1e5 | 2 | 1.487 | 1.35 | 0.154 |
| layer_a_n_1e5 | 4 | 1.868 | 2.14 | 0.201 |
| layer_a_n_1e5 | 8 | 2.720 | 2.94 | 0.325 |
| b2_n_1e6_1000_intervals | 1 | 3.098 | 0.32 | 1.845 |
| b2_n_1e6_1000_intervals | 2 | 3.123 | 0.64 | 1.778 |
| b2_n_1e6_1000_intervals | 4 | 3.385 | 1.18 | 1.911 |
| b2_n_1e6_1000_intervals | 8 | 3.803 | 2.10 | 2.045 |

## 8. Projections (NOT measured)

PROJECTION, not measured: linear extrapolation of the medians above, this machine, one process.

- Layer A: 1.058e-06 s per nucleus, hence about 1058 s for 1e9 nuclei (memory would not fit in one piece; Layer A is chunked over nuclei but keeps every event).
- B2: 2.519e-03 s per interval at 1e6 nuclei, hence about 2519 s for 1e6 intervals.
