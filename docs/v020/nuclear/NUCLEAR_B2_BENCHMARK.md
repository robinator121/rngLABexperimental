# RNGLab v0.20 nuclear module, milestone B2: single-machine benchmarks

**Status: 2026-11-19. Measured on one workstation; none of this is an HPC measurement.** Produced by `research/v020-nuclear-b2-bench/benchmark.py`, run once at commit `d91cb45` (`results/benchmark.json`, started 2026-11-19T17:48:50Z, 581 s in total). The file records `git_dirty: true`: the run started two seconds after the numerics script had written its output, which was not yet committed. No tracked file differed from `d91cb45`; the twelve module hashes recorded in `benchmark.json` equal the files of that commit. Every number below comes from that file or from `results/tables.md`, which `report_tables.py` writes from it. Projections are labelled as such and kept apart from measurements.

## 1. How it was measured

Every configuration ran in fresh Python processes, so peak resident memory, CPU time and wall time belong to that configuration alone. Times are medians of three untraced processes, whose samples are kept; Python-level allocation peaks come from a fourth process under `tracemalloc`, whose times are discarded. Repeats must reproduce every field that the configuration fixes (words, bytes, parts, words per interval; B2's output size; Layer A's conserved population), or the benchmark stops. Random numbers come from the engineering namespace `rnglab-v020-nuclear-bench`. Scheme S2 of the Layer A protocol, started from Xa-400, horizon 12.

**Machine.** AMD Ryzen 7 7435HS, 16 logical CPUs, 15.3 GiB RAM, Linux 7.2.14 (Arch), Python 3.14.9, NumPy 2.5.3, SciPy 1.18.1. The interpreter with the imported libraries occupies 64 MiB (import-only baseline). **The workstation was in ordinary use during the run:** other applications were active and took several cores at times. `benchmark.json` does not record the system load, so the timings include an unknown amount of that interference; the spread of the three samples per configuration shows its size.

## 2. Layer A, for reference (S2 from Xa-400, horizon 12)

| Nuclei | Wall s (median) | Events | Events per s | Output MiB | Python peak MiB | Peak RSS MiB |
|---:|---:|---:|---:|---:|---:|---:|
| 10⁵ | 0.093 | 309,195 | 3.3 × 10⁶ | 12.7 | 40.3 | 104 |
| 10⁶ | 1.06 | 3,093,913 | 2.9 × 10⁶ | 126.9 | 365.6 | 434 |
| 2 × 10⁶ | 2.47 | 6,187,478 | 2.5 × 10⁶ | 253.7 | 731.5 | 797 |

Layer A keeps every event, so its memory grows with the population: about 366 MiB of Python allocations per million nuclei over this horizon.

## 3. B2 against the population (S2, six intervals, plan construction included)

| Nuclei | Parts | Words per interval | Wall s (median) | Plan s | Python peak MiB |
|---:|---:|---:|---:|---:|---:|
| 10³ | 1 | 25 | 0.041 | 0.029 | 0.17 |
| 10⁶ | 1 | 25 | 0.051 | 0.026 | 0.18 |
| 10¹² | 1 | 25 | 0.106 | 0.028 | 0.15 |
| 2⁵⁰ | 1 | 25 | 0.321 | 0.025 | 0.15 |
| 2⁵³ | 8 | 200 | 0.393 | 0.025 | 0.15 |
| 2⁵⁶ | 64 | 1,600 | 1.23 | 0.030 | 0.19 |
| 2⁶⁰ | 1,024 | 25,600 | 22.8 | 0.029 | 0.92 |

Memory does not grow with the population (peak RSS 84–85 MiB throughout). Time grows slowly within one part (the binomial distribution function costs more at large n) and with the number of parts above 2⁵⁰, reaching 3.8 s per interval at 2⁶⁰ (1,024 parts).

## 4. B2 against the number of intervals (10⁶ nuclei, uniform grid on [0, 12])

| Intervals | Wall s, outputs kept | Wall s, not kept | ms per interval (not kept) | Python peak MiB, kept / not kept |
|---:|---:|---:|---:|---:|
| 10 | 0.067 | 0.072 | 3.9 | 0.21 / 0.20 |
| 100 | 0.283 | 0.310 | 2.7 | 0.44 / 0.24 |
| 1,000 | 2.09 | 1.94 | 1.9 | 2.67 / 0.37 |
| 10,000 | 14.0 | 13.2 | 1.3 | 24.5 / 1.58 |

The time per interval falls as fixed costs are spread over more intervals, to 1.3 ms at 10,000 intervals. Kept outputs take 5.7 MiB at 10,000 intervals; without them the Python peak stays at 1.6 MiB.

## 5. Checkpoints, plan construction and threads

- **Checkpoints** (2,000 intervals, outputs not kept). Without checkpoints the run took 3.49 s, with one every 100 intervals 3.20 s and with one every interval 3.23 s. Building a checkpoint costs nothing measurable here. The run without checkpoints was the slowest of the three, which shows the size of the background variation. The callback discards each checkpoint, so writing checkpoints to disk was not measured.
- **Plan construction** (200 intervals). A uniform grid needs one set of interval laws: 0.033 s. A grid whose 200 lengths all differ needs 200 sets: 0.96 s, about 5 ms per distinct length.
- **Threads** (1,000 intervals). With the default BLAS/OpenMP settings the run took 1.67 s of wall time and 2.78 s of CPU; with one thread, 2.09 s and 2.06 s. More threads save about a fifth of the wall time at about a third more CPU time, so one thread per process is the efficient setting for many independent tasks.

## 6. Independent processes started together (one thread each)

| Workload | Processes | Runs per second | Median task time, s |
|---|---:|---:|---:|
| Layer A, 10⁵ nuclei | 1 / 2 / 4 / 8 | 0.71 / 1.35 / 2.14 / 2.94 | 0.105 / 0.154 / 0.201 / 0.325 |
| B2, 10⁶ nuclei, 1,000 intervals | 1 / 2 / 4 / 8 | 0.32 / 0.64 / 1.18 / 2.10 | 1.85 / 1.78 / 1.91 / 2.05 |

Throughput includes interpreter start-up and imports, about 1.3 s per process here. B2 tasks slow by about a tenth from one to eight processes; Layer A tasks, which allocate far more memory (40 MiB against under 1 MiB per task), by a factor of about three.

## 7. Projections (not measured)

Linear extrapolation of the medians above, for this machine and one process: Layer A at 1.06 × 10⁻⁶ s per nucleus, about 1,058 s for 10⁹ nuclei (whose events would not fit in this machine's memory in one piece); B2 at 2.5 × 10⁻³ s per interval at 10⁶ nuclei, about 2,519 s for 10⁶ intervals. The B2 figure is the median over 10 to 10,000 intervals, so it overstates the cost of long runs, which measured 1.3 ms per interval.

## 8. What this does not show

- One workstation, one operating system, one Python build. Cluster nodes differ in CPU, memory bandwidth, file system and Python (Komondor has no Python 3.14; a container is needed).
- Wall times of short runs (milliseconds) are dominated by fixed costs and vary between repeats; the samples are in `results/benchmark.json`.
- The projections assume linear scaling of medians measured here. They are not measurements, and nothing here estimates HPC CPU-hours or storage; that needs profiling on the target nodes (Phase G).
