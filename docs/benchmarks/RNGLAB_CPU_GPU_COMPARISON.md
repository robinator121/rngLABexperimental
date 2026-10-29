# RNGLab CPU/GPU comparison — actual laptop measurements

The same harmonic FP64 supplied-input task was executed on the Ryzen 7 7435 HS
compiled CPU backend and RTX 4070 Laptop GPU. Each history has 73728 attempted
updates,8192 burn-in and 65536 retained rows. Fixed namespace 93901/rho +.025,
identical clean initialization and supplied uniforms remove RNG confounding.
This namespace is a benchmark fixture, not independent scientific replication.

Three fixed timed repetitions followed full-shape warmup, alternating backend
order. Times include host preparation, transfers, synchronized physics kernel
and host statistical diagnostics; RNG generation is recorded separately.

|Batch|CPU pipeline median, s|GPU pipeline median, s|GPU/CPU time|CPU updates/s|GPU updates/s|
|---:|---:|---:|---:|---:|---:|
|1|0.024366|0.182409|7.49×|3,025,884|404,192|
|2|0.047873|0.228347|4.77×|3,080,140|645,753|
|4|0.068831|0.272870|3.96×|4,284,573|1,080,780|

**No GPU speedup was observed at these tiny batches.** More chains increased
GPU aggregate throughput, but the compiled CPU remained faster. Only 1–4 serial
chains were active, deliberately leaving ample exporter headroom. This is not
an optimized GPU throughput comparison, a desktop-versus-laptop comparison or
a claim that the GPU will always be slower. Low occupancy, serial chain
latency, FP64 arithmetic and global-state accesses are candidate explanations;
the experiment does not separately identify their contributions.

GPU median timing breakdown:

|Batch|Kernel, s|H2D, s|D2H, s|Host preparation, s|Host statistics, s|
|---:|---:|---:|---:|---:|---:|
|1|0.166303|0.001379|0.000787|0.000540|0.013856|
|2|0.197463|0.001966|0.004640|0.002427|0.021143|
|4|0.215519|0.002515|0.007638|0.003279|0.043268|

The first GPU warmup included 0.586852 seconds in its synchronized kernel/JIT
interval versus the subsequent batch 1 median 0.166303 seconds. The difference
is not a pure JIT estimate because clock/startup effects also occur. All raw
warmup and repetition timings are preserved. Transfers were much smaller than
kernel time in these measurements; eliminating transfers alone would not make
the tiny-batch GPU faster than the compiled CPU.

Telemetry showed total desktop+test VRAM approximately 1.4–1.5 GiB and benchmark
process RSS below 470 MiB. GPU sensor utilization samples can lag and indicate
active time rather than useful SM occupancy; they do not contradict the tiny
grid warning. One-shot CPU percentage samples are uninformative for subsecond
runs. No calibrated power integration was performed, so throughput per watt
and energy efficiency remain unmeasured.

The included time-to-SE extrapolations are unstable pilot-variance projections
from at most 4 independent chains, not reliable measured time-to-precision.
They must not be promoted to a hardware ranking. A useful later benchmark,
after preservation is secure, would prospectively test larger well-populated
batches with matched CPU threading and enough independent chains to estimate
variance. No such scaling was performed here.

Cross-platform scientific parity succeeded:33 actual CUDA tests and the full
clock fixture passed unchanged FP64 tolerances; all exercised acceptance/event
decisions agreed. The fresh 32-identity pilot reproduced expected centroid-shift
directions, with wide probability uncertainty. Neither parity nor these tiny
timings establish a new physical discovery or independent publication-grade
confirmation. No fair VM speedup is reported because the VM campaign includes
different batching, condition reuse, detector policies and parallelism.
