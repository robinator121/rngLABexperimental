# RNGLab v0.15 measured-source pipeline

> [!WARNING]
> **Production Cryptography Disclaimer:** The conditioning mechanisms (SHA-256 block hashing, ChaCha20 extraction, von Neumann debiasing) implemented in RNGLab are strictly for academic demonstration, stochastic simulation benchmarking, and research treatments. They are **not certified or audited for production cryptographic systems or key generation**.
> 
> Passing statistical test batteries (such as NIST SP 800-22) on cryptographically conditioned output demonstrates only the pseudorandom diffusion of the hash function, **not** the physical entropy of the raw entropy source. True physical entropy must be evaluated directly on raw, unconditioned noise using NIST SP 800-90B min-entropy estimators and continuous health testing.

This release is an additive measured-data research lane. It preserves rnglab through
rnglab14 and the existing rnglab12 command. Historical numerical and standards claims
are not inherited as validation of the new measurements.

## Recovered design and evidence

The physical source pipeline incorporates experimental hardware measurements
including Arduino sensor acquisition and muon detector waveform streams.
The architecture establishes raw data ingest, cryptographic provenance tracking,
separate source conditioning stages, finite supplied-input Monte Carlo testing,
and NIST SP 800-90B entropy evaluation.

The implemented reconstruction is:

1. Explicit non-sensitive physical-source declaration and an explicit file schema.
2. Hash original bytes, make a content-addressed read-only snapshot, verify both.
3. Parse in acquisition order, preserving resets, missing values and malformed-data
   failures. Unknown units stay unknown; no silent unwrapping, sorting or repair.
4. Characterize waveform shape/rails/trigger alignment/recovery/baseline correlation,
   timestamp lattice/order/gaps, and per-sensor missingness/range.
5. Define each experimental symbol stream separately. Current detector experiments
   use one sample-0 ADC LSB per event, peak LSB, timestamp LSB, positive-gap LSB,
   and the positive-gap LSB after an explicit exploratory integer divisor.
   Dividing by an observed lattice is not evidence of entropy. Sensor experiments
   use one register LSB per available reading; NA omission can bridge time gaps.
6. Run continuous binary RCT/APT diagnostics on the complete symbol stream. H=0.5
   is an explicit exploratory assumption, not an inferred or validated lower bound.
   State survives input chunks, APT uses exact binomial tails, and failures latch.
7. Apply an IID-only simultaneous Hoeffding marginal bound and held-out context
   predictors. The predictor tables are frozen on the first half; the second half
   is scored sequentially. These screens do not validate arbitrary non-IID entropy.
8. Export one binary symbol per byte for the official external IID/non-IID suite.
   The external adapter records executable/input hashes, invocation, exit code,
   stdout/stderr and timeouts. Fewer than one million symbols or no executable is
   a concrete blocker. Restart assessment remains separate and is required.
9. Produce finite raw-packed, Von Neumann, and 512-to-256-bit SHA-256 research
   treatments. Carry all tail/drop/output counts. SHA uses domain separation; it
   is a conditioning experiment with no unconditional extraction guarantee. Health
   alarms remain visible even if hashing produces visually uniform output.
10. Compare raw, deliberately degraded, and conditioned research streams using
    finite file segments in pi, VMC, Ising and PIMC. No padding, cyclic replay or
    DRBG expansion disguises insufficient physical input. Keep input offsets,
    exact consumed-byte hashes, separate proposal/acceptance/initialization roles,
    and paired auxiliary seeds. VMC includes x² because optimum local energy is
    identically 0.5 even when the sampled distribution is wrong.
11. Preserve source alarms separately from application tolerance failures. The
    inherited forecasting/RUS lane remains available; applying it to physical
    data needs calibrated trajectories and independent future runs. No RUS claim
    is made from this single acquisition.

## Entropy and extraction boundary

Every measured stream currently receives **zero validated entropy credit**. A
balanced marginal or a hashing result is not a conditional min-entropy bound.
Official non-IID/restart assessment, a noise-source model, operating-condition
coverage, and independent review are missing. No secure output bit rate is assigned.
The Toeplitz reference implements universal hashing with the explicit sufficient
budget m <= floor(k - 2 log2(1/epsilon)). It refuses overdrawn budgets. The caller
must establish conditional block min-entropy and an independent uniform seed; this
function cannot validate either premise and is not used to grant entropy here.
A public fixed file is replayable and cannot supply fresh private entropy.

A production conditioned source or reseeded DRBG must remain blocked until this
boundary is satisfied. Its control plane must latch health failures before emission,
retain/quarantine affected blocks, refuse entropy-budget overdraw, document seed
independence/reseed limits, and demonstrate validated security bounds.

## Formats and streaming limits

The inventory lists every regular file, byte size, permissions, owner IDs, hash and
CSV header. It never treats a filename as proof of origin. Waveforms require a
configured length and two provisional metadata fields. Every integer/range/row
length is checked; errors stop the run with a failed registry. CSV rejects duplicate
headers and malformed row widths. Literal NA is retained. Event CSV requires explicit
channel and integer timestamp columns; count-window data is not event timing.
The explicit rng_bytes adapter preserves a regular binary file with caller-supplied
provenance; it does not identify ESP32 output from its byte distribution.

Raw copy/hash I/O uses 1 MiB chunks; line input is capped at 128 KiB. Features retain
one waveform at a time and at most 16 channels. Per-sample moments use arrays sized
by the configured waveform length. A bounded 2M-symbol prefix is used for the
predictability screen; full health diagnostics and treatment output are streamed.
Trend summaries cap at 4096 blocks of 512 events; longer acquisitions must split
into registered runs. A finite run refuses existing output directories. Partial
failures retain their evidence and failed registry, and are never silently resumed.
Raw capture leaves .partial files after transport failure; do not append to them.
Close and register a new segment with an explicit continuity relationship.

Serial capture is optional (pyserial), records every received byte including startup
messages, and has duration/byte/time-out bounds. It does not promise lossless
transport or infer a device clock from host reception. Opening UNO ports may reset
hardware. The previous standalone acquisition package at commit 55b2714 includes
framed MCU/serial and firmware support; it stays separate under its existing AGPL
license. This release does not copy that package or change ownership/licensing.
No live hardware was opened or flashed for this task.

## Experimental inference

Synthetic sensitivity protocol is saved before outcomes. Two independent OS-byte
cohorts use distinct role/identity seed namespaces and archived raw reference bytes.
The four primary tolerance rules and model horizons are frozen. Exact paired
sign/McNemar tests compare discordant failures to OS; Holm controls the complete
28-comparison confirmatory family. Outcomes are short-horizon sensitivity endpoints,
not equilibrium-convergence evidence. Inherited within-chain standard errors are
not used to claim significance. Ising uses exact L=4 enumeration; PIMC uses an exact
finite-bead harmonic covariance reference to avoid calling discretization bias an
RNG defect. Application runs use the same 8-bit midpoint uniform grid for all arms.

Physical records and all their treatments are already observed: physical results
are **exploratory only**. Adjacent file segments from one acquisition are not
independent physical replications. Unequal finite output lengths limit the paired
trial set; report counts and abstain when insufficient. Future confirmatory
physical analysis must use a new acquisition, frozen source mapping, hardware and
run-level held-out splits, no training/selection leakage, power planning and
correction across every prespecified channel/observable/treatment.

Authoritative references (checked 2026-10-19):
- [NIST SP 800-90B and listed errata](https://csrc.nist.gov/pubs/sp/800/90/b/final)
- [Official assessment software and byte-symbol interface](https://github.com/usnistgov/SP800-90B_EntropyAssessment/blob/master/README.md)
- [ESP32 hardware RNG prerequisites](https://docs.espressif.com/projects/esp-idf/en/latest/esp32/api-reference/system/random.html)
