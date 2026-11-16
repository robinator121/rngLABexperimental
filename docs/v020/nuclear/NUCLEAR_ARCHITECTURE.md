# RNGLab v0.20 nuclear module: architecture and scientific assessment

**Status: draft for the author's review, 2026-11-16.** Branch `feature/criticality-rng-sensitivity`, created from `feature/v020-stage-d-readiness` at `6b75205`; not pushed. Prepared for the author. Nothing here changes a frozen protocol, the Stage D pilot, its files or its namespaces.

## 1. Questions, scope and relation to the central hypothesis

- **Primary question (author's brief).** How do degraded random-number sources influence statistical reliability, uncertainty calibration and observable output distributions in stochastic nuclear-physics simulations?
- **Secondary question.** Are RNG defects more detectable near a critical transition than in noncritical regimes? This is treated as a hypothesis to falsify; no direction is assumed.
- **Scope.** Computational reliability only. Every model is synthetic and non-operational: decay schemes with Z ≥ 140 (outside the chart of known nuclides), abstract branching processes indexed by a mean offspring number, and reduced population kinetics with synthetic precursor groups. There are no geometries, inventories, control mechanisms, facility parameters or operating procedures. Nothing in the module is reactor-physics, neutron-transport or reactor-design software, and passing tests do not make it validated physics software.
- **Relation to H1–H3.** The module does not redefine them. It adds application workloads whose healthy output laws are known exactly, which is what application-level detection (H1) needs: a distortion can then be measured without the sampler-failure confound that Phase II found in healthy PIMC. Whether it counts as roadmap capability C10 (a second application) or as a separate supporting study is the author's decision (Section 9). H2 and H3 are not addressed.

## 2. What the repository provides (inspection of 2026-11-16)

| Item | Where | Role in this module |
|---|---|---|
| `RandomSource.read_bytes`: one audited boundary that counts bytes and SHA-256-hashes them | `src/rnglab14/sources/base.py`, re-exported by `rnglab.sources` | every random number passes through it |
| Healthy generators (PCG64, xoshiro256++, MT19937, ChaCha, OS entropy) | `rnglab14.sources.reference` | xoshiro256++ is the second healthy control. PCG64, xoshiro256++ and ChaCha are pure Python (microseconds per word, too slow for campaigns); MT19937 uses Python's C implementation; OS entropy cannot be replayed |
| `namespace_seed`; `SeedSequence` → PCG64 convention | `rnglab.seeds` | seed label of every stream |
| Bit-level fault adapter (bias, Markov ρ on bits, stuck bit, periodic coupling) | `rnglab14.sources.faults.DegradedFaultSource` | existing defect classes for Phase E. Its ρ path loops per bit in Python |
| Word-level faults SHIFT, LAG1, DUP with a time profile | `research/v019-d-pilot/dlib.py` (`FaultyWords`) | definitions to port in Phase E. Research code: it accepts only single 8-byte reads and is indexed by sweep |
| Holm, Bonferroni, Benjamini–Hochberg | `rnglab.stats` | test families |
| Test suite on this branch before any change: **718 passed, 2 skipped** in 193.7 s | `pytest tests` | regression baseline (same counts as `docs/v020/evidence/pytest_v020_confirmed.txt`) |

Constraints carried over: freeze before generate; fresh seed namespaces; every number from a script; negative results with equal weight; no test samples in a reserved namespace; CPU only. Python 3.14 runs locally, but Komondor has none, so a container is needed there.

## 3. Architecture

```text
RandomSource  (healthy generator | fault adapter wrapping one | later: recorded physical stream)
  └─ WordStream: 8-byte big-endian words → uniforms (word >> 11)·2⁻⁵³; counts words, never buffers
       ├─ Layer A  decay           event engine (vectorised) + pure-Python reference engine   [implemented]
       │                           aggregated multinomial mode                                [planned, B2]
       ├─ Layer B  branching and reduced kinetics                                             [planned]
       └─ Layer C  detector response: physical emissions → detector records                   [planned]
analytic references (draw nothing): bateman.py; planned: PGF fixed points, moment equations
checks.py (exact-law tests, invariants) → frozen studies research/v020-nuclear-*/ → manifests, results
```

Design rules:

1. **One random-number boundary.** All randomness goes through `RandomSource.read_bytes`. Roles (physics, detector, analysis) get separate seed labels. A fault adapter is a `RandomSource` that wraps another one, so it acts on exactly the words the physics consumes. A uniform here equals `RandomSource.random_float` on the same bytes.
2. **Inversion sampling with a fixed word layout.** In Layer A each decay step owns three words (waiting time, channel, photon), and each nucleus owns K = 3·D words, D being the longest decay path. Consumption therefore does not depend on the outcome, chunked and unchunked runs agree bit for bit, a run can be replayed exactly, and every stream position maps to one (nucleus, step, purpose), so a fault's footprint can be traced. The cost is that unused slots are read and discarded; reports give words read and events separately.
3. **Physical events are separate from measurement.** Layers A and B emit physical event records with identifiers and lineage. Layer C only transforms them. The generator that drives a simulation is never treated as the physical randomness a detector measures. A future recorded stream enters as a `RandomSource` with its own provenance and no assumed uniformity or independence.
4. **Two engines where that is cheap.** A vectorised engine for production and a pure-Python reference engine read the same stream. Their agreement is a software check, not physics validation.
5. **Validation against analytic laws.** For independent nuclei the count laws are exact (binomial, multinomial), and the state probabilities come from two independent analytic routes: the matrix exponential of the decay generator and the closed-form Bateman solution summed over decay paths.
6. **Manifests.** Each run records the code (Git commit, dirty flag, SHA-256 of every module file), the environment, the stream (source identity, seed label, first word, words and bytes read, digest of the consumed bytes), the configuration and the output digests (`rnglab/nuclear/manifest.py`).

## 4. Layers

### Layer A: radioactive decay (implemented in this milestone)

- **Data model** (`decay.py`). Synthetic nuclides with Z, A, isomer level, half-life (abstract time units) and channels. Modes: α, β⁻, β⁺, electron capture, isomeric transition, β-delayed neutron emission. The scheme validator rejects unnormalised branching (tolerance 10⁻¹², never renormalised), a (ΔZ, ΔA) that contradicts the mode, isomeric transitions that do not lower the level, duplicate (Z, A, level), unknown daughters and cycles.
- **Events.** Each event carries an identifier fixed by the layout (nucleus·D + step), the nucleus, the step, a timestamp, parent and daughter, channel, mode, a photon flag, a neutron flag and its predecessor (lineage). Emissions (α, e⁻, e⁺, γ, X-ray, n) are derived from events and sorted in time. Photons carry the decay's timestamp and are emitted with the channel's yield.
- **Neutrons.** Only the β-delayed-neutron channel emits one. `DecayScheme.precursors()` lists the delayed-neutron precursors with their neutron emission probability per decay. The aggregate precursor groups of Layer B are a separate abstraction and will be checked against this event-level one.
- **Not modelled.** Internal conversion and Auger electrons (a photon yield below 1 stands in for them), bremsstrahlung, neutrinos, annihilation photons, recoil, energies, geometry, environment dependence. Waiting times come from 53-bit uniforms, so they are capped at 53·ln 2/λ (about 36.7 mean lives); the sampled law differs from the exponential by at most 2⁻⁵³ in its distribution function.
- **Macroscopic populations (planned, milestone B2).** For independent nuclei the joint law of what happens in a time bin is a multinomial over decay paths, with exact path probabilities from the Bateman path terms. Sampling it by binomial inversion keeps one uniform per binomial and a fixed consumption, so the aggregated mode is exact up to the numerical accuracy of the binomial distribution function. It will be checked against the event engine and the analytic laws.

### Layer B: branching processes and reduced kinetics (planned)

- **B-1, individual histories.** Galton–Watson and continuous-time Markov branching processes, with offspring drawn by inversion. The regime is set by the mean offspring number m: subcritical m < 1, critical m = 1, supercritical m > 1. Exact or classical references: extinction probability as the smallest root of f(s) = s; E Z_n = mⁿ and Var Z_n = σ² m^(n−1)(mⁿ − 1)/(m − 1), equal to nσ² at m = 1; P(Z_n > 0) by iterating the generating function; at criticality, Kolmogorov's P(Z_n > 0) ~ 2/(σ² n) and Yaglom's exponential limit; closed-form extinction-time law for the linear birth–death process.
- **B-2, reduced kinetics.** A jump process on the neutron count and synthetic precursor-group counts, with prompt births and losses and precursor decay, simulated exactly (Gillespie). All rates are linear in the state, so the first and second moments obey closed linear equations: exact mean and covariance come from a matrix exponential, and the deterministic mean-field model is the reference. The regime descriptor is the synthetic system's multiplication factor or reactivity. No spatial transport and no real delayed-neutron data.
- The two families are built independently (individual histories against an aggregate jump process) and are cross-checked where they describe the same law.

### Layer C: synthetic detector response (planned)

Detection efficiency (thinning), background (independent Poisson superposition), timing jitter, finite timestamp resolution, non-paralysable and paralysable dead time (analytic count-rate relations m = n/(1 + nτ) and m = n·e^(−nτ)), event loss, saturation flags, count and interarrival distributions, and correlated sources (branching chains, with the variance-to-mean excess as an observable). It has its own RNG role. An ingestion interface for external recordings (chunked reading, offsets, checksums, provenance) will assume neither uniformity nor independence. No BME data is assumed. The muon-detector export already examined in `docs/physical/` is not a decay stream and carries no entropy claim.

## 5. RNG defects and controls (Phase E plan)

- **Classes, mapped to existing definitions.** Bias (SHIFT of uniforms; byte-shift bias), autocorrelation (LAG1 Gaussian copula; Markov bits), repeated subsequences (DUP; replay of a window), reduced effective entropy (words restricted to 2ᵇ values; stuck low bits), stuck bits. Severity and time profile (constant, step, ramp) are indexed by **word position**, so that "RNG consumption before detection" is a direct measurement. The definitions are ported into a package module and tested against the Stage D and `DegradedFaultSource` definitions on shared inputs.
- **Predicted effects.** With inversion and a fixed layout, a uniform shift δ moves each waiting-time distribution function and the first channel threshold by about δ. The first-order effect of each defect on each observable can be computed before generation (a fault-to-parameter map), so each campaign can state a prediction and test it.
- **Pairing.** Common-random-number pairs (the same base stream with and without the adapter) are used only where a design declares them. Otherwise seed labels are disjoint, and checks such as D5 of the Layer A study (distinct seeds and stream digests across replicates) guard against accidental sharing.

## 6. The criticality question: design considerations (Phase F)

- Near m = 1 both sensitivity and noise grow. Var Z_n grows linearly in n at m = 1 and geometrically above it. The extinction probability has infinite slope at m = 1⁺ (q ≈ 1 − 2(m − 1)/σ²), so a defect that shifts the effective m by Δm changes q by an amount that depends sharply on the regime. A defect that only changes correlations may instead be masked by the larger intrinsic variance. Either direction is possible.
- The **matching basis** must be declared. Near criticality histories are longer, so a comparison per history confounds RNG consumption. Proposal: primary basis = RNG bytes consumed (the unit of H1 and H2), secondary = histories.
- Endpoints to fix before generation: bias of the mean population and of the extinction probability by generation g; variance ratio; coverage of a declared interval estimator; distance between output distributions; detection latency in bytes for a sequential monitor at matched false-alarm rates. Each family has healthy negative controls, positive controls, Holm correction and a power analysis based on the fault-to-parameter map.
- Null hypothesis to try to reject, per defect class and severity: detection power at matched bytes in the critical regime is not greater than in both noncritical regimes. Regime × defect maps are reported whichever way they fall.

## 7. HPC and storage (Phase G plan; nothing has been measured on an HPC system)

- **Work unit.** One work unit is a (cell, replicate block), run as one task of a SLURM job array. Streams are regenerated from seed labels and never stored. Each task writes compact output packed into few files, because Komondor's Lustre quotas limit inodes (300k in scratch, 1M in project). A completion journal makes interrupted arrays restartable, aggregation is deterministic (sorted unit keys), and memory is bounded by the chunk size: Layer A holds chunk × K × 8 bytes of words plus the chunk's events.
- **Storage classes.** (1) Raw physical RNG recordings: immutable, checksummed, read in chunks with offsets. (2) Synthetic streams: never stored. (3) Intermediate simulation data: per task, deleted after aggregation unless a retention rule keeps it. (4) Final results: small JSON or NPZ files.
- **How the estimate will be made.** Peak storage and CPU-hours will be estimated from measured per-unit costs multiplied by the declared campaign matrix, after profiling. The local Layer A timings in the verification results are single-machine measurements, not HPC benchmarks.

## 8. Status at this milestone

| Component | Implemented | Tested | Independently verified | Experimentally validated |
|---|---|---|---|---|
| RNG bridge (`streams.py`) | yes | unit tests | conventions checked against `RandomSource.random_float` | n/a |
| Layer A event engines | yes | unit tests | second independent review done before the freeze; frozen study against analytic laws (`research/v020-nuclear-a-verify`) pending | no: synthetic schemes have no experimental counterpart |
| Layer A analytic references | yes | unit tests | two independent routes agree | n/a |
| Layer A aggregated mode | no (B2) | — | — | — |
| Layer B, Layer C | no | — | — | — |
| Fault injection, campaigns | no (existing adapters identified) | — | — | — |
| HPC execution and estimates | design only | — | — | — |
| Criticality question | hypothetical | — | — | not yet investigated |

## 9. Decisions for the author

1. Does the module count as C10 under H1, or as a separate supporting study?
2. A name and a separate pre-registration for the criticality question.
3. The matching basis for regime comparisons (proposal: RNG bytes).
4. Whether the Phase E adapters must reproduce the Stage D definitions bit for bit or become new, versioned definitions.
5. Merge order. This branch contains the unreviewed Stage D readiness commits it was created from; it should be merged after them, or rebased if they change.
6. The container approach for Komondor.

## 10. Restart instructions

In a checkout of branch `feature/criticality-rng-sensitivity`, with the environment of the README:

```bash
PYTHONPATH=src .venv/bin/python -m pytest tests/v020 -q     # Layer A unit tests
PYTHONPATH=src .venv/bin/python -m pytest tests -q          # whole suite
```

Next milestone: B2 (aggregated decay mode), then Phase C (Layer B), each with its own frozen verification study in a fresh namespace.
