# Phase II: harmonic PIMC uncertainty and random-channel faults

Completed 30 October 2026. Reference repository revision: `dd96adba23b8100aba9147d71fffc63932fc6b0d`.

## What the experiments established

The original healthy-source undercoverage was verified from 786,432 saved observations. Independent reference calculations found no finite-M target, energy-formula or interval-arithmetic explanation. Undercoverage persisted from exact stationary starts. Slow centroid dynamics dominated the variability of CV energy averages, while short-batch variance estimates missed much of that variability.

A new, locally frozen confirmation used M=16, beta=2, proposal half-width 0.4, 8,192 retained observations, exact Gaussian starts and CV BM512 intervals with 16 batches and a Student-t15 critical value:

| Update | Covered /1,024 | Coverage | Exact binomial 95% interval |
|---|---:|---:|---:|
| Local Metropolis | 811 | 79.20% | 76.58–81.65% |
| Local + centroid refresh | 964 | 94.14% | 92.52–95.50% |
| Independent Gaussian paths | 969 | 94.63% | 93.07–95.93% |

Both required healthy controls passed the frozen coverage, variance-ratio and bias tolerances. This is validation of the specified harmonic workload and interval, not a universal 95% coverage guarantee. A separate analytically solvable stationary dependence control reproduced short-batch variance underestimation with a correct target and healthy random inputs.

![Independent confirmation](experiment/figures/confirmation_primary_endpoints.png)

Only after that healthy gate passed, acceptance-input bias and proposal/acceptance coupling were injected into otherwise matched runs. Healthy coverage was **487/512**, versus **0/512** and **13/512**. These are controlled adapter effects, not discovered intrinsic PCG64 defects. A marginal source-mean check missed coupling, the joint source check detected it, and the application-shape monitor missed acceptance bias. No general application-monitor superiority was established.

The RUS counterexample has an identical observed history through the alarm but different future failure laws. It demonstrates the need for an explicit continuation model; no successful general forecasting claim is made.

## Read the evidence

- [PIMC root-cause investigation](reports/PIMC_root_cause_report.md)
- [Controlled random-input fault study](reports/Controlled_RNG_fault_report.md)
- [Research assessment and one next milestone](reports/Research_advancement_and_next_milestone.md)
- [Exact dependence control](analytic-dependence/analytic_dependence_report.md)
- [RUS counterexample and protocol](rus/RUS_identifiability_and_protocol.md)
- [Original baseline recount](baseline/baseline_verification_report.md)
- [Optional centroid update and its regression checks](environment/centroid_update_patch.md)
- [Public reproduction instructions](REPRODUCE.md)
- [Data availability and historical-record conventions](DATA_AVAILABILITY.md)
- [Public file manifest](MANIFEST.json) and [checksums](SHA256SUMS)

The repository includes the tested optional centroid update in `src/rnglab14/models`, disabled by default. The coverage campaign used the independently checked standalone implementation preserved here. Equivalent coverage validation through the package API remains the explicitly proposed integration milestone.

## What is distributed

Reports, unchanged frozen generation code, regression checks, protocols, seed manifests, analysis tables, replica endpoints, figures, dependency pins and recorded verification evidence are included. Seven lightweight Phase I artifacts are in [the historical pilot directory](../phase1-20261030/README.md).

The **19 large Phase II NPZ archives are not in Git**. Their expected sizes and hashes are recorded, but no hosted download is promised. The public replay preparer constructs a new workspace without private inputs or a bundled source tarball. A full replay generates approximately 2 GiB of traces. See [data availability](DATA_AVAILABILITY.md).

The unchanged repository passed 238 tests; the optional patch passed 257 core tests; the research checks passed 35 tests; the full local evidence audit passed 132 checks. Those are historical executed results with distinct scopes. Publication did not rerun the numerical campaigns or claim a new CI pass. The new public preparer was exercised without running simulations. The historical full-data audit requires omitted inputs and cannot complete from this lightweight checkout.

Generation scripts and frozen protocols retain their original bytes, including the pre-generation timestamp wording. Reports have publication context added, and local path prefixes in public documentation/log copies are normalized. [The transformation record](provenance/publication_transformations.json) maps those copies to their original hashes. Historical component manifests describe the complete local experiment; use this directory's new public manifest for the files actually distributed.

Phase II computations, code, checks and report drafting were produced with Codex assistance during the research session. The original experimental evidence and the limitations of each claim are documented explicitly; independent supervisory review remains appropriate for thesis or publication use.
