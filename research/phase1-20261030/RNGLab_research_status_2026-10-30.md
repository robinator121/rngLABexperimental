# RNGLab research reconnaissance and status

**Prepared:** 30 October 2026  
**Repository audited:** [`rngLABexperimental`, main at `dd96adba`](https://github.com/robinator121/rngLABexperimental/commit/dd96adba23b8100aba9147d71fffc63932fc6b0d)  
**Scope:** current-state audit, literature and novelty check, new PIMC pilot, physical-data inventory, and research roadmap.

## Executive findings

RNGLab has grown into a broad research prototype: it contains RNG implementations and fault injectors, Monte Carlo and Markov-chain models, physical acquisition and entropy-source code, monitoring and forecasting experiments, and financial-risk workloads. The most mature current line is the v0.17 CALIBRATE work, which separates healthy behavior, detection, prediction and decision usefulness more carefully than earlier versions.

Three findings should guide the next phase:

1. **The highest-value scientific gap is still healthy uncertainty calibration.** The project's own v0.18 priorities place this ahead of more monitoring or forecasting. My independent 64-chain scalar PIMC replay found only **20/64, 38/64 and 48/64** coverage for nominal 95% CV intervals at batch sizes 32, 128 and 512. Those counts are close to the project's documented **21/64, 39/64 and 51/64**. The finite-bead reference and interval quantiles passed independent numerical checks. This is exploratory evidence, not confirmation.
2. **The public evidence package is incomplete.** Current v0.17 documents refer to calibration JSON, analysis files and full NPZ traces that are absent from the checked-in main tree. The latest public GitHub Actions run is marked failed in both Python jobs, but the API exposes no job steps or usable logs, so the failure's cause is unknown. The local full test suite could not start because this runtime lacks pytest, NumPy and SciPy.
3. **The detector exports are useful measurements, not yet a validated random source.** I read the full primary file and verified its SHA-256 against the progress report. It has 117,962 records with 1,023 waveform codes each, a constant first field, an 80-unit timestamp lattice, and many ceiling codes. A second export has a different row width (64 waveform values rather than 1,023). Neither file establishes independent random bits, detector identity, calibrated timestamps or NIST entropy compliance.

The defensible BSc route is to finish the healthy PIMC uncertainty baseline, then add one prespecified RNG fault and one application-level reliability endpoint. Broad claims about a universal detector, physical entropy generation, reliable remaining-useful-samples forecasts or general HPC acceleration are not supported yet.

## What I inspected

- Current public GitHub `main` at commit `dd96adba23b8100aba9147d71fffc63932fc6b0d`; local clone matched that commit.
- Source packages, experiment scripts, tests, README, current scientific reports, validation notes, reproducibility notes, benchmark documents and the latest CI run.
- The eight-page progress report *RNGLab - Véletlenforrások és Monte Carlo-becslések megbízhatósága* and its original Xournal file.
- Both detector exports under `~/rnglab-data`, using the repository's standard-library export inspector. The primary file's SHA-256 was independently checked; the source files were not edited.
- A separate standard-library PIMC replay in this output folder. No files in the checked-out RNGLab repository or its main branch were changed.

Evidence labels below distinguish **directly rechecked** facts from results **documented by the project but not rerun here**. A documented value is not treated as independently reproduced merely because it appears in a report.

## Current-state audit

### Architecture and maturity

| Area | What is present | Evidence status |
|---|---|---|
| Randomness sources | Seeded PCG/MT/Xoshiro/ChaCha-style sources, file-backed sources, channel abstractions and controlled bias/dependence/repetition faults | Implemented research code; no claim here that each implementation is statistically or cryptographically certified |
| Physical acquisition | Arduino/serial acquisition components, offline digitized-waveform inspection, entropy/health-test modules described as NIST-inspired | Pipeline and software exist; physical source model and NIST validation remain unestablished |
| Applications | Classical integration, Ising Metropolis/Wolff, VMC and PIMC harmonic oscillator, plus financial risk examples | Multiple runnable prototypes; scientific acceptance varies by workload and saved experiment |
| Reliability and monitoring | Source/application monitors, threshold calibration, multiplicity control, healthy/fault cohorts and several prediction/action policies | v0.17 is substantially more disciplined, but its own reports show failed or scoped gates |
| Forecasting | Fixed-horizon incident risk and a preliminary remaining-useful-samples workflow | A narrowly evaluated classifier exists; alarm-conditioned RUS is not demonstrated |
| GPU/HPC | Laptop GPU benchmarking and separate L4 audit materials; the University of Debrecen also operates the Komondor HPC system | Real infrastructure exists, but an allocation for this project is not confirmed; current work does not need it |

The repository is versioned by successive research generations (`rnglab11` through `rnglab17`) plus finance and physical-acquisition areas. This preserves history but makes it easy to confuse a legacy model, a current implementation and a saved result from an earlier version. In particular, the PIMC family has both a user-facing `rnglab14` model and the vectorized `rnglab16` campaign implementation, with `rnglab17` calibration helpers. The thesis package should identify one supported path and clearly label archived outputs.

### Reported findings: supported, limited and still unverified

The main project sources are its [v0.17 scientific report](https://github.com/robinator121/rngLABexperimental/blob/dd96adba23b8100aba9147d71fffc63932fc6b0d/docs/calibrate/V017_SCIENTIFIC_REPORT.md), [PIMC validation note](https://github.com/robinator121/rngLABexperimental/blob/dd96adba23b8100aba9147d71fffc63932fc6b0d/docs/calibrate/V017_VMC_PIMC_VALIDATION.md), [predictive-reliability report](https://github.com/robinator121/rngLABexperimental/blob/dd96adba23b8100aba9147d71fffc63932fc6b0d/docs/calibrate/V017_PREDICTIVE_RELIABILITY.md), and [physical detector assessment](https://github.com/robinator121/rngLABexperimental/blob/dd96adba23b8100aba9147d71fffc63932fc6b0d/docs/physical/MUON_DETECTOR_ASSESSMENT.md). The values below are documented in those files; they were not all rerun in this session.

| Finding | What the evidence supports | What it does not support |
|---|---|---|
| Repeated/periodic streams can return plausible point estimates with poor interval coverage | The progress report documents a specific controlled repetition design with near-equivalent means/variances and nominal 95% coverage between 1.11% and 34.38% | A universal property of periodic sources, or evidence that every healthy estimator is calibrated |
| VMC zero-variance blind spot | At the harmonic ground-state construction α=0.5, local energy is constant while (x^2) can reveal a sampling problem; the report documents 4/128 healthy and 45/64 faulty detections for one calibrated rule | A generally effective fault detector or an application-independent threshold |
| PIMC intervals under-cover | The progress report and v0.17 notes show poor healthy-chain coverage for short blocks; a separate long-chain CV probe improved from 21/64 to 51/64 as batch size increased from 32 to 512 | That longer chains alone fix the problem, or that a particular cause (burn-in, mixing or interval method) has been isolated |
| Combined monitoring and forecasts | v0.17 reports a restricted L3/Tc test Brier score of 0.12078 versus a 0.24206 prevalence baseline | General superiority: the OOD score gain is inconclusive, the full no-mitigation decision cost is worse than always-stop, and the warnings were not compared at matched false-positive rates |
| RUS | The report describes a fixed-time, non-alarm-triggered analysis with 24 overoptimistic positive bounds out of 234, and 18 abstentions among 252 eligible paths | Alarm-conditioned forecasting or a general RUS guarantee |
| Healthy calibration | The v0.17 report says 10/12 primary cells met its stated simultaneous upper-error ceiling; L4/Tc and PIMC16 failed with 38/512 events each | A verified 5% conditional false-alarm guarantee or a fully calibrated multi-model deployment |

The project itself withdraws universal early-warning claims and states that the current platform is not a generally validated reliability product. That caution is appropriate and should remain prominent.

### Scientific and reproducibility risks

**PIMC references.** The current PIMC source contains the corrected thermodynamic and centroid-virial estimators and an exact finite-bead reference. At M=16 and beta=2, the finite-bead energy is 0.655474, while the continuum value is 0.656518. Sampling uncertainty and discretization error must stay separate. The current code uses the finite-bead reference for its uncertainty summaries and also reports a continuum bias field. I found no new formula error in the estimator code I inspected; this is not a full line-by-line audit of every model. The historical integrated v0.14 JSON contains large PIMC discrepancies and should be regenerated before it is presented as evidence for the corrected implementation.

**Missing current outputs.** The current v0.17 reports refer to `calibration_analysis.json`, `forecast_analysis.json`, `physics-controls/results.json`, trace NPZ files and other artifacts. Those files are not in the current main tree. Older v0.2 synthetic and fault-injection outputs exist locally, but they are not a substitute for the current v0.17 traces. I found no current PIMC trace bundle in the accessible Codex workspace.

**Tests and CI.** The [README](https://github.com/robinator121/rngLABexperimental/blob/dd96adba23b8100aba9147d71fffc63932fc6b0d/README.md) advertises 238 passing tests and shows a green test badge. The [latest workflow run](https://github.com/robinator121/rngLABexperimental/actions/runs/40931415718) is marked `failure` for both Python 3.11 and 3.12 jobs; GitHub's job summaries contain no steps, and the log endpoint was unavailable. This therefore does not show whether an assertion failed or a runner failed to start. The local command `PYTHONPATH=src python -m pytest` stopped before collection because `pytest` is absent; a separate import check also found no NumPy or SciPy. No dependencies were installed in this session. The number of passing tests in the README is not verified here.

**Benchmark claims.** The repository has one [report](https://github.com/robinator121/rngLABexperimental/blob/dd96adba23b8100aba9147d71fffc63932fc6b0d/docs/benchmarks/RNGLAB_CPU_GPU_COMPARISON.md) describing an RTX 4070 laptop comparison at tiny batch sizes where the compiled CPU was faster, and another [L4 audit](https://github.com/robinator121/rngLABexperimental/blob/dd96adba23b8100aba9147d71fffc63932fc6b0d/docs/benchmarks/GPU_L4_AUDIT_SUMMARY.md) with different timing scopes. The README links the laptop report under an NVIDIA L4 claim. These should be separated and linked to their exact run manifests. Komondor is a real Debrecen HPC resource, but that does not imply a user allocation or a speedup for a particular RNGLab workload.

**Physical measurements.** The primary export's full integrity and layout counts are verified by the attached inventory. The detailed pulse and timing conclusions below come from the project assessment report, not from a second full waveform-feature analysis in this session.

| File | Directly verified | Important limit |
|---|---|---|
| [`meas_data_2026-10-19 14.05.19.txt` inventory](../phase2-20261030/DATA_AVAILABILITY.md#omitted-large-and-local-artifacts) | 606,632,259 bytes; SHA-256 `feba450d95b074d3b0fff50df8357a682e90ea62614a5948e3e4e106149b4e10`; 117,962 well-formed rows, 1,025 integer fields, 1,023 assumed waveform codes, all CRLF; field 0 always zero; timestamp gaps and values share GCD 80; waveform range 0–65,535 | No header, units, sample-rate calibration, exporter schema or confirmed mapping from triggers to muons. The assessment report describes 91.96% of records reaching the ceiling and near-fixed pulse alignment. |
| [`meas_data_2026-10-19 14.30.18.txt` inventory](../phase2-20261030/DATA_AVAILABILITY.md#omitted-large-and-local-artifacts) | 35,928,158 bytes; SHA-256 `b03a780d326d35c54a8916774209226ca5a44f6410459fe782391ce84d7a6964`; 103,336 rows, 66 integers per row, 64 assumed waveform codes, all CRLF; 544,145 of 6,613,504 codes equal 65,535 | Different width suggests a changed acquisition/export protocol. It cannot be pooled with the first file until that change is explained. |

The larger export is event-triggered, highly aligned and heavily clipped according to the project assessment. It is not a free-running list of independent noise symbols. Its raw timestamp low bits are fixed by the 80-unit lattice. A balanced bit histogram does not provide a min-entropy estimate. A one-file recording is also replayable, so it cannot be an ongoing source of fresh entropy. SP 800-90B expects at least one million sequential source samples for the main assessment and separate restart data; if one detector event were the declared source symbol, this file would be too small. Treating each of its 120.7 million ADC codes as a symbol would require a justified source/dependence model and byte mapping, and would not make those codes independent. No 90B entropy estimate was run here.

## Literature and novelty assessment

The core framing is scientifically worthwhile, but several broad formulations already have direct predecessors. The potentially new work is the *specific controlled evidence and operational protocol*, not the general idea of testing an RNG inside an application.

### What is already established

**Statistical source tests are not certification.** [TestU01](https://doi.org/10.1145/1268776.1268777) is a mature empirical test library for generators. NIST [SP 800-22](https://csrc.nist.gov/pubs/sp/800/22/r1/upd1/final) is a statistical test suite, and [NIST's 2022 revision decision](https://www.nist.gov/news-events/news/2022/04/decision-revise-nist-sp-800-22-rev-1a) explicitly says its use should be clarified and rejected for assessing cryptographic RNGs. For an entropy source, [SP 800-90B](https://csrc.nist.gov/pubs/sp/800/90/b/final) is the relevant source-design, entropy-estimation and health-testing framework; [SP 800-90A](https://csrc.nist.gov/pubs/sp/800/90/a/r1/final) specifies DRBGs, and final [SP 800-90C](https://csrc.nist.gov/pubs/sp/800/90/c/final) (September 2025) specifies constructions combining sources and DRBGs. The 90B publication page currently notes two identified errata. Keep the physical noise source, any conditioning function and any DRBG as separate pipeline stages with separate evidence. A passing SP 800-22 battery cannot establish source entropy, cryptographic unpredictability or downstream correctness.

**Application-based RNG testing is not new.** [Sim and Nitschke (1993)](https://pubmed.ncbi.nlm.nih.gov/8470994/) compared generators inside a radiation-transport Monte Carlo application. [Coddington (1996)](https://surface.syr.edu/npac/34/) tested generators through Ising simulations with Metropolis, Swendsen-Wang and Wolff updates. RNGLab should not claim novelty for the broad proposition that application behavior can reveal RNG problems.

**MCMC variance estimation has an established theory.** [Geyer's (1992)](https://doi.org/10.1214/ss/1177011137) initial-positive/monotone sequence and later [batch-means and spectral-variance work](https://doi.org/10.1214/09-AOS735) require chain assumptions and sufficient dependence information. A Geyer ESS, standard error or batch interval is diagnostic; it is not a certificate that a chain is stationary or has mixed. The project's own uncertainty module says this explicitly. [Murray and Elliott (2012)](https://arxiv.org/abs/1204.3187) also show that some MCMC algorithms can be modified to tolerate a dependent random stream while preserving the target invariant distribution. This rules out a blanket claim that any dependent stream necessarily corrupts every MCMC procedure.

**Cosmic-radiation entropy generation already exists.** [Kutschera et al. (2023)](https://doi.org/10.3390/e25060854) report an MRNG prototype using cosmic-ray detections. A newer [arXiv preprint](https://arxiv.org/abs/2608.24729) submitted 25 August 2026 reports extracting entropy from muon-track directions over six months. It is a preprint, not evidence that the RNGLab waveform exports are suitable as an entropy source. The data representation, detector, entropy assessment and application question differ.

**Change detection and censored prediction are established methods.** [Page's 1954 continuous inspection scheme](https://doi.org/10.1093/biomet/41.1-2.100) is an origin for cumulative-sum monitoring. [Conformalized survival analysis](https://doi.org/10.1093/jrsssb/qkac004) and [Mondrian predictive systems for censored data](https://proceedings.mlr.press/v204/bostrom23a.html) already provide calibrated time-to-event bounds or distributions under stated assumptions. Applying such tools to simulation degradation is an application question; the existence of an alarm-conditioned RUS target is not by itself a new prediction method.

**Parallel streams are a known requirement.** [L'Ecuyer and colleagues](https://doi.org/10.1016/j.matcom.2016.05.005) discuss how to construct and manage many streams for parallel computers and GPUs. HPC work should use deliberately designed stream splitting or counter-based generators, record stream identities, and test for cross-worker overlap; assigning nearby integer seeds is not a sufficient methods section.

### Defensible contribution candidates

1. **Best BSc contribution:** a reproducible benchmark of how one specified RNG defect changes uncertainty coverage in a validated exact-reference PIMC workload, after the healthy baseline itself has been calibrated. Report point-estimate bias, interval coverage, source diagnostics and application diagnostics separately.
2. **Follow-on monitoring contribution:** determine whether an application monitor detects that defect earlier or more usefully than a source monitor at a matched whole-run false-alarm rate, using held-out independent trajectories and explicit misses.
3. **Physical-source contribution:** characterize a defined raw symbol stream from a documented detector and acquisition chain, estimate entropy under a justified 90B route, then measure effects on a downstream simulation. This needs new settings-rich, restart-rich data; current waveform files alone are not enough.
4. **Later RUS contribution:** estimate an alarm-conditioned distribution or conservative bound for useful remaining simulation work under a specified future degradation scenario, with trajectory-level train/test splits, censoring, calibration and abstention. Distinguish raw samples from effective samples.
5. **Recovery contribution:** compare source reseed, PRNG state reset, chain reinitialization and full simulation restart. Reseeding changes future draws; it does not erase biased samples already accumulated or automatically put an existing MCMC state back into the target distribution.

These are **potential** contributions. A literature search cannot prove novelty. Before any paper claim, compare against the nearest papers above and do a broader database search with the final question and fault definition.

## Experiment completed: PIMC pilot

I ran [the scalar pilot](pimc_uncertainty_pilot.py), with [methodology](PIMC_uncertainty_pilot_methodology.md), [machine-readable results](pimc-uncertainty-pilot-results.json), [chain and seed summary](pimc-uncertainty-pilot-chain-summary.csv), [complete compressed traces](../phase2-20261030/DATA_AVAILABILITY.md#omitted-large-and-local-artifacts) and [coverage plot](pimc-uncertainty-pilot-coverage.svg). It completed in 25.7 seconds using the standard library.

The reference check returned finite-bead energy 0.6554741523, continuum energy 0.6565176427, and a 1.9e-11 difference from the independent partition-function derivative. Student-t critical-value checks passed at all degrees of freedom used. The long CV window yielded 31.3%, 59.4% and 75.0% coverage at the three batch sizes; the corresponding Wilson intervals were 21.2–43.4%, 47.1–70.5% and 63.2–84.0%. A second run returned identical coverage values and byte-identical trace, chain-summary and figure files.

This is a useful independent *qualitative* reproduction of the documented pattern, not a final validation. It suggests the first research question should separate initialization, retained length and block length; the current pilot changes burn-in and retained length together. It also reinforces that healthy-source MCMC dependence can defeat an uncertainty estimator, so RNG-source faults must not be blamed before this baseline is sound.

### Execution and tests

- The full repository test suite did not run locally: `pytest`, NumPy and SciPy are missing from the current environment, so module collection was not reached.
- The public latest CI run is marked failed for both Python versions, but the returned job records have empty step lists and the logs were unavailable. The failure cause is unresolved.
- The pilot itself ran and its internal exact-reference and t-quantile checks passed. It did not exercise the repo package or replace its test suite.
- I did not install dependencies or edit the main repository.

## Academic path and collaborators

### BSc thesis

The current topic title and ID in your brief are a good umbrella, but the submission should focus on one falsifiable claim. Proposed title:

> **Effects of controlled random-number defects on confidence-interval coverage in a finite-bead PIMC benchmark**

Start with a healthy PIMC uncertainty baseline; use one exact finite-bead reference, one primary estimator and one prespecified source fault. Include one source-level diagnostic and one application-level endpoint. Add a second model only after the first passes the held-out criteria. This keeps the BSc tied to RNG reliability while correcting the project's present healthy-control weakness.

### TDK

The Informatics Faculty's [official TDK page](https://www.ik.unideb.hu/smartcampus/tdk/) lists the 2026 conference for **5 November 2026**, with submission deadline **26 October 2026 at 12:00**. As of this report the deadline has passed. If you submitted, prepare the result story and the limits above for the presentation. If not, ask the coordinator when the next call opens; I found no next-cycle date on the current page.

### People and infrastructure

| Opportunity | Confirmed public information | Still unconfirmed |
|---|---|---|
| Dr. Tamás Herendi | University profile lists [PRNGs, number theory, cryptography, parallel algorithms and FPGA programming](https://inf.unideb.hu/en/dr-tamas-herendi-associate-professor); a strong fit for the RNG and reproducibility core | Your allocation to Komondor, desired thesis scope and access to particular hardware |
| Dr. Anett Rácz | University profile lists her in [Applied Mathematics and Probability](https://tudoster.unideb.hu/hu/szerzok/1431); recent work includes simulation analysis of multi-objective cutting problems and optimization | The meeting and possible Polish/Czech/Italian links are from your brief; I did not verify those arrangements externally |
| Dr. Gselmann Eszter | Your brief says she reviewed an earlier summary | You said the core questions were outside her current research expertise, so I have not treated this as an active collaboration |
| BME Nuclear Techniques Institute | Official research page lists [gamma/alpha spectrometry, neutron activation analysis and Monte Carlo calibration](https://www.reak.bme.hu/en/research/radioanalytics-and-chemistry-and-radiochemistry.html), a plausible methods collaboration | Detector availability, supervision, facility permission and any radiation-safety approval |
| Debrecen Komondor HPC | Official HPC page describes [CPU, GPU, AI and big-data partitions](https://hpc.unideb.hu/komondor), including A100 GPU nodes | A personal account or allocation is not established. Ask before planning the thesis around it. |

The direct thesis-topic page for ID 25682 did not load in the web viewer. I therefore treat the title and ID as supplied by you; current faculty pages verify the supervisor's research fit, not your individual registration record.

## Revised roadmap

### 1. Immediate: make the evidence reconstructible

- Freeze the exact commit, environment lock, data schema and configuration used for each primary result.
- Restore or regenerate the v0.17 JSON/NPZ files referred to by the reports. Publish trace files or stable, versioned artifact links with hashes and a regeneration command.
- Resolve the CI job-start failure and get a successful run on both supported Python versions before repeating broad claims. Align the README test badge and GPU label with what the actual runs show.
- Keep old v0.14 results marked as historical; regenerate the integrated result after the corrected PIMC estimator and finite-bead target are in place.

### 2. Main science: confirm healthy PIMC uncertainty

- Freeze (M=16,\beta=2,\delta=0.4), finite-bead reference and a primary endpoint before a confirmation run.
- Use an exact finite-bead IID/Gaussian path control and the local Metropolis control. Keep the two source/model conditions separate.
- Run enough independent chains to estimate coverage precisely (512 is a practical first target, with more if the chosen coverage margin requires it). Store one full trace per independent seed.
- Change initialization, burn-in and retained length one at a time; ensure each batch-size cell has adequate batches. Keep batch-32/128/512 comparisons secondary until a primary estimator is chosen.
- Use a calibration/development set and a new locked confirmation namespace. Report binomial intervals for coverage, point-estimate error and abstentions. Do not tune thresholds on the confirmation cohort.

### 3. Then: one controlled source fault and fair monitoring comparison

- Select one fault tied to a documented prior, such as periodic block reuse or a specified dependence model.
- Compare healthy and faulty streams through identical application code and workload settings. Log the source stream and application trajectory separately.
- Calibrate source and application monitors on development trajectories; compare them at matched whole-run false-positive rates on independent held-out trajectories. Report misses and lead time, not only a combined detection score.

### 4. Physical source and independent-system comparison

- Before deriving bits, obtain the digitizer model, native resolution, scaling, sample interval, trigger threshold/holdoff, timestamp unit, exporter version, dropped-event counters, detector identity and operating conditions.
- Acquire separate free-running baseline, source-control and restart runs with raw data and settings. Explain why the 1,023-sample and 64-sample formats differ.
- Define the source symbol and entropy model before estimating entropy. If a record is one symbol, the main file has only 117,962 samples; if an ADC code is a symbol, justify the dependence model and byte mapping. Collect restart data. Run the official 90B tool only after those inputs are defensible.
- If there are two detectors, characterize each separately. Feed each through the same documented conditioner/generator and workload protocol. Compare across-run error, coverage, convergence and reliability distributions; do not expect identical raw streams or outputs.

### 5. Forecasting, recovery and HPC

- Keep alarm-conditioned RUS for after a stable failure definition exists. Define forecast time, future degradation regime, event horizon, censoring, positive-bound calibration and abstention before fitting. Compare with simple hazard/survival baselines and a static policy.
- For recovery, distinguish source reseeding from resetting the sampler state and from restarting the whole simulation; evaluate residual bias in already accumulated estimates.
- The current 64-chain pilot is small enough for local compute. Use Komondor only when independent-chain calibration or parallel-stream experiments exceed local resources. Make random streams deterministic across workers and use a pre-specified CPU baseline before any GPU-performance claim.

### 6. Publication and longer study

The project is not ready for a general reliability or entropy-source paper. A narrow application-methods paper could become plausible after a locked held-out PIMC/Ising study, public traces, a passing reproduction path and external statistical review. A physical-source paper requires a separate instrument/entropy contribution. MSc/PhD work could connect validated physical sources, parallel RNG streams and application-level calibration across more than one workload after the BSc question is settled.

## Action list

| Owner | Next action |
|---|---|
| You + Dr. Herendi | Agree on the one-sentence BSc claim, primary PIMC endpoint, fault family and whether the existing v0.17 evidence can be recovered from archived traces. Ask about Komondor eligibility and how to obtain a student allocation. |
| You + Dr. Rácz | If the Monte Carlo decision-support collaboration proceeds, choose one real decision objective and cost function; keep it as a second application after the PIMC baseline is fixed. |
| You + BME contact, if approved | Request instrument/export settings and discuss a supervised, settings-rich restart acquisition. Treat facility access as pending until explicitly approved. |
| You + TDK coordinator | If the November submission was not made, ask for the next submission cycle and dates. If it was made, prepare a talk focused on the narrow evidence and limitations. |
| Research software | Restore trace artifacts, lock dependencies, resolve CI startup and publish manifests/hashes before using the repository as a thesis reproduction package. |

## Sources

### RNG standards and empirical testing

- NIST, [SP 800-22 Rev. 1a](https://csrc.nist.gov/pubs/sp/800/22/r1/upd1/final), 2010; and [NIST's 2022 decision to revise it](https://www.nist.gov/news-events/news/2022/04/decision-revise-nist-sp-800-22-rev-1a).
- NIST, [SP 800-90A Rev. 1: DRBG mechanisms](https://csrc.nist.gov/pubs/sp/800/90/a/r1/final), 2015.
- Turan et al., [SP 800-90B: entropy sources](https://csrc.nist.gov/pubs/sp/800/90/b/final), 2018 (page notes two identified errata).
- Barker et al., [SP 800-90C: RBG constructions](https://csrc.nist.gov/pubs/sp/800/90/c/final), final publication, 25 September 2025.
- L'Ecuyer & Simard, [TestU01](https://doi.org/10.1145/1268776.1268777), *ACM Transactions on Mathematical Software* 33(4), 2007.

### Application behavior and MCMC uncertainty

- Sim & Nitschke, [Testing random number generators for Monte Carlo applications](https://pubmed.ncbi.nlm.nih.gov/8470994/), 1993.
- Coddington, [Tests of random number generators using Ising model simulations](https://surface.syr.edu/npac/34/), 1996 technical report.
- Geyer, [Practical Markov Chain Monte Carlo](https://doi.org/10.1214/ss/1177011137), *Statistical Science* 7(4), 1992.
- Flegal & Jones, [Batch means and spectral variance estimators in MCMC](https://doi.org/10.1214/09-AOS735), *The Annals of Statistics* 38(2), 2010.
- Murray & Elliott, [Driving MCMC with a dependent random stream](https://arxiv.org/abs/1204.3187), 2012.
- Page, [Continuous inspection schemes](https://doi.org/10.1093/biomet/41.1-2.100), *Biometrika* 41(1–2), 1954.

### Physical sources, prediction, parallel RNG and institutions

- Kutschera et al., [MRNG: Accessing cosmic radiation as an entropy source](https://doi.org/10.3390/e25060854), *Entropy* 25(6), 2023.
- Samuel, [Cosmic muon arrival directions as a source of entropy](https://arxiv.org/abs/2608.24729), arXiv preprint, submitted 25 August 2026.
- Candès, Lei & Ren, [Conformalized survival analysis](https://doi.org/10.1093/jrsssb/qkac004), *JRSS B* 85(1), 2023.
- Boström, Linusson & Vesterberg, [Mondrian predictive systems for censored data](https://proceedings.mlr.press/v204/bostrom23a.html), PMLR 204, 2023.
- L'Ecuyer, Munger, Oreshkin & Simard, [Random numbers for parallel computers](https://doi.org/10.1016/j.matcom.2016.05.005), *Mathematics and Computers in Simulation* 135, 2017.
- University of Debrecen, [Dr. Tamás Herendi profile](https://inf.unideb.hu/en/dr-tamas-herendi-associate-professor), [Rácz Anett profile](https://tudoster.unideb.hu/hu/szerzok/1431), and [current Faculty Informatics TDK page](https://www.ik.unideb.hu/smartcampus/tdk/).
- University of Debrecen, [Komondor HPC specifications](https://hpc.unideb.hu/komondor); BME Institute of Nuclear Techniques, [radioanalytics and radiochemistry](https://www.reak.bme.hu/en/research/radioanalytics-and-chemistry-and-radiochemistry.html).

For editable BibTeX, see [RNGLab_literature.bib](RNGLab_literature.bib).
