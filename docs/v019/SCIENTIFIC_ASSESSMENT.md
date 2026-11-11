# RNGLab after Stage C: which original contribution is defensible?

**Status: assessment for the author's decision, 2026-11-11.** It changes no protocol and no claim. It rests on results already in the repository, on the audit and the surrogate study of v0.19 (`IMPLEMENTATION_AUDIT.md`, `FORECASTING_METHODS.md`), and on a literature check whose limits are in Section 5.

## 0. Recommendation in brief

- **Direction A (reliable uncertainty under several time scales) is the stronger candidate for an original result today.** RNGLab already holds a clean, exactly referenced failure: a batch-length rule that was calibrated on AR(1) accepts every refreshed path-integral run and then covers 0.8828. The phenomenon behind it is known in lattice field theory. What is not established is whether the standard remedies restore calibration on this benchmark, and where exactly rules of this kind fail. One experiment on stored traces would settle most of that.
- **Direction B (degradation-aware prognosis) has no positive result yet.** Its solid results are negative or conditional: H1-C was not supported in any of four families; a forecast is impossible without a declared continuation, with a computable loss; and on a surrogate the forecasters work only inside the declared family. Its practical relevance is the weakest point: no field evidence was found that a generator degrades along a ramp.
- **Proposed framing, which keeps the central question as the author wrote it.** The question has two halves: how input defects propagate, and how a sampler's own uncertainty failure is told apart from a source fault. Direction A is the second half. It is the control without which no anomaly can be attributed to the source. It need not replace H1 to H3; it can be the part of the programme that carries a new result while H2 stays conditional on the pilot.
- **Whether A leads the TDK is the author's decision.** The roadmap treats estimator calibration as supporting work. This assessment only reports where the evidence is strongest.

## 1. Direction A: uncertainty estimation under complex stochastic dependence

**What RNGLab has shown** (`docs/v018/RESULTS_SUMMARY.md`; `research/v018-api-validation/REPORT.md` §4.1):
- On AR(1) and Ising the rule "batch at least 10·τ̂, τ̂ from Geyer's initial sequence" never accepted a run and then mis-covered (21 of 21 evaluable gates).
- On the refreshed path-integral sampler the same rule accepted all 4,096 runs and covered 3616/4096 = 0.8828 (CP95 [0.8726, 0.8925]). Geyer's τ̂ had median 3.1. The batch-means variance ratio was 0.601 at batch 32, 0.875 at 128 and 0.971 at 512. Fixed batches of 512 covered 0.939.
- The stated reading: a fast component plus a slow shape mode of small weight, missed by the truncation and by short batches. This is a hypothesis; no fit was made.

| Question | Assessment |
|---|---|
| **Prior art** | The mechanism is known. In lattice QCD, slow modes that couple weakly to an observable make window-based error estimates too small, and a tail built from an estimate of the slowest time scale is the usual remedy ([SSV11], [Wol04], [Ram19]). In statistics, batch-means and spectral estimators are consistent only when the batch grows with the run ([FJ10], [Dam94], [GW91]); batch-size selection is an active topic ([SS95], [LVF22]); blocking looks for a plateau [FP89]; initial-sequence truncation assumes reversibility [Gey92]. Reviews warn that τ_exp can far exceed τ_int [Sok97] |
| **Smallest defensible new claim** | For processes with two time scales, rules that set the batch as a multiple of an initial-sequence τ̂ have a failure region that can be written down in the weight and the ratio of the scales; an exactly solvable path-integral benchmark lies inside it; and a stated set of established methods does, or does not, restore nominal coverage there at the same run length. The claim is a characterisation and a benchmark. It is not "a new estimator" unless one is shown to beat the established ones |
| **Evidence required** | (1) The two-scale model fitted to the measured variance-ratio curve, with the failure region derived and checked by simulation. (2) Coverage of each baseline on the stored healthy traces against the exact energy. (3) A prospective confirmation on fresh seeds under a frozen protocol. (4) At least one more system: Ising is in the repository |
| **Independent baselines** | Fixed-batch t intervals; blocking plateau [FP89]; the Γ-method with and without a tail ([Wol04], [SSV11]); batch-size selection [LVF22]; Geyer's estimators [Gey92]; overlapping batch means and spectral variance [FJ10] |
| **Principal confounders** | Residual initialisation effects (the exact start removes them); the sweep is not reversible, which Geyer's theorem assumes; one parameter point (16 beads, β = 2); a harmonic system, where the centroid refresh is exact; the rule's constant was chosen on AR(1) |
| **Practical relevance** | High and general. Every Markov chain Monte Carlo result carries an error bar; path-integral practice relies on blocking |
| **What would falsify it** | An established method already gives calibrated intervals on this benchmark at this run length. The contribution is then a replication with an exact reference, which is still reportable but not new. Or the two-scale model does not fit the measured curve |
| **TDK suitability** | Good. The data exist (Phase 2: three conditions of 4,096 runs; Stage C: 4,096 healthy and 999 calibration runs with traces). The analysis is cheap. The result is informative either way |
| **Before a journal** | A second sampler with a reference that does not rely on harmonicity; the derivation written as a proposition with its assumptions; a database search forward from [SSV11], [Wol04] and [FJ10] |

## 2. Direction B: degradation-aware Monte Carlo prognosis

**What RNGLab has shown:**
- **H1-C not supported** in SHIFT, LAG1, DUP and GFSR at matched false alarms: the battery alarmed first or alone in every cell (`research/v018-h1-detection/REPORT.md`). A post-hoc analysis shows that for LAG1 this rested on the test aimed at it.
- **Identifiability:** two continuations with the same history give any history-based forecaster an excess Brier score of at least E[(p₁ − p₂)²]/4 in one of them. The size is 0.2488 in the Phase II counterexample (`FORECASTING_METHODS.md` §1).
- **Methods on a surrogate:** under the declared ramp a posterior forecaster is calibrated and has skill 0.228 [0.189, 0.268] over an alarm-time baseline; under two other families it is miscalibrated and worse than that baseline (skill −0.256 and −0.218).
- **The Stage D pilot is frozen and not run.** Nothing is known yet about forecasting in the real application.

| Question | Assessment |
|---|---|
| **Prior art** | Remaining-useful-life estimation is mature ([SWHZ11], [GLLR05], [LLGLL18]) with its own evaluation metrics [SCSSG10]. Health tests and online quality control act at the source ([SP800-90B], [GF25]). Application tests of generators are classical ([FLW92], [Cod96], [OS04], [BM07], [Hel98]). Runtime detection of silent data corruption is an adjacent field ([Sni14], [DC16]). No work was found that forecasts the validity of a Monte Carlo estimate under a degrading generator |
| **Smallest defensible new claim** | In one application, for synthetic faults and a declared ramp, forecasts of interval validity issued at a calibrated battery alarm are calibrated and beat an alarm-time baseline by a stated margin; outside the declared family they are not; and the loss of any history-based forecaster under an undeclared continuation is bounded below by a computable quantity. A laboratory characterisation, stated as such |
| **Evidence required** | The pilot's criteria passed; then a campaign with held-out rates and onsets and unseen families; a check that the sampler's batch means behave as the validity model assumes; an independent recount |
| **Independent baselines** | Base rate; alarm time; persistence; the oracle in its doubly noncentral form |
| **Principal confounders** | The faults are synthetic and the family is told to the forecaster. Two battery tests are aimed at two of the faults. The alarm time alone carries much of the information under a ramp. The alarm selects the statistic at its threshold. Short issuance intervals under-cover in healthy runs. One application at one parameter point |
| **Practical relevance** | Uncertain. Field data on one entropy source showed entropy rising with age [WSMG20]. No measured ramp-like degradation of a generator in use was found. Without it the declared family is a convenience |
| **What would falsify it** | Pilot criterion 2 or 4 fails; skill disappears on held-out rates; the batch-mean model fails in the sampler |
| **TDK suitability** | As the central positive claim: risky, since it depends on a pilot that may well be negative. As a reported programme with negative and conditional results: sound |
| **Before a journal** | A positive campaign; a second application; evidence that some real source degrades in a way the declared family covers, or a reframing as validity forecasting under input drift in general |

## 3. Side by side

| | Direction A | Direction B |
|---|---|---|
| Result in hand | an exactly referenced failure of a calibrated rule | negative H1-C; an impossibility bound; methods on a surrogate |
| Closest prior art | strong and close (lattice QCD, output analysis) | strong on methods, absent on the application |
| Risk of "already known" | **high** for the phenomenon, lower for the characterisation | low for the application, but relevance is questioned |
| Cost of the decisive experiment | hours of analysis on stored traces, plus one fresh confirmation | the frozen pilot (about nine hours of computing), then a campaign |
| Outcome if negative | a replication with an exact reference | "not feasible as designed", already foreseen by the stop rule |

## 4. Does RNGLab have a defensible original contribution today?

- **For a TDK: yes, of a specific kind.** A pre-registered, independently recounted negative result on application-aware detection at matched false alarms, and a documented transfer failure of an uncertainty rule with an exact reference. Neither is claimed as new beyond "not found by these searches".
- **For a journal: not yet.** No positive result is established. The identifiability bound is elementary. The forecasting results are on a surrogate.
- **Closest to a new positive result:** the Direction A characterisation, if the established remedies turn out not to restore calibration on the benchmark, or do so only at a stated cost.

## 5. Limits of this assessment

- Literature: three web searches on 2026-11-11 and Crossref checks of author, title, venue, year, volume and pages (`docs/v019/evidence/`). **Full texts were not read.** Statements about a paper's content rest on its abstract, on a search summary, or on secondary description. "Not found" means not found by these searches.
- Three of the first-pass entries had a wrong detail (a DOI that belonged to the neighbouring article, a publication year, an author order); they were corrected against Crossref (`refs_check_first_pass.json` against `refs_check.json`). [BA68] is flagged there only because Crossref's title carries markup; its metadata are in `crossref_lookup.json`.
- The mechanism in Section 1 is a hypothesis until the two-scale model is fitted.

## References (new in v0.19; Crossref-checked)

- [BM07] H. Bauke, S. Mertens. Random numbers for large-scale distributed Monte Carlo simulations. *Phys. Rev. E* 75, 066701, 2007. doi:10.1103/PhysRevE.75.066701
- [Dam94] H. Damerdji. Strong consistency of the variance estimator in steady-state simulation output analysis. *Math. Oper. Res.* 19(2), 494–512, 1994. doi:10.1287/moor.19.2.494
- [DC16] S. Di, F. Cappello. Adaptive impact-driven detection of silent data corruption for HPC applications. *IEEE Trans. Parallel Distrib. Syst.* 27(10), 2809–2823, 2016. doi:10.1109/TPDS.2016.2517639
- [FP89] H. Flyvbjerg, H. G. Petersen. Error estimates on averages of correlated data. *J. Chem. Phys.* 91, 461–466, 1989. doi:10.1063/1.457480
- [GF25] C. Gerolimetto Fabrello, L. Rossi, M. Witek, D. Trombetta, M. Baszczyk, P. Dorosz, W. Kucewicz, M. Caccia. Entropy measurement and online quality control of bit streams by a true random bit generator. *Front. Comput. Sci.* 7, 2025. doi:10.3389/fcomp.2025.1642566 (first names not checked)
- [GW91] P. W. Glynn, W. Whitt. Estimating the asymptotic variance with batch means. *Oper. Res. Lett.* 10(8), 431–435, 1991. doi:10.1016/0167-6377(91)90019-L
- [Hel98] P. Hellekalek. Good random number generators are (not so) easy to find. *Math. Comput. Simul.* 46(5–6), 485–505, 1998. doi:10.1016/S0378-4754(98)00078-0
- [LLGLL18] Y. Lei, N. Li, L. Guo, N. Li, T. Yan, J. Lin. Machinery health prognostics: a systematic review from data acquisition to RUL prediction. *Mech. Syst. Signal Process.* 104, 799–834, 2018. doi:10.1016/j.ymssp.2017.11.016 (first names not checked)
- [LVF22] Y. Liu, D. Vats, J. M. Flegal. Batch size selection for variance estimators in MCMC. *Methodol. Comput. Appl. Probab.* 24(1), 65–93, 2022 (online 2021). doi:10.1007/s11009-020-09841-7
- [OS04] G. Ossola, A. D. Sokal. Systematic errors due to linear congruential random-number generators with the Swendsen–Wang algorithm: a warning. *Phys. Rev. E* 70, 027701, 2004. doi:10.1103/PhysRevE.70.027701
- [Ram19] A. Ramos. Automatic differentiation for error analysis of Monte Carlo data. *Comput. Phys. Commun.* 238, 19–35, 2019. doi:10.1016/j.cpc.2018.12.020
- [Sni14] M. Snir et al. Addressing failures in exascale computing. *Int. J. High Perform. Comput. Appl.* 28(2), 129–173, 2014. doi:10.1177/1094342014522573
- [Sok97] A. D. Sokal. Monte Carlo methods in statistical mechanics: foundations and new algorithms. In *Functional Integration*, NATO ASI Series, 131–192, 1997. doi:10.1007/978-1-4899-0319-8_6
- [SS95] W. T. Song, B. W. Schmeiser. Optimal mean-squared-error batch sizes. *Manage. Sci.* 41(1), 110–123, 1995. doi:10.1287/mnsc.41.1.110
- [SSV11] S. Schaefer, R. Sommer, F. Virotta. Critical slowing down and error analysis in lattice QCD simulations. *Nucl. Phys. B* 845, 93–119, 2011. doi:10.1016/j.nuclphysb.2010.11.020
- [Wol04] U. Wolff. Monte Carlo errors with less errors. *Comput. Phys. Commun.* 156(2), 143–153, 2004. doi:10.1016/S0010-4655(03)00467-3
- Already in the programme: [FJ10], [Gey92], [FLW92], [Cod96], [SWHZ11], [GLLR05], [SCSSG10], [SP800-90B], [WSMG20] (`docs/v018/RELATED_WORK.md`, `STAGE_DE_REVIEW.md`).
