# Related work and novelty review (v0.18 Stage B, widened Phase 3.5)

**Status: adopted by the author, 2026-11-06 ("Adopt GO").** Decision: GO for Stage C, with the H1 claim narrowed as in Section 1. Stages D and E are conditional GO, to be revisited after Stage C (revisited in `STAGE_DE_REVIEW.md`, decided 2026-11-08). The review was carried out by Claude Code (Opus 5.5) using web searches on 2026-11-06. Every reference in the list was checked against a publisher or indexing record for its authors, title, venue, year and DOI. Full texts were **not** read unless an entry says so. Statements about a paper's content therefore rest on its abstract and on bibliographic or secondary summaries. Section 6 sets out the limits.

Scope, as adopted in `ROADMAP.md` (decision 2): the whole RNGLab sequence (degradation, detection, lifetime forecast, intervention) across five areas:
1. application-based RNG testing;
2. online RNG health tests;
3. MCMC diagnostics;
4. prognostics and remaining-useful-life (RUL) estimation;
5. restart and intervention in Monte Carlo.

## 1. Bottom line: stop/go per stage

| Stage | Claim RNGLab could make | Closest prior work found | Status of the claim | Recommendation |
|---|---|---|---|---|
| C (H1) | Application observables detect RNG degradation while conventional tests still pass | Application tests catching generators that pass standard tests: [FLW92], [VAK94], [Cod94], [Cod96], [SMC03], [CLK11]. Online health tests at the source: [SP800-90B], [AIS31], [VRV12], [RYMV16] | **Not new as an observation**; established since 1992. **Not found:** a calibrated comparison in which the application monitors and a conventional battery run online on the same consumed stream, at matched whole-run false-alarm rates, with latencies in sweeps and bytes and their uncertainty, and with healthy sampler failure separated out | **GO**, with the claim narrowed to that calibrated, latency-matched comparison. Prior work must be cited as establishing the phenomenon |
| D (H2) | Remaining useful lifetime of a Monte Carlo application under RNG degradation, with prediction intervals | RUL methodology and its evaluation: [SWHZ11], [SCSSG10]. Field ageing data for an entropy source: [WSMG20] | **No RNG or entropy-source RUL work found.** The methods are mature in prognostics. Ageing data show entropy can rise rather than fall, so a degradation trajectory cannot be assumed | **GO, conditional on C**, with a declared continuation model (RNGLab's identifiability result) and the prognostics evaluation metrics of [SCSSG10]. The novelty is modest: an established methodology applied to a new domain, plus the identifiability boundary |
| E (H3) | Targeted RNG interventions extend the application's lifetime, measured in matched branched experiments | Software rejuvenation [HKKF95]; PRNG recovery after compromise through reseeding [DPRVW13]; restart strategies [LSZ93]; coupled chains [JOA20] for the matched design | **No measured lifetime extension of a Monte Carlo application through an RNG intervention found.** The intervention idea is old in software ageing and in PRNG design | **Conditional GO after C and D.** The novelty lies in the measurement, not in the idea of intervening. This stage needs the most engineering |

**Overall: GO for Stage C**, with two conditions:
- **Claim wording.** Wherever H1 is stated, the claim is the calibrated, matched-false-alarm, latency comparison with sampler failure controlled. It is never "application tests find what standard tests miss", which [FLW92], [Cod94] and [VAK94] already show.
- **A non-trivial conventional arm.** Section 4 discusses this. Comparing against tests built only to catch catastrophic failure would be a straw man.

## 2. Area by area

### 2.1 Application-based RNG testing

- **The phenomenon is established.** With the Wolff cluster algorithm, "high quality" generators (R250 among them) gave systematic errors in Ising simulations with exact answers [FLW92].
  - [VAK94] proposed three "physical" tests: Wolff autocorrelation times, random walks and n-blocks. They traced such errors to short-range correlations in shift-register generators.
  - [Cod94] reported that many widely used generators that pass standard statistical tests fail Ising Monte Carlo tests; [Cod96] extended these tests (Metropolis, Swendsen–Wang, Wolff).
  - [MB04] attributed the cluster-algorithm failures to the generators' low conditional entropy.
  - [CLK11] found that a standard LCG changed the physical results of liquid and peptide simulations.
  - [SMC03] combined statistical and application-based tests for parallel generators (SPRNG).
- **Batteries.** TestU01 [LS07] supplies the SmallCrush, Crush and BigCrush batteries. Lattice structure in LCGs is classical [Mar68].
- **What these works do not do**, judging from their abstracts:
  - They test a *fixed generator* offline against known answers. None monitors a *degrading* source during a run.
  - None calibrates the application test and the conventional battery to a common false-alarm rate, or compares detection latency.
  - None separates a sampler's own uncertainty failure from an RNG fault. RNGLab Phases 1–2 build exactly that control.
- **One search found no head-to-head sensitivity study** of application tests against BigCrush. One documentation source notes that multi-stream application-style tests and long-stream batteries can each miss what the other finds (Section 6).

### 2.2 Online RNG health tests

- **NIST SP 800-90B** [SP800-90B] specifies continuous health tests for entropy sources. The repetition count test is aimed at "stuck" outputs, and the adaptive proportion test at a loss of entropy. Both are designed around a chosen false-positive probability.
- **AIS 31** [AIS31] distinguishes total-failure, start-up and online tests, the last tailored to a stochastic model of the source. I verified the 2011 reference document only through secondary sources and did not read it.
- **Hardware work** puts statistical tests on chip:
  - on-the-fly NIST tests [VRV12];
  - "canary numbers" [RYMV16]. These are a deliberately weaker copy of the noise that gives an early warning before the main output is compromised. In spirit this is the closest work to RNGLab's "detect before failure", but it operates at the entropy source, not at the application.
- **Ageing.** In a two-year field study of SRAM sources (16 boards), min-entropy *rose* on average [WSMG20]. Other sources found in the search argue that ageing can also reduce entropy; those were not verified for this list.
- **Implications for RNGLab:**
  - Health tests target catastrophic or gross failures by design. That they pass under a subtle defect is expected, not a finding.
  - H1 therefore needs a stronger conventional arm as well (Section 4).
  - H2 cannot assume that degradation is monotone.

### 2.3 MCMC diagnostics

- **Established tools:** multiple-chain R̂ ([GR92], refined in [VGSCB21]); initial-sequence estimators [Gey92]; batch-means and spectral long-run variance with consistency conditions [FJ10]. [Roy20] reviews these diagnostics and recent stopping rules.
- **The link to RNG quality:** [CDO11] proves consistency of MCMC driven by completely uniformly distributed sequences. As a by-product it covers the usual practice of using pseudo-random numbers.
- **These diagnostics target convergence and mixing, not RNG faults.** RNGLab uses them to rule out sampler failure before attributing an anomaly to the source. That is a role the reviewed literature does not take up.
- **Phase 2's finding is consistent with this literature's caution.** The AR(1)-derived 10·τ̂ batch rule failed on the refreshed PIMC sampler. The results of [FJ10] are asymptotic, and [Roy20] stresses the limits of diagnostics.

### 2.4 Prognostics and remaining useful life

- **Statistical, data-driven RUL estimation is mature.** It works from indirect condition data through degradation models and first-passage times [SWHZ11]. Prognostics has its own evaluation metrics, including ones that use the forecast's uncertainty, for offline evaluation against known end-of-life [SCSSG10].
- **No RUL or prognostics work on RNGs or entropy sources was found.** One extended search returned only mechanical and material RUL.
- **Implications for RNGLab:**
  - Stage D should adopt established RUL baselines and the metrics of [SCSSG10], not invent its own.
  - The identifiability result (`research/phase2-20261030/rus/`) is the condition that makes any forecast meaningful.

### 2.5 Restart, intervention and matched designs

- **Software ageing and rejuvenation** [HKKF95]: periodic restarts that delay failure in long-running systems. This is the closest conceptual precedent for H3's "intervention extends lifetime".
- **Universal restart strategies for Las Vegas algorithms** [LSZ93].
- **Recovery through reseeding:** robustness of PRNGs with input, including recovery after state compromise by fresh entropy [DPRVW13]. This is directly relevant to the H3 interventions "independent reseed" and "source switch". Its security notion differs from application lifetime, however.
- **Matched designs:** coupled Markov chains [JOA20] supply machinery for the identical-prefix branched design that H3 needs.
- **Sequential change detection** with calibrated alarm thresholds starts with CUSUM [Pag54]. It is relevant to the latency component of H1 (capability C5).
- **What was not found:** a matched experiment measuring how far an RNG intervention extends a Monte Carlo application's time to failure.

## 3. What RNGLab must not claim

- To be the first to show that application tests catch generators that pass standard tests ([FLW92], [Cod94], [VAK94]).
- That conventional online health tests "fail" in general. SP 800-90B tests target gross failures by design [SP800-90B].
- That entropy sources degrade monotonically with age ([WSMG20] reports the opposite for SRAM).
- That the centroid refresh or the virial estimator is new (the Phase II milestone report covers this prior work).
- Any statement about real generators that rests on synthetic adapter faults alone (`ROADMAP.md` §8).

## 4. Consequences for the Stage C (H1) protocol, to be decided when it is drafted

- **Conventional arm.** Use the in-repository SP 800-90B RCT/APT and the SP 800-22 subset, both calibrated to the same whole-run false-alarm rate as the application monitors. Adding TestU01 [LS07] (SmallCrush or Crush on windows of the consumed stream, version pinned) would answer the obvious referee objection that the comparison was against weak tests. This is open decision 3 in `ROADMAP.md`.
- **Real weak generators (C3).** Prefer generators whose application failures are documented, such as R250-type shift registers [FLW92], [VAK94] and small LCGs with lattice structure [Mar68]. They anchor H1 to known external evidence.
- **Sequential monitors (C5).** CUSUM-type rules [Pag54], with thresholds calibrated on healthy runs.
- **Reporting.** Report the regions where the battery detects first with the same prominence as the "pass but broken" regions.

## 5. Bibliographic notes on terms used

- "Conventional tests pass" always names the battery, stream window, sample size and false-alarm rate (`ROADMAP.md` §8).
- "Early warning" in [RYMV16] means detection before the generator's *output* is compromised, at the source. RNGLab's sense is detection before the *application* result fails. The two are not the same claim.

## 6. Limitations of this review

- **Search method.** Web search only, on 2026-11-06, in English. No Scopus, Web of Science or Google Scholar citation-graph sweep was made. Recall is therefore not exhaustive, and "not found" means "not found by these searches", not "does not exist".
- **Full texts.** Mostly not read. Content statements rest on abstracts and secondary summaries.
- **Adjacent field not reviewed:** silent-data-corruption detection in HPC, i.e. runtime statistical detection of hardware faults from application data. One search found HPC work on iterative and stencil codes but no Monte Carlo-specific method. It deserves a check before a paper is submitted.
- **Secondary verification only:** AIS 31 [AIS31] was verified only through secondary sources (slides and summaries).
- **Recommended before Stage F (paper):** a citation-graph search forward from [FLW92], [Cod94], [VAK94] and [RYMV16]. Its purpose is to catch any later calibrated comparison this search missed.

## References (verified 2026-11-06: authors, title, venue, year, DOI against publisher or indexing records)

- [AIS31] W. Killmann, W. Schindler. *A Proposal for: Functionality Classes for Random Number Generators*, Version 2.0, 18 September 2011. BSI (mathematical-technical reference for AIS 20/31). Verified via secondary sources (BSI/Schindler slides); full text not read.
- [CDO11] S. Chen, J. Dick, A. B. Owen. Consistency of Markov chain quasi-Monte Carlo on continuous state spaces. *Ann. Statist.* 39(2), 673–701, 2011. doi:10.1214/10-AOS831
- [CLK11] T. H. Click, A. Liu, G. A. Kaminski. Quality of random number generators significantly affects results of Monte Carlo simulations for organic and biological systems. *J. Comput. Chem.* 32(3), 513–524, 2011. doi:10.1002/jcc.21638
- [Cod94] P. D. Coddington. Analysis of random number generators using Monte Carlo simulation. *Int. J. Mod. Phys. C* 5(3), 1994. doi:10.1142/S0129183194000726 (volume, issue and DOI from the RePEc record; page range not verified). Preprint arXiv:cond-mat/9309017.
- [Cod96] P. D. Coddington. Tests of random number generators using Ising model simulations. *Int. J. Mod. Phys. C* 7(3), 295–303, 1996. doi:10.1142/S0129183196000235
- [DPRVW13] Y. Dodis, D. Pointcheval, S. Ruhault, D. Vergnaud, D. Wichs. Security analysis of pseudo-random number generators with input: /dev/random is not robust. *ACM CCS 2013*, 647–658. doi:10.1145/2508859.2516653
- [FJ10] J. M. Flegal, G. L. Jones. Batch means and spectral variance estimators in Markov chain Monte Carlo. *Ann. Statist.* 38(2), 1034–1070, 2010. doi:10.1214/09-AOS735 (checked in Phase 1)
- [FLW92] A. M. Ferrenberg, D. P. Landau, Y. J. Wong. Monte Carlo simulations: hidden errors from "good" random number generators. *Phys. Rev. Lett.* 69, 3382–3384, 1992. doi:10.1103/PhysRevLett.69.3382
- [Gey92] C. J. Geyer. Practical Markov chain Monte Carlo. *Statist. Sci.* 7(4), 473–483, 1992. doi:10.1214/ss/1177011137 (checked in Phase 1)
- [GR92] A. Gelman, D. B. Rubin. Inference from iterative simulation using multiple sequences. *Statist. Sci.* 7(4), 457–472, 1992. doi:10.1214/ss/1177011136
- [HKKF95] Y. Huang, C. Kintala, N. Kolettis, N. D. Fulton. Software rejuvenation: analysis, module and applications. *Proc. FTCS-25*, Pasadena, 381–390, 1995.
- [JOA20] P. E. Jacob, J. O'Leary, Y. F. Atchadé. Unbiased Markov chain Monte Carlo methods with couplings. *J. R. Stat. Soc. B* 82(3), 543–600, 2020. doi:10.1111/rssb.12336
- [LS07] P. L'Ecuyer, R. Simard. TestU01: a C library for empirical testing of random number generators. *ACM Trans. Math. Softw.* 33(4), Article 22, 2007. doi:10.1145/1268776.1268777
- [LSZ93] M. Luby, A. Sinclair, D. Zuckerman. Optimal speedup of Las Vegas algorithms. *Inf. Process. Lett.* 47(4), 173–180, 1993. doi:10.1016/0020-0190(93)90029-9
- [Mar68] G. Marsaglia. Random numbers fall mainly in the planes. *PNAS* 61(1), 25–28, 1968. doi:10.1073/pnas.61.1.25
- [MB04] S. Mertens, H. Bauke. Entropy of pseudo-random-number generators. *Phys. Rev. E* 69, 055702(R), 2004. doi:10.1103/PhysRevE.69.055702
- [Pag54] E. S. Page. Continuous inspection schemes. *Biometrika* 41(1/2), 100–115, 1954. doi:10.1093/biomet/41.1-2.100
- [Roy20] V. Roy. Convergence diagnostics for Markov chain Monte Carlo. *Annu. Rev. Stat. Appl.* 7, 387–412, 2020. (Preprint arXiv:1909.11827; the published page was not opened.)
- [RYMV16] V. Rožić, B. Yang, N. Mentens, I. Verbauwhede. Canary numbers: design for light-weight online testability of true random number generators. IACR ePrint 2016/386.
- [SCSSG10] A. Saxena, J. Celaya, B. Saha, S. Saha, K. Goebel. Metrics for offline evaluation of prognostic performance. *Int. J. Progn. Health Manag.* 1(1), 1–20, 2010. (The author order differs between the DOAJ and NASA records; NASA's order used.)
- [SMC03] A. Srinivasan, M. Mascagni, D. Ceperley. Testing parallel random number generators. *Parallel Comput.* 29(1), 69–94, 2003. doi:10.1016/S0167-8191(02)00163-1
- [SP800-90B] M. S. Turan, E. Barker, J. Kelsey, K. McKay, M. Baish, M. Boyle. *Recommendation for the Entropy Sources Used for Random Bit Generation*. NIST SP 800-90B, January 2018. doi:10.6028/NIST.SP.800-90B
- [SWHZ11] X.-S. Si, W. Wang, C.-H. Hu, D.-H. Zhou. Remaining useful life estimation – a review on the statistical data driven approaches. *Eur. J. Oper. Res.* 213(1), 1–14, 2011. doi:10.1016/j.ejor.2010.11.018
- [VAK94] I. Vattulainen, T. Ala-Nissila, K. Kankaala. Physical tests for random numbers in simulations. *Phys. Rev. Lett.* 73, 2513, 1994. doi:10.1103/PhysRevLett.73.2513
- [VGSCB21] A. Vehtari, A. Gelman, D. Simpson, B. Carpenter, P.-C. Bürkner. Rank-normalization, folding, and localization: an improved R̂. *Bayesian Anal.* 16(2), 667–718, 2021. doi:10.1214/20-BA1221 (checked in Phase 1)
- [VRV12] F. Veljković, V. Rožić, I. Verbauwhede. Low-cost implementations of on-the-fly tests for random number generators. *DATE 2012*. doi:10.1109/DATE.2012.6176635
- [WSMG20] R. Wang, G. Selimis, R. Maes, S. Goossens. Long-term continuous assessment of SRAM PUF and source of random numbers. *DATE 2020*. arXiv:2007.15909
