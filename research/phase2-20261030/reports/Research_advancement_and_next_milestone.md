> Publication note: this is the completed local experiment report. Statements about files not yet pushed or an isolated patch describe that run before publication. The repository now includes the optional patch. Full raw traces are not distributed; see [data availability](../DATA_AVAILABILITY.md).

# Research advancement and next executable milestone

30 October 2026. This assessment uses the completed Phase II computational evidence. The integration milestone at the end is **proposed and unexecuted**.

## Recommendation

Choose **Direction A, uncertainty calibration under Markov-chain dependence**, as the primary near-term BSc and TDK contribution. Use the controlled random-channel study as a bounded Direction B extension. Direction C should remain an identifiability and protocol result until a defensible model of future degradation is supplied.

A suitable focused thesis question is: **How can sampler-induced uncertainty failure be separated from random-source faults in finite-bead harmonic PIMC?** This keeps a direct connection to the official random-number-generator topic while building the argument on a validated healthy baseline. A possible working title is *Uncertainty calibration and random-channel fault diagnosis in harmonic path-integral Monte Carlo*.

The strongest result is a causal diagnosis with explicit limits. The finite-M reference and energy formulas survived independent checks. Exact stationary starts still exhibited severe local-Metropolis undercoverage. A target-preserving centroid refresh reduced the dominant variance component and passed a frozen independent healthy-control gate. The analytically solvable dependence control reproduces the mechanism with known finite-n variance and shows why a mean MC-variance ratio alone does not determine coverage. This is enough to support a coherent thesis investigation; it is not a new PIMC estimator or universal convergence theorem. See the [root-cause report](PIMC_root_cause_report.md) and [analytic control](../analytic-dependence/analytic_dependence_report.md).

## What changed scientifically

Phase I established reproducible undercoverage without identifying its cause. Phase II distinguishes several explanations: finite-M discretization and formula errors do not explain the recorded failure; initialization can introduce a real transient but is not necessary; equilibrium dependence, the centroid contribution, and inadequate finite-run variance estimation are directly implicated by targeted controls. Remaining uncertainty concerns the full distribution and long-run covariance of the systematic Metropolis estimator, including the joint behavior of mean errors and MCSE estimates.

The newly frozen healthy confirmation covers 811/1024 trajectories with local Metropolis, 964/1024 with centroid refresh, and 969/1024 with independent Gaussian paths under the same primary CV BM512 rule. The two required healthy controls meet all locally frozen coverage, variance-ratio and bias tolerances. This supports calibration in the specific tested domain: M=16, β=2, δ=0.4, exact equilibrium starts, 8192 retained sweeps and BM512. It does not transfer automatically to other models, starts, run lengths or uncertainty methods.

The controlled adapter study now provides a separate Direction B result. Under a newly validated healthy control, acceptance bias and proposal–acceptance coupling produce additional energy and coverage failures. The application shape monitor responds very differently to the two faults. The marginal source-mean monitor detects acceptance bias but often misses channel coupling; the explicit cross-channel monitor detects coupling. That is evidence about those specified diagnostics and channel contracts. It does not establish that application monitoring is generally superior to source diagnostics, nor does a two-fault study establish broad fault-family generalization. The [controlled fault report](Controlled_RNG_fault_report.md) gives the paired effects, uncertainty, calibration and failure timing.

The RUS work provides a useful negative result: identical observed histories can have different future degradation paths and hence different alarm-conditioned remaining lifetimes. A forecast requires a declared continuation model or an explicit set of allowed future paths. More features cannot recover a future branch absent from the observed information. The completed counterexample and protocol are in [RUS identifiability and protocol](../rus/RUS_identifiability_and_protocol.md).

## Comparison of research directions

| Criterion | Direction A: uncertainty calibration | Direction B: application-aware fault diagnosis | Direction C: alarm-conditioned RUS |
|---|---|---|---|
| Scientific novelty. | The underlying estimators and centroid updates are established. A carefully bounded diagnosis and independently validated implementation can be a reproducibility contribution. | Application-based RNG testing is established. The specific channel contracts, paired attribution and diagnostic disagreement may support a narrow contribution. | The alarm-conditioned target is ambitious, but novelty and achievable prediction guarantees remain to be established against prior survival and reliability work. |
| Existing evidence. | Independent references, stationary controls, an intervention, a frozen confirmation and an exact dependence control provide the strongest evidence. | The completed comparison establishes two strong adapter effects after a healthy gate. It does not yet establish performance on weak or unseen faults. | A reproducible identifiability counterexample is supported. Successful general forecasting is not established. |
| Experimental feasibility. | A local CPU suffices. The next unresolved step is a bounded integration validation through the actual repository API. | A local CPU suffices for the implemented channel faults. Broader fault strength and family generalization would need another fixed study. | A small model is feasible locally, but a defensible forecasting study needs explicit future-path families, censoring and alarm issuance rules. |
| Statistical defensibility. | The present study separates exploration and new validation, uses whole trajectories as units, and reports coverage and variance uncertainty. | Pairing and whole-run healthy calibration support the measured comparisons. Actual false-alarm rates differ and must remain visible. | A forecasting claim is not defensible without future-path assumptions and a held-out alarm-conditioned evaluation. |
| Technical complexity. | Complexity is moderate and most mathematical and implementation risks are now identified. | Complexity is moderate to high because source consumption, channel dependence, fault onset and monitor calibration must all agree. | Complexity is high because prediction, censoring, selection at alarm time, abstention and identifiability interact. |
| BSc suitability. | This is the strongest primary thesis scope because it has a precise question and completed controlled evidence. | This is a useful extension that preserves the connection to random-number sources. It should remain limited to the demonstrated mechanisms. | This is a weaker primary thesis choice at present. A clearly presented impossibility result can still be included. |
| TDK potential. | A closed causal investigation with independent reproduction and a reviewable patch is a credible presentation. | The contrasting diagnostic responses can strengthen the presentation when their limits and whole-run false alarms are shown. | A rigorous identifiability argument is valuable, but it should not be presented as a successful predictor. |
| Publication potential. | A methods-validation or reproducibility study may be possible after independent API validation. Established centroid sampling alone is not a publication novelty claim. | A focused contribution may emerge if channel-contract failures and diagnostic performance survive independent replication and stronger controls. Two large faults alone are a limited foundation. | Publication would need either a new formal result or a genuinely validated forecasting method under explicit assumptions. Neither follows from ambition alone. |
| HPC and physical equipment. | Neither HPC nor a detector is required for the proposed milestone. | Neither is required for the current computational faults. Physical-source claims require a separate acquisition and entropy-assessment protocol. | Neither resolves identifiability. Hardware access should follow a validated computational question. |

These judgments concern suitability and evidence, not a guarantee of thesis acceptance, TDK placement or publication.

## Closest prior work and the bounded contribution

Flegal and Jones establish batch-means and spectral-variance consistency under stated assumptions and examine finite-sample behavior. RNGLab does not discover that short batches can underestimate long-run variance; it diagnoses and quantifies that failure in its own saved and new experiments, with explicit stationary and intervention controls ([Flegal and Jones, 2010](https://doi.org/10.1214/09-AOS735)).

Yamamoto derives virial estimators by fluctuation-coordinate scaling, and Moustafa and Schultz study harmonic-coordinate collective sampling, including the centroid. RNGLab's optional harmonic heat bath implements a known conditional Gaussian idea. Its evidence contributes implementation, source-accounting and uncertainty validation for the identified pipeline, rather than a new sampling principle ([Yamamoto, 2005](https://arxiv.org/abs/physics/0505109), [Moustafa and Schultz, 2024](https://doi.org/10.1021/acs.jctc.4c00522)).

Ferrenberg, Landau and Wong show simulation errors arising from subtle RNG correlations in a cluster algorithm; Coddington explicitly uses Ising simulations as RNG tests. Therefore the broad statement that a generator can pass some diagnostics yet fail in an application is established prior research. A defensible RNGLab distinction is the controlled separation of healthy sampler uncertainty failure from two explicit random-channel contract violations, with paired effects and separately calibrated marginal, joint-channel and application monitors ([Ferrenberg, Landau and Wong, 1992](https://doi.org/10.1103/PhysRevLett.69.3382), [Coddington, 1996](https://surface.syr.edu/npac/34/)).

The present evidence supports statements about the tested harmonic workload and injected adapters. It supports neither universal RNG certification nor universal application-monitor superiority. The physical detector inventories establish file identity and structure; their entropy properties remain a separate research question.

## What to present to Dr. Herendi

Present the closed argument in this order:

1. Show the independently verified original undercoverage counts and the correct finite-M target.
2. Show that failure persists from exact equilibrium starts, then show the centroid variance decomposition and target-preserving intervention.
3. Show the frozen confirmation with coverage, variance-ratio and bias intervals, including the difference between standalone numerical validation and the repository API patch.
4. Show the two channel-fault effects and the diagnostic disagreement, together with each monitor's actual healthy whole-run false-alarm rate.
5. Propose the single integration milestone below as the completion criterion for the thesis's validated computational core. Keep the detector extension and unconstrained RUS prediction outside that milestone.

The scientific claim to defend is that RNGLab now separates a demonstrated healthy-sampler uncertainty failure from additional controlled channel faults in one exact-reference workload. The next question is whether that validated behavior survives integration through the package interface used by future research.

## Proposed next milestone: independent validation through the patched repository API

**Status: proposed and unexecuted.** The Phase II coverage campaign used the independently checked standalone NumPy/Numba implementation. The optional upstream patch has mathematical and regression validation, but has not yet been the implementation used in an independent R=1024 coverage campaign. That distinction is the narrow remaining integration question.

The current patched `run_simulation` API always replaces the initial path with uniform coordinates and calls uncertainty analysis with its default batch sizes 32 and 128. Setting `model.path` before calling it would not produce an exact stationary start. Thus the milestone needs a minimal, reviewed API extension before its protocol can be frozen:

- Add an explicit optional exact initial-path argument, validate shape and finiteness before sampling, copy its contents, and record its provenance. Keep default initialization and default source consumption unchanged when this argument is absent.
- Expose the uncertainty batch sizes or an audited raw-trace output so BM512 uses all 8192 retained observations through the actual repository uncertainty machinery. Preserve full precision in saved estimates and traces.
- Add regression cases for these explicit options, default behavior, source counters and recorded source digests. Run the complete suite in the pinned environment. Do not infer coverage from regression success.

After those changes pass review and tests, freeze a new locally recorded protocol and code/patch hashes. Use M=16, β=2, δ=0.4, zero burn-in from independently generated exact finite-M Gaussian starts, and n=8192 retained sweeps. Run R=1024 whole-trajectory replicas for the local sampler and for the local sampler with the explicit harmonic centroid source. Generate a fresh R=1024 IID Gaussian-path control using the package reference sampler and apply the same CV BM512 analysis. The primary interval is CV BM512 with 16 batches and t₁₅. TD and other variance methods remain secondary.

Use a previously unused namespace such as `rnglab-phase2-upstream-validation-v1/<condition>/<role>/<replica>` with separate initialization, proposal, acceptance, centroid and IID-reference roles. If local/refreshed conditions deliberately share seeds for pairing, record that design. All replicas remain independent units under the declared healthy PRNG assumption. Keep this campaign separate from both Phase II holdout and fault streams.

The upstream source interface consumes acceptance words only for uphill proposals. The standalone study consumed one acceptance word for every proposal. The upstream centroid helper uses an audited 16-byte Box–Muller draw, whereas the standalone generator uses its own normal transform. These differences prevent bitwise matching across the two implementations. Record source bytes, channel manifests and digests; validate distributional targets, numerical formulas and aggregate coverage rather than expecting identical trajectories.

Predefine the same primary healthy gates as Phase II. For **both** the refreshed sampler and the IID control, require the Clopper–Pearson 95% coverage interval wholly inside [0.92,0.98], the replica bootstrap 95% interval for estimated/empirical variance wholly inside [0.75,1.25], and the across-replica mean-bias 95% interval wholly inside [−0.01,0.01]. Retain a separately reported paired refreshed-minus-local coverage comparison with the same prespecified improvement margin; failure of any required gate is a failed integration validation, not a reason to choose another batch on those data.

A small separate smoke namespace may be used to measure runtime and check artifact plumbing before freezing the final protocol; those trajectories must not count toward the final campaign. Preserve the entire final traces, per-replica endpoint table, independent seed manifest, software lock, code/patch hashes and checksums. Independently recount all primary coverage events from the saved traces.

The milestone is complete when the actual patched package API meets the frozen gates and the separate recount agrees, or when a discrepancy is localized and documented with a minimal reproducer. It needs local CPU time and the existing pinned environment. It does not depend on HPC allocation, new instruments, physical entropy claims or a larger forecasting system.
