> Publication note: this is the completed local experiment report. Statements about files not yet pushed or an isolated patch describe that run before publication. The repository now includes the optional patch. Full raw traces are not distributed; see [data availability](../DATA_AVAILABILITY.md).

# RNGLab Phase II: mathematical reference and method review

Prepared 2026-10-30. This is an independent analytic review of the primitive-action harmonic model in the isolated RNGLab snapshot at commit `dd96adba23b8100aba9147d71fffc63932fc6b0d`. It supplies mathematical checks and interpretation for the new experiments; it does not claim the experiments have a particular outcome.

## 1. The finite-M target and reference are correct

In units m = ℏ = ω = 1, with periodic bead indices, the target is

\[
 \pi_M(x)\propto e^{-S_M(x)},\qquad
 S_M(x)=\frac{M}{2\beta}\sum_{j=0}^{M-1}(x_{j+1}-x_j)^2+
             \frac{\beta}{2M}\sum_{j=0}^{M-1}x_j^2.
\]

Let D be the periodic first-difference matrix and L = DᵀD. Its Gaussian precision is

\[
 Q=\frac{\beta}{M}I+\frac{M}{\beta}L,
 \qquad q_k=\frac{\beta}{M}+\frac{4M}{\beta}\sin^2(\pi k/M).
\]

All qₖ are positive for β > 0. The exact path distribution is N(0,Q⁻¹). An independent real orthogonal diagonalization followed by normals with variances 1/qₖ, or a Cholesky factor of Q⁻¹, is a valid IID reference sampler under the assumed normal source. A complex FFT must enforce conjugate symmetry and the correct real-mode normalization; drawing arbitrary complex independent coefficients and taking a real part is generally wrong.

Including the primitive free-particle normalization gives

\[
 Z_M=\left[2\sinh\{M\operatorname{asinh}(\beta/(2M))\}\right]^{-1},
 \qquad E_M=-\partial_\beta\log Z_M
 =\frac{1}{2\sqrt{1+(\beta/(2M))^2}}
  \coth\{M\operatorname{asinh}(\beta/(2M))\}.
\]

At M = 16 and β = 2,

| Quantity | Value |
|---|---:|
| Finite-M energy Eₘ | 0.6554741522981734 |
| Continuum energy ½ coth(β/2) | 0.6565176427496656 |
| Discretization bias Eₘ − E∞ | −0.0010434904514922 |

The current `rnglab14/models/pimc_harmonic.py` computes this finite-M reference. Its uncertainty intervals use the finite-M value. Its separate `PIMCResult.energy_bias` uses the continuum value, so these two bias quantities should be named distinctly in experimental tables. That naming issue does not explain coverage against the finite-M reference.

The supplied `harmonic_reference_checks.py` independently verifies the closed form by eigenmode sums, the numerical derivative of log Z, and Q Q⁻¹ = I on 18 combinations of M and β, including M = 1 and 2 where careless periodic-edge handling is easy to get wrong. These deterministic checks were executed and passed; the saved values are in `harmonic_reference_values.json`.

## 2. Both energy formulas are correct; they probe different dynamics

The thermodynamic (TD) estimator follows from differentiating both the primitive normalization and action:

\[
 E_{TD}(x)=\frac{M}{2\beta}-\frac{M}{2\beta^2}\sum_j(x_{j+1}-x_j)^2
                  +\frac{1}{2M}\sum_jx_j^2.
\]

The centroid virial (CV) estimator is

\[
 E_{CV}(x)=\frac{1}{2\beta}+\frac{1}{M}\sum_jx_j^2-\frac12\bar x^2.
\]

For a direct finite-dimensional proof, write x = c1 + y, with c = x̄ and Σyⱼ = 0. Then

\[
 S_M=\frac\beta2c^2+\frac12 y^TQy,
 \quad c\sim N(0,1/\beta),\quad c\perp y,
 \quad r_g^2=M^{-1}\sum_jy_j^2,
 \quad E_{CV}=\frac{1}{2\beta}+\frac12c^2+r_g^2.
\]

Integration by parts in the M−1 dimensional shape subspace gives
\(\mathbb{E}[2S_{spring}+\beta r_g^2]=M-1\). Substitution in the TD expectation gives exactly the CV expectation. The latter also follows directly from the precision eigenvalues. This establishes the equality for the finite-M distribution, without appealing to the continuum virial theorem. Yamamoto's fluctuation-coordinate derivation supplies the wider path-integral context for centroid virial estimators ([Yamamoto, 2005](https://arxiv.org/abs/physics/0505109)).

For a zero-mean Gaussian x with covariance Σ and a symmetric quadratic observable a+xᵀAx,
\(\mathbb{E}[a+x^TAx]=a+\operatorname{tr}(A\Sigma)\) and
\(\operatorname{Var}(x^TAx)=2\operatorname{tr}[(A\Sigma)^2]\).
For TD, the eigenmode coefficient is
\(a_k=1/(2M)-M[4\sin^2(\pi k/M)]/(2\beta^2)\).
For CV it is b₀ = 1/(2M), bₖ = 1/M for k ≠ 0. Thus the exact marginal moments at M = 16, β = 2 are:

| Quantity | Exact value |
|---|---:|
| E[TD] | 0.6554741522981731 |
| E[CV] | 0.6554741522981733 |
| Var(TD) | 1.854157292369273 |
| Var(CV) | 0.13463144466744617 |
| Cov(TD,CV) | 0.056894368518359564 |
| E[c²] | 0.5 |
| E[r_g²] | 0.1554741522981732 |
| Var(½c²) | 0.125 |
| Var(r_g²) | 0.00963144466744617 |
| Centroid fraction of marginal CV variance | 92.8460660202846% |

The large centroid fraction makes centroid dynamics a strong mechanism candidate for CV uncertainty failure. This is a prediction from the model, not yet a result about the sampler. TD's larger raw variance can coexist with a shorter integrated autocorrelation time: smaller marginal variance alone does not imply better Monte Carlo efficiency.

The marginal CV observable is a weighted sum of independent squared normals plus a constant. It is not itself normal (skewness ≈ 2.56388). The IID mean of 8192 observations has skewness ≈ 0.02833, which suggests a normal approximation is reasonable but does not guarantee exact 95% t coverage. The actual IID coverage control remains necessary.

## 3. Exact centroid refresh is a valid, sharply targeted intervention

After a local sweep, draw c′ independently from N(0,1/β) and shift each bead by c′ − x̄. The springs and y do not change. Because the target factorizes as π(c,y) = π(c)π(y), this is an exact Gibbs update of the centroid and preserves π. The identity

\[
 S_M(x+(c'-c)1)-S_M(x)=\frac\beta2[(c')^2-c^2]
\]

was checked deterministically in the supplied reference script. It requires no Metropolis rejection for the harmonic model. Replacing the centroid with a normal draw of variance M/β would be wrong: that is the variance of the normalized Fourier coefficient √M c, not of c.

An exact Gaussian start followed by target-preserving local sweeps is stationary from the first recorded state under the ideal PRNG model. Therefore persistence of undercoverage under exact starts rejects initialization transient as a necessary cause. If centroid refresh then corrects the variance ratio and coverage, it provides much stronger mechanism evidence than a correlation between acceptance rate and coverage. It does not validate local Metropolis for other models or prove every variance estimator adequate.

The 92.85% number concerns same-time variance. Long-run variance of the CV mean includes cross-lag centroid-shape covariances, since local Metropolis decisions can dynamically couple c and y despite their independence at one equilibrium time. Across independent replica means, compute

\[
 \operatorname{Var}(\bar E_{CV})=
  \operatorname{Var}(\overline{c^2/2})+
  \operatorname{Var}(\overline{r_g^2})+
  2\operatorname{Cov}(\overline{c^2/2},\overline{r_g^2}).
\]

Even after exact centroid refresh, a refreshed c at time t can influence shape at t+1. Do not automatically drop cross-lag terms. Moustafa and Schultz explicitly study harmonic-coordinate collective sampling, normal modes, and special treatment of the external centroid mode; these are prior methods, so a centroid or normal-mode improvement is not a new algorithmic invention ([Moustafa and Schultz, 2024](https://arxiv.org/html/2404.12551v3), [published article](https://doi.org/10.1021/acs.jctc.4c00522)).

## 4. What batch means and Geyer do, and where their assumptions matter

For a stationary scalar trace with autocovariance γₖ, the finite-length variance is

\[
 \operatorname{Var}(\bar X_n)=\frac1n\left[\gamma_0+
    2\sum_{k=1}^{n-1}(1-k/n)\gamma_k\right].
\]

Divide n = ab observations into a batch means Bᵢ of length b. The usual estimate is s²_B/a. If the batch means were IID normal, the t interval with a−1 degrees of freedom would be exact. Short, mutually dependent batches do not satisfy that condition. In a stationary process with Cₕ = Cov(B₀,Bₕ), the exact bias is

\[
 \mathbb{E}[s_B^2/a]-\operatorname{Var}(\bar X_n)
 =-\frac{2}{a(a-1)}\sum_{h=1}^{a-1}(a-h)C_h.
\]

Consequently, positive correlations between batches drive the variance estimate downward even without burn-in bias. Increasing the number of short batches cannot by itself make them independent. As n grows with b fixed, the estimator approaches the tapered, finite-lag variance sum
\(\gamma_0+2\sum_{k=1}^{b-1}(1-k/b)\gamma_k\), divided by n, rather than automatically estimating the full long-run variance. Consistency arguments require suitable growth of batch length and batch count and further dependence/moment conditions. Overlapping batch means can reduce variability but does not remove inadequacy of a short batch length. Flegal and Jones establish sufficient conditions for batch means and spectral variance estimators; their results are asymptotic, not guarantees for a particular 8192-step run ([Flegal and Jones, 2010](https://doi.org/10.1214/09-AOS735), [author preprint](https://arxiv.org/abs/0811.1729)).

**Systematic scan caveat.** Let Kⱼ be the reversible one-bead Metropolis kernel. The full sweep P = K₁⋯Kₘ preserves π, but P* = Kₘ⋯K₁, so the full-sweep chain is generally nonreversible. Geyer's positive/monotone adjacent-autocovariance-pair theorem assumes reversibility. It cannot be invoked merely because each bead update is a Metropolis step ([Geyer, 1992](https://doi.org/10.1214/ss/1177011137); [author's method documentation](https://www.stat.umn.edu/geyer/mcmc/library/mcmc/html/initseq.html)). This qualification does not itself explain the Phase I failure, because Phase I reported batch means, not Geyer intervals.

Use IPS/IMS on systematic-sweep traces as a sensitivity analysis with an explicit assumption caveat. A 50/50 mixture of complete forward and reverse sweeps is reversible, as is a uniform random-scan kernel (or powers of it). A uniformly selected full permutation is also reversible when each permutation and its reverse have equal probability. Composing a reversible sweep with a separate centroid refresh still need not be reversible; a symmetric composition or a mixture with its adjoint is needed to claim that property. A Bartlett/HAC estimate or batch means need not impose Geyer's specific reversibility shape restriction, but their own truncation, stationarity and large-run assumptions still matter.

An MCMC CLT needs conditions beyond preservation of the target. Finite analytic Gaussian moments establish observable integrability, while drift/mixing properties of the actual transition kernel determine a CLT and its variance; an observed ESS value does not prove those conditions. Jones reviews these distinct sufficient conditions ([Jones, 2004](https://arxiv.org/abs/math/0409112)).

## 5. What the controlled experiment can establish

The most informative comparisons have exact stationary starts and independent new seeds. IID Gaussian paths apply the same CV/TD and interval code while removing time dependence. Exact-start local Metropolis retains equilibrium dependence while removing initialization transient. Exact-start local Metropolis plus centroid refresh targets one mode without changing the equilibrium law. Fixed retained length with different initialization and warmup should be a separate comparison, so warmup and sampling effort are not confounded.

Report mean bias with a standard error from independent replica means, empirical variance of those means, mean estimated MC variance, their ratio, and binomial intervals for coverage. Pooling within-chain observations is not an independent replication count. For the variance ratio, the denominator is itself estimated and uncertainty can be substantial. Replica bootstrap intervals can assess that uncertainty. Coverage indicators across lengths of the same chain are paired; comparisons should retain that pairing or avoid an independence claim.

With 512 independent replicas, the binomial standard error at 95% coverage is about 0.96 percentage point and the approximate two-sided 95% precision is ±1.89 percentage points. This can distinguish severe 75% undercoverage from 95%, but cannot certify an exact nominal target. A predefined acceptance criterion and its uncertainty interval should replace declaring success because 95% happens to lie in one confidence interval. Exploratory selection of the best batch size on the same final replicas would bias a claimed validation.

Separate acceptance and proposal channels for fault comparisons. The Phase I scalar pilot consumes an acceptance word on every bead proposal, while upstream consumes it only for uphill moves. Both have the same healthy target law, but their stream histories differ, which matters for periodic reuse and fault timing. Equal arithmetic workloads alone do not guarantee equal consumed stream positions. Do not infer stochastic independence from unique SHA-derived seeds; it is a modeling assumption of the healthy PRNG experiment.

## 6. Scope of possible research claims

Correct finite-M formulas, low marginal CV variance, and 100% acceptance exact harmonic sampling are established mathematics. Established MCMC theory already explains why short batches can underestimate long-run variance. RNGLab can contribute a carefully controlled, reproducible causal diagnosis in its own pipeline, including the empirical size of the problem, which modes explain it, and how a frozen correction behaves on independent holdout chains. That is a defensible validation contribution; it is not a new general theorem about all PIMC or RNGs.

A healthy PRNG's local-sampler failure cannot support RNG degradation claims. If the repaired healthy workload is calibrated, a specified proposal-only or acceptance-only source fault can be compared against it. Corrupting normal inputs to the exact Gaussian sampler instead studies failure of that transform's source assumptions; it should not be described as a Metropolis acceptance fault. Fault study conclusions should name the adapter, channel, strength, onset, measured diagnostic, and application endpoint. A monitor needs separate healthy calibration and final whole-run false-alarm evaluation.

The exact harmonic benchmark is useful for isolating mechanisms but is too easy to justify a general forecast claim. Alarm-conditioned Remaining Useful Samples also requires an assumption on future degradation. Two fault processes can have identical observable histories through alarm time and then respectively recover or worsen. Their RUS distributions differ despite identical available features. Without a future-path model or a robust worst-case set of paths, RUS is not identifiable from the observed prefix alone.

## Files and execution

- `harmonic_reference_checks.py`: independent standard-library deterministic formulas and assertions; executed successfully on 18 model settings.
- `harmonic_reference_values.json`: saved exact moments, eigenvalues, and validation residuals.
- Inspected source: `src/rnglab14/models/pimc_harmonic.py`, `src/rnglab12/uncertainty.py`, and the Phase I `outputs/pimc_uncertainty_pilot.py`.
- This review does not change the root experiment implementation or upstream repository.
