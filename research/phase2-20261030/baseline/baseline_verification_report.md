> Publication note: this is the completed local experiment report. Statements about files not yet pushed or an isolated patch describe that run before publication. The repository now includes the optional patch. Full raw traces are not distributed; see [data availability](../DATA_AVAILABILITY.md).

# Independent verification of the preserved Phase-I PIMC baseline

Date: 30 October 2026. Experiment: `phase2-baseline-independent-verification-v1`.

All 12 saved coverage counts are correct. The checker read 786,432 rows from the original compressed traces, without importing or executing the original simulator. It independently recomputed every batch interval, all chain means, median reported standard errors, Wilson coverage intervals and seed derivations. No original artifact was modified.

## Recovered configuration

The dimensionless harmonic-oscillator primitive-discretization action is

\[
S(x)=\frac{M}{2\beta}\sum_{j=0}^{M-1}(x_{j+1}-x_j)^2+
\frac{\beta}{2M}\sum_{j=0}^{M-1}x_j^2,\qquad x_M=x_0.
\]

Here M=16, β=2, and the proposal half-width is 0.4. Each sweep proposes **every bead once**, in fixed index order, with `x'_j=x_j+0.4(2u-1)`. Each proposal draws an additional uniform acceptance variate from the same replica's `random.Random` stream and accepts with probability min(1,exp(−ΔS)). Drawing the acceptance variate even for automatically accepted proposals makes the random stream schedule explicit. Symmetric local proposals and the action difference implement ordinary Metropolis updates. Fixed-order compositions preserve the invariant distribution, although the complete sweep need not be reversible.

Each of 64 replicas starts from independent uniform bead coordinates in [−0.5,0.5) and has 12,288 recorded sweeps. The short window contains sweeps 1,025–4,096 (3,072 retained); the long window contains sweeps 4,097–12,288 (8,192 retained). The window comparisons share trajectories and change both warmup and retained length.

The seed namespace is `rnglab-pimc-pilot-20261030-scalar-v1`. The seed for replica r is the first 128 bits of SHA-256 over `namespace/{r:03d}`, interpreted as an integer. All 64 derived seeds are distinct, match the CSV manifest and have distinct full TD/CV trajectories. This checks duplication and seed bookkeeping; it does not establish mathematical independence of PRNG streams.

The estimators are

\[
E_{TD}(x)=\frac{M}{2\beta}-\frac{M}{2\beta^2}\sum_j(x_{j+1}-x_j)^2+\frac 12\overline{x^2},
\]

\[
E_{CV}(x)=\frac{1}{2\beta}+\overline{x^2}-\frac 12\bar x^2.
\]

The center of a batch interval is the mean of a=⌊n/b⌋ non-overlapping block means. Its estimated variance is `sample_variance(block_means)/a`; its half-width is `t_(a−1,.975) × sqrt(estimated_variance)`. All window lengths are divisible by 32, 128 and 512, so no observations are discarded. Coverage is inclusive at both endpoints. The pilot permits four batches; the short b=512 row has six, below the general RNGLab summarizer's eight-batch minimum. Its label and interpretation must preserve that distinction.

## Independent Gaussian reference and interval arithmetic

Write S=½xᵀAx. The circulant precision eigenvalues are

\[
\lambda_k=\beta/M+(4M/\beta)\sin^2(\pi k/M),\quad k=0,\ldots,M-1.
\]

Thus E[mean(x²)]=M⁻¹Σ λₖ⁻¹ and E[bar(x)²]=1/β. Substitution into CV gives E[CV]=M⁻¹Σ λₖ⁻¹. The expected spring sum is Σ[4 sin²(πk/M)/λₖ]; substitution into TD produces the same value. The checker obtains **0.6554741522981732**, agreeing with the pilot's closed form to floating-point precision. This derivation uses the Gaussian precision eigenvalues rather than differentiating the pilot's partition-function expression.

The continuum target is 0.6565176427496657; the 0.00104349 finite-M difference is a discretization effect. It is not a Monte Carlo error and cannot explain the severe observed interval undercoverage.

The checker derives the Student-t critical values by direct integration of the normalized t density, using adaptive Simpson quadrature and bisection. This is independent of the pilot's incomplete-beta implementation. All six required quantiles agree to better than 10⁻⁹. Wilson intervals were also recomputed independently.

## Verified coverage and variance estimates

The variance ratio is the mean within-chain estimated variance of the mean divided by the empirical variance of the 64 chain means. It is descriptive and itself uncertain; it is not treated as a fixed true ratio.

| Window | Observable | Batch | Covered | Coverage | Wilson 95% interval | Estimated / empirical variance |
|---|---|---:|---:|---:|---:|---:|
| short | TD | 32 | 42/64 | 65.6% | 53.4–76.1% | 0.175 |
| short | TD | 128 | 49/64 | 76.6% | 64.9–85.3% | 0.330 |
| short | TD | 512 | 58/64 | 90.6% | 81.0–95.6% | 0.588 |
| short | CV | 32 | 24/64 | 37.5% | 26.7–49.7% | 0.063 |
| short | CV | 128 | 37/64 | 57.8% | 45.6–69.1% | 0.219 |
| short | CV | 512 | 53/64 | 82.8% | 71.8–90.1% | 0.572 |
| long | TD | 32 | 40/64 | 62.5% | 50.3–73.3% | 0.191 |
| long | TD | 128 | 51/64 | 79.7% | 68.3–87.7% | 0.347 |
| long | TD | 512 | 57/64 | 89.1% | 79.1–94.6% | 0.600 |
| long | CV | 32 | 20/64 | 31.3% | 21.2–43.4% | 0.049 |
| long | CV | 128 | 38/64 | 59.4% | 47.1–70.5% | 0.169 |
| long | CV | 512 | 48/64 | 75.0% | 63.2–84.0% | 0.429 |

The long CV mean is 0.652845976, bias −0.002628177. The empirical variance of chain means is 0.008754462; the mean reported variance is 0.000431840 at b=32, 0.001482912 at b=128 and 0.003754337 at b=512. The empirical standard deviation of chain means is much greater than these typical estimated standard errors. The mean bias is small relative to the across-replica uncertainty, so these data do not establish a persistent mean bias. They do establish that these particular reported uncertainty estimates are too small for the observed between-replica variability.

The mean acceptance rate is 0.695621 (range 0.693380–0.697805). High acceptance does not establish rapid mixing of the observable. Mean empirical long-window CV autocorrelations are 0.9884 at lag 1, 0.7645 at lag 32, 0.4228 at lag 128 and 0.0375 at lag 512. These are estimated after demeaning each finite trajectory and can underestimate very slow dependence. They support investigating dependence; they do not by themselves rule out initialization effects or establish an integrated autocorrelation time.

Among 2,016 pairs of long CV traces, the mean synchronized Pearson correlation is 0.00208, RMS 0.14189, and range −0.3264–0.5969. Large individual correlations are possible with strongly correlated time series and many pairs. No IID-based significance threshold was applied. This is a descriptive check, not evidence proving stream independence or a finding that the streams are dependent.

## Corrections and limitations

The pilot JSON wording “one bead per sweep” is inaccurate: the actual implementation updates each of 16 beads once per sweep. The methodology note describes this correctly. The four-batch minimum is an explicit standalone-pilot choice, not the general RNGLab eight-batch rule. Neither issue changes the long-window counts.

No reference, energy-formula or interval-arithmetic defect was found. These checks do not verify the original simulator's full trajectory against another implementation, establish equilibrium, or provide a causal identification of the undercoverage mechanism. Independent Gaussian sampling, exact stationary starts, alternative transition kernels and newly seeded controlled experiments are needed for that step. Replaying the same data or seeds does not count as independent replication.

Run `python verify_baseline.py` from any directory. The checker resolves the preserved artifacts relative to itself. `baseline_verification.json` includes the original artifact hashes, recovered settings, seeds, numerical checks, ACF summaries and all variance ratios. `recomputed_intervals.csv` contains the 768 individually recomputed intervals. The originals remain unchanged.
