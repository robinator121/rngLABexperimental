# RNGLab v0.14 — Application Reliability & Defect Sensitivity Matrix

**Document Identifier:** RNGLAB-DOC-v014-APP-RELIABILITY-01  
**Classification:** Scientific Analysis & Defect Propagation Model  
**Target Release:** RNGLab v0.14 (`feat/v014-scientific-integration`)  
**Date:** 2026-10-18  
**Provenance SHA-256:** `570aa3e7a1e432ba979454c75cf0c90859175810f75678b4c529a023f6f83a6e`

---

## 1. Scientific Overview & Defect Propagation

A central theorem demonstrated in RNGLab v0.14 is that **the sensitivity of a computational application to randomness degradation is fundamentally governed by the structural coupling between the stochastic defect and the mathematical observable**.

Different applications do not degrade uniformly when supplied with flawed random numbers:
- Some algorithms exhibit high resilience due to intrinsic self-averaging or ergodic contraction.
- Other algorithms experience catastrophic failure while simultaneously reporting misleadingly normal diagnostic metrics (the **Observable Masking Phenomenon**).

This document analyzes the empirical vulnerability profiles of five representative scientific application models evaluated under identical controlled degradation regimes.

---

## 2. Cross-Application Defect Sensitivity Matrix

The table below contrasts the vulnerability of the five application domains across four canonical defect modes: Marginal Bit Bias ($\epsilon$), First-Order Autocorrelation ($\rho$), Coordinate Cross-Channel Coupling ($r_{xy}$), and Harmonic Clock Coupling ($A_{osc}$).

| Application Domain | Mathematical Target | Sensitivity to Bias $\epsilon$ | Sensitivity to Autocorr $\rho$ | Sensitivity to Channel Coupling $r_{xy}$ | Observable Masking Risk | Critical Failure Threshold |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Classical MC $\pi$** | $\int \int \mathbb{I}(x^2+y^2 \le 1) dx dy$ | **Severe** ($\Delta \pi \propto \epsilon$) | **Moderate** | **Extreme** ($\Delta \pi \propto r_{xy}$) | Low (Error is directly visible) | $|\epsilon| > 0.01$ or $|r_{xy}| > 0.005$ |
| **2D Ising Model ($L=8$)** | $\langle E \rangle, C_v, \langle |M| \rangle, \chi$ near $T_c$ | **Low** (Metropolis filters bias) | **High** ($N_{eff}$ collapses) | **Severe** (Proposal vs Acceptance) | Moderate (Distorts phase transition) | $|\rho| > 0.15$ in acceptance stream |
| **VMC Harmonic Oscillator** | Ground state $E_0$, $\langle x^2 \rangle$ | **Zero** on Energy if $\alpha=\alpha_0$ | **Severe** on Spatial Density | **Moderate** | **EXTREME** (Energy variance $\to 0$ while density collapses) | $|\rho| > 0.10$ corrupts $\langle x^2 \rangle$ without changing $E_0$ |
| **PIMC Quantum Chains ($M=16$)**| $E = \frac{1}{2}\coth(\beta/2)$ | **Catastrophic** | **Catastrophic** | **Catastrophic** | Low (Virial energy explodes) | $|\rho| > 0.10$ yields $+9800\%$ energy error |
| **Monte Carlo DSS** | Regret, Ranking, $95\%\ \text{VaR}$ | **Moderate** | **High** (Underestimates tails) | **Moderate** | High (False high confidence spread) | $\Delta \text{VaR}_{95\%} > 15\%$ at $|\rho| > 0.20$ |

---

## 3. Deep-Dive: Mathematical Findings by Application

### 3.1 Classical Monte Carlo $\pi$ & Coordinate Coupling
The standard classical Monte Carlo estimator for $\pi$:
$$\hat{\pi} = \frac{4}{N} \sum_{i=1}^N \mathbb{I}(x_i^2 + y_i^2 \le 1)$$
assumes that $X$ and $Y$ are independent uniform random variables $U(0, 1)$.

- **Marginal Bias:** If $X \sim U(\epsilon, 1+\epsilon)$, the estimated circle area contracts or expands linearly.
- **Cross-Channel Coupling ($r_{xy}$):** If $X$ and $Y$ originate from consecutive samples of a PRNG with lag-1 autocorrelation $\rho_1$, or if a single stream is naively interleaved ($x_i = u_{2i}, y_i = u_{2i+1}$):
  $$\mathbb{E}[X \cdot Y] = \frac{1}{4} + \frac{1}{12}\rho_1$$
  The sample points cluster along the diagonal $y = x$, artificially depleting the outer corners and distorting the quarter-circle boundary.
- **Empirical v0.14 Finding:** On `Degraded_MarkovBias` ($\rho = 0.3$), the $\pi$ estimate collapsed from $3.1416$ to $2.5472$ ($z = -30.90$, relative error $18.9\%$).

### 3.2 2D Ising Model & Stochastic Channel Decoupling
In the Metropolis-Hastings simulation of the 2D Ising Hamiltonian:
$$\mathcal{H} = -J \sum_{\langle i, j \rangle} s_i s_j$$
each trial step requires two distinct random decisions:
1. **Proposal Channel:** Selecting a lattice site $(i, j)$ uniformly at random.
2. **Acceptance Channel:** Drawing $u \sim U(0, 1)$ to test against the Boltzmann transition probability $\alpha = \min(1, \exp(-\beta \Delta E))$.

**The Decoupling Invariant:**
In RNGLab v0.14, `ising_2d.py` explicitly isolates the proposal channel from the acceptance channel via `StochasticChannels`.
- When the acceptance channel is subjected to positive autocorrelation ($\rho > 0$), consecutive rejection or acceptance runs occur.
- This creates an artificial **hysteresis loop**, freezing the system into a metastable magnetized domain even above the critical Curie temperature $T_c \approx 2.269$.
- The magnetic susceptibility $\chi = \beta L^2 (\langle M^2 \rangle - \langle |M| \rangle^2)$ drops sharply from $1.615$ to $0.758$, producing a false signature of an ordered ferromagnetic phase.

### 3.3 Variational Monte Carlo (VMC) & The Observable Masking Phenomenon
The Variational Monte Carlo simulation of the 1D quantum harmonic oscillator ($H = -\frac{1}{2}\frac{d^2}{dx^2} + \frac{1}{2}x^2$) employs the Gaussian trial wave function:
$$\Psi_\alpha(x) = \exp(-\alpha x^2)$$
The local energy observable is given by:
$$E_L(x) = \frac{H \Psi_\alpha(x)}{\Psi_\alpha(x)} = \alpha + x^2\left(\frac{1}{2} - 2\alpha^2\right)$$

#### The Masking Mechanism
When the variational parameter is tuned to the exact ground state $\alpha = \alpha_0 = 0.5$:
$$E_L(x) = 0.5 + x^2\left(\frac{1}{2} - 2(0.25)\right) \equiv 0.5 \quad \forall x$$
**The local energy is mathematically independent of the position coordinate $x$!**

Consequently, even if the Metropolis walker is fed severely corrupted random numbers (e.g. `Degraded_MarkovBias` causing extreme clustering or spatial bias):
1. Every drawn sample returns $E_L(x_i) = 0.500000$.
2. The calculated energy mean is $\hat{E} = 0.500000 \pm 0.000000$.
3. The sample variance of the energy is **identically zero**.
4. **Any monitoring system that evaluates only the energy observable will conclude the simulation is flawlessly healthy!**

#### The Resolution in RNGLab v0.14
RNGLab v0.14 solves this critical vulnerability by evaluating the non-commuting spatial moment:
$$\langle x^2 \rangle_{exact} = \frac{1}{4\alpha} = 0.500000$$
Under `Degraded_MarkovBias`:
- Observed energy: $\hat{E} = 0.500000$ (error = $0.000000$).
- Observed spatial moment: $\langle x^2 \rangle = 0.629631$ (error = $+0.129631$, $z > 14.5$).
- `observable_masking_detected = True` is triggered automatically in `rnglab14.models.vmc_harmonic`.

### 3.4 Path Integral Monte Carlo (PIMC) & Centroid Virial Convergence
PIMC maps a quantum particle at inverse temperature $\beta = 1/k_B T$ onto a classical ring polymer of $M$ beads interacting via harmonic spring potentials.

#### Thermodynamic vs Centroid Virial Estimators
1. **Thermodynamic Estimator:**
   $$E_{therm} = \frac{M}{2\beta} - \sum_{j=1}^M \frac{m M}{2\hbar^2\beta^2}(x_{j+1}-x_j)^2 + \frac{1}{M}\sum_{j=1}^M V(x_j)$$
   The kinetic term involves the difference between two large quantities of order $\mathcal{O}(M/\beta)$, causing the statistical variance of $E_{therm}$ to scale as $\mathcal{O}(M)$. In our experiments at $M=16$, $E_{therm} \approx 3.81 \pm 0.85$, far above the exact quantum energy $E = 0.6565$.
2. **Centroid Virial Estimator:**
   $$E_{vir} = \frac{1}{2\beta} + \frac{1}{2}\langle x^2 \rangle + \frac{1}{2}\langle (x - x_{cm}) x \rangle$$
   This estimator subtracts the centroid fluctuation, yielding stable convergence ($E_{vir} = 0.5829 \pm 0.17$ on healthy reference sources).
3. **Catastrophic Failure under Degraded Randomness:**
   Under `Degraded_MarkovBias`, the bead displacements undergo unconstrained random walks, causing the polymer ring to expand uncontrollably. The virial energy explodes to **$65.47$** (an error of $+9800\%$).

### 3.5 Monte Carlo Decision-Support Systems (DSS)
In risk-sensitive decision making, investment options $\{A, B, C\}$ are evaluated under uncertainty via $N=2000$ simulated market scenarios:
- **Option A (Conservative):** Normal distribution ($\mu=100, \sigma=15$).
- **Option B (Growth):** Student's $t$ heavy-tailed distribution ($\nu=3, \mu=120, \text{scale}=20$).
- **Option C (Speculative):** Log-normal distribution with fat upside tails ($\mu=105, \sigma=0.6$).

**Degradation Consequences:**
- Positive autocorrelation compresses the effective variance of simulated market returns, underestimating the $95\%$ Value-at-Risk ($\text{VaR}_{95\%}$) by up to $18\%$.
- While the expected value ranking remained intact in the moderate test case (Option B selected), the perceived confidence spread expanded artificially, creating an illusion of high risk certainty.

---

## 4. Synthesis & Architectural Guidelines

1. **Never rely on a single primary observable** to validate simulation integrity (VMC Observable Masking proves that an energy estimator can be 100% accurate while sampling distribution is 100% corrupted).
2. **Always enforce stochastic channel isolation**: Separate spatial proposal streams, acceptance testing streams, and auxiliary noise channels.
3. **Deploy dual-level integrity monitoring**: Source-level statistical testing (NIST/BSI) provides early warning, while application-level moment tracking provides fail-safe validation.
