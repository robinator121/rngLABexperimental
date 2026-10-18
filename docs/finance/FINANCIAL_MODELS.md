# RNGLab-Finance — Mathematical & Algorithmic Model Specifications

**Document Identifier:** SPEC-FIN-v01-01  
**Classification:** Quantitative Financial Mathematics & Reference Benchmarks  
**Version:** v0.1-PoC (`feat/finance-poc`)  
**Date:** 2026-10-18  

---

## 1. Overview & Scope

RNGLab-Finance establishes a benchmark suite of quantitative financial models spanning derivative valuation, market risk quantiles, and coherent tail risk measures.

Every model is anchored by an **independent, closed-form analytical reference ground truth** to eliminate model ambiguity and isolate computational randomness defects.

---

## 2. Model A: European Option Pricing under Geometric Brownian Motion

### 2.1 Continuous-Time Dynamics
Under the risk-neutral probability measure $\mathbb{Q}$, the asset price $S_t$ follows the stochastic differential equation (SDE):
$$dS_t = r S_t dt + \sigma S_t dW_t$$
where $r$ is the continuously compounded risk-free interest rate, $\sigma > 0$ is the constant volatility, and $W_t$ is a standard Brownian motion.

By Itô's Lemma, the exact terminal price at maturity $T$ is given analytically by:
$$S_T = S_0 \exp\left( \left(r - \frac{1}{2}\sigma^2\right)T + \sigma \sqrt{T} Z \right), \quad Z \sim \mathcal{N}(0, 1)$$

### 2.2 Analytical Black-Scholes Reference Ground Truth
The exact price of a European Call option $C(S_0, K, r, \sigma, T)$ is:
$$C_{BS} = S_0 \Phi(d_1) - K e^{-r T} \Phi(d_2)$$
and European Put option:
$$P_{BS} = K e^{-r T} \Phi(-d_2) - S_0 \Phi(-d_1)$$
where:
$$d_1 = \frac{\ln(S_0 / K) + (r + \frac{1}{2}\sigma^2)T}{\sigma \sqrt{T}}, \quad d_2 = d_1 - \sigma \sqrt{T}$$
and $\Phi(\cdot)$ denotes the standard normal cumulative distribution function (CDF).

### 2.3 Economic Invariants for Integrity Monitoring
1. **The Risk-Neutral Martingale Invariant:**
   $$\mathbb{E}^{\mathbb{Q}}\left[ e^{-r T} S_T \right] = S_0$$
   *Diagnostic Invariant Test:*
   $$\hat{S}_{disc} = \frac{1}{N} \sum_{i=1}^N e^{-r T} S_T^{(i)}, \quad Z_{martingale} = \frac{\hat{S}_{disc} - S_0}{\hat{\sigma}_{disc} / \sqrt{N}}$$
   If $|Z_{martingale}| > 3.0$, the simulation fails the martingale condition ($p < 0.0027$).
2. **Put-Call Parity Identity:**
   $$C - P = S_0 - K e^{-r T}$$
   Any asymmetric defect in positive vs negative innovations $Z$ produces an immediate discrepancy between Call and Put valuations.

### 2.4 Variance Reduction via Antithetic Variates
For every normal innovation $Z_i \sim \mathcal{N}(0, 1)$, an antithetic innovation $\tilde{Z}_i = -Z_i$ is evaluated:
$$S_T^{(1)} = S_0 \exp(\mu_d + \sigma_d Z_i), \quad S_T^{(2)} = S_0 \exp(\mu_d - \sigma_d Z_i)$$
The antithetic estimator $\hat{C}_{anti} = \frac{1}{2}(C(S_T^{(1)}) + C(S_T^{(2)}))$ guarantees an exact sample mean of $\bar{Z} = 0$, eliminating first-order innovation bias and significantly reducing Monte Carlo standard error.

---

## 3. Model B: Multi-Asset Portfolio Value at Risk (VaR)

### 3.1 Mathematical Specification
Consider a portfolio of $K$ assets with initial value $V_0 = \sum_{k=1}^K w_k V_0$, where $w = [w_1, \dots, w_K]^T$ ($\sum w_k = 1$). Asset returns $R = [R_1, \dots, R_K]^T$ over horizon $\Delta t$ follow a multivariate normal distribution:
$$R \sim \mathcal{N}(\mu \Delta t, \Sigma \Delta t)$$
where $\mu$ is the annual expected drift vector and $\Sigma$ is the $K \times K$ positive-definite covariance matrix.

### 3.2 Correlated Simulation via Cholesky Factorization
Let $\Sigma = L L^T$ be the lower-triangular Cholesky factor. Correlated return vectors are generated from independent standard normal draws $Z \sim \mathcal{N}(0, I_K)$:
$$R = \mu \Delta t + \sqrt{\Delta t} L Z$$
The portfolio return is $R_p = w^T R$, and the dollar portfolio loss is:
$$L = - V_0 R_p$$

### 3.3 Analytical Gaussian VaR Benchmark
For a linear portfolio, the distribution of loss $L$ is Gaussian:
$$\mu_L = - V_0 (w^T \mu \Delta t), \quad \sigma_L = V_0 \sqrt{w^T \Sigma w \Delta t}$$
The exact Value at Risk at confidence level $\alpha \in (0, 1)$ is:
$$\text{VaR}_\alpha = \mu_L + \Phi^{-1}(\alpha) \sigma_L$$

### 3.4 Kupiec Proportion of Failures (POF) Likelihood Ratio Test
Given $N$ simulated loss scenarios, the count of exceedances beyond the theoretical $\text{VaR}_\alpha$ is $X = \sum_{i=1}^N \mathbb{I}(L_i > \text{VaR}_\alpha)$. Under the null hypothesis that the simulation is correctly calibrated:
$$X \sim \text{Binomial}(N, 1 - \alpha)$$
The Kupiec Likelihood Ratio statistic is:
$$LR_{POF} = -2 \ln \left[ \frac{(1 - \alpha)^{N - X} \alpha^X}{\left(1 - \frac{X}{N}\right)^{N - X} \left(\frac{X}{N}\right)^X} \right] \sim \chi^2(1)$$
If $LR_{POF} > 3.841$, the risk model is rejected at the $\alpha = 0.05$ significance level.

---

## 4. Model C: Expected Shortfall (ES / CVaR) under Basel III / FRTB

### 4.1 Definition & Coherent Risk Properties
Under the Basel III Fundamental Review of the Trading Book (FRTB), Expected Shortfall replaced 99% VaR as the primary regulatory market risk metric.

Expected Shortfall at confidence level $\alpha$ is defined as the expected loss conditional on the loss exceeding the Value at Risk quantile:
$$\text{ES}_\alpha = \mathbb{E}\left[ L \;\middle|\; L \ge \text{VaR}_\alpha \right] = \frac{1}{1 - \alpha} \int_\alpha^1 \text{VaR}_u du$$
Unlike VaR, Expected Shortfall is a **coherent risk measure** satisfying sub-additivity:
$$\text{ES}_\alpha(L_1 + L_2) \le \text{ES}_\alpha(L_1) + \text{ES}_\alpha(L_2)$$

### 4.2 Analytical Gaussian Reference Ground Truth
For normally distributed loss $L \sim \mathcal{N}(\mu_L, \sigma_L^2)$:
$$\text{ES}_\alpha = \mu_L + \sigma_L \frac{\phi(\Phi^{-1}(\alpha))}{1 - \alpha}$$
where $\phi(\cdot)$ is the standard normal probability density function.

At the FRTB standard $\alpha = 0.975$:
$$\frac{\phi(\Phi^{-1}(0.975))}{1 - 0.975} = \frac{\phi(1.95996)}{0.025} \approx 2.3378$$

### 4.3 Monte Carlo Tail Estimator
$$\widehat{\text{ES}}_\alpha = \frac{1}{N_{tail}} \sum_{L_i \ge \widehat{\text{VaR}}_\alpha} L_i, \quad N_{tail} = \sum_{i=1}^N \mathbb{I}(L_i \ge \widehat{\text{VaR}}_\alpha)$$
Because Expected Shortfall is estimated exclusively from the $(1 - \alpha)$ tail sample (e.g. 1,000 out of 40,000 paths at $\alpha=0.975$), **it is exponentially more vulnerable to PRNG tail corruption, periodic spikes, and buffer clipping than mean option pricing.**

---

## 5. Candidate Advanced Models for Future Commercial Expansion

| Candidate Model | Mathematical Formulation | Primary Sensitivity to RNG Defects | Commercial Target |
| :--- | :--- | :--- | :--- |
| **Heston Stochastic Volatility** | $dS_t = r S_t dt + \sqrt{v_t} S_t dW_t^{(1)}$, $dv_t = \kappa(\theta - v_t)dt + \xi \sqrt{v_t} dW_t^{(2)}$ with $dW_1 dW_2 = \rho dt$ | Flawed PRNG coupling alters correlation $\rho$, causing smile/skew collapse | Exotic derivatives & FX desks |
| **Path-Dependent Asian Options** | Payoff based on arithmetic average $\bar{S} = \frac{1}{M}\sum_{m=1}^M S_{t_m}$ | Serial autocorrelation in PRNG artificially expands or contracts $\text{Var}(\bar{S})$ | Commodity & energy trading |
| **Credit Risk / Default Copulas** | Latent factor $X_i = \sqrt{\rho} Y + \sqrt{1-\rho} Z_i$; default if $X_i \le \Phi^{-1}(PD_i)$ | Seed collisions cause simultaneous false defaults, distorting tranches | CDO / credit portfolio risk |
| **Quasi-Monte Carlo (Sobol Sequences)** | Low-discrepancy deterministic sequences achieving $O((\ln N)^d / N)$ error rate | Distinguishing intentional low discrepancy from true RNG corruption | High-frequency pricing engines |
