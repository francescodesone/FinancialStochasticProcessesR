# 📈 FinancialStochasticProcessesR

**FinancialStochasticProcessR** is an R package designed to simulate and visualize a variety of **stochastic processes** commonly used in finance, econometrics, and time-series modeling.

The package provides flexible tools to generate, analyze, and plot random processes driven by different probability distributions — including normal, Student’s *t*, generalized error, and skewed distributions.

---

## 🚀 Features

- Simulation of **multiplicative random walks** with customizable drift and distributions  
- Simulation of **ARMA(p, q)** processes with flexible innovation distributions  
- Simulation of **GARCH(p, q)** volatility processes  
- Support for **heavy-tailed** and **skewed** distributions (GED, skewed-t, skewed-normal, etc.)  
- Simulation of **interest rate models** (Vasicek, CIR, Hull–White, Black–Karasinski, CKLS, CEV)  
- Parameter estimation functions with **bootstrap confidence intervals** (GBM, Vasicek, CIR)  
- Easy visualization of simulated paths with automatic coloring options  
- Unified plotting style with mean trajectory highlighted  

---

## 🧩 Processes Included

The package currently supports the following stochastic processes:

1. **Random Walk (RW)**  
   - Additive or multiplicative form  
   - Custom drift and volatility  
   - Distributions: Normal, Student’s *t*, Skewed Normal, Skewed Student’s *t*, GED, Skewed GED  

2. **ARMA(p, q) Processes**  
   - Flexible autoregressive and moving-average dynamics  
   - Customizable innovation distributions Normal, Student’s *t*, Skewed Normal, Skewed Student’s *t*, GED, Skewed GED 

3. **GARCH(p, q) Processes**  
   - Conditional heteroskedasticity models for volatility clustering  
   - Simulation of returns with time-varying variance  

4. **Geometric Brownian Motion (GBM)**  
   - Classic model for asset prices  
   - Simulation via `GBM.sim()`  
   - Parameter estimation via `GBM.estimate.params()` with optional bootstrap  

5. **Vasicek Model**  
   - Mean-reverting Ornstein–Uhlenbeck process  
   - Simulation via `Vasicek.sim()`  
   - Parameter estimation via `Vasicek.estimate.params()`  

6. **Cox–Ingersoll–Ross (CIR) Model**  
   - Mean-reverting square-root diffusion  
   - Simulation via `CIR.sim()`  
   - Parameter estimation via `CIR.estimate.params()`  

7. **Black–Karasinski Model**  
   - Lognormal short-rate model  
   - Functions:  
     - `BlackKarasinski.multi()` (linear approximation)  
     - `BlackKarasinskiLog.multi()` (lognormal version ensuring positivity)  

8. **CKLS Model (Chan–Karolyi–Longstaff–Sanders)**  
   - Generalized diffusion model  
   - Special cases: Vasicek, CIR, GBM  
   - Function: `CKLS.sim()`  

9. **CEV Model (Constant Elasticity of Variance)**  
    - Generalization of GBM with elasticity parameter γ  
    - Function: `CEVProcess.single()`  
    - Supports theoretical expectation and variance for γ = 1 (GBM case)  

---

## 📊 Parameter Estimation Functions

- **GBM.estimate.params()**  
  - Estimates drift (μ) and volatility (σ) from price series  
  - Supports block bootstrap confidence intervals  

- **Vasicek.estimate.params()**  
  - Estimates mean-reversion speed (a), long-term mean (b), and volatility (σ)  
  - Supports block bootstrap confidence intervals  

- **CIR.estimate.params()**  
  - Estimates CIR parameters (a, b, σ) using transformation method  
  - Supports block bootstrap confidence intervals  

All estimation functions print results with **colored console headings** using the `cli` package.

---

## 📘 Examples

```r
# GBM parameter estimation
set.seed(123)
prices <- cumprod(100 * exp(rnorm(252, mean = 0.0003, sd = 0.01)))
GBM.estimate.params(price.series = prices, dt = 1/252, name = "Simulated Asset")

# Vasicek simulation
set.seed(42)
vas <- Vasicek.sim(initial.rate = 0.02, a = 0.1, b = 0.03,
                   sigma_r = 0.01, T = 1, N = 100)

# CIR simulation
set.seed(42)
cir <- CIR.sim(initial.rate = 0.02, a = 0.1, b = 0.03,
               sigma_r = 0.01, T = 1, N = 100)

# Black–Karasinski lognormal simulation
set.seed(42)
bk_log <- BlackKarasinskiLog.multi(r0 = 0.03, n = 1000, dt = 1/252,
                                   theta = seq(0.03, 0.05, length.out = 1000),
                                   sigma = seq(0.02, 0.04, length.out = 1000),
                                   alpha = rep(0.1, 1000), num.paths = 20)

# CKLS simulation
set.seed(42)
ckls <- CKLS.sim(initial.rate = 0.02, alpha = 0.5, beta = 0.5,
                 sigma = 0.1, gamma = 1.5, T = 1, N = 50)

# CEV simulation (GBM case)
set.seed(42)
cev <- CEVProcess.single(mu = 0.1, sigma = 0.3, gamma = 1,
                         x0 = 1, T = 1, n = 200, N = 100)

# ARMA process simulation
set.seed(42)
arma_paths <- ARMA.sim(p = 2, q = 1, n = 500, N = 20,
                       dist = "normal", plot = TRUE)

# GARCH process simulation
set.seed(42)
garch_paths <- GARCH.sim(p = 1, q = 1, n = 1000, N = 50,
                         omega = 0.01, alpha = 0.1, beta = 0.85,
                         plot = TRUE)
