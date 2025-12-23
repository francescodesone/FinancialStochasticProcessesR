# 📈 FinancialStochasticProcessesR

**FinancialStochasticProcessR** is an R package designed to simulate and visualize a variety of **stochastic processes** commonly used in finance, econometrics, and time-series modeling.

The package provides flexible tools to generate, analyze, and plot random processes driven by different probability distributions — including normal, Student’s *t*, generalized error, and skewed distributions.

---

## 🚀 Features

- Simulation of **random walks** with customizable drift and distributions  
- Simulation of **ARMA(p, q)** processes with flexible innovation distributions  
- Simulation of **GARCH(p, q)** volatility processes  
- Support for **heavy-tailed** and **skewed** distributions (GED, skewed-t, skewed-normal, etc.)  
- Simulation of **interest rate models** (Vasicek, CIR, Black–Karasinski, CKLS)  
- Parameter estimation functions with **bootstrap confidence intervals** (GBM, Vasicek, CIR, Jump Diffusion)  
- Easy visualization of simulated paths with automatic coloring options 

---

## 🧩 Processes Included

The package currently supports the following stochastic processes:

1. **Random Walk (RW)**  
   - Additive or multiplicative form  
   - Custom drift and volatility  
   - Distributions: Normal, Student’s *t*, Skewed Normal, Skewed Student’s *t*, GED, Skewed GED  

2. **ARMA(p, q) Processes**  
   - Flexible autoregressive and moving-average dynamics  
   - Distributions: Normal, Student’s *t*, Skewed Normal, Skewed Student’s *t*, GED, Skewed GED

3. **GARCH(p, q) and APARCH Processes**  
   - Conditional heteroskedasticity models for volatility clustering  
   - Simulation of returns with time-varying variance
   - Distributions: Normal, Student’s *t*, Skewed Normal, Skewed Student’s *t*, GED, Skewed GED

4. **Geometric Brownian Motion (GBM)**  
   - Classic model for asset prices  
   - Simulation via `GBM.sim()`  
   - Parameter estimation via `GBM.estimate.params()` with optional bootstrap  

5. **Jump Diffusion Model**  
   - Jump Diffusion Process, combining Poisson Jumps and Brownian Motion 
   - Simulation via `JumpDiffusion.sim()`  
   - Parameter estimation via `JD.estimate.params()`
     
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
 
- **JD.estimate.params()**  
  - Estimates the mean (μ), volatility (σ), Jump size (g), and Jump frequency (phi)
  - Supports block bootstrap confidence intervals  

- **CIR.estimate.params()**  
  - Estimates CIR parameters (a, b, σ) using transformation method  
  - Supports block bootstrap confidence intervals  

All estimation functions print results with **colored console headings** using the `cli` package.


---

## 🔄 Continuous Improvement

**FinancialStochasticProcessesR** is an evolving project.  
The package will be **updated and improved over time**, with new models, enhanced parameter estimation methods, and expanded visualization tools.  

Future updates will include:
- Additional stochastic processes and volatility models  
- Extended parameter estimation functions 
- Vignettes and tutorials for practical applications  
- Performance optimizations and expanded plotting options  

---

## 🙌 Final Notes

This package has been built to provide a **comprehensive toolkit** for simulating, visualizing, and estimating parameters of stochastic processes in finance and econometrics.  

All functions are designed to be **flexible, easy to use, and freely available** for experimentation, research, and learning purposes.  

Feel free to explore, adapt, and extend the functions to suit your own projects — whether academic, professional, or personal.  

---

## 👋 Greetings

Thank you for using **FinancialStochasticProcessesR**!  
We hope this package helps you deepen your understanding of stochastic modeling and empowers you to run your own simulations with confidence.  

Enjoy experimenting with the models, and remember:  
**all functions can be used freely, customized, and integrated into your workflows.**

## ⚠️ Disclaimer This package is intended **solely for educational, research, and illustrative purposes**. It provides tools to simulate and estimate stochastic processes commonly used in finance and econometrics. It is **not a financial advisory tool** and should not be used for trading, investment decisions, or professional risk management without proper validation. Users are responsible for verifying the accuracy and suitability of the models for their own applications.
