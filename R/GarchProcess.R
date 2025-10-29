#------------------------------------------------------------
# Function: sim_garch()
# Purpose:  Simulate time series data from an ARMA-GARCH or
#           APARCH model under various conditional distributions
#           (normal, GED, t-Student, skewed variants).
#------------------------------------------------------------

#' Simulate ARMA-GARCH/APARCH Time Series
#'
#' This function simulates one or more realizations of an ARMA-GARCH or
#' APARCH (Asymmetric Power ARCH) process with user-defined parameters
#' for mean and volatility dynamics. It supports several conditional
#' distributions including normal, GED, Student-t, and their skewed versions.
# This function is based on fGarch Package
#' @param n Integer. Number of observations to simulate (after burn-in).
#' @param n.start Integer. Number of initial "burn-in" observations discarded to allow the process to stabilize.
#' @param ar Numeric vector. AR (autoregressive) coefficients for the conditional mean equation.
#' @param ma Numeric vector. MA (moving average) coefficients for the conditional mean equation.
#' @param omega Numeric. Constant term in the variance equation.
#' @param alpha Numeric vector. ARCH parameters controlling reaction to past shocks.
#' @param beta Numeric vector. GARCH parameters controlling persistence of volatility.
#' @param gamma Numeric vector. Asymmetric (leverage) parameters for APARCH specification.
#'              If omitted (default `numeric(0)`), symmetry is assumed.
#' @param delta Numeric. Power parameter for APARCH model (default = 2 for standard GARCH).
#' @param mu Numeric. Constant term in the mean equation.
#' @param dist Character. Conditional distribution of innovations:
#'             one of `"norm"`, `"ged"`, `"std"`, `"snorm"`, `"sged"`, `"sstd"`.
#' @param nu Numeric. Shape parameter for heavy-tailed distributions
#'            (used in GED and Student-t).
#' @param skew Numeric. Skewness parameter for skewed distributions
#'             (used in "snorm", "sged", "sstd").
#' @param seed Optional integer. Random seed for reproducibility.
#' @param n_sim Integer. Number of independent simulated series to generate.
#' @param extended Logical. If `TRUE`, returns a list containing simulated
#'                 values, volatilities, residuals, and innovations.
#'                 If `FALSE` (default), only simulated series `y` are returned.
#' @param plot Logical. If `TRUE`, plots the simulated series using `matplot()`.
#'
#' @return If `extended = FALSE`, returns a numeric matrix of dimension
#'         \code{n × n_sim} containing the simulated series.
#'         If `extended = TRUE`, returns a list with components:
#'         \itemize{
#'           \item \code{y} — simulated series
#'           \item \code{sigma} — conditional standard deviations
#'           \item \code{eps} — conditional residuals
#'           \item \code{z} — standardized innovations
#'         }
#'
#' @details
#' The model simulated is of the general ARMA–APARCH form:
#' \deqn{
#' y_t = \mu + \sum \phi_i y_{t-i} + \sum \theta_j \epsilon_{t-j} + \epsilon_t
#' }
#' with
#' \deqn{
#' \epsilon_t = h_t^{1/\delta} z_t, \quad
#' h_t = \omega + \sum \alpha_i (|\epsilon_{t-i}| - \gamma_i \epsilon_{t-i})^{\delta}
#'       + \sum \beta_j h_{t-j}
#' }
#' where \(z_t\) are i.i.d. draws from the chosen distribution.
#'
#' @examples
#' set.seed(123)
#'
#' # Basic GARCH(1,1) with normal innovations
#' sim_garch(n = 200, n.start = 50, dist = "norm", n_sim = 2)
#'
#' # GED (Generalized Error Distribution)
#' sim_garch(n = 200, n.start = 50, dist = "ged", nu = 5, n_sim = 2)
#'
#' # Student-t with 3 degrees of freedom
#' sim_garch(n = 200, n.start = 50, dist = "std", nu = 3, n_sim = 2)
#'
#' # Skewed normal
#' sim_garch(n = 200, n.start = 50, dist = "snorm", skew = 1.5, n_sim = 2)
#'
#' # Skewed GED
#' sim_garch(n = 200, n.start = 50, dist = "sged", nu = 5, skew = 1.5, n_sim = 2)
#'
#' # Skewed Student-t
#' sim_garch(n = 200, n.start = 50, dist = "sstd", nu = 3, skew = 1.5, n_sim = 2)
#'
#' @export
sim_garch <- function(n = 100,
                      n.start = 100,
                      ar = numeric(0),
                      ma = numeric(0),
                      omega = 0.1,
                      alpha = 0.1,
                      beta = 0.8,
                      gamma = numeric(0),    # leverage effect vector (APARCH)
                      delta = 2,             # power term for volatility
                      mu = 0,
                      dist = "norm",
                      nu = 5,                # shape parameter for GED/t etc.
                      skew = 1.5,            # skewness param
                      seed = NULL,
                      n_sim = 1,
                      extended = FALSE,
                      plot = TRUE) {

  if (!is.null(seed)) set.seed(seed)

  order.ar <- length(ar)
  order.ma <- length(ma)
  order.alpha <- length(alpha)
  order.beta <- length(beta)
  order.gamma <- length(gamma)

  stopifnot(order.alpha == order.gamma || order.gamma == 0)

  sims_y <- matrix(NA, nrow = n, ncol = n_sim)
  sims_sigma <- matrix(NA, nrow = n, ncol = n_sim)
  sims_eps <- matrix(NA, nrow = n, ncol = n_sim)
  sims_z <- matrix(NA, nrow = n, ncol = n_sim)

  for (sim_i in 1:n_sim) {
    n_total <- n + n.start

    z <- switch(dist,
                "norm"  = rnorm(n_total),
                "ged"   = rged(n_total, nu = nu),
                "std"   = rstd(n_total, nu = nu),
                "snorm" = rsnorm(n_total, xi = skew),
                "sged"  = rsged(n_total, nu = nu, xi = skew),
                "sstd"  = rsstd(n_total, nu = nu, xi = skew),
                stop("Unsupported distribution"))

    h <- rep(NA, n_total)
    eps <- rep(NA, n_total)
    y <- rep(NA, n_total)

    h[1] <- omega / (1 - sum(alpha) - sum(beta))
    eps[1] <- h[1]^(1/delta) * z[1]
    y[1] <- mu + eps[1]

    for (t in 2:n_total) {
      arch_term <- 0
      if (order.alpha > 0) {
        lags <- eps[(t - 1):(t - order.alpha)]
        gamma_vec <- if (order.gamma == 0) rep(0, order.alpha) else gamma
        arch_term <- sum(alpha * (abs(lags) - gamma_vec * lags)^delta)
      }

      garch_term <- if (order.beta > 0) sum(beta * h[(t - 1):(t - order.beta)]) else 0

      h[t] <- omega + arch_term + garch_term

      eps[t] <- h[t]^(1/delta) * z[t]

      ar_part <- if (order.ar > 0 && t > order.ar) sum(ar * y[(t - 1):(t - order.ar)]) else 0
      ma_part <- if (order.ma > 0 && t > order.ma) sum(ma * eps[(t - 1):(t - order.ma)]) else 0

      y[t] <- mu + ar_part + ma_part + eps[t]
    }

    y <- y[(n.start + 1):n_total]
    h <- h[(n.start + 1):n_total]
    eps <- eps[(n.start + 1):n_total]
    z <- z[(n.start + 1):n_total]

    sims_y[, sim_i] <- y
    sims_sigma[, sim_i] <- h^(1/delta)
    sims_eps[, sim_i] <- eps
    sims_z[, sim_i] <- z
  }

  colnames(sims_y) <- paste0("Sim", 1:n_sim)
  colnames(sims_sigma) <- paste0("Sim", 1:n_sim)
  colnames(sims_eps) <- paste0("Sim", 1:n_sim)
  colnames(sims_z) <- paste0("Sim", 1:n_sim)

  if (plot) {
    matplot(1:n, sims_y, type = "l", lty = 1, col = rainbow(n_sim),
            main = paste0("ARMA-GARCH Simulation (dist = ", dist, ")"),
            xlab = "Time", ylab = "y")
  }

  if (extended) {
    return(list(y = sims_y,
                sigma = sims_sigma,
                eps = sims_eps,
                z = sims_z))
  } else {
    return(sims_y)
  }
}


