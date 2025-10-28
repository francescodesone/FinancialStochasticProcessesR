#------------------------------------------------------------
# Function: sim_rw()
# Purpose: Simulate a multiplicative random walk with
#          optional drift and various distributions.
#------------------------------------------------------------
#' Simulate Multiplicative Random Walks with Custom Distributions
#'
#' @param n Integer. Number of time steps per simulation (default = 100).
#' @param sigma Numeric. Standard deviation of the innovations (default = 0.01).
#' @param S0 Numeric. Initial value of the random walk (default = 1).
#' @param n_sim Integer. Number of simulated paths (default = 1).
#' @param dist Character. Type of distribution for innovations.
#'   Options: `"norm"`, `"ged"`, `"sged"`, `"snorm"`, `"std"`, `"sstd"`.
#'   (default = `"norm"`).
#' @param df Numeric. Degrees of freedom for Student’s t or skewed t distributions (default = 5).
#' @param xi Numeric. Skewness parameter for skewed distributions (default = 1.5).
#' @param drift Numeric. Additive drift applied to each shock (default = 0).
#' @param seed Integer or NULL. Random seed for reproducibility (default = NULL).
#' @param plot Logical. If TRUE, plot the simulated random walks (default = TRUE).
#' @param use.rainbow Logical. If TRUE, use rainbow colors for multiple paths (default = TRUE).
#' @param fixed.color Character. Color used when `use.rainbow = FALSE` (default = "darkred").
#'
#' @return
#' A numeric matrix with \eqn{n + 1} rows (including initial value) and `n_sim` columns.
#' Each column represents one simulated random walk path.
#'
#' @export
#'
#' @examples
#' # Example 1: Default normal distribution with 3 simulated paths
#' set.seed(123)
#' sim_rw(n = 100, n_sim = 3, S0 = 100, drift = 0.002)
#'
#' # Example 2: Student's t-distribution with 5 degrees of freedom
#' sim_rw(n = 100, n_sim = 5, dist = "std", df = 5)
#'
#' # Example 3: Skewed Student's t-distribution
#' sim_rw(n = 100, n_sim = 10, dist = "sstd", xi = 10, df = 3, S0 = 10)
sim_rw <- function(n = 100,                   # number of time steps
                   sigma = 0.01,              # standard deviation of returns
                   S0 = 1,                    # initial value
                   n_sim = 1,                 # number of simulations
                   dist = "norm",             # distribution type
                   df = 5,                    # degrees of freedom (t / std)
                   xi = 1.5,                  # skewness parameter
                   drift = 0,                 # optional additive drift
                   seed = NULL,               # random seed
                   plot = TRUE,               # whether to plot the result
                   use.rainbow = TRUE,        # use rainbow colors for multiple sims
                   fixed.color = "darkred") { # fallback color if not using rainbow

  #-------------------------------
  # Set seed if provided
  #-------------------------------
  if (!is.null(seed)) set.seed(seed)

  n_total <- n * n_sim

  #-------------------------------
  # Generate shocks (innovations)
  #-------------------------------
  shocks <- switch(dist,
                   "norm"  = matrix(rnorm(n_total, mean = 0, sd = sigma), nrow = n),
                   "ged"   = matrix(rged(n_total, mean = 0, sd = sigma, nu = 2), nrow = n),
                   "sged"  = matrix(rsged(n_total, mean = 0, sd = sigma, nu = 2, xi = xi), nrow = n),
                   "snorm" = matrix(rsnorm(n_total, mean = 0, sd = sigma, xi = xi), nrow = n),
                   "std"   = matrix(rstd(n_total, mean = 0, sd = sigma, nu = df), nrow = n),
                   "sstd"  = matrix(rsstd(n_total, mean = 0, sd = sigma, nu = df, xi = xi), nrow = n),
                   stop("Unsupported distribution type. Choose one of: 'norm', 'ged', 'sged', 'snorm', 'std', 'sstd'.")
  )

  #-------------------------------
  # Add drift term
  #-------------------------------
  shocks <- shocks + drift

  #-------------------------------
  # Compute multiplicative random walk
  #-------------------------------
  rw <- apply(shocks + 1, 2, cumprod)

  # Prepend initial value S0
  rw <- rbind(rep(S0, n_sim), rw * S0)

  colnames(rw) <- paste0("Sim", 1:n_sim)

  #-------------------------------
  # Plot simulated paths
  #-------------------------------
  if (plot) {
    colors <- if (use.rainbow) {
      rainbow(n_sim, alpha = 0.6)
    } else {
      rep(adjustcolor(fixed.color, alpha.f = 0.6), n_sim)
    }

    matplot(rw, type = "l", lty = 1, col = colors,
            main = paste("Random Walk (dist =", dist, ", drift =", drift, ")"),
            xlab = "Time", ylab = "Value")
  }

  return(rw)
}

#------------------------------------------------------------
# Function: sim_arma()
# Purpose: Simulate ARMA(p, q) processes with different
#          innovation distributions.
#------------------------------------------------------------
#' Simulate ARMA(p, q) Processes with Flexible Innovation Distributions
#'
#' This function generates one or multiple simulations of an
#' Autoregressive Moving Average (ARMA) process with optional intercept,
#' and innovations drawn from a variety of distributions such as
#' normal, Student’s t, skewed t, GED, and their skewed counterparts.
#'
#' This is particularly useful in financial econometrics where return series
#' often exhibit skewness and heavy tails.
#'
#' @param n Integer. Length of the simulated series (default = 100).
#' @param ar Numeric vector. Autoregressive coefficients (e.g., `c(0.5, -0.2)` for AR(2)).
#' @param ma Numeric vector. Moving average coefficients (e.g., `0.3` for MA(1)).
#' @param intercept Numeric. Constant term added at each step (default = 0).
#' @param sd Numeric. Standard deviation of the innovations (default = 1).
#' @param dist Character. Type of innovation distribution.
#'   Options: `"norm"`, `"ged"`, `"sged"`, `"snorm"`, `"std"`, `"sstd"`.
#'   Default is `"norm"`.
#' @param df Numeric. Degrees of freedom for t or skewed t distributions (default = 5).
#' @param xi Numeric. Skewness parameter for skewed distributions (default = 1.5).
#' @param nu Numeric. Shape parameter for GED and skewed GED distributions (default = 2).
#' @param n_sim Integer. Number of simulated ARMA series (default = 1).
#' @param seed Integer or NULL. Random seed for reproducibility (default = NULL).
#' @param plot Logical. If TRUE, plot the simulated time series (default = TRUE).
#'
#' @return
#' A numeric matrix with `n` rows and `n_sim` columns, where each column
#' represents one simulated ARMA series.
#'
#'
#' @export
#'
#' @examples
#' # Example 1: ARMA(1,1) with skewed Student's t innovations
#' set.seed(123)
#' sim_arma(n = 200, ar = 0.5, ma = -0.3, dist = "sstd", xi = 1.8, df = 10, n_sim = 3)
#'
#' # Example 2: ARMA(2,1) with GED innovations
#' sim_arma(n = 250, ar = c(0.3, -0.2), ma = 0.5, dist = "ged", nu = 1.3, n_sim = 3)
#'
#' # Example 3: Simple AR(1) Gaussian process
#' sim_arma(n = 100, ar = 0.8, dist = "norm", n_sim = 5)
sim_arma <- function(n = 100,
                     ar = numeric(0),
                     ma = numeric(0),
                     intercept = 0,
                     sd = 1,
                     dist = "norm",
                     df = 5,
                     xi = 1.5,
                     nu = 2,
                     n_sim = 1,
                     seed = NULL,
                     plot = TRUE) {

  #-------------------------------
  # Set seed if provided
  #-------------------------------
  if (!is.null(seed)) set.seed(seed)

  p <- length(ar)
  q <- length(ma)

  #-------------------------------
  # Helper: simulate one ARMA series
  #-------------------------------
  sim_single <- function() {
    # Generate innovations based on chosen distribution
    innovations <- switch(dist,
                          "norm"  = rnorm(n, mean = 0, sd = sd),
                          "ged"   = rged(n, mean = 0, sd = sd, nu = nu),
                          "sged"  = rsged(n, mean = 0, sd = sd, nu = nu, xi = xi),
                          "snorm" = rsnorm(n, mean = 0, sd = sd, xi = xi),
                          "std"   = rstd(n, mean = 0, sd = sd, nu = df),
                          "sstd"  = rsstd(n, mean = 0, sd = sd, nu = df, xi = xi),
                          stop("Unsupported distribution: use one of 'norm', 'ged', 'sged', 'snorm', 'std', 'sstd'")
    )

    # Initialize output
    x <- numeric(n)

    # Recursive computation
    for (t in 1:n) {
      ar_part <- if (p > 0 && t > p) sum(ar * x[(t-1):(t-p)]) else 0
      ma_part <- if (q > 0 && t > q) sum(ma * innovations[(t-1):(t-q)]) else 0
      x[t] <- intercept + ar_part + innovations[t] + ma_part
    }

    return(x)
  }

  #-------------------------------
  # Replicate simulation for n_sim paths
  #-------------------------------
  sims <- replicate(n_sim, sim_single())
  colnames(sims) <- paste0("Sim", 1:n_sim)
  rownames(sims) <- 1:n

  #-------------------------------
  # Plot the simulated series
  #-------------------------------
  if (plot) {
    matplot(1:n, sims, type = "l", lty = 1, col = rainbow(n_sim),
            xlab = "Time", ylab = "Value",
            main = paste0("ARMA(", length(ar), ",", length(ma), ") - dist = ", dist))
  }

  return(sims)
}
