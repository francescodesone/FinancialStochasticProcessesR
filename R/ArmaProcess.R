#------------------------------------------------------------
# Function: sim_rw()
# Purpose: Simulate a multiplicative random walk with
#          optional drift and various distributions.
#------------------------------------------------------------
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

  if (!is.null(seed)) set.seed(seed)

  n_total <- n * n_sim

  # Generate shocks (mean = 0, sd = sigma)
  shocks <- switch(dist,
                   "norm"  = matrix(rnorm(n_total, mean = 0, sd = sigma), nrow = n),
                   "ged"   = matrix(rged(n_total, mean = 0, sd = sigma, nu = 2), nrow = n),
                   "sged"  = matrix(rsged(n_total, mean = 0, sd = sigma, nu = 2, xi = xi), nrow = n),
                   "snorm" = matrix(rsnorm(n_total, mean = 0, sd = sigma, xi = xi), nrow = n),
                   "std"   = matrix(rstd(n_total, mean = 0, sd = sigma, nu = df), nrow = n),
                   "sstd"  = matrix(rsstd(n_total, mean = 0, sd = sigma, nu = df, xi = xi), nrow = n),
                   stop("Unsupported distribution")
  )

  # Add drift to shocks
  shocks <- shocks + drift

  # Compute the multiplicative random walk (1 + return)
  rw <- apply(shocks + 1, 2, cumprod)

  # Prepend the initial value S0
  rw <- rbind(rep(S0, n_sim), rw * S0)

  colnames(rw) <- paste0("Sim", 1:n_sim)

  # Plot the simulated paths
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
# Examples
#------------------------------------------------------------

# Example 1: Default normal distribution with 3 simulations
sim1 <- sim_rw(n = 100, n_sim = 3, S0 = 100, drift = 0.002)
# → Generates 3 multiplicative random walk paths with normal shocks (mean = 0, sd = 0.01)

# Example 2: Student’s t-distribution with 5 degrees of freedom, 100 simulations
sim2 <- sim_rw(n = 100, n_sim = 100, dist = "std", df = 5, use.rainbow = TRUE)

# Example 3: Skewed Student’s t-distribution with xi = 10 and df = 3, 100 simulations
sim3 <- sim_rw(n = 100, n_sim = 100, dist = "sstd", xi = 10, df = 3, S0 = 10)



#------------------------------------------------------------
# Function: sim_arma()
# Purpose: Simulate ARMA(p, q) processes with different
#          innovation distributions.
#------------------------------------------------------------
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

  if (!is.null(seed)) set.seed(seed)

  p <- length(ar)
  q <- length(ma)

  sim_single <- function() {
    # Generate innovations based on the chosen distribution
    innovations <- switch(dist,
                          "norm"  = rnorm(n, mean = 0, sd = sd),
                          "ged"   = rged(n, mean = 0, sd = sd, nu = nu),
                          "sged"  = rsged(n, mean = 0, sd = sd, nu = nu, xi = xi),
                          "snorm" = rsnorm(n, mean = 0, sd = sd, xi = xi),
                          "std"   = rstd(n, mean = 0, sd = sd, nu = df),
                          "sstd"  = rsstd(n, mean = 0, sd = sd, nu = df, xi = xi),
                          stop("Unsupported distribution: use one of 'norm', 'ged', 'sged', 'snorm', 'std', 'sstd'")
    )

    x <- numeric(n)
    for (t in 1:n) {
      ar_part <- if (p > 0 && t > p) sum(ar * x[(t-1):(t-p)]) else 0
      ma_part <- if (q > 0 && t > q) sum(ma * innovations[(t-1):(t-q)]) else 0
      x[t] <- intercept + ar_part + innovations[t] + ma_part
    }

    return(x)
  }

  sims <- replicate(n_sim, sim_single())
  colnames(sims) <- paste0("Sim", 1:n_sim)
  rownames(sims) <- 1:n

  # Plot the simulated ARMA series
  if (plot) {
    matplot(1:n, sims, type = "l", lty = 1, col = rainbow(n_sim),
            xlab = "Time", ylab = "Value",
            main = paste0("ARMA(", length(ar), ",", length(ma), ") - dist = ", dist))
  }

  return(sims)
}


#------------------------------------------------------------
# Examples
#------------------------------------------------------------

# Example 4: ARMA(1,1) with skewed Student’s t-distribution
sim_arma(n = 200, ar = 0.5, ma = -0.3, dist = "sstd", xi = 1.8, df = 10, n_sim = 5)

# Example 5: ARMA(2,1) with GED innovations
sim_arma(n = 250, ar = c(0.3, -0.2), ma = 0.5, dist = "ged", nu = 1.3, n_sim = 3)

