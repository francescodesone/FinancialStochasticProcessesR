#------------------------------------------------------------
# Function: GBM.sim()
# Purpose: Simulate correlated Geometric Brownian Motion (GBM)
#          for two assets with specified drift, volatility,
#          and correlation.
#------------------------------------------------------------
GBM.sim <- function(inV, r, sigma.par, rho.par,
                    expiry, num.time.steps, num.paths,
                    plot = TRUE, use.rainbow = TRUE, fixed.colors = c("blue", "darkgreen")) {

  if (length(inV) != 2 || length(r) != 2 || length(sigma.par) != 2) {
    stop("inV, r, and sigma.par must all be vectors of length 2.")
  }

  dt <- expiry / num.time.steps

  # Generate standard Brownian increments
  B.X <- matrix(rnorm(num.time.steps * num.paths), nrow = num.paths)
  B.Z <- matrix(rnorm(num.time.steps * num.paths), nrow = num.paths)

  # Create correlated Brownian increments
  xi.X <- B.X
  xi.Z <- rho.par * B.X + sqrt(1 - rho.par^2) * B.Z

  # Log-price increments
  dX <- r[1] * dt + sigma.par[1] * xi.X * sqrt(dt)
  dZ <- r[2] * dt + sigma.par[2] * xi.Z * sqrt(dt)

  # Log-price paths
  logX <- t(apply(cbind(rep(log(inV[1]), num.paths), dX), 1, cumsum))
  logZ <- t(apply(cbind(rep(log(inV[2]), num.paths), dZ), 1, cumsum))

  # Price paths
  X <- exp(logX)
  Z <- exp(logZ)

  # Theoretical mean paths
  time.axis <- seq(0, expiry, length.out = num.time.steps + 1)
  mean_X <- inV[1] * exp(r[1] * time.axis)
  mean_Z <- inV[2] * exp(r[2] * time.axis)

  # Plot simulated paths if requested
  if (plot) {
    old.par <- par(mfrow = c(2, 1))
    on.exit(par(old.par))

    cols_X <- if (use.rainbow) rainbow(num.paths, alpha = 0.5)
    else rep(adjustcolor(fixed.colors[1], alpha.f = 0.5), num.paths)
    matplot(time.axis, t(X), type = "l", col = cols_X, lty = 1,
            xlab = "Time", ylab = "X price", main = "Simulated GBM paths: X")
    lines(time.axis, mean_X, col = "red", lwd = 2)
    legend("topright", legend = "Mean", col = "red", lty = 1, lwd = 2, bty = "n")

    cols_Z <- if (use.rainbow) rainbow(num.paths, alpha = 0.5)
    else rep(adjustcolor(fixed.colors[2], alpha.f = 0.5), num.paths)
    matplot(time.axis, t(Z), type = "l", col = cols_Z, lty = 1,
            xlab = "Time", ylab = "Z price", main = "Simulated GBM paths: Z")
    lines(time.axis, mean_Z, col = "red", lwd = 2)
    legend("topright", legend = "Mean", col = "red", lty = 1, lwd = 2, bty = "n")
  }

  return(list(X = X, Z = Z))
}

#------------------------------------------------------------
# Function: ABM.sim()
# Purpose: Simulate Arithmetic Brownian Motion (ABM) paths
#          with constant drift and volatility.
#------------------------------------------------------------
ABM.sim <- function(inV, r, sigma.par,
                    expiry, num.time.steps, num.paths,
                    plot = TRUE, use.rainbow = TRUE, fixed.color = "darkred") {

  if (length(inV) != 1 || length(r) != 1 || length(sigma.par) != 1) {
    stop("inV, r, and sigma.par must be scalars (length 1).")
  }

  dt <- expiry / num.time.steps

  dX <- r * dt + sigma.par * matrix(rnorm(num.time.steps * num.paths),
                                    nrow = num.paths) * sqrt(dt)

  X <- t(apply(cbind(rep(inV, num.paths), dX), 1, cumsum))

  time.axis <- seq(0, expiry, length.out = num.time.steps + 1)
  mean_path <- inV + r * time.axis

  if (plot) {
    cols <- if (use.rainbow) rainbow(num.paths, alpha = 0.5)
    else rep(adjustcolor(fixed.color, alpha.f = 0.5), num.paths)
    matplot(time.axis, t(X), type = "l", col = cols, lty = 1,
            xlab = "Time", ylab = "Value", main = "Simulated ABM Paths")
    lines(time.axis, mean_path, col = "red", lwd = 2)
    legend("topright", legend = "Mean", col = "red", lty = 1, lwd = 2, bty = "n")
  }

  return(X)
}

#------------------------------------------------------------
# Function: BB.sim()
# Purpose: Simulate Brownian Bridge paths starting at x0
#          and ending at y over [t0, T].
#------------------------------------------------------------
BB.sim <- function(N = 1000, M = 1, x0 = 0, y = 0,
                   t0 = 0, T = 1, Dt = NULL,
                   plot = TRUE, use.rainbow = TRUE, fixed.color = "steelblue") {

  if (!is.numeric(x0) || !is.numeric(y)) stop("'x0' and 'y' must be numeric.")
  if (!is.numeric(t0) || !is.numeric(T) || T <= t0) stop("Invalid times: ensure that 0 <= t0 < T.")
  if (!is.numeric(N) || N <= 1 || N != floor(N)) stop("'N' must be a positive integer > 1.")
  if (!is.numeric(M) || M <= 0 || M != floor(M)) stop("'M' must be a positive integer.")

  if (is.null(Dt)) Dt <- (T - t0) / N

  time.axis <- seq(t0, T, by = Dt)
  if (length(time.axis) != N + 1) {
    time.axis <- seq(t0, T, length.out = N + 1)
    Dt <- diff(time.axis[1:2])
  }

  simulate_bridge <- function() {
    w <- c(0, cumsum(rnorm(N, mean = 0, sd = sqrt(Dt))))
    x <- x0 + w - ((time.axis - t0) / (T - t0)) * (w[N + 1] - y + x0)
    return(x)
  }

  paths <- replicate(M, simulate_bridge())
  colnames(paths) <- paste0("Path_", 1:M)
  mean_path <- x0 + ((time.axis - t0) / (T - t0)) * (y - x0)

  if (plot) {
    colors <- if (use.rainbow) rainbow(M, alpha = 0.8)
    else rep(adjustcolor(fixed.color, alpha.f = 0.6), M)
    old.par <- par(no.readonly = TRUE)
    on.exit(par(old.par))
    layout(matrix(1, ncol = 1))
    matplot(time.axis, paths, type = "l", lty = 1, col = colors,
            xlab = "Time", ylab = "Value", main = "Simulated Brownian Bridges")
    lines(time.axis, mean_path, col = "red", lwd = 2)
    legend("topright", legend = "Mean", col = "red", lty = 1, lwd = 2, bty = "n")
  }

  return(list(time = time.axis, paths = paths, mean = mean_path))
}

#------------------------------------------------------------
# Function: GBM.single()
# Purpose: Simulate a single-asset Geometric Brownian Motion (GBM)
#------------------------------------------------------------
GBM.single <- function(S0 = 100, mu = 0.05, sigma = 0.2,
                       expiry = 1, num.time.steps = 252,
                       num.paths = 10, plot = TRUE,
                       use.rainbow = TRUE, fixed.color = "blue") {

  dt <- expiry / num.time.steps
  epsilon <- matrix(rnorm(num.time.steps * num.paths), nrow = num.paths)
  dS <- (mu - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * epsilon
  logS <- t(apply(cbind(rep(log(S0), num.paths), dS), 1, cumsum))
  S <- exp(logS)

  time.axis <- seq(0, expiry, length.out = num.time.steps + 1)
  mean_path <- S0 * exp(mu * time.axis)

  if (plot) {
    cols <- if (use.rainbow) rainbow(num.paths, alpha = 0.5)
    else rep(adjustcolor(fixed.color, alpha.f = 0.5), num.paths)
    matplot(time.axis, t(S), type = "l", col = cols, lty = 1,
            xlab = "Time", ylab = "Asset Price", main = "Simulated GBM Paths")
    lines(time.axis, mean_path, col = "red", lwd = 2)
    legend("topright", legend = "Mean", col = "red", lty = 1, lwd = 2, bty = "n")
  }

  return(S)
}

#------------------------------------------------------------
# Example simulations
#------------------------------------------------------------

set.seed(42)
bb <- BB.sim(N = 500, M = 100, x0 = 0, y = 1, T = 1,
             plot = TRUE, use.rainbow = FALSE, fixed.color = "darkred")

set.seed(42)
ABM.sim(inV = 100, r = 0.03, sigma.par = 0.1,
        expiry = 1, num.time.steps = 252, num.paths = 100,
        plot = TRUE, use.rainbow = FALSE, fixed.color = "darkred")

set.seed(123)
sim <- GBM.sim(inV = c(100, 50), r = c(0.03, 0.04),
               sigma.par = c(0.2, 0.25), rho.par = 0.5,
               expiry = 1, num.time.steps = 250, num.paths = 20,
               use.rainbow = FALSE)

set.seed(123)
gbm_single_example <- GBM.single(
  S0 = 100,
  mu = 0.05,
  sigma = 0.2,
  expiry = 1,
  num.time.steps = 252,
  num.paths = 10,
  plot = TRUE,
  use.rainbow = FALSE,
  fixed.color = "blue"
)
