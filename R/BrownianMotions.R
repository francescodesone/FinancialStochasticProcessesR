#------------------------------------------------------------
# Function: GBM.sim()
# Purpose: Simulate correlated Geometric Brownian Motion (GBM)
#          for two assets with specified drift, volatility,
#          and correlation.
#------------------------------------------------------------
#' Simulate Correlated Geometric Brownian Motions (GBM)
#'
#' @param inV Numeric vector of length 2. Initial asset prices \eqn{S_1(0)} and \eqn{S_2(0)}.
#' @param r Numeric vector of length 2. Drift rates (expected returns) for each asset.
#' @param sigma.par Numeric vector of length 2. Volatility parameters for the two assets.
#' @param rho.par Numeric scalar in [-1, 1]. Correlation coefficient between the two Brownian motions.
#' @param expiry Numeric scalar. Time horizon of the simulation (e.g., in years).
#' @param num.time.steps Integer. Number of time steps in each simulated path.
#' @param num.paths Integer. Number of simulated paths.
#' @param plot Logical. If TRUE (default), plots the simulated GBM paths for both assets.
#' @param use.rainbow Logical. If TRUE, uses rainbow colors for each simulated path;
#' otherwise, uses the fixed colors provided in \code{fixed.colors}.
#' @param fixed.colors Character vector of length 2. Fixed colors for the two assets if \code{use.rainbow = FALSE}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{X}{Matrix of simulated paths for the first asset (\eqn{S_1(t)}).}
#'   \item{Z}{Matrix of simulated paths for the second asset (\eqn{S_2(t)}).}
#' }
#'
#' @details
#' The function generates correlated normal increments using the correlation parameter \eqn{\rho}.
#' Each asset price evolves according to its drift and volatility parameters.
#' The correlation is introduced through a linear combination of two independent standard normal variables:
#' \deqn{\xi_Z = \rho \xi_X + \sqrt{1 - \rho^2} Z}
#' ensuring the desired correlation between the two Wiener processes.
#'
#' @examples
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
#' Simulate Arithmetic Brownian Motion (ABM) Paths
#'
#' @param inV Numeric scalar. Initial value of the process \eqn{X(0)}.
#' @param r Numeric scalar. Drift parameter \eqn{\mu}, representing the expected rate of change.
#' @param sigma.par Numeric scalar. Volatility parameter \eqn{\sigma}, controlling the diffusion magnitude.
#' @param expiry Numeric scalar. Total simulation time horizon (e.g., in years).
#' @param num.time.steps Integer. Number of discrete time steps in each simulated path.
#' @param num.paths Integer. Number of independent paths to simulate.
#' @param plot Logical. If TRUE (default), produces a plot of the simulated paths.
#' @param use.rainbow Logical. If TRUE, colors each simulated path differently using rainbow colors;
#' otherwise, all paths use the color specified in \code{fixed.color}.
#' @param fixed.color Character string. Default color for all paths when \code{use.rainbow = FALSE}.
#'
#' @return A matrix where each row corresponds to a simulated path and each column to a time step.
#' Columns are ordered chronologically from \eqn{t = 0} to \eqn{t = T}.
#'
#' @details
#' The ABM model assumes additive noise, meaning that the increments are normally distributed
#' and independent across time:
#' \deqn{X_{t+\Delta t} - X_t = \mu \, \Delta t + \sigma \sqrt{\Delta t} Z, \quad Z \sim N(0,1).}
#' The simulation is performed via Euler–Maruyama discretization, and the expected value evolves linearly as:
#' \deqn{E[X_t] = X_0 + \mu t.}
#'
#' @examples
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
#' Simulate Brownian Bridge (BB) Paths
#'
#' This function simulates one or more realizations of a **Brownian Bridge** process
#' over a finite time interval \eqn{[t_0, T]}, starting at value \eqn{x_0} and
#' conditioned to end at value \eqn{y} at time \eqn{T}.
#'
#' A Brownian Bridge can be seen as a standard Brownian Motion \eqn{W_t}
#' constrained to return to a specified endpoint, with the stochastic differential form:
#' \deqn{dX_t = \frac{y - X_t}{T - t} dt + dW_t,}
#' where \eqn{W_t} denotes a standard Wiener process.
#'
#' @param N Integer. Number of time steps for each simulated path.
#' @param M Integer. Number of Brownian Bridge paths to simulate.
#' @param x0 Numeric. Starting value of the process at time \eqn{t_0}.
#' @param y Numeric. Terminal value of the process at time \eqn{T}.
#' @param t0 Numeric. Initial time of the process (default = 0).
#' @param T Numeric. Terminal time of the process (default = 1).
#' @param Dt Numeric (optional). Time step size; if \code{NULL}, it is computed as \eqn{(T - t_0)/N}.
#' @param plot Logical. If TRUE (default), plots the simulated Brownian Bridge paths.
#' @param use.rainbow Logical. If TRUE (default), assigns a rainbow color palette to each path;
#' otherwise, all paths are drawn using \code{fixed.color}.
#' @param fixed.color Character string. Color to use for all paths if \code{use.rainbow = FALSE}.
#'
#' @return A list with the following components:
#' \item{time}{A numeric vector of equally spaced time points between \eqn{t_0} and \eqn{T}.}
#' \item{paths}{A matrix of simulated paths, with each column corresponding to a single realization.}
#' \item{mean}{A numeric vector representing the theoretical expected value of the Brownian Bridge at each time step.}
#'
#' @details
#' The Brownian Bridge is constructed from a standard Brownian Motion by conditioning on
#' the terminal value \eqn{X_T = y}. Given \eqn{W_t}, the bridge can be expressed as:
#' \deqn{X_t = x_0 + W_t - \frac{t - t_0}{T - t_0}(W_T - y + x_0).}
#'
#' This implementation follows this construction directly, ensuring that each path
#' starts at \eqn{x_0} and ends at \eqn{y}. The expected value of the process evolves linearly as:
#' \deqn{E[X_t] = x_0 + \frac{t - t_0}{T - t_0}(y - x_0),}
#' while the variance decreases as the process approaches \eqn{T}, reflecting the constraint at the endpoint.
#'
#' @examples
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
#' Simulate Single-Asset Geometric Brownian Motion (GBM)
#'
#' This function simulates one or more sample paths of a **Geometric Brownian Motion (GBM)**
#' for a single financial asset. GBM is one of the fundamental models used to describe
#' the stochastic evolution of asset prices in continuous time under the assumption of
#' constant drift and volatility.
#'
#' The process follows the stochastic differential equation (SDE):
#' \deqn{dS_t = \mu S_t \, dt + \sigma S_t \, dW_t,}
#' where:
#' \itemize{
#'   \item \eqn{S_t} is the asset price at time \eqn{t},
#'   \item \eqn{\mu} is the drift (expected return),
#'   \item \eqn{\sigma} is the volatility,
#'   \item and \eqn{W_t} is a standard Wiener process.
#' }
#'
#' The analytical solution of the SDE is:
#' \deqn{S_t = S_0 \exp \left[ \left( \mu - \frac{1}{2}\sigma^2 \right)t + \sigma W_t \right].}
#'
#' @param S0 Numeric. Initial asset price (default 100).
#' @param mu Numeric. Drift coefficient (expected return) of the asset (default 0.05).
#' @param sigma Numeric. Volatility of the asset (default 0.2).
#' @param expiry Numeric. Time horizon for simulation in years (default 1).
#' @param num.time.steps Integer. Number of discrete time steps (default 252).
#' @param num.paths Integer. Number of simulated paths (default 10).
#' @param plot Logical. If TRUE, plots the simulated paths (default TRUE).
#' @param use.rainbow Logical. If TRUE, uses rainbow colors for paths; otherwise uses \code{fixed.color} (default TRUE).
#' @param fixed.color Character. Fixed color to use when \code{use.rainbow = FALSE} (default "blue").
#'
#' @return
#' A numeric matrix of simulated asset prices with dimensions:
#' \itemize{
#'   \item Rows = number of simulated paths (\code{num.paths}),
#'   \item Columns = number of time points (\code{num.time.steps + 1}).
#' }
#' Each row represents a simulated price trajectory of the asset.
#'
#' @details
#' The simulation is based on the discretized version of the GBM model under the **Euler–Maruyama scheme**:
#' \deqn{S_{t+\Delta t} = S_t \exp \left[ \left( \mu - \frac{1}{2}\sigma^2 \right)\Delta t + \sigma \sqrt{\Delta t} Z \right],}
#' where \eqn{Z \sim N(0,1)} are independent standard normal random variables.
#'
#' The function optionally plots all simulated paths along with the **theoretical mean path**
#' \eqn{E[S_t] = S_0 e^{\mu t}} shown in red.
#'
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
