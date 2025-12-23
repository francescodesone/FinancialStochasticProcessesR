#------------------------------------------------------------
# Function: CEVProcess.single()
# Purpose: Simulate multiple paths of the Constant Elasticity
#          of Variance (CEV) process using Euler–Maruyama.
#------------------------------------------------------------
#' Simulate Constant Elasticity of Variance (CEV) Process
#'
#' This function simulates \code{N} trajectories of the
#' **Constant Elasticity of Variance (CEV)** stochastic process
#' using Euler–Maruyama discretization.
#'
#' The CEV model generalizes Geometric Brownian Motion (GBM) by
#' allowing volatility to depend on the level of the process:
#' \deqn{
#' dX_t = \mu X_t\,dt + \sigma X_t^\gamma\, dW_t,
#' }
#' where:
#' \itemize{
#'   \item \eqn{\mu} is the drift coefficient,
#'   \item \eqn{\sigma} is the volatility scale,
#'   \item \eqn{\gamma} is the elasticity parameter,
#'   \item \eqn{W_t} is a Wiener process.
#' }
#'
#' Special cases:
#' \itemize{
#'   \item \eqn{\gamma = 1}: Geometric Brownian Motion (GBM),
#'   \item \eqn{\gamma = 0.5}: Square-root diffusion (CIR-like),
#'   \item \eqn{\gamma > 1}: Volatility grows faster than linearly.
#' }
#'
#' @param mu Numeric scalar. Drift coefficient.
#' @param sigma Numeric scalar. Volatility scale.
#' @param gamma Numeric scalar. Elasticity parameter.
#' @param x0 Numeric scalar. Initial value \eqn{X(0)}.
#' @param T Numeric scalar. Total time horizon of the simulation.
#' @param n Integer. Number of time steps.
#' @param N Integer. Number of simulated paths.
#' @param plot Logical. If TRUE (default), plots the simulated paths.
#' @param envelope Logical. If TRUE, adds theoretical expectation and
#'   confidence bands (only valid for \eqn{\gamma = 1}).
#' @param use.rainbow Logical. If TRUE, uses rainbow colors for each path;
#' otherwise uses a fixed semi-transparent color.
#' @param fixed.color Character. Color used when \code{use.rainbow = FALSE}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{paths}: Matrix of simulated paths (\code{N × (n+1)}).
#'   \item \code{times}: Time grid.
#'   \item \code{expectation}: Theoretical expectation (only for \eqn{\gamma = 1}).
#'   \item \code{variance}: Theoretical variance (only for \eqn{\gamma = 1}).
#'   \item \code{lower}, \code{upper}: Confidence bands (only for \eqn{\gamma = 1}).
#' }
#'
#' @details
#' The simulation uses Euler–Maruyama discretization:
#' \deqn{
#' X_{t+dt} = X_t + \mu X_t\,dt + \sigma X_t^\gamma \sqrt{dt}\,Z_t,
#' }
#' with \eqn{Z_t \sim N(0,1)}.
#'
#' For \eqn{\gamma = 1} (GBM), closed-form expressions for expectation
#' and variance are available:
#' \deqn{E[X_t] = x_0 e^{\mu t},}
#' \deqn{Var[X_t] = x_0^2 e^{2\mu t}\left(e^{\sigma^2 t} - 1\right).}
#'
#' @examples
#' #------------------------------------------------------------
#' # Example simulations — CEV process
#' #------------------------------------------------------------
#'
#' # Case 1: Generic gamma (no theoretical envelope)
#' set.seed(42)
#' CEVProcess.single(mu = 0.5, sigma = 0.2, gamma = 1.25,
#'                   x0 = 1, T = 1, n = 200, N = 20)
#'
#' # Case 2: Gamma = 1 (GBM, theoretical envelope active)
#' set.seed(42)
#' CEVProcess.single(mu = 0.1, sigma = 0.3, gamma = 1,
#'                   x0 = 1, T = 1, n = 200, N = 100)
#'
#' @export
CEVProcess.single <- function(mu = 0.5,
                              sigma = 0.1,
                              gamma = 1.5,
                              x0 = 1.0,
                              T = 1.0,
                              n = 100,
                              N = 10,
                              plot = TRUE,
                              use.rainbow = TRUE,
                              fixed.color = "blue") {

  dt <- T / n
  time.grid <- seq(0, T, length.out = n + 1)

  # Matrice dei cammini
  paths <- matrix(NA, nrow = N, ncol = n + 1)
  paths[, 1] <- x0

  # Simulazione con schema di Eulero-Maruyama
  for (i in 1:N) {
    for (j in 1:n) {
      X_prev <- paths[i, j]
      dW <- rnorm(1, mean = 0, sd = sqrt(dt))
      X_next <- X_prev + mu * X_prev * dt + sigma * X_prev^gamma * dW
      paths[i, j + 1] <- max(X_next, 0)  # Evita valori negativi
    }
  }

  # Expectation teorica solo per gamma = 1 (GBM)
  if (gamma == 1) {
    expectation <- x0 * exp(mu * time.grid)
  } else {
    expectation <- rep(NA, n + 1)
  }

  # Plotting
  if (plot) {
    cols <- if (use.rainbow) {
      rainbow(N)
    } else {
      rep(adjustcolor(fixed.color, alpha.f = 0.6), N)
    }

    matplot(time.grid, t(paths), type = "l", lty = 1, col = cols,
            xlab = "Time", ylab = "X(t)",
            main = paste0("CEV Process: μ=", mu,
                          ", σ=", sigma, ", γ=", gamma))

    # Linea della media teorica (solo per gamma = 1)
    if (gamma == 1) {
      lines(time.grid, expectation, col = "black", lwd = 2, lty = 2)
      legend("topleft", legend = "Mean",
             col = "black", lty = 2, lwd = 2, bty = "n")
    }
  }

  # Output
  return(list(paths = paths,
              times = time.grid,
              expectation = expectation))
}












#------------------------------------------------------------
# Example simulations — CEV Process
#------------------------------------------------------------

set.seed(42)

# Case 1: Gamma generic
cev_paths <- CEVProcess.single(
  mu = 0.5,
  sigma = 0.2,
  gamma = 1.25,
  x0 = 1,
  T = 1,
  n = 200,
  N = 20,
  plot = TRUE,
  use.rainbow = FALSE,
  fixed.color = "darkblue"
)

# Case 2: Gamma = 1 (GBM)
cev_paths_gbm <- CEVProcess.single(
  mu = 0.1,
  sigma = 0.3,
  gamma = 1,
  x0 = 1,
  T = 1,
  n = 200,
  N = 100,
  plot = TRUE,
  use.rainbow = TRUE
)





