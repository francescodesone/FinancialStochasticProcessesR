#------------------------------------------------------------
# Function: CKLS.sim()
# Purpose: Simulate multiple paths of the CKLS stochastic
#          interest-rate model using Euler–Maruyama.
#------------------------------------------------------------
#' Simulate CKLS Stochastic Process (Multiple Paths)
#'
#' This function simulates \code{N} trajectories of the
#' **Chan–Karolyi–Longstaff–Sanders (CKLS)** stochastic process,
#' a flexible model that generalizes Vasicek and CIR dynamics.
#'
#' The CKLS model is defined by the stochastic differential equation:
#' \deqn{
#' dX_t = (\alpha + \beta X_t)\,dt + \sigma X_t^\gamma\, dW_t,
#' }
#' where:
#' \itemize{
#'   \item \eqn{\alpha} controls the constant drift component,
#'   \item \eqn{\beta} controls the linear drift (mean-reversion or growth),
#'   \item \eqn{\sigma} is the volatility scale,
#'   \item \eqn{\gamma} determines how volatility depends on the level of the process,
#'   \item \eqn{W_t} is a Wiener process.
#' }
#'
#' Special cases include:
#' \itemize{
#'   \item Vasicek model: \eqn{\gamma = 0},
#'   \item CIR model: \eqn{\gamma = 1/2},
#'   \item Geometric Brownian Motion: \eqn{\gamma = 1}.
#' }
#'
#' @param initial.rate Numeric scalar. Initial value \eqn{X(0)}.
#' @param alpha Numeric scalar. Constant drift parameter.
#' @param beta Numeric scalar. Linear drift parameter.
#' @param sigma Numeric scalar. Volatility scale.
#' @param gamma Numeric scalar. Elasticity of volatility.
#' @param dt Numeric scalar. Time step size (default: 1/250).
#' @param T Numeric scalar. Total time horizon of the simulation.
#' @param N Integer. Number of simulated paths.
#' @param plot Logical. If TRUE (default), plots the simulated CKLS paths.
#' @param use.rainbow Logical. If TRUE, uses rainbow colors for each path;
#' otherwise uses a fixed semi-transparent color.
#' @param fixed.color Character. Color used when \code{use.rainbow = FALSE}.
#'
#' @return A matrix of dimension \code{num_steps × N} containing the simulated
#' CKLS paths. Each column corresponds to one trajectory.
#'
#' @details
#' The simulation uses the Euler–Maruyama discretization:
#' \deqn{
#' X_{t+dt} = X_t + (\alpha + \beta X_t)\,dt
#'            + \sigma X_t^\gamma \sqrt{dt}\,Z_t,
#' }
#' where \eqn{Z_t \sim N(0,1)}.
#'
#' Note that the CKLS model does **not** guarantee positivity unless
#' \eqn{\gamma \ge 1/2} and certain parameter conditions hold.
#'
#' If \code{plot = TRUE}, the function displays:
#' \itemize{
#'   \item all simulated CKLS paths,
#'   \item the average path across simulations.
#' }
#'
#' @examples
#' #------------------------------------------------------------
#' # Example simulations — CKLS model
#' #------------------------------------------------------------
#'
#' set.seed(42)
#'
#' ckls_paths <- CKLS.sim(
#'   initial.rate = 0.02,
#'   alpha = 0.5,
#'   beta = 0.5,
#'   sigma = 0.1,
#'   gamma = 1.5,
#'   dt = 1/250,
#'   T = 1,
#'   N = 50,
#'   plot = TRUE,
#'   use.rainbow = FALSE,
#'   fixed.color = "blue"
#' )
#'
CKLS.sim <- function(initial.rate = 0.02,
                     alpha = 0.5,
                     beta = 0.5,
                     sigma = 0.1,
                     gamma = 1.5,
                     dt = 1 / 250,
                     T = 1,
                     N = 100,
                     plot = TRUE,
                     use.rainbow = FALSE,
                     fixed.color = "blue") {

  # Numero di step temporali
  num_steps <- T / dt
  time.axis <- seq(0, T, length.out = num_steps)

  # Preallocazione matrici
  dW <- matrix(rnorm(num_steps * N, 0, sqrt(dt)), nrow = num_steps, ncol = N)
  X <- matrix(NA, nrow = num_steps, ncol = N)
  X[1, ] <- initial.rate

  # Simulazione con schema di Eulero-Maruyama
  for (i in 2:num_steps) {
    X_prev <- X[i - 1, ]
    drift <- (alpha + beta * X_prev) * dt
    diffusion <- sigma * (X_prev^gamma) * dW[i - 1, ]
    X[i, ] <- X_prev + drift + diffusion
  }

  # Plot
  if (plot) {
    cols <- if (use.rainbow) {
      rainbow(N)
    } else {
      rep(adjustcolor(fixed.color, alpha.f = 0.5), N)
    }

    matplot(time.axis, X, type = "l", col = cols, lty = 1,
            xlab = "Time (years)", ylab = "Process Value",
            main = "Simulated CKLS Process Paths")

    # Linea media
    matlines(time.axis, rowMeans(X), lwd = 2, col = "black", lty = 1)
    legend("topright", legend = "Mean", col = "black", lwd = 2, lty = 1, bty = "n")
  }

  return(invisible(X))
}



#------------------------------------------------------------
# Example simulations
#------------------------------------------------------------

set.seed(123)

paths <- CKLS.sim(
  initial.rate = 0.02,
  alpha = 0.01,
  beta = 0.05,
  sigma = 0.1,
  gamma = 0,
  dt = 1 / 252,
  T = 2,
  N = 30,
  plot = TRUE,
  use.rainbow = TRUE
)
