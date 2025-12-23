


#------------------------------------------------------------
# Function: Vasicek.sim()
# Purpose: Simulate multiple paths of the Vasicek
#          mean-reverting short-rate model.
#------------------------------------------------------------
#' Simulate Vasicek Interest Rate Paths
#'
#' @description
#' This function simulates \eqn{N} paths of the Vasicek short-rate model,
#' a mean-reverting stochastic process widely used in interest rate modelling.
#'
#' The Vasicek dynamics are given by:
#' \deqn{dr_t = a(b - r_t)\,dt + \sigma_r\,dW_t}
#' where:
#' \itemize{
#'   \item \eqn{a} is the speed of mean reversion,
#'   \item \eqn{b} is the long-term mean level,
#'   \item \eqn{\sigma_r} is the volatility of the short rate,
#'   \item \eqn{W_t} is a Wiener process.
#' }
#'
#' @param initial.rate Numeric scalar. Initial short rate \eqn{r(0)}.
#' @param a Numeric scalar. Speed of mean reversion.
#' @param b Numeric scalar. Long-term mean level of the short rate.
#' @param sigma_r Numeric scalar. Volatility parameter of the Vasicek model.
#' @param dt Numeric scalar. Time step size (default: 1/250).
#' @param T Numeric scalar. Total time horizon of the simulation (in years).
#' @param N Integer. Number of simulated paths.
#' @param plot Logical. If TRUE (default), plots the simulated short-rate paths.
#' @param use.rainbow Logical. If TRUE, uses rainbow colors for each path;
#' otherwise uses a fixed semi-transparent color.
#' @param fixed.color Character. Color used for all paths when
#' \code{use.rainbow = FALSE}.
#'
#' @return A matrix of dimension \code{num_steps × N} containing the simulated
#' short-rate paths. Each column corresponds to one simulated trajectory.
#'
#' @details
#' The function generates standard normal increments scaled by \eqn{\sqrt{dt}}
#' to approximate Brownian motion.
#' The short rate is updated using the Euler–Maruyama discretization:
#' \deqn{
#' r_{t+dt} = r_t + a(b - r_t)\,dt + \sigma_r \sqrt{dt}\,Z_t
#' }
#' where \eqn{Z_t \sim N(0,1)}.
#'
#' If \code{plot = TRUE}, the function produces:
#' \itemize{
#'   \item all simulated paths,
#'   \item the average path across simulations.
#' }
#'
#' @examples
#' # Simulate 50 Vasicek paths over 1 year
#' out <- Vasicek.sim(initial.rate = 0.02, a = 0.15, b = 0.03,
#'                    sigma_r = 0.01, T = 1, N = 50)
#'
#' # Disable plotting
#' out <- Vasicek.sim(plot = FALSE)
#'
Vasicek.sim <- function(initial.rate = 0.02,
                        a = 0.1,
                        b = 0.03,
                        sigma_r = 0.01,
                        dt = 1 / 250,
                        T = 1,
                        N = 100,
                        plot = TRUE,
                        use.rainbow = FALSE,
                        fixed.color = "gray") {

  # Compute number of steps
  num_steps <- T / dt
  time.axis <- seq(0, T, length.out = num_steps)

  # Preallocate matrices
  dW <- array(rnorm(num_steps * N, 0, sqrt(dt)), dim = c(num_steps, N))
  r <- array(NA, dim = c(num_steps, N))
  r[1, ] <- initial.rate

  # Simulate Vasicek process
  for (i in 2:num_steps) {
    r[i, ] <- r[i - 1, ] + a * (b - r[i - 1, ]) * dt + sigma_r * dW[i - 1, ]
  }

  # Plot
  if (plot) {
    cols <- if (use.rainbow) {
      rainbow(N)
    } else {
      rep(adjustcolor(fixed.color, alpha.f = 0.5), N)
    }

    matplot(time.axis, r, type = 'l', col = cols, lty = 1,
            xlab = "Time (years)", ylab = "Interest Rate",
            main = "Simulated Vasicek Interest Rate Paths")

    # Add mean line
    matlines(time.axis, rowMeans(r), lwd = 2, col = "black", lty = 1)
    legend("topright", legend = c("Mean"), col = "black", lwd = 2, lty = 1, bty = "n")
  }

  return(invisible(r))
}



#------------------------------------------------------------
# Function: CIR.sim()
# Purpose: Simulate multiple paths of the Cox–Ingersoll–Ross
#          (CIR) mean‑reverting short‑rate model, ensuring
#          positivity of interest rates.
#------------------------------------------------------------
#' Simulate Cox–Ingersoll–Ross (CIR) Interest Rate Paths
#'
#' @description
#' This function simulates \eqn{N} paths of the Cox–Ingersoll–Ross (CIR)
#' short‑rate model, a mean‑reverting and strictly positive stochastic
#' process widely used in interest‑rate modelling.
#'
#' The CIR dynamics follow:
#' \deqn{
#' dr_t = a(b - r_t)\,dt + \sigma_r \sqrt{r_t}\,dW_t
#' }
#' where:
#' \itemize{
#'   \item \eqn{a} is the speed of mean reversion,
#'   \item \eqn{b} is the long‑term mean level,
#'   \item \eqn{\sigma_r} is the volatility,
#'   \item \eqn{W_t} is a Wiener process.
#' }
#'
#' @param initial.rate Numeric scalar. Initial short rate \eqn{r(0)}.
#' @param a Numeric scalar. Speed of mean reversion.
#' @param b Numeric scalar. Long‑term mean level of the short rate.
#' @param sigma_r Numeric scalar. Volatility parameter of the CIR model.
#' @param dt Numeric scalar. Time step size (default: 1/250).
#' @param T Numeric scalar. Total time horizon of the simulation (in years).
#' @param N Integer. Number of simulated paths.
#' @param plot Logical. If TRUE (default), plots the simulated CIR paths.
#' @param use.rainbow Logical. If TRUE, uses rainbow colors for each path;
#' otherwise uses a fixed semi‑transparent color.
#' @param fixed.color Character. Color used for all paths when
#' \code{use.rainbow = FALSE}.
#'
#' @return A matrix of dimension \code{num_steps × N} containing the simulated
#' CIR short‑rate paths. Each column corresponds to one simulated trajectory.
#'
#' @details
#' The CIR model ensures non‑negative interest rates due to the
#' \eqn{\sqrt{r_t}} diffusion term.
#' The Euler–Maruyama discretization used here is:
#' \deqn{
#' r_{t+dt} = r_t + a(b - r_t)\,dt + \sigma_r \sqrt{r_t}\sqrt{dt}\,Z_t
#' }
#' with \eqn{Z_t \sim N(0,1)}.
#'
#' To avoid numerical issues (e.g., negative values inside the square root),
#' the implementation applies:
#' \itemize{
#'   \item \code{pmax(r, 0)} before taking the square root,
#'   \item a final \code{pmax()} to enforce positivity.
#' }
#'
#' If \code{plot = TRUE}, the function displays:
#' \itemize{
#'   \item all simulated CIR paths,
#'   \item the average path across simulations.
#' }
#'
#' @examples
#' # Simulate 100 CIR paths over 1 year
#' cir <- CIR.sim(initial.rate = 0.02, a = 0.2, b = 0.04,
#'                sigma_r = 0.015, T = 1, N = 100)
#'
#' # Disable plotting
#' cir <- CIR.sim(plot = FALSE)
#'
CIR.sim <- function(initial.rate = 0.02,
                    a = 0.1,
                    b = 0.03,
                    sigma_r = 0.01,
                    dt = 1 / 250,
                    T = 1,
                    N = 100,
                    plot = TRUE,
                    use.rainbow = FALSE,
                    fixed.color = "gray") {

  # Numero di step
  num_steps <- T / dt
  time.axis <- seq(0, T, length.out = num_steps)

  # Genera incrementi Wiener
  dW <- array(rnorm(num_steps * N, 0, sqrt(dt)), dim = c(num_steps, N))
  r <- array(NA, dim = c(num_steps, N))
  r[1, ] <- initial.rate

  # Simulazione CIR
  for (i in 2:num_steps) {
    r_prev_sqrt <- sqrt(pmax(r[i - 1, ], 0))
    r[i, ] <- r[i - 1, ] +
      a * (b - r[i - 1, ]) * dt +
      sigma_r * r_prev_sqrt * dW[i - 1, ]
    r[i, ] <- pmax(r[i, ], 0)
  }

  # Plot
  if (plot) {
    cols <- if (use.rainbow) rainbow(N)
    else rep(adjustcolor(fixed.color, alpha.f = 0.5), N)

    matplot(time.axis, r, type = 'l', col = cols, lty = 1,
            xlab = "Time (years)", ylab = "Interest Rate",
            main = "Simulated CIR Interest Rate Paths")

    matlines(time.axis, rowMeans(r), lwd = 2, col = "black", lty = 1)
    legend("topright", legend = c("Mean"), col = "black",
           lwd = 2, lty = 1, bty = "n")
  }

  return(invisible(r))
}















#------------------------------------------------------------
# Example simulations
#------------------------------------------------------------

set.seed(42)

vas <- Vasicek.sim(
  initial.rate = 0.02,
  a = 0.1,
  b = 0.03,
  sigma_r = 0.01,
  dt = 1/250,
  T = 1,
  N = 200,
  plot = TRUE,
  use.rainbow = FALSE,
  fixed.color = "darkblue"
)


set.seed(42)

cir <- CIR.sim(
  initial.rate = 0.02,
  a = 0.1,
  b = 0.03,
  sigma_r = 0.01,
  dt = 1/250,
  T = 1,
  N = 300,
  plot = TRUE,
  use.rainbow = FALSE,
  fixed.color = "darkgreen"
)















