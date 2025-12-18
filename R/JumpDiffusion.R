#------------------------------------------------------------
# Function: PoissonJump()
# Purpose: Simulate multiple paths of a Poisson jump process
#          with specified jump size, intensity, and horizon.
#------------------------------------------------------------
#' Simulate Poisson Jump Process Paths
#'
#' @param X0 Numeric scalar. Initial value of the process at time \eqn{t = 0}.
#' @param g Numeric scalar. Jump size (magnitude of each upward jump).
#' @param phi Numeric scalar. Intensity (rate) of the Poisson process, i.e. expected number of jumps per unit time.
#' @param T Numeric scalar. Time horizon of the simulation (e.g., in years).
#' @param dt Numeric scalar. Time step size for discretization.
#' @param N Integer. Number of simulated paths.
#' @param plot Logical. If TRUE (default), plots the simulated Poisson jump paths.
#' @param use.rainbow Logical. If TRUE (default), uses rainbow colors for each simulated path;
#' otherwise, uses the fixed color provided in \code{fixed.color}.
#' @param fixed.color Character scalar. Fixed color for all paths if \code{use.rainbow = FALSE}.
#'
#' @return A numeric matrix of dimension \code{(time_steps + 1) x N}, where each column
#' represents one simulated path of the Poisson jump process.
#'
#' @details
#' The function simulates increments of a Poisson process with intensity \eqn{\phi}.
#' At each time step of length \eqn{dt}, the number of jumps is drawn from a Poisson distribution
#' with mean \eqn{\phi \cdot dt}. Each jump increases the process by \eqn{g}.
#' Paths are built cumulatively over the time horizon \eqn{T}.
#'
#' If plotting is enabled, all simulated paths are displayed along with the mean path
#' (shown in thick black). Colors can be customized via rainbow or fixed options.
#'
#' @examples
#' # Simulate 50 paths of a Poisson jump process
#' set.seed(123)
#' X <- PoissonJump(X0 = 0, g = 1, phi = 2, T = 2,
#'                  dt = 1/100, N = 50, plot = TRUE,
#'                  use.rainbow = FALSE, fixed.color = "red")
#'
#' # Inspect the first few rows of the simulated matrix
#' head(X[, 1:5])
#'
#' # Compute mean path without plotting
#' mean_path <- rowMeans(PoissonJump(plot = FALSE))
#'
#' @export
PoissonJump <- function(X0 = 0, g = 1, phi = 2, T = 5, dt = 1/250,
                        N = 100, plot = TRUE, use.rainbow = TRUE,
                        fixed.color = "blue") {

  time_steps <- T / dt
  time_axis <- seq(0, T, length.out = time_steps + 1)

  # Pre-allocate matrix for jump process paths
  X <- matrix(NA, nrow = time_steps + 1, ncol = N)
  X[1, ] <- X0

  # Generate increments of Poisson process (number of jumps per dt)
  dPi <- matrix(rpois(time_steps * N, lambda = phi * dt), nrow = time_steps, ncol = N)

  # Build paths cumulatively
  for (i in 2:(time_steps + 1)) {
    X[i, ] <- X[i - 1, ] + g * dPi[i - 1, ]
  }

  # Plot
  if (plot) {
    cols <- if (use.rainbow) {
      rainbow(N)
    } else {
      rep(adjustcolor(fixed.color, alpha.f = 0.5), N)
    }

    matplot(time_axis, X, type = "l", col = cols, lty = 1,
            xlab = "Time (years)", ylab = "Process Value",
            main = "Simulated Poisson Jump Paths")

    # Plot mean path thicker in black
    mean_path <- rowMeans(X)
    matlines(time_axis, mean_path, col = "black", lwd = 2)
  }

  return(X)
}

# Example usage:
set.seed(42)
PoissonJump(X0 = 0, g = 1, phi = 2, T = 5, dt = 1/250, N = 100)




#------------------------------------------------------------
# Function: JumpDiffusion.sim()
# Purpose: Simulate multiple paths of a Jump-Diffusion process
#          combining Brownian motion with Poisson jumps.
#------------------------------------------------------------
#' Simulate Jump-Diffusion Process Paths
#'
#' @param X0 Numeric scalar. Initial value of the process at time \eqn{t = 0}.
#' @param g Numeric scalar. Jump size (relative impact of each jump, e.g. -0.2 for a 20% drop).
#' @param phi Numeric scalar. Intensity (rate) of the Poisson process, i.e. expected number of jumps per unit time.
#' @param mu Numeric scalar. Drift parameter of the diffusion component.
#' @param sigma Numeric scalar. Volatility parameter of the diffusion component.
#' @param T Numeric scalar. Time horizon of the simulation (e.g., in years).
#' @param dt Numeric scalar. Time step size for discretization.
#' @param N Integer. Number of simulated paths.
#' @param plot Logical. If TRUE (default), plots the simulated jump-diffusion paths.
#' @param use.rainbow Logical. If TRUE (default), uses rainbow colors for each simulated path;
#' otherwise, uses the fixed color provided in \code{fixed.color}.
#' @param fixed.color Character scalar. Fixed color for all paths if \code{use.rainbow = FALSE}.
#'
#' @return A list with two components:
#' \describe{
#'   \item{time}{Numeric vector of time points.}
#'   \item{paths}{Numeric matrix of dimension \code{n.steps x N}, where each column
#'   represents one simulated path of the jump-diffusion process.}
#' }
#'
#' @details
#' The jump-diffusion process combines a continuous diffusion (Brownian motion) with
#' discrete jumps governed by a Poisson process. The dynamics are given by:
#' \deqn{dX_t = X_t \left[ (\mu - \phi g) dt + \sigma dW_t + g d\Pi_t \right]}
#' where:
#' \itemize{
#'   \item \eqn{\mu} is the drift rate,
#'   \item \eqn{\sigma} is the volatility,
#'   \item \eqn{dW_t} are increments of a Wiener process,
#'   \item \eqn{d\Pi_t} are increments of a Poisson process with intensity \eqn{\phi},
#'   \item \eqn{g} is the relative jump size.
#' }
#'
#' The term \eqn{(\mu - \phi g)} adjusts the drift to account for the expected jump impact.
#'
#' If plotting is enabled, all simulated paths are displayed along with the mean path
#' (shown in thick black). Colors can be customized via rainbow or fixed options.
#'
#' @examples
#' # Simulate 50 paths of a jump-diffusion process
#' set.seed(123)
#' sim <- JumpDiffusion.sim(X0 = 1, g = -0.2, phi = 1,
#'                          mu = 0.05, sigma = 0.2,
#'                          T = 2, dt = 1/100, N = 50,
#'                          plot = TRUE, use.rainbow = FALSE,
#'                          fixed.color = "darkgreen")
#'
#' # Inspect the first few simulated paths
#' head(sim$paths[, 1:5])
#'
#' # Compute mean path without plotting
#' mean_path <- rowMeans(JumpDiffusion.sim(plot = FALSE)$paths)
#'
#' @export
JumpDiffusion.sim <- function(X0 = 1, g = -0.2, phi = 1,
                              mu = 0.07, sigma = 0.2,
                              T = 5, dt = 1/250, N = 100,
                              plot = TRUE,
                              use.rainbow = TRUE,
                              fixed.color = "gray") {

  n.steps <- T / dt
  time.axis <- seq(0, T, length.out = n.steps)

  # Inizializza matrici
  X <- matrix(NA, nrow = n.steps, ncol = N)
  X[1, ] <- X0
  dPi <- matrix(rpois(n.steps * N, phi * dt), nrow = n.steps)
  dW <- matrix(rnorm(n.steps * N, mean = 0, sd = sqrt(dt)), nrow = n.steps)

  # Simulazione dei path
  for (i in 2:n.steps) {
    X[i, ] <- X[i - 1, ] +
      X[i - 1, ] * (mu - phi * g) * dt +
      X[i - 1, ] * sigma * dW[i - 1, ] +
      X[i - 1, ] * g * dPi[i - 1, ]
  }

  # Plot
  if (plot) {
    cols <- if (use.rainbow) rainbow(N) else rep(adjustcolor(fixed.color, alpha.f = 0.5), N)

    matplot(time.axis, X, type = "l", col = cols, lty = 1,
            xlab = "Time", ylab = "Value", main = "Simulated Jump Diffusion Paths")
    matlines(time.axis, rowMeans(X), col = "black", lwd = 2)
    legend("topleft", legend = "Mean", col = "black", lwd = 2, bty = "n")
  }

  return(list(time = time.axis, paths = X))
}




#------------------------------------------------------------
# Example simulations
#------------------------------------------------------------

set.seed(123)
sim_poisson <- PoissonJump(
  X0 = 0,
  g = 1,
  phi = 2,
  T = 2,
  dt = 1/250,
  N = 20,
  plot = TRUE,
  use.rainbow = FALSE,
  fixed.color = "steelblue"
)

# Access simulated data
head(sim_poisson[, 1:5])

# Plot a single trajectory
plot(seq(0, 2, length.out = nrow(sim_poisson)),
     sim_poisson[, 1],
     type = "l",
     main = "First simulated Poisson jump path",
     ylab = "Process value",
     xlab = "Time")


set.seed(123)
sim_jump <- JumpDiffusion.sim(
  X0 = 100,
  g = -0.2,
  phi = 1,
  mu = 0.07,
  sigma = 0.2,
  T = 2,
  dt = 1/252,
  N = 20,
  plot = TRUE,
  use.rainbow = FALSE,
  fixed.color = "darkorange"
)

# Access simulated data
head(sim_jump$paths)

plot(sim_jump$time, sim_jump$paths[, 1],
     type = "l",
     main = "First simulated trajectory",
     ylab = "Value",
     xlab = "Time")




