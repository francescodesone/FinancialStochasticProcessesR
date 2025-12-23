#------------------------------------------------------------
# Function: BlackKarasinski.multi()
# Purpose: Simulate multiple paths of the Black–Karasinski
#          short-rate model with time‑varying parameters.
#------------------------------------------------------------
#' Simulate Black–Karasinski Short-Rate Model (Multiple Paths)
#'
#' @description
#' This function simulates \code{num.paths} trajectories of the
#' Black–Karasinski (BK) short-rate model in its *non‑lognormal*
#' Euler‑discretized form with time‑varying parameters.
#'
#' The continuous‑time BK dynamics (in log form) are:
#' \deqn{
#' d(\ln r_t) = \left(\theta(t) - \alpha(t)\ln r_t\right) dt
#'              + \sigma(t)\, dW_t
#' }
#'
#' Qui viene utilizzata una versione lineare approssimata:
#' \deqn{
#' dr_t = \left(\theta(t) - \alpha(t) r_t\right) dt
#'        + \sigma(t)\sqrt{dt}\, Z_t
#' }
#' con \eqn{Z_t \sim N(0,1)}.
#'
#' @param r0 Numeric scalar. Initial short rate \eqn{r(0)}.
#' @param n Integer. Number of time steps in each simulated path.
#' @param dt Numeric scalar. Time step size (default: 1/252).
#' @param theta Numeric vector of length \code{n}. Time‑varying drift component.
#' @param sigma Numeric vector of length \code{n}. Time‑varying volatility.
#' @param alpha Numeric vector of length \code{n}. Mean‑reversion speed.
#' @param num.paths Integer. Number of simulated paths.
#' @param plot Logical. If TRUE (default), plots the simulated BK paths.
#' @param use.rainbow Logical. If TRUE, uses rainbow colors for each path;
#' otherwise uses a fixed semi‑transparent color.
#' @param fixed.color Character. Color used when \code{use.rainbow = FALSE}.
#'
#' @return A matrix of dimension \code{num.paths × n} containing the simulated
#' Black–Karasinski short‑rate paths. Each row corresponds to one trajectory.
#'
#' @details
#' Il modello Black–Karasinski è un modello di tasso a media‑reversione
#' che garantisce positività del tasso quando implementato nella sua forma
#' lognormale.
#'
#' In questa implementazione viene utilizzata una discretizzazione diretta
#' del tasso (non lognormale), utile per simulazioni veloci con parametri
#' time‑dependent:
#'
#' \deqn{
#' r_{t+dt} = r_t + (\theta_t - \alpha_t r_t) dt
#'            + \sigma_t \sqrt{dt}\, Z_t
#' }
#'
#' Se \code{plot = TRUE}, la funzione mostra:
#' \itemize{
#'   \item tutte le traiettorie simulate,
#'   \item colori personalizzabili,
#'   \item asse temporale coerente con \code{dt}.
#' }
#'
#' @examples
#' #------------------------------------------------------------
#' # Example simulations
#' #------------------------------------------------------------
#'
#' set.seed(42)
#'
#' bk_paths <- BlackKarasinski.multi(
#'   r0 = 0.03,
#'   n = 1000,
#'   dt = 1 / 252,
#'   theta = seq(0.03, 0.05, length.out = 1000),
#'   sigma = seq(0.02, 0.04, length.out = 1000),
#'   alpha = rep(0.1, 1000),
#'   num.paths = 20,
#'   plot = TRUE,
#'   use.rainbow = FALSE
#' )
#'
BlackKarasinski.multi <- function(r0 = 0.05,
                                  n = 1000,
                                  dt = 1 / 252,
                                  theta = seq(0.05, 0.06, length.out = n),
                                  sigma = seq(0.03, 0.06, length.out = n),
                                  alpha = rep(0.1, n),
                                  num.paths = 10,
                                  plot = TRUE,
                                  use.rainbow = TRUE,
                                  fixed.color = "darkgreen") {

  # Inizializza matrice dei path
  r.paths <- matrix(0, nrow = num.paths, ncol = n)
  r.paths[, 1] <- r0

  # Simulazione multipla
  for (path in 1:num.paths) {
    for (i in 2:n) {
      dr <- theta[i - 1] - alpha[i - 1] * r.paths[path, i - 1]
      r.paths[path, i] <- r.paths[path, i - 1] + dr * dt +
        sigma[i - 1] * sqrt(dt) * rnorm(1)
    }
  }

  # Plot
  if (plot) {
    time.axis <- seq(0, (n - 1) * dt, by = dt)
    colors <- if (use.rainbow) {
      rainbow(num.paths)
    } else {
      rep(adjustcolor(fixed.color, alpha.f = 0.5), num.paths)
    }

    matplot(time.axis, t(r.paths), type = "l", col = colors, lty = 1,
            xlab = "Time", ylab = "Short Rate",
            main = paste("Black-Karasinski Model (", num.paths, " paths)", sep = ""))
  }

  return(r.paths)
}




#------------------------------------------------------------
# Function: BlackKarasinskiLog.multi()
# Purpose: Simulate multiple lognormal Black–Karasinski
#          short-rate paths with time‑varying parameters.
#------------------------------------------------------------
#' Simulate Lognormal Black–Karasinski Short-Rate Model (Multiple Paths)
#'
#' @description
#' This function simulates \code{num.paths} trajectories of the
#' lognormal Black–Karasinski (BK) short-rate model with
#' time‑dependent parameters \code{theta}, \code{alpha}, and \code{sigma}.
#'
#' The Black–Karasinski model is defined on the logarithm of the short rate:
#' \deqn{
#' d(\ln r_t) = \left(\theta(t) - \alpha(t)\ln r_t\right) dt
#'              + \sigma(t)\, dW_t
#' }
#'
#' The simulation uses Euler discretization on \eqn{\ln r_t}, ensuring
#' strictly positive interest rates:
#' \deqn{
#' \ln r_{t+dt} = \ln r_t
#'               + (\theta_t - \alpha_t \ln r_t) dt
#'               + \sigma_t \sqrt{dt}\, Z_t
#' }
#' \deqn{
#' r_{t+dt} = \exp(\ln r_{t+dt})
#' }
#'
#' @param r0 Numeric scalar. Initial short rate \eqn{r(0)}.
#' @param n Integer. Number of time steps in each simulated path.
#' @param dt Numeric scalar. Time step size (default: 1/252).
#' @param theta Numeric vector of length \code{n}. Time‑varying drift component.
#' @param sigma Numeric vector of length \code{n}. Time‑varying volatility.
#' @param alpha Numeric vector of length \code{n}. Mean‑reversion speed.
#' @param num.paths Integer. Number of simulated paths.
#' @param plot Logical. If TRUE (default), plots the simulated BK lognormal paths.
#' @param use.rainbow Logical. If TRUE, uses rainbow colors for each path;
#' otherwise uses a fixed semi‑transparent color.
#' @param fixed.color Character. Color used when \code{use.rainbow = FALSE}.
#'
#' @return A matrix of dimension \code{num.paths × n} containing the simulated
#' lognormal Black–Karasinski short‑rate paths. Each row corresponds to one trajectory.
#'
#' @details
#' Il modello Black–Karasinski lognormale garantisce tassi sempre positivi
#' grazie alla simulazione diretta del processo su \eqn{\ln r_t}.
#'
#' La discretizzazione utilizzata è:
#' \deqn{
#' \ln r_{t+dt} = \ln r_t
#'               + (\theta_t - \alpha_t \ln r_t) dt
#'               + \sigma_t \sqrt{dt}\, Z_t
#' }
#' con \eqn{Z_t \sim N(0,1)}.
#'
#' Se \code{plot = TRUE}, la funzione mostra tutte le traiettorie simulate
#' con colori personalizzabili.
#'
#' @examples
#' #------------------------------------------------------------
#' # Example simulations
#' #------------------------------------------------------------
#'
#' set.seed(42)
#'
#' bk_log_paths <- BlackKarasinskiLog.multi(
#'   r0 = 0.03,
#'   n = 1000,
#'   dt = 1 / 252,
#'   theta = seq(0.03, 0.05, length.out = 1000),
#'   sigma = seq(0.02, 0.04, length.out = 1000),
#'   alpha = rep(0.1, 1000),
#'   num.paths = 20,
#'   plot = TRUE,
#'   use.rainbow = FALSE
#' )
#'
BlackKarasinskiLog.multi <- function(r0 = 0.05,
                                     n = 1000,
                                     dt = 1 / 252,
                                     theta = seq(0.05, 0.06, length.out = n),
                                     sigma = seq(0.03, 0.06, length.out = n),
                                     alpha = rep(0.1, n),
                                     num.paths = 10,
                                     plot = TRUE,
                                     use.rainbow = TRUE,
                                     fixed.color = "darkred") {

  # Inizializza matrice dei path
  r.paths <- matrix(0, nrow = num.paths, ncol = n)
  r.paths[, 1] <- r0

  # Simulazione lognormale
  for (path in 1:num.paths) {
    for (i in 2:n) {
      ri_1_ln <- log(r.paths[path, i - 1])
      dri_ln <- theta[i - 1] - alpha[i - 1] * ri_1_ln
      ri_ln <- ri_1_ln + dri_ln * dt + sigma[i - 1] * sqrt(dt) * rnorm(1)
      r.paths[path, i] <- exp(ri_ln)
    }
  }

  # Plot
  if (plot) {
    time.axis <- seq(0, (n - 1) * dt, by = dt)
    colors <- if (use.rainbow) {
      rainbow(num.paths)
    } else {
      rep(adjustcolor(fixed.color, alpha.f = 0.5), num.paths)
    }

    matplot(time.axis, t(r.paths), type = "l", col = colors, lty = 1,
            xlab = "Time", ylab = "Short Rate",
            main = paste("Black-Karasinski Lognormal (", num.paths, " paths)", sep = ""))
  }

  return(r.paths)
}






#------------------------------------------------------------
# Example simulations — Black–Karasinski (linear version)
#------------------------------------------------------------

set.seed(42)

bk_paths <- BlackKarasinski.multi(
  r0 = 0.03,
  n = 1000,
  dt = 1 / 252,
  theta = seq(0.03, 0.05, length.out = 1000),
  sigma = seq(0.02, 0.04, length.out = 1000),
  alpha = rep(0.1, 1000),
  num.paths = 20,
  plot = TRUE,
  use.rainbow = FALSE,
  fixed.color = "darkgreen"
)




#------------------------------------------------------------
# Example simulations — Black–Karasinski (lognormal version)
#------------------------------------------------------------

set.seed(42)

bk_log_paths <- BlackKarasinskiLog.multi(
  r0 = 0.03,
  n = 1000,
  dt = 1 / 252,
  theta = seq(0.03, 0.05, length.out = 1000),
  sigma = seq(0.02, 0.04, length.out = 1000),
  alpha = rep(0.1, 1000),
  num.paths = 20,
  plot = TRUE,
  use.rainbow = FALSE,
  fixed.color = "darkred"
)


















