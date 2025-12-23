

#------------------------------------------------------------
# Required Packages
#------------------------------------------------------------
# This function requires the following R packages:
#
# 1. cli   : For colored console output and styled headings.
# 2. boot  : For block bootstrap estimation of Vasicek parameters.
#
# Make sure both packages are installed:
# install.packages(c("cli", "boot"))
#------------------------------------------------------------
library(cli)
library(boot)

#------------------------------------------------------------
# Function: print_centered_heading()
# Purpose: Print a centered heading with decorative lines
#------------------------------------------------------------
print_centered_heading <- function(text) {
  width <- getOption("width", 80)
  text_len <- nchar(text)
  padding <- max(0, floor((width - text_len) / 2))
  line <- paste0(strrep("-", text_len + 4))

  cat(strrep(" ", padding), line, "\n", sep = "")
  cat(strrep(" ", padding), "-- ", text, " --\n", sep = "")
  cat(strrep(" ", padding), line, "\n\n", sep = "")
}

#------------------------------------------------------------
# Function: Vasicek.estimate.params()
# Purpose: Estimate Vasicek model parameters (a, b, sigma)
#          from an observed short-rate time series, with
#          optional block bootstrap confidence intervals.
#------------------------------------------------------------
#' Estimate Parameters of the Vasicek Interest Rate Model
#'
#' This function estimates the parameters of the **Vasicek mean‑reverting
#' short‑rate model** from an observed time series of interest rates.
#'
#' The Vasicek model follows the stochastic differential equation:
#' \deqn{
#' dr_t = a(b - r_t)\,dt + \sigma\,dW_t,
#' }
#' where:
#' \itemize{
#'   \item \eqn{a} is the speed of mean reversion,
#'   \item \eqn{b} is the long‑term mean level,
#'   \item \eqn{\sigma} is the volatility,
#'   \item \eqn{W_t} is a Wiener process.
#' }
#'
#' Parameters are estimated using **ordinary least squares (OLS)** on the
#' discrete‑time approximation:
#' \deqn{
#' r_{t+dt} = \beta_0 + \beta_1 r_t + \epsilon_t,
#' }
#' where:
#' \deqn{
#' a = \frac{1 - \beta_1}{dt}, \qquad
#' b = \frac{\beta_0}{a\,dt}.
#' }
#'
#' @param rate_series Numeric vector of observed short rates.
#'   Must contain at least two values.
#' @param dt Numeric. Time step between observations.
#'   If not provided, defaults to \code{1/250} with a warning.
#' @param name Character. Name used in printed output headings.
#' @param bootstrap Logical. If \code{TRUE}, performs block bootstrap
#'   estimation using \pkg{boot}.
#' @param R Integer. Number of bootstrap replications (default \code{1000}).
#' @param block_length Integer. Block length for moving block bootstrap.
#'
#' @details
#' Volatility is estimated from first differences:
#' \deqn{
#' \hat{\sigma} = \frac{\mathrm{sd}(r_{t+1} - r_t)}{\sqrt{dt}}.
#' }
#'
#' When \code{bootstrap = TRUE}, the function uses
#' \code{boot::tsboot()} with fixed block resampling to compute
#' empirical confidence intervals for \eqn{a}, \eqn{b}, and \eqn{\sigma}.
#'
#' Results are printed using colored console output via the \pkg{cli} package.
#'
#' @return
#' A named list containing:
#' \itemize{
#'   \item \code{a}: Estimated mean‑reversion speed.
#'   \item \code{b}: Estimated long‑term mean.
#'   \item \code{sigma}: Estimated volatility.
#'   \item \code{dt}: Time step used.
#'   \item \code{bootstrap} (optional): Bootstrap result object.
#' }
#'
#' @examples
#' #------------------------------------------------------------
#' # Example: Estimate Vasicek parameters from synthetic data
#' #------------------------------------------------------------
#' set.seed(123)
#' r <- 0.03 + arima.sim(list(ar = 0.98), n = 500, sd = 0.001)
#'
#' Vasicek.estimate.params(
#'   rate_series = r,
#'   dt = 1/250,
#'   name = "Synthetic Short Rate"
#' )
#'
#' #------------------------------------------------------------
#' # Example with bootstrap (may take longer)
#' #------------------------------------------------------------
#' \dontrun{
#' Vasicek.estimate.params(
#'   rate_series = r,
#'   dt = 1/250,
#'   bootstrap = TRUE,
#'   R = 500,
#'   block_length = 20
#' )
#' }
#'
#' @export
Vasicek.estimate.params <- function(rate_series, dt = NULL, name = "Interest Rate",
                                    bootstrap = FALSE, R = 1000, block_length = 20) {

  if (!is.numeric(rate_series) || length(rate_series) < 2) {
    stop("Input must be a numeric vector of interest rates with at least two values.")
  }

  if (is.null(dt)) {
    cli::cli_alert_warning("dt (time step) not specified. Assuming daily data, setting dt = 1/250.")
    dt <- 1 / 250
  }

  #------------------------------------------------------------
  # Internal bootstrap statistic function
  #------------------------------------------------------------
  stat_fun <- function(series, indices) {
    sample_series <- series[indices]
    dr <- diff(sample_series)
    sigma_hat <- sd(dr) / sqrt(dt)

    regr <- lm(tail(sample_series, -1) ~ head(sample_series, -1))
    beta0 <- coef(regr)[1]
    beta1 <- coef(regr)[2]

    a_hat <- (1 - beta1) / dt
    b_hat <- beta0 / (a_hat * dt)

    c(a = a_hat, b = b_hat, sigma = sigma_hat)
  }

  #------------------------------------------------------------
  # Parameter estimation on original series
  #------------------------------------------------------------
  dr <- diff(rate_series)
  sigma_hat <- sd(dr) / sqrt(dt)

  regr <- lm(tail(rate_series, -1) ~ head(rate_series, -1))
  beta0 <- coef(regr)[1]
  beta1 <- coef(regr)[2]

  a_hat <- (1 - beta1) / dt
  b_hat <- beta0 / (a_hat * dt)

  #------------------------------------------------------------
  # Print formatted heading
  #------------------------------------------------------------
  heading_text <- paste0("Vasicek Parameter Estimates for ", cli::style_bold(name))
  print_centered_heading(heading_text)

  cat(cli::col_green("✓ Mean-reversion speed (a): "), sprintf("%10.5f\n", a_hat))
  cat(cli::col_blue("✓ Long-term mean (b):        "), sprintf("%10.5f\n", b_hat))
  cat(cli::col_red("✓ Volatility (sigma):        "), sprintf("%10.5f\n", sigma_hat))
  cat(cli::col_silver(paste0("Assumed dt = ", dt, " (", round(1/dt), " obs/year)")), "\n\n")

  #------------------------------------------------------------
  # Optional bootstrap
  #------------------------------------------------------------
  if (bootstrap) {
    cli::cli_alert_info("Bootstrap enabled: running block bootstrap with fixed block length.")

    boot_res <- tsboot(
      tseries = rate_series,
      statistic = stat_fun,
      R = R,
      l = block_length,
      sim = "fixed"
    )

    a_boot <- boot_res$t[, 1]
    b_boot <- boot_res$t[, 2]
    sigma_boot <- boot_res$t[, 3]

    ci_a <- quantile(a_boot, probs = c(0.025, 0.975))
    ci_b <- quantile(b_boot, probs = c(0.025, 0.975))
    ci_sigma <- quantile(sigma_boot, probs = c(0.025, 0.975))

    cat(cli::col_blue("Bootstrap 95% CI for a:       "),
        sprintf("[%10.5f, %10.5f]\n", ci_a[1], ci_a[2]))
    cat(cli::col_blue("Bootstrap 95% CI for b:       "),
        sprintf("[%10.5f, %10.5f]\n", ci_b[1], ci_b[2]))
    cat(cli::col_blue("Bootstrap 95% CI for sigma:   "),
        sprintf("[%10.5f, %10.5f]\n", ci_sigma[1], ci_sigma[2]))
    cat("\n")
  }

  params <- list(a = a_hat, b = b_hat, sigma = sigma_hat, dt = dt)
  if (bootstrap) params$bootstrap <- boot_res

  return(invisible(params))
}



#------------------------------------------------------------
# Required Packages
#------------------------------------------------------------
# This function requires the following R packages:
#
# 1. cli   : For colored console output and styled headings.
# 2. boot  : For block bootstrap estimation of CIR parameters.
#
# Make sure both packages are installed:
# install.packages(c("cli", "boot"))
#------------------------------------------------------------
library(cli)
library(boot)

#------------------------------------------------------------
# Function: print_centered_heading()
# Purpose: Print a centered heading with decorative lines
#------------------------------------------------------------
print_centered_heading <- function(text) {
  width <- getOption("width", 80)
  text_len <- nchar(text)
  padding <- max(0, floor((width - text_len) / 2))
  line <- paste0(strrep("-", text_len + 4))

  cat(strrep(" ", padding), line, "\n", sep = "")
  cat(strrep(" ", padding), "-- ", text, " --\n", sep = "")
  cat(strrep(" ", padding), line, "\n\n", sep = "")
}

#------------------------------------------------------------
# Function: CIR.estimate.params()
# Purpose: Estimate CIR model parameters (a, b, sigma)
#          from an observed short-rate time series, with
#          optional block bootstrap confidence intervals.
#------------------------------------------------------------
#' Estimate Parameters of the Cox–Ingersoll–Ross (CIR) Interest Rate Model
#'
#' This function estimates the parameters of the **CIR mean‑reverting,
#' strictly positive short‑rate model** from an observed time series of
#' interest rates.
#'
#' The CIR model follows the stochastic differential equation:
#' \deqn{
#' dr_t = a(b - r_t)\,dt + \sigma\sqrt{r_t}\,dW_t,
#' }
#' where:
#' \itemize{
#'   \item \eqn{a} is the speed of mean reversion,
#'   \item \eqn{b} is the long‑term mean level,
#'   \item \eqn{\sigma} is the volatility,
#'   \item \eqn{W_t} is a Wiener process.
#' }
#'
#' To linearize the model, the transformation
#' \deqn{
#' y_t = 2\sqrt{r_t}
#' }
#' is used, yielding the approximate discrete‑time regression:
#' \deqn{
#' y_{t+dt} = \beta_1 y_t + \beta_2 y_t^{-1} + \epsilon_t,
#' }
#' from which:
#' \deqn{
#' a = \frac{2(1 - \beta_1)}{dt}, \qquad
#' b = \frac{\beta_2 + 0.5\sigma^2 dt}{4(1 - \beta_1)}.
#' }
#'
#' @param rate_series Numeric vector of observed short rates.
#'   Must contain at least two values and at least three positive values.
#' @param dt Numeric. Time step between observations.
#'   If not provided, defaults to \code{1/250} with a warning.
#' @param name Character. Name used in printed output headings.
#' @param bootstrap Logical. If \code{TRUE}, performs block bootstrap
#'   estimation using \pkg{boot}.
#' @param R Integer. Number of bootstrap replications (default \code{1000}).
#' @param block_length Integer. Block length for moving block bootstrap.
#'
#' @details
#' Volatility is estimated from the transformed increments:
#' \deqn{
#' \hat{\sigma} = \frac{\mathrm{sd}(y_{t+1} - y_t)}{\sqrt{dt}}.
#' }
#'
#' When \code{bootstrap = TRUE}, the function uses
#' \code{boot::tsboot()} with fixed block resampling to compute
#' empirical confidence intervals for \eqn{a}, \eqn{b}, and \eqn{\sigma}.
#'
#' Results are printed using colored console output via the \pkg{cli} package.
#'
#' @return
#' A named list containing:
#' \itemize{
#'   \item \code{a}: Estimated mean‑reversion speed.
#'   \item \code{b}: Estimated long‑term mean.
#'   \item \code{sigma}: Estimated volatility.
#'   \item \code{dt}: Time step used.
#'   \item \code{bootstrap} (optional): Bootstrap result object.
#' }
#'
#' @examples
#' #------------------------------------------------------------
#' # Example: Estimate CIR parameters from synthetic data
#' #------------------------------------------------------------
#' set.seed(123)
#' r <- abs(arima.sim(list(ar = 0.97), n = 500, sd = 0.002))  # ensure positivity
#'
#' CIR.estimate.params(
#'   rate_series = r,
#'   dt = 1/250,
#'   name = "Synthetic CIR Short Rate"
#' )
#'
#' #------------------------------------------------------------
#' # Example with bootstrap (may take longer)
#' #------------------------------------------------------------
#' \dontrun{
#' CIR.estimate.params(
#'   rate_series = r,
#'   dt = 1/250,
#'   bootstrap = TRUE,
#'   R = 500,
#'   block_length = 20
#' )
#' }
#'
#' @export
CIR.estimate.params <- function(rate_series, dt = NULL, name = "Interest Rate",
                                bootstrap = FALSE, R = 1000, block_length = 20) {

  if (!is.numeric(rate_series) || length(rate_series) < 2) {
    stop("Input must be a numeric vector of interest rates with at least two values.")
  }

  if (is.null(dt)) {
    cli::cli_alert_warning("dt (time step) not specified. Assuming daily data, setting dt = 1/250.")
    dt <- 1 / 250
  }

  #------------------------------------------------------------
  # Filter positive rates and apply CIR transformation
  #------------------------------------------------------------
  rate_positive <- rate_series[rate_series > 0]
  if (length(rate_positive) < 3) {
    stop("Too few positive interest rates for CIR estimation.")
  }

  y <- 2 * sqrt(rate_positive)

  #------------------------------------------------------------
  # Internal bootstrap statistic function
  #------------------------------------------------------------
  stat_fun <- function(series, indices) {
    sample_series <- series[indices]
    sample_positive <- sample_series[sample_series > 0]
    y_s <- 2 * sqrt(sample_positive)

    dy <- diff(y_s)
    sigma_hat <- sd(dy) / sqrt(dt)

    Y <- tail(y_s, -1)
    X1 <- head(y_s, -1)
    X2 <- head(y_s, -1)^(-1)

    regr <- lm(Y ~ X1 + X2 - 1)
    beta0 <- coef(regr)[1]
    beta1 <- coef(regr)[2]

    a_hat <- 2 * (1 - beta0) / dt
    b_hat <- (beta1 + 0.5 * sigma_hat^2 * dt) / (4 * (1 - beta0))

    c(a = a_hat, b = b_hat, sigma = sigma_hat)
  }

  #------------------------------------------------------------
  # Parameter estimation on original series
  #------------------------------------------------------------
  dy <- diff(y)
  sigma_hat <- sd(dy) / sqrt(dt)

  Y <- tail(y, -1)
  X1 <- head(y, -1)
  X2 <- head(y, -1)^(-1)

  regr <- lm(Y ~ X1 + X2 - 1)
  beta0 <- coef(regr)[1]
  beta1 <- coef(regr)[2]

  a_hat <- 2 * (1 - beta0) / dt
  b_hat <- (beta1 + 0.5 * sigma_hat^2 * dt) / (4 * (1 - beta0))

  #------------------------------------------------------------
  # Print formatted heading
  #------------------------------------------------------------
  heading_text <- paste0("CIR Parameter Estimates for ", cli::style_bold(name))
  print_centered_heading(heading_text)

  cat(cli::col_green("✓ Mean-reversion speed (a): "), sprintf("%10.5f\n", a_hat))
  cat(cli::col_blue("✓ Long-term mean (b):        "), sprintf("%10.5f\n", b_hat))
  cat(cli::col_red("✓ Volatility (sigma):        "), sprintf("%10.5f\n", sigma_hat))
  cat(cli::col_silver(paste0("Assumed dt = ", dt, " (", round(1/dt), " obs/year)")), "\n\n")

  #------------------------------------------------------------
  # Optional bootstrap
  #------------------------------------------------------------
  if (bootstrap) {
    cli::cli_alert_info("Bootstrap enabled: running block bootstrap with fixed block length.")

    boot_res <- tsboot(
      tseries = rate_series,
      statistic = stat_fun,
      R = R,
      l = block_length,
      sim = "fixed"
    )

    a_boot <- boot_res$t[, 1]
    b_boot <- boot_res$t[, 2]
    sigma_boot <- boot_res$t[, 3]

    ci_a <- quantile(a_boot, probs = c(0.025, 0.975))
    ci_b <- quantile(b_boot, probs = c(0.025, 0.975))
    ci_sigma <- quantile(sigma_boot, probs = c(0.025, 0.975))

    cat(cli::col_blue("Bootstrap 95% CI for a:       "),
        sprintf("[%10.5f, %10.5f]\n", ci_a[1], ci_a[2]))
    cat(cli::col_blue("Bootstrap 95% CI for b:       "),
        sprintf("[%10.5f, %10.5f]\n", ci_b[1], ci_b[2]))
    cat(cli::col_blue("Bootstrap 95% CI for sigma:   "),
        sprintf("[%10.5f, %10.5f]\n", ci_sigma[1], ci_sigma[2]))
    cat("\n")
  }

  params <- list(a = a_hat, b = b_hat, sigma = sigma_hat, dt = dt)
  if (bootstrap) params$bootstrap <- boot_res

  return(invisible(params))
}










#------------------------------------------------------------
# Example simulations — Vasicek
#------------------------------------------------------------

set.seed(42)
r <- numeric(251)
r[1] <- 0.03
a_true <- 0.15
b_true <- 0.05
sigma_true <- 0.01
dt <- 1/250

for (t in 2:251) {
  r[t] <- r[t-1] + a_true*(b_true - r[t-1])*dt + sigma_true*sqrt(dt)*rnorm(1)
}

# Stima senza bootstrap
Vasicek.estimate.params(r, dt = dt, name = "Simulated Rate")

# Stima con bootstrap
Vasicek.estimate.params(r, dt = dt, name = "Simulated Rate",
                        bootstrap = TRUE, R = 500, block_length = 20)







#------------------------------------------------------------
# Example simulations — CIR
#------------------------------------------------------------

set.seed(42)
r <- numeric(251)
r[1] <- 0.03
a_true <- 0.15
b_true <- 0.05
sigma_true <- 0.01
dt <- 1/250

for (t in 2:251) {
  r[t] <- r[t-1] + a_true*(b_true - r[t-1])*dt + sigma_true*sqrt(r[t-1])*sqrt(dt)*rnorm(1)
  # Garantiamo positività (con max)
  r[t] <- max(r[t], 0.0001)
}

# Stima senza bootstrap
CIR.estimate.params(r, dt = dt, name = "Simulated CIR Rate")

# Stima con bootstrap
CIR.estimate.params(r, dt = dt, name = "Simulated CIR Rate",
                    bootstrap = TRUE, R = 500, block_length = 20)

























