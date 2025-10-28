#------------------------------------------------------------
# Required Packages
#------------------------------------------------------------
# This script requires the following R packages:
#
# 1. cli   : Provides functions for beautiful and user-friendly
#            command line interfaces. Used here for:
#            - Colored alerts and messages
#            - Bold/colored headings in console output
#            Reference: https://cran.r-project.org/package=cli
#
# 2. boot  : Implements bootstrap resampling techniques for
#            estimating distributions, confidence intervals,
#            and other statistics. Used here for:
#            - Block bootstrap of GBM parameter estimates
#            Reference: https://cran.r-project.org/package=boot
#
# Make sure both packages are installed before running this script:
# install.packages(c("cli", "boot"))
#------------------------------------------------------------
library(cli)
library(boot)

#------------------------------------------------------------
# Function: center_text()
# Purpose: Center a text string according to console width
#------------------------------------------------------------
center_text <- function(text) {
  width <- getOption("width", 80)  # default 80 if not set
  padding <- max(0, floor((width - nchar(text)) / 2))
  paste0(strrep(" ", padding), text)
}

#------------------------------------------------------------
# Function: print_centered_heading()
# Purpose: Print a heading with decorative lines centered
#------------------------------------------------------------
print_centered_heading <- function(text) {
  width <- getOption("width", 80)
  text_len <- nchar(text)
  padding <- max(0, floor((width - text_len) / 2))
  line <- paste0(strrep("-", text_len + 4))  # +4 for margin

  cat(strrep(" ", padding), line, "\n", sep = "")
  cat(strrep(" ", padding), "-- ", text, " --\n", sep = "")
  cat(strrep(" ", padding), line, "\n\n", sep = "")
}

#------------------------------------------------------------
# Function: GBM.estimate.params()
# Purpose: Estimate GBM parameters (drift mu and volatility sigma)
#          from a price series, with optional block bootstrap.
#          Requires 'boot' and 'cli' packages.
#------------------------------------------------------------
#' Estimate Drift and Volatility Parameters of a Geometric Brownian Motion (GBM)
#'
#' This function estimates the drift (\eqn{\mu}) and volatility (\eqn{\sigma}) parameters
#' of a **Geometric Brownian Motion (GBM)** from a given price time series.
#' Optionally, it can compute **bootstrap confidence intervals** using a block bootstrap
#' procedure (via the \pkg{boot} package) to account for time-series dependence.
#'
#' The GBM model assumes that asset prices follow the stochastic differential equation:
#' \deqn{dS_t = \mu S_t \, dt + \sigma S_t \, dW_t,}
#' where \eqn{W_t} is a standard Wiener process.
#'
#' @param price.series Numeric vector. Asset price time series (must contain at least two values).
#' @param dt Numeric. Time step between observations (e.g., \code{1/250} for daily data).
#'   If not provided, it defaults to \code{1/250} with a warning message.
#' @param name Character. Optional asset name used in printed output headers (default: \code{"Asset"}).
#' @param bootstrap Logical. If \code{TRUE}, performs a block bootstrap estimation of parameters (default \code{FALSE}).
#' @param R Integer. Number of bootstrap replications (default \code{1000}).
#' @param block_length Integer. Length of each block for block bootstrap (default \code{20}).
#'
#' @details
#' The estimation is based on **log-returns** of the price series:
#' \deqn{r_t = \ln\left(\frac{S_t}{S_{t-1}}\right).}
#'
#' The parameters are estimated using the following relationships:
#' \deqn{\hat{\sigma} = \sqrt{\frac{\mathrm{Var}(r_t)}{dt}},}
#' \deqn{\hat{\mu} = \frac{\mathrm{E}[r_t]}{dt} + \frac{1}{2}\hat{\sigma}^2.}
#'
#' When \code{bootstrap = TRUE}, the function applies a **moving block bootstrap**
#' using \code{boot::tsboot()} to resample the series and estimate empirical confidence intervals
#' for \eqn{\mu} and \eqn{\sigma}.
#'
#' The function prints results in color using the \pkg{cli} package and returns a list
#' containing point estimates and, if requested, bootstrap results.
#'
#' @return
#' A named list containing:
#' \itemize{
#'   \item \code{mu}: Estimated drift parameter (\eqn{\mu}).
#'   \item \code{sigma}: Estimated volatility parameter (\eqn{\sigma}).
#'   \item \code{dt}: Time step used for estimation.
#'   \item \code{bootstrap} (optional): Full bootstrap result object (if \code{bootstrap = TRUE}).
#' }
#'
#' @note
#' Requires the \pkg{boot} and \pkg{cli} packages for bootstrap computation and formatted console output.
#'
#' @examples
#' # Example 1: Estimate GBM parameters for a simple synthetic price series
#' set.seed(123)
#' prices <- cumprod(100 * exp(rnorm(252, mean = 0.0003, sd = 0.01)))
#' GBM.estimate.params(price.series = prices, dt = 1/252, name = "Simulated Asset")
#'
#' # Example 2: Estimate parameters with automatic dt and display warning
#' GBM.estimate.params(price.series = prices)
#'
#' # Example 3: Run with block bootstrap and 500 replications
#' \dontrun{
#' GBM.estimate.params(price.series = prices, bootstrap = TRUE, R = 500, block_length = 10)
#' }
#'
#' @export
GBM.estimate.params <- function(price.series,   # numeric vector of prices
                                dt = NULL,      # time step (e.g., 1/250 for daily)
                                name = NULL,    # asset name for display
                                bootstrap = FALSE, # whether to perform block bootstrap
                                R = 1000,       # number of bootstrap replicates
                                block_length = 20) { # block length for bootstrap

  # Default name
  if (is.null(name) || name == "") {
    name <- "Asset"
  }

  # Input validation
  if (!is.numeric(price.series) || length(price.series) < 2) {
    stop("Input must be a numeric vector of prices with at least two values.")
  }

  # Handle dt
  if (is.null(dt)) {
    cli::cli_alert_warning("dt not specified. Assuming daily data, dt = 1/250.")
    dt <- 1/250
  } else if (dt <= 0 || dt > 1) {
    cli::cli_alert_warning("dt seems unusual. Check time step.")
  }

  # Compute log-returns
  dlnS <- diff(log(price.series))

  # Function to compute mu and sigma for bootstrap samples
  stat_fun <- function(series, indices) {
    sample_prices <- series[indices]
    sample_dlnS <- diff(log(sample_prices))
    sigma_hat <- sqrt(var(sample_dlnS) / dt)
    mu_hat <- mean(sample_dlnS) / dt + 0.5 * sigma_hat^2
    c(mu = mu_hat, sigma = sigma_hat)
  }

  # Estimate parameters on original series
  sigma <- sqrt(var(dlnS) / dt)
  mu <- mean(dlnS) / dt + 0.5 * sigma^2

  # Print heading
  heading_text <- paste0("GBM Parameter Estimates for ", cli::style_bold(name))
  print_centered_heading(heading_text)

  # Display estimates
  cat(cli::col_green("✓ Drift (mu):        "), sprintf("%10.5f  (per year)", mu), "\n")
  cat(cli::col_red("✓ Volatility (sigma):"), sprintf("%10.5f  (per √year)", sigma), "\n")
  cat(cli::col_silver(paste0("Assumed dt = ", dt, " (", round(1/dt), " obs/year)")), "\n\n")

  # Optional bootstrap
  if (bootstrap) {
    cli::cli_alert_warning("Bootstrap enabled using tsboot package")

    boot_res <- tsboot(tseries = price.series,
                       statistic = stat_fun,
                       R = R,
                       l = block_length,
                       sim = "fixed")

    mu_boot <- boot_res$t[,1]
    sigma_boot <- boot_res$t[,2]

    ci_mu <- quantile(mu_boot, probs = c(0.025, 0.975))
    med_mu <- median(mu_boot)

    ci_sigma <- quantile(sigma_boot, probs = c(0.025, 0.975))
    med_sigma <- median(sigma_boot)

    cat(cli::col_blue("Bootstrap 95% CI for mu:     "), sprintf("[%8.5f, %8.5f]\n", ci_mu[1], ci_mu[2]))
    cat(cli::col_green("Median (mu):                 "), sprintf("%8.5f\n", med_mu))
    cat(cli::col_blue("Bootstrap 95% CI for sigma:  "), sprintf("[%8.5f, %8.5f]\n", ci_sigma[1], ci_sigma[2]))
    cat(cli::col_red("Median (sigma):              "), sprintf("%8.5f\n", med_sigma))
    cat("\n")
  }

  # Return parameters (and bootstrap results if requested)
  params <- list(mu = mu, sigma = sigma, dt = dt)
  if (bootstrap) params$bootstrap <- boot_res

  return(invisible(params))
}

#------------------------------------------------------------
# Examples
#------------------------------------------------------------

# Example 1: Simple estimation from simulated price series
set.seed(42)
SP <- cumprod(1 + rnorm(251, mean = 0.002, sd = 0.01))
params <- GBM.estimate.params(SP)

# Example 2: Estimation with bootstrap and named asset
params_boot <- GBM.estimate.params(
  SP,
  dt = 1/250,
  name = "Simulated S&P 500",
  bootstrap = TRUE,
  R = 1000,
  block_length = 20
)









