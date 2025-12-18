
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
#            - Block bootstrap of Jump-Diffusion parameter estimates
#            Reference: https://cran.r-project.org/package=boot
#
# Make sure both packages are installed before running this script:
# install.packages(c("cli", "boot"))
#------------------------------------------------------------
library(cli)
library(boot)

#------------------------------------------------------------
# Function: JD.estimate.params()
# Purpose: Estimate Jump-Diffusion parameters (mu, sigma, g, phi)
#          from a price series, with optional block bootstrap.
#------------------------------------------------------------
#' Estimate Parameters of a Jump-Diffusion Model
#'
#' This function estimates the drift (\eqn{\mu}), volatility (\eqn{\sigma}),
#' jump size (\eqn{g}), and jump intensity (\eqn{\phi}) of a **Jump-Diffusion model**
#' from a given asset price series. Optionally, it can compute **bootstrap confidence intervals**
#' using a block bootstrap procedure (via the \pkg{boot} package).
#'
#' The Jump-Diffusion model assumes:
#' \deqn{dX_t = X_t \left[(\mu - \phi g) dt + \sigma dW_t + g d\Pi_t \right]}
#' where \eqn{W_t} is a Wiener process and \eqn{\Pi_t} a Poisson process with intensity \eqn{\phi}.
#'
#' @param price.series Numeric vector. Asset price time series (must contain at least two values).
#' @param dt Numeric. Time step between observations (e.g., \code{1/250} for daily data).
#'   If not provided, defaults to \code{1/250} with a warning message.
#' @param name Character. Optional asset name used in printed output headers (default: \code{"Asset"}).
#' @param bootstrap Logical. If \code{TRUE}, performs a block bootstrap estimation of parameters (default \code{FALSE}).
#' @param R Integer. Number of bootstrap replications (default \code{1000}).
#' @param block_length Integer. Length of each block for block bootstrap (default \code{20}).
#'
#' @return
#' A named list containing:
#' \itemize{
#'   \item \code{mu}: Estimated drift parameter.
#'   \item \code{sigma}: Estimated volatility parameter.
#'   \item \code{g}: Estimated jump size.
#'   \item \code{phi}: Estimated jump intensity.
#'   \item \code{dt}: Time step used for estimation.
#'   \item \code{bootstrap} (optional): Full bootstrap result object (if \code{bootstrap = TRUE}).
#' }
#'
#' @examples
#' set.seed(123)
#' sim <- JumpDiffusion.sim(X0 = 100, g = -0.2, phi = 1,
#'                          mu = 0.07, sigma = 0.2,
#'                          T = 2, dt = 1/252, N = 1, plot = FALSE)
#' prices <- sim$paths[,1]
#' JD.estimate.params(price.series = prices, dt = 1/252, name = "Simulated JD")
#'
#' @export
JD.estimate.params <- function(price.series, dt = NULL, name = NULL,
                               bootstrap = FALSE, R = 1000, block_length = 20) {
  if (is.null(name) || name == "") {
    name <- "Asset"
  }

  if (!is.numeric(price.series) || length(price.series) < 2) {
    stop("Input must be a numeric vector of prices with at least two values.")
  }

  if (is.null(dt)) {
    cli::cli_alert_warning("dt (time step) not specified. Assuming daily data, setting dt = 1/250.")
    dt <- 1/250
  } else if (dt <= 0 || dt > 1) {
    cli::cli_alert_warning("dt value seems unusual. Please ensure it correctly represents the time step.")
  }

  dlnS <- diff(log(price.series))
  M1 <- mean(dlnS)
  M2 <- var(dlnS)
  M3 <- mean(dlnS^3)
  M4 <- mean(dlnS^4)

  g <- exp(M4 / M3) - 1
  phi <- (1 / dt) * M3^4 / M4^3
  sigma_term <- M2 / dt - phi * log(1 + g)^2
  sigma <- ifelse(sigma_term > 0, sqrt(sigma_term), NA_real_)
  mu <- if (!is.na(sigma)) M1 / dt + phi * g + 0.5 * sigma^2 - phi * log(1 + g) else NA_real_

  heading_text <- paste0("Jump Diffusion Parameter Estimates for ", cli::style_bold(name))
  print_centered_heading(heading_text)

  cat(cli::col_green("✓ Drift (mu):        "), sprintf("%10.5f\n", mu))
  cat(cli::col_red("✓ Volatility (sigma):"), sprintf("%10.5f\n", sigma))
  cat(cli::col_cyan("✓ Jump size (g):     "), sprintf("%10.5f\n", g))
  cat(cli::col_magenta("✓ Jump freq (phi):  "), sprintf("%10.5f\n", phi))
  cat(cli::col_silver(paste0("Assumed dt = ", dt, " (i.e. ", round(1/dt), " obs/year)")), "\n\n")

  if (bootstrap) {
    cli::cli_alert_warning("Bootstrap enabled: dt assumed fixed during bootstrapping.")

    stat_fun <- function(series, indices) {
      sample_prices <- series[indices]
      dlnS <- diff(log(sample_prices))
      M1 <- mean(dlnS)
      M2 <- var(dlnS)
      M3 <- mean(dlnS^3)
      M4 <- mean(dlnS^4)

      g <- exp(M4 / M3) - 1
      phi <- (1 / dt) * M3^4 / M4^3
      sigma_term <- M2 / dt - phi * log(1 + g)^2

      if (sigma_term <= 0 || any(!is.finite(c(g, phi, sigma_term)))) {
        return(rep(NA, 4))
      }

      sigma <- sqrt(sigma_term)
      mu <- M1 / dt + phi * g + 0.5 * sigma^2 - phi * log(1 + g)
      c(mu = mu, sigma = sigma, g = g, phi = phi)
    }

    boot_res <- tsboot(tseries = price.series,
                       statistic = stat_fun,
                       R = R,
                       l = block_length,
                       sim = "fixed")

    boot_vals <- na.omit(as.data.frame(boot_res$t))
    colnames(boot_vals) <- c("mu", "sigma", "g", "phi")

    if (nrow(boot_vals) < 30) {
      cli::cli_alert_danger("Too few valid bootstrap samples. Cannot compute confidence intervals reliably.")
    } else {
      summaries <- apply(boot_vals, 2, function(x) c(median = median(x), quantile(x, c(0.025, 0.975))))
      cat(cli::col_blue("Bootstrap 95% CI and Medians:\n"))
      cat(sprintf("  mu   : [%8.5f, %8.5f], median = %8.5f\n", summaries[2, "mu"], summaries[3, "mu"], summaries[1, "mu"]))
      cat(sprintf("  sigma: [%8.5f, %8.5f], median = %8.5f\n", summaries[2, "sigma"], summaries[3, "sigma"], summaries[1, "sigma"]))
      cat(sprintf("  g    : [%8.5f, %8.5f], median = %8.5f\n", summaries[2, "g"], summaries[3, "g"], summaries[1, "g"]))
      cat(sprintf("  phi  : [%8.5f, %8.5f], median = %8.5f\n\n", summaries[2, "phi"], summaries[3, "phi"], summaries[1, "phi"]))
    }
  }

  out <- list(mu = mu, sigma = sigma, g = g, phi = phi, dt = dt)
  if (bootstrap) {
    out$bootstrap <- boot_res
  }

  invisible(out)
}





#------------------------------------------------------------
# Example parameter estimation for Jump-Diffusion
#------------------------------------------------------------

set.seed(123)

# Simulate a Jump-Diffusion price series
sim_jd <- JumpDiffusion.sim(
  X0 = 100,
  g = -0.2,
  phi = 1,
  mu = 0.07,
  sigma = 0.2,
  T = 2,
  dt = 1/252,
  N = 1,
  plot = FALSE
)

# Extract the simulated price series
price_series <- sim_jd$paths[, 1]

#------------------------------------------------------------
# Estimate parameters without bootstrap
#------------------------------------------------------------
params_jd <- JD.estimate.params(
  price.series = price_series,
  dt = 1/252,
  name = "Simulated JD"
)

#------------------------------------------------------------
# Estimate parameters with bootstrap
#------------------------------------------------------------
params_jd_boot <- JD.estimate.params(
  price.series = price_series,
  dt = 1/252,
  name = "Simulated JD",
  bootstrap = TRUE,
  R = 500,
  block_length = 20
)

