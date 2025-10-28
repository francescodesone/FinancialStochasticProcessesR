

# ============================================================
# File: distributions.R
# Author: francescodesone
# Package: FinancialStochasticProcessR
# License: MIT
#
# Description:
#   This file contains custom implementations of random number
#   generators for several financial distributions, including:
#     - Generalized Error Distribution (GED)
#     - Skewed GED (SGED)
#     - Skew-Normal Distribution
#     - Student’s t-Distribution
#     - Skewed Student’s t-Distribution
#
#   These functions are conceptually based on the methods and
#   implementations available in the R package **fGarch**,
#   but are rewritten here for educational and customization
#   purposes, maintaining compatibility with financial modeling
#   workflows in R.
#
# Dependencies:
#   - Base R functions (gamma, rgamma, runif, rt, etc.)
#   - Optionally compatible with the 'fGarch' package
#
# Reference:
#   The fGarch package by Diethelm Würtz et al.
#   (https://cran.r-project.org/package=fGarch)
#
# ============================================================
# ============================================================
# Function: rged()
# Purpose: Generate random numbers from a Generalized Error Distribution (GED)
# ============================================================
#' Generate random numbers from a Generalized Error Distribution (GED)
#'
#' @param n Number of observations to generate.
#' @param mean Mean of the distribution. Default is 0.
#' @param sd Standard deviation of the distribution. Default is 1.
#' @param nu Shape parameter of the GED (nu > 0). Default is 2 (normal distribution).
#'
#' @return A numeric vector of length `n` with GED random values.
#' @export
#' @examples
#' rged(5)
rged <- function(n, mean = 0, sd = 1, nu = 2) {
  if (nu <= 0) stop("nu must be positive")
  if (sd <= 0) stop("sd must be positive")
  lambda <- sqrt(2^(-2 / nu) * gamma(1 / nu) / gamma(3 / nu))
  r <- rgamma(n, shape = 1 / nu)
  z <- lambda * (2 * r)^(1 / nu) * sign(runif(n) - 0.5)
  result <- z * sd + mean
  return(result)
}

# ============================================================
# Function: rsged()
# Purpose: Generate random numbers from a Skewed Generalized Error Distribution (SGED)
# ============================================================
#' Generate random numbers from a Skewed Generalized Error Distribution (SGED)
#'
#' @param n Number of observations to generate.
#' @param mean Mean of the distribution. Default is 0.
#' @param sd Standard deviation of the distribution. Default is 1.
#' @param nu Shape parameter of the underlying GED (nu > 0). Default is 2.
#' @param xi Skewness parameter (xi > 0). Default is 1.5.
#'
#' @return A numeric vector of length `n` with SGED random values.
#' @export
#' @examples
#' rsged(5, xi = 2)
rsged <- function(n, mean = 0, sd = 1, nu = 2, xi = 1.5) {
  if (nu <= 0) stop("nu must be positive")
  if (sd <= 0) stop("sd must be positive")
  if (xi <= 0) stop("xi must be positive")
  lambda <- sqrt(2^(-2 / nu) * gamma(1 / nu) / gamma(3 / nu))
  g  <- nu / (lambda * (2^(1 + 1 / nu)) * gamma(1 / nu))
  m1 <- 2^(1 / nu) * lambda * gamma(2 / nu) / gamma(1 / nu)
  mu <- m1 * (xi - 1 / xi)
  sigma <- sqrt((1 - m1^2) * (xi^2 + 1 / xi^2) + 2 * m1^2 - 1)
  weight <- xi / (xi + 1 / xi)
  z <- runif(n, -weight, 1 - weight)
  Xi <- xi^sign(z)
  rged_inner <- function(n, nu) {
    lambda <- sqrt(2^(-2 / nu) * gamma(1 / nu) / gamma(3 / nu))
    r <- rgamma(n, shape = 1 / nu)
    z <- lambda * (2 * r)^(1 / nu) * sign(runif(n) - 0.5)
    return(z)
  }
  raw <- -abs(rged_inner(n, nu = nu)) / Xi * sign(z)
  standardized <- (raw - mu) / sigma
  result <- standardized * sd + mean
  return(result)
}

# ============================================================
# Function: rsnorm()
# Purpose: Generate random numbers from a Skew-Normal distribution
# ============================================================
#' Generate random numbers from a Skew-Normal distribution
#'
#' @param n Number of observations to generate.
#' @param mean Mean of the distribution. Default is 0.
#' @param sd Standard deviation of the distribution. Default is 1.
#' @param xi Skewness parameter (xi > 0). Default is 1.5.
#'
#' @return A numeric vector of length `n` with Skew-Normal random values.
#' @export
#' @examples
#' rsnorm(5, xi = 2)
rsnorm <- function(n, mean = 0, sd = 1, xi = 1.5) {
  weight <- xi / (xi + 1 / xi)
  z <- runif(n, -weight, 1 - weight)
  Xi <- xi^sign(z)
  random <- -abs(rnorm(n)) / Xi * sign(z)
  m1 <- 2 / sqrt(2 * pi)
  mu <- m1 * (xi - 1 / xi)
  sigma <- sqrt((1 - m1^2) * (xi^2 + 1 / xi^2) + 2 * m1^2 - 1)
  standardized <- (random - mu) / sigma
  result <- standardized * sd + mean
  return(result)
}

# ============================================================
# Function: rstd()
# Purpose: Generate random numbers from a Student’s t-distribution
# ============================================================
#' Generate random numbers from a Student's t-distribution
#'
#' @param n Number of observations to generate.
#' @param mean Mean of the distribution. Default is 0.
#' @param sd Standard deviation of the distribution. Default is 1.
#' @param nu Degrees of freedom (nu > 2). Default is 5.
#'
#' @return A numeric vector of length `n` with Student’s t random values.
#' @export
#' @examples
#' rstd(5, nu = 10)
rstd <- function(n, mean = 0, sd = 1, nu = 5) {
  if (nu <= 2) stop("nu must be greater than 2 for finite variance.")
  s <- sqrt(nu / (nu - 2))
  result <- rt(n = n, df = nu) * sd / s + mean
  return(result)
}

# ============================================================
# Function: rsstd()
# Purpose: Generate random numbers from a Skewed Student’s t-distribution
# ============================================================
#' Generate random numbers from a Skewed Student's t-distribution
#'
#' @param n Number of observations to generate.
#' @param mean Mean of the distribution. Default is 0.
#' @param sd Standard deviation of the distribution. Default is 1.
#' @param nu Degrees of freedom (nu > 2). Default is 5.
#' @param xi Skewness parameter (xi > 0). Default is 1.5.
#'
#' @return A numeric vector of length `n` with Skewed Student’s t random values.
#' @export
#' @examples
#' rsstd(5, nu = 10, xi = 2)
rsstd <- function(n, mean = 0, sd = 1, nu = 5, xi = 1.5) {
  if (nu <= 2) stop("nu must be greater than 2 for finite variance.")
  if (xi <= 0) stop("xi must be positive.")
  beta_fn <- function(a, b) exp(lgamma(a) + lgamma(b) - lgamma(a + b))
  weight <- xi / (xi + 1 / xi)
  z <- runif(n, -weight, 1 - weight)
  Xi <- xi^sign(z)
  random <- -abs(rstd(n, mean = 0, sd = 1, nu = nu)) / Xi * sign(z)
  m1 <- 2 * sqrt(nu - 2) / ((nu - 1) * beta_fn(1/2, nu/2))
  mu <- m1 * (xi - 1 / xi)
  sigma <- sqrt((1 - m1^2) * (xi^2 + 1 / xi^2) + 2 * m1^2 - 1)
  standardized <- (random - mu) / sigma
  result <- standardized * sd + mean
  return(result)
}















