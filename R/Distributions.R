

# ============================================================
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
rged <- function(n, mean = 0, sd = 1, nu = 2) {
  # Validate input
  if (nu <= 0) stop("nu must be positive")
  if (sd <= 0) stop("sd must be positive")

  # Compute scale factor (lambda)
  lambda <- sqrt(2^(-2 / nu) * gamma(1 / nu) / gamma(3 / nu))

  # Generate random values from a Gamma distribution
  r <- rgamma(n, shape = 1 / nu)

  # Apply transformation to obtain GED
  z <- lambda * (2 * r)^(1 / nu) * sign(runif(n) - 0.5)

  # Scale and shift
  result <- z * sd + mean

  return(result)
}

# ============================================================
# Function: rsged()
# Purpose: Generate random numbers from a Skewed Generalized Error Distribution (SGED)
# ============================================================
rsged <- function(n, mean = 0, sd = 1, nu = 2, xi = 1.5) {
  # Validate parameters
  if (nu <= 0) stop("nu must be positive")
  if (sd <= 0) stop("sd must be positive")
  if (xi <= 0) stop("xi (skewness parameter) must be positive")

  # Compute auxiliary parameters
  lambda <- sqrt(2^(-2 / nu) * gamma(1 / nu) / gamma(3 / nu))
  g  <- nu / (lambda * (2^(1 + 1 / nu)) * gamma(1 / nu))
  m1 <- 2^(1 / nu) * lambda * gamma(2 / nu) / gamma(1 / nu)
  mu <- m1 * (xi - 1 / xi)
  sigma <- sqrt((1 - m1^2) * (xi^2 + 1 / xi^2) + 2 * m1^2 - 1)

  # Generate standardized SGED random variables
  weight <- xi / (xi + 1 / xi)
  z <- runif(n, -weight, 1 - weight)
  Xi <- xi^sign(z)

  # Nested GED generator
  rged <- function(n, nu) {
    lambda <- sqrt(2^(-2 / nu) * gamma(1 / nu) / gamma(3 / nu))
    r <- rgamma(n, shape = 1 / nu)
    z <- lambda * (2 * r)^(1 / nu) * sign(runif(n) - 0.5)
    return(z)
  }

  raw <- -abs(rged(n, nu = nu)) / Xi * sign(z)
  standardized <- (raw - mu) / sigma

  # Scale and shift
  result <- standardized * sd + mean
  return(result)
}

# ============================================================
# Function: rsnorm()
# Purpose: Generate random numbers from a Skew-Normal distribution
# ============================================================
rsnorm <- function(n, mean = 0, sd = 1, xi = 1.5) {
  # Arguments:
  #   n    - number of observations
  #   mean - desired mean
  #   sd   - desired standard deviation
  #   xi   - skewness parameter (xi > 0)
  #
  # Returns:
  #   A vector of n random values from a Skew-Normal distribution

  # Step 1: Generate standardized skew-normal random values
  weight <- xi / (xi + 1 / xi)
  z <- runif(n, -weight, 1 - weight)
  Xi <- xi^sign(z)
  random <- -abs(rnorm(n)) / Xi * sign(z)

  # Step 2: Compute scale and location parameters
  m1 <- 2 / sqrt(2 * pi)
  mu <- m1 * (xi - 1 / xi)
  sigma <- sqrt((1 - m1^2) * (xi^2 + 1 / xi^2) + 2 * m1^2 - 1)

  # Standardize and rescale
  standardized <- (random - mu) / sigma
  result <- standardized * sd + mean

  return(result)
}

# ============================================================
# Function: rstd()
# Purpose: Generate random numbers from a Student’s t-distribution
# ============================================================
rstd <- function(n, mean = 0, sd = 1, nu = 5) {
  # Arguments:
  #   n    - number of observations
  #   mean - desired mean
  #   sd   - desired standard deviation
  #   nu   - degrees of freedom (> 2)
  #
  # Returns:
  #   A vector of n random values from a Student’s t-distribution

  if (nu <= 2) stop("nu must be greater than 2 for finite variance.")

  s <- sqrt(nu / (nu - 2))
  result <- rt(n = n, df = nu) * sd / s + mean

  return(result)
}

# ============================================================
# Function: rsstd()
# Purpose: Generate random numbers from a Skewed Student’s t-distribution
# ============================================================
rsstd <- function(n, mean = 0, sd = 1, nu = 5, xi = 1.5) {
  # Arguments:
  #   n    - number of observations
  #   mean - desired mean
  #   sd   - desired standard deviation
  #   nu   - degrees of freedom (> 2)
  #   xi   - skewness parameter (xi > 0)
  #
  # Returns:
  #   A vector of n random values from a Skewed Student’s t-distribution

  if (nu <= 2) stop("nu must be greater than 2 for finite variance.")
  if (xi <= 0) stop("xi must be positive.")

  # Beta function (numerically stable version)
  beta_fn <- function(a, b) exp(lgamma(a) + lgamma(b) - lgamma(a + b))

  # Step 1: Skew transformation
  weight <- xi / (xi + 1 / xi)
  z <- runif(n, -weight, 1 - weight)
  Xi <- xi^sign(z)
  random <- -abs(rstd(n, mean = 0, sd = 1, nu = nu)) / Xi * sign(z)

  # Step 2: Correction parameters
  m1 <- 2 * sqrt(nu - 2) / ((nu - 1) * beta_fn(1/2, nu/2))
  mu <- m1 * (xi - 1 / xi)
  sigma <- sqrt((1 - m1^2) * (xi^2 + 1 / xi^2) + 2 * m1^2 - 1)

  # Standardize and rescale
  standardized <- (random - mu) / sigma
  result <- standardized * sd + mean

  return(result)
}




















