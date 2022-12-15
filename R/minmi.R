#' Estimate Extinction Times using the MINMI Estimator
#'
#' This function estimates extinction times using the MINMI procedure
#'
#' @param W Numeric vector of fossil ages.
#' @param sd Numeric vector of measurement error uncertainty associated with each fossil.
#' @param K Numeric upper bound for fossil ages.
#' @param alpha Numeric between 0 and 1. Used to find 100(1-alpha)\% confidence intervals. Defaults to 0.05 (95\% confidence intervals)
#' @param B Optional numeric greater than 1 specifying the number of Monte Carlo samples to use.
#' @param .B_init Optional numeric greater than 1 specifying the number of Monte Carlo samples to use in the pilot estimates. Defaults to 500.
#' @param .max_var Optional numeric greater than 0 specifying the maximum Monte Carlo error variance we aim to have associated with our MINMI estimates.
#'
#' @return A list with estimates for the lower end point of the 100(1-alpha)\% confidence interval, point estimate, upper end point, and a list containing the B's used for each.
#' @export
#' @importFrom stats dnorm pnorm runif var
minmi <- function (W, sd, K, alpha = 0.05, B = NULL, .B_init = 500, .max_var = NULL) {
  result <- list(lower=NULL, point = NULL, upper = NULL, B = list(lower=B, point=B, upper=B))
  n <- length(W)
  m <- min(W)
  dating.sd = mean(sd)
  q <- c(lower = alpha/2, point = 0.5, upper = 1-alpha/2)

  flag.delta_model = (all(sd) == 0)

  if (flag.delta_model) {
    for (i in 1:length(q)) {
      result[[i]] <- estimate_quantile.minmi(K = K, W = W, eps.mean = 0, eps.sigma = dating.sd, q = q[i])
    }
    return(result)
  }
  else {
    # Set the maximum variance
    max_var <- 0.2 * (dating.sd) ^ 2
    if (!is.null(.max_var)) {
      max_var <- .max_var
    }

    # Calculate number of monte carlo samples to use
    if (is.null(B)) {
      u.init <- runif(.B_init, 0, 1)
      for (i in 1:length(q)) {
        result$B[[i]] <- find_optimal_B(max_var = max_var, K = K, m = m, n = n, u = u.init, eps.mean = 0, eps.sigma = dating.sd, q = q[i])
      }
    }

    # Generate Monte Carlo Samples
    B.max <- max(unlist(result$B))
    mc.samples <- runif(B.max, min = 0, max = 1)

    # Calculate estimates
    for (i in 1:length(q)) {
      result[[i]] <- estimate_quantile.minmi(K = K, W = W, u = mc.samples[1:result$B[[i]]], eps.mean = 0, eps.sigma = dating.sd, q = q[i])
    }
    return(result)
  }
}

#' Estimate a specific quantile using the MINMI method.
#'
#' @param K Numeric upper bound for fossil ages.
#' @param W Numeric vector of fossil ages.
#' @param u Optional numeric vector of values sampled from a uniform distribution from 0 to 1 and with same length as W.
#' @param eps.mean Optional numeric mean of the measurement error distribution. Defaults to 0.
#' @param eps.sigma Optional numeric standard error of the measurement error distribution. Defaults to 0.
#' @param q Numeric between 0 and 1 specifying the quantile we want to estimate.
#'
#' @return A numeric indicating the MINMI estimate of the qth quantile of our sample statistic.
#' @export
estimate_quantile.minmi <- function (K, W, u=NULL, eps.mean = 0, eps.sigma = 0, q) {
  n <- length(W)
  m <- min(W)
  # No Measurement Error case
  theta_q.hat <- K - q ^ (-1 / n) * (K - m)
  if (eps.mean != 0 || eps.sigma != 0) {
    # Measurement Error case
    newton.res <- pracma::newtonRaphson(
      fun = function(theta) estimating_eqn(theta, q, K, u, m, n, eps.mean, eps.sigma),
      x0 = theta_q.hat
    )
    theta_q.hat <- newton.res$root
  }
  return(theta_q.hat)
}

#' Function for Finding the Number of Monte Carlo Samples for MINMI
#'
#' @param max_var Numeric greater than 0 specifying the maximum Monte Carlo error variance we aim to have associated with our MINMI estimates.
#' @param K Numeric. Upper bound for fossil record.
#' @param m Numeric. Sample minimum from fossil record.
#' @param n Number of fossils.
#' @param u Numeric vector of values sampled from a uniform distribution from 0 to 1 and with same length as W.
#' @param eps.mean Numeric mean of the measurement error distribution. Defaults to 0.
#' @param eps.sigma Numeric standard error of the measurement error distribution. Defaults to 0.
#' @param q Numeric between 0 and 1 specifying the quantile we want to estimate.
#'
#' @return A number.
#' @export
find_optimal_B <- function (max_var, K, m, n, u, eps.mean, eps.sigma, q) {
  # Initial estimate of theta_q using no measurement error case
  theta_q.hat.init <- K - q ^ (-1 / n) * (K - m)

  # pdfs and CDF evaluations (for convenience)
  f_eps.K <- dnorm(K - theta_q.hat.init, eps.mean, eps.sigma)
  F_eps.K <- pnorm(K - theta_q.hat.init, eps.mean, eps.sigma)
  # Estimate sigma^2_var
  B <- length(u)
  e <- uniform_to_tnorm(u, eps.mean, eps.sigma, a = -Inf, b = m - theta_q.hat.init)
  sample_var.psi_hat <- var((m - e - theta_q.hat.init) / (K - e - theta_q.hat.init))

  # Estimate \hat\psi and \hat\psi prime
  psi_hat <- mean((m - e - theta_q.hat.init) / (K - e - theta_q.hat.init))
  psi_hat.prime <- -mean((K - m) / (K - e - theta_q.hat.init) ^ 2)

  optimal_B <- 1 / max_var * sample_var.psi_hat * (f_eps.K / F_eps.K * psi_hat + psi_hat.prime) ^ (-2)

  return(max(ceiling(optimal_B), 100))
}

#' Estimating equation for Newton-Raphson's Algorithm
#'
#' @param theta Number. Estimate for theta.
#' @param q Numeric between 0 and 1 specifying the quantile we want to estimate.
#' @param K Numeric. Upper bound for fossil record.
#' @param u Numeric vector of values sampled from a uniform distribution from 0 to 1 and with same length as W.
#' @param m Numeric. Sample minimum from fossil record.
#' @param n Number of fossils.
#' @param eps.mean Numeric mean of the measurement error distribution. Defaults to 0.
#' @param eps.sigma Numeric standard error of the measurement error distribution. Defaults to 0.
#'
#' @return A number
estimating_eqn <- function (theta, q, K, u, m, n, eps.mean, eps.sigma) {
  F.eps.m <- pnorm(m - theta, mean = eps.mean, sd = eps.sigma)
  F.eps.K <- pnorm(K - theta, mean = eps.mean, sd = eps.sigma)
  psi.hat <- estimate_psi(u = u, mean = eps.mean, sd = eps.sigma, a = -Inf, b = m - theta, K = K, m = m, theta = theta)
  return(1 - F.eps.m / F.eps.K * psi.hat - q ^ (1 / n))
}

# Helper function
estimate_psi <- function (u, mean, sd, a, b, K, m, theta) {
  e <- uniform_to_tnorm(u, mean, sd, a, b)
  psi.hat <- mean((m - e - theta) / (K - e - theta))
  return(psi.hat)
}

# Helper function
uniform_to_tnorm <- function (u, mean, sd, a, b) {
  mc.samples <- extraDistr::qtnorm(p = u, mean = mean, sd = sd, a = a, b = b)
  return(mc.samples)
}
