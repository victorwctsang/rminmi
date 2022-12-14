#' Estimate Extinction Times using the MINMI Estimator
#'
#' This function estimates extinction times using the MINMI procedure
#'
#' @param ages Numeric vector of fossil ages.
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
minmi <- function (ages, sd, K, alpha = 0.05, B = 100, .B_init = 500, .max_var = NULL) {
  B = 100 # TODO
  result <- list(lower=NULL, point = NULL, upper = NULL, B = list(lower=B, point=B, upper=B))
  n <- length(ages)
  m <- min(ages)

  q <- c(lower = alpha/2, point = 0.5, upper = 1-alpha/2)

  flag.delta_model = (all(sd == 0))

  if (flag.delta_model) {
    for (i in 1:length(q)) {
      result[[i]] <- estimate_quantile.minmi(K = K, W = ages, u=NULL, eps.sigma = sd, q = q[i])
    }
    return(result)
  }
  else {
    # Calculate number of monte carlo samples to use
    # TODO
    # Set the maximum variance
    # max_var <- 0.2 * (sd) ^ 2
    # if (!is.null(.max_var)) {
    #   max_var <- .max_var
    # }
    # if (is.null(B)) {
    #   u.init <- runif(.B_init, 0, 1)
    #   for (i in 1:length(q)) {
    #     result$B[[i]] <- find_optimal_B(max_var = max_var, K = K, m = m, n = n, u = u.init, eps.sigma = dating.sd, q = q[i])
    #   }
    # }

    # Generate Monte Carlo Samples
    B.max <- max(unlist(result$B))
    mc.samples <- matrix(runif(n*B.max, min = 0, max = 1), ncol=B.max)

    # Calculate estimates
    for (i in 1:length(q)) {
      result[[i]] <- estimate_quantile.minmi(K = K, W = ages, u = mc.samples[, 1:result$B[[i]]], eps.sigma = sd, q = q[i])
    }
    return(result)
  }
}

#' Estimate a specific quantile using the MINMI method.
#'
#' @param K Numeric upper bound for fossil ages.
#' @param W Numeric vector of fossil ages.
#' @param u Optional Matrix of numbers sampled from a uniform distribution from 0 to 1. Rows correspond to fossils and Columns correspond to Monte Carlo samples
#' @param eps.sigma Optional vector of standard errors for each fossil. Defaults to 0.
#' @param q Numeric between 0 and 1 specifying the quantile we want to estimate.
#'
#' @return A numeric indicating the MINMI estimate of the qth quantile of our sample statistic.
#' @export
estimate_quantile.minmi <- function (K, W, u=NULL, eps.sigma = 0, q) {
  n <- length(W)
  m <- min(W)
  # No Measurement Error case
  theta_q.hat <- K - q ^ (-1 / n) * (K - m)
  if (!all(eps.sigma == 0)) {
    # Measurement Error case
    newton.res <- pracma::newtonRaphson(
      fun = function(theta) estimating_eqn(theta, q, K, u, m, n, eps.sigma),
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
#' @param u Matrix of numbers sampled from a uniform distribution from 0 to 1. Rows correspond to fossils and Columns correspond to Monte Carlo samples
#' @param eps.sigma Vector of standard errors for each fossil. Length n.
#' @param q Numeric between 0 and 1 specifying the quantile we want to estimate.
#'
#' @return Vector of numbers for each fossil.
find_optimal_B <- function (max_var, K, m, n, u, eps.sigma, q) {
  # TODO
  # Initial estimate of theta_q using no measurement error case
  theta_q.hat.init <- K - q ^ (-1 / n) * (K - m)

  # pdfs and CDF evaluations (for convenience)
  f_eps.K <- dnorm(K - theta_q.hat.init, mean=0, sd=eps.sigma)
  F_eps.K <- pnorm(K - theta_q.hat.init, mean=0, sd=eps.sigma)
  # Estimate sigma^2_var
  B <- length(u)
  e <- uniform_to_tnorm(u, eps.sigma, a = -Inf, b = m - theta_q.hat.init)
  sample_var.psi_hat <- apply((m - e - theta_q.hat.init) / (K - e - theta_q.hat.init), MARGIN=1, FUN=var)

  # Estimate \hat\psi and \hat\psi prime
  psi_hat <- apply((m - e - theta_q.hat.init) / (K - e - theta_q.hat.init), MARGIN=1, FUN=mean)
  psi_hat.prime <- - apply((K - m) / (K - e - theta_q.hat.init) ^ 2, MARGIN=1, FUN=mean)

  optimal_B <- 1 / max_var * sample_var.psi_hat * (f_eps.K / F_eps.K * psi_hat + psi_hat.prime) ^ (-2)

  return(max(ceiling(optimal_B), 100))
}

#' Estimating equation for Newton-Raphson's Algorithm
#'
#' @param theta Number. Estimate for theta.
#' @param q Numeric between 0 and 1 specifying the quantile we want to estimate.
#' @param K Numeric. Upper bound for fossil record.
#' @param u Matrix of numbers sampled from a uniform distribution from 0 to 1. Rows correspond to fossils and Columns correspond to Monte Carlo samples
#' @param m Numeric. Sample minimum from fossil record.
#' @param n Number of fossils.
#' @param eps.sigma Vector of standard errors for each fossil. Length n.
#'
#' @return A number
estimating_eqn <- function (theta, q, K, u, m, n, eps.sigma) {
  F.eps.m <- pnorm(m - theta, mean = 0, sd = eps.sigma)
  F.eps.K <- pnorm(K - theta, mean = 0, sd = eps.sigma)
  e <- uniform_to_tnorm(u, eps.sigma, a = -Inf, b = m - theta)
  psi.hat <- apply((m - e - theta) / (K - e - theta), MARGIN = 1, FUN=mean)
  return(sum(log(1 - F.eps.m / F.eps.K * psi.hat)) - log(q))
}

#' Helper function to transform samples from a uniform dist to a truncated normal
#'
#' @param u Matrix of numbers. Rows = fossils; columns = MC samples
#' @param eps.sigma Vector of numbers. Standard deviations for each fossil.
#' @param a Lower bound for truncated normal
#' @param b Upper bound for truncated normal
#'
#' @return Matrix of numbers. Same dimensions as u
uniform_to_tnorm <- function (u, eps.sigma, a, b) {
  mc.samples <- matrix(extraDistr::qtnorm(p = u, mean = 0, sd = eps.sigma, a = a, b = b), ncol=ncol(u))
  return(mc.samples)
}
