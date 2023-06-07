#' Estimate Extinction Times using the MINMI Estimator
#'
#' A generic function for estimating extinction times using the MINMI procedure
#'
#' @param ages Numeric vector of fossil ages.
#' @param sd Numeric vector of measurement error uncertainty associated with each fossil.
#' @param K Numeric upper bound for fossil ages.
#' @param alpha Numeric between 0 and 1. Used to find 100(1-alpha)\% confidence intervals. Defaults to 0.05 (95\% confidence intervals)
#' @param B Optional numeric greater than 1 specifying the number of Monte Carlo samples to use.
#' @param .B_init Optional numeric greater than 1 specifying the number of Monte Carlo samples to use in the pilot estimates. Defaults to 500.
#' @param A Optional numeric greater than 0 specifying the maximum Monte Carlo error variance we aim to have associated with our MINMI estimates.
#'
#' @details Note that providing a value for `B` means that the value for `A` is made redundant, as `A` is used to find `B` but not the other way around.
#' @returns minmi() returns a list with estimates for the lower end point of the 100(1-alpha)\% confidence interval, point estimate, upper end point, and a list containing the B's used for each.
#' @export
#' @examples
#' ages = runif(20, 10000, 25000)
#' sd = runif(20, 50, 100)
#'
#' minmi(ages=ages, sd=sd, K=22000, alpha=0.05)
#' minmi(ages=ages, sd=sd, K=22000, alpha=0.05, B = 100)
#' minmi(ages=ages, sd=sd, K=22000, alpha=0.05, A = 1000)
#' @importFrom stats dnorm pnorm runif var
minmi <- function (ages, sd, K, alpha = 0.05, B = NULL, A = NULL, .B_init = 500) {
  result <- list(lower=NULL, point = NULL, upper = NULL, B = list(lower=B, point=B, upper=B))
  n <- length(ages)
  m <- min(ages)

  q <- c(lower = alpha/2, point = 0.5, upper = 1-alpha/2)

  flag.delta_model <- (all(sd == 0))

  if (flag.delta_model) {
    for (i in 1:length(q)) {
      result[[i]] <- estimate_quantile.minmi(K = K, W = ages, u=NULL, eps.sigma = sd, q = q[i])
    }
    return(result)
  }
  else {
    # Set A (our target MCE Variation)
    if (is.null(A)) {
      A <- 0.2 * min(sd^2)
      message(sprintf('No value for A provided, using 20%% of the minimum variance instead (A = %.4f).\n\n', A))
    }

    # Choose B
    if (is.null(B)) {
      u.init <- matrix(runif(n*.B_init, 0, 1), ncol=.B_init)
      for (i in 1:length(q)) {
        result$B[[i]] <- choose_B(A=A, K=K, m=m, n=n, u=u.init, eps.sigma=sd, q=q[i])
      }
    }

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
  if (!all(eps.sigma == 0) & !is.null(u)) {
    # Measurement Error case
    newton.res <- pracma::newtonRaphson(
      fun = function(.th) estimating_eqn(.th, q, K, u, m, n, eps.sigma),
      x0 = theta_q.hat
    )
    theta_q.hat <- newton.res$root
  }
  return(theta_q.hat)
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
