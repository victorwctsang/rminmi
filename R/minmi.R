#' Confidence Interval for Extinction Time using the MINMI Estimator
#'
#' Estimates a confidence interval for extinction time using the MINMI procedure, taking into account sampling error (the fact that
#' the most recent fossil date is not necessarily the most recent time that the species was extant) and measurement error (error dating fossils).
#'
#' @param ages Numeric vector of fossil ages.
#' @param sd Numeric vector of measurement error standard deviations for each fossil (listed in the same order as they appear in `ages`.
#' @param K Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis.
#' @param alpha Numeric between 0 and 1. Used to find 100(1-alpha)\% confidence intervals. Defaults to 0.05 (95\% confidence intervals)
#' @param B Optional numeric greater than 1 specifying the number of Monte Carlo samples to use.
#' @param .B_init Optional numeric greater than 1 specifying the number of Monte Carlo samples to use in the pilot estimates. Defaults to 500.
#' @param A Optional numeric greater than 0 specifying the maximum Monte Carlo error variance we aim to have associated with our MINMI estimates.
#'
#' @details 
#' The MINMI procedure involves assuming:
#' - That measurement error for each fossil is normally distributed around the provided point estimate of fossil age, with the provided standard deviation
#' - That fossil dates are uniformly distributed over the interval of allowable dates (which goes from estimated extinction time up until `K` minus measurement error)
#' We then estimate extinction time by inversion of the sample minimum, that is, we find the estimate of extinction time \eqn{\theta}{t} at quantile level `q`
#' such that the probability of seeing a sampling minimum less than the smallest fossil date observed in the sample `ages` is equal to `q`. This function returns thee values:
#' a point estimator for extinction time (solving at `q=0.5`), and upper and lower limits that give us an `alpha`-level confidence interval.
#' 
#' When there is measurement error (that is, when `sd` is not a vector of zeros), Monte Carlo estimation is used, sampling a set `B` of fossil datasets.
#' For each fossil, we first simulate a measurement error value `w` by sampling from a normal distributions centered on zero with standard deviation determined by the relevant entry in `sd`.
#' Then we simulate fossil dates by sampling at random from a uniform distribution over [\eqn{\theta}{t},'K'-'w']. The number of Monte Carlo datasets `B`
#' can be controlled in two ways: `B` can be specified directly in your function call (otherwise its default value is 500), or you can set `A`, the
#' desired Monte Carlo variance the final extinction time estimates should have. If both `A` and `B` are specified then `A` is ignored. The bigger `B` is, or the smaller `A` is,
#' the less Monte Carlo error there will be (and the longer this code will take to run, but it is usually pretty fast)/
#' 
#' @returns minmi() returns a list with estimates for the lower end point of the 100(1-alpha)\% confidence interval, point estimate, upper end point, and a list containing the B's used for each.
#' @export
#' @examples
#' ages = runif(20, 10000, 25000) #simulating some random data
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
      message(sprintf('No value for A provided, using 20%% of the minimum variance instead (A = %.4f)', A))
    }

    # Choose B
    if (is.null(B)) {
      u.init <- matrix(runif(n*.B_init, 0, 1), ncol=.B_init)
      for (i in 1:length(q)) {
        result$B[[i]] <- choose_B(A=A, K=K, m=m, n=n, u=u.init, eps.sigma=sd, q=q[i]) # TODO
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
  if (!all(eps.sigma == 0)) {
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
