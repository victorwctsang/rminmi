#' Confidence Interval for Extinction Time using the MINMI Estimator
#'
#' Estimates a confidence interval for extinction time using the MINMI procedure, which finds a range of values for extinction time that
#' are plausible considering the date of the most recently observed fossil, accounting for sampling error (the fact that
#' the most recent fossil date is not necessarily the most recent time that the species was extant) and measurement error (error dating fossils).
#'
#' @param ages Numeric vector of fossil ages, with smaller values being more recent fossil ages.
#' @param sd Numeric vector of measurement error standard deviations for each fossil (listed in the same order as they appear in \code{ages}).
#' @param K Numeric upper bound for fossil ages - how old fossils can be before they are ignored, for the purpose of this analysis. A sensible choice of \code{B} is
#' close to the age of the oldest fossil.
#' @param alpha Numeric between 0 and 1. Used to find 100(1-\code{alpha})\% confidence intervals. Defaults to 0.05 (95\% confidence intervals)
#' @param q Numeric vector of values between 0 and 1, specifying the quantile at which we want to solve for extinction time. Defaults to \code{c(alpha/2,0.5,1-alpha/2)},
#' which gives the limits of a 100(1-\code{alpha})\% confidence interval and a point estimate obtained by solving at 0.5. If \code{q} is specified it overrides any input for \code{alpha}.
#' @param B Optional numeric greater than 1 specifying the number of Monte Carlo samples to use in quantile estimation.
#' @param A Optional numeric greater than 0 specifying the maximum Monte Carlo standard error we aim to have associated with our MINMI estimates. Defaults to
#' 10\% of the smallest non-zero value of \code{sd}.
#' @param .B_init In choosing the required sample size \code{B} to keep Monte Carlo error variance below \code{A}, the error variance must first be approximated.
#' This is done using a starting value for extinction time (assuming no measurement error) from \code{.B_init} samples
#'
#' @details 
#' The MINMI procedure involves assuming:
#' \itemize{
#' \item That measurement error for each fossil is normally distributed around the provided point estimate of fossil age, with the provided standard deviation
#' \item That fossil dates are uniformly distributed over the interval of allowable dates (which goes from estimated extinction time up until \code{K} minus measurement error)
#' }
#' We then estimate extinction time by inversion of the sample minimum, that is, we find the estimate of extinction time \eqn{\theta}{t} at quantile level \code{q}
#' such that the probability of seeing a sample minimum less than the most recent fossil date observed in the sample \code{ages} is equal to \code{q}. This function returns thee values:
#' a point estimator for extinction time (solving at \code{q=0.5}), and upper and lower limits that give us an \code{alpha}-level confidence interval.
#' 
#' It is assumed that \code{ages} has been specified with smaller values representing more recent fossils, for example, \code{ages} could be specified in years before present.
#' 
#' When there is measurement error (that is, when \code{sd} is not a vector of zeros), Monte Carlo estimation is used to find quantiles, from a set \code{B}
#' random samples. The number of Monte Carlo samples \code{B} that are used can be controlled in two ways: \code{B} can be specified directly in your function 
#' call, or you can set \code{A}, the desired Monte Carlo variance the final extinction time estimates should have. If neither is specified, \code{A} defaults to
#' a tenth of the smallest non-zero value in \code{sd}, and the required number of Monte Carlo samples to achieve this is approximated using the \link{choose_B} function. 
#' If both \code{A} and \code{B} are specified then \code{A} is ignored. The bigger \code{B} is, or the smaller \code{A} is,
#' the less Monte Carlo error there will be (and the longer this code will take to run, but it is usually pretty fast).
#' 
#' @return minmi() returns a list with estimates for the lower end point of the 100(1-\code{alpha})\% confidence interval, point estimate, upper end point, and a list containing the values of \code{B} used for each.
#' @return An object of class "minmi" with the following components:
#'
#'
#'  \item{theta}{ a vector of estimated extinction times at each of a set of quantiles specified in \code{q}.}
#'  \item{B}{ a vector containing the number of Monte Carlo samples used to estimate each \code{theta}.}
#'  \item{q}{ the vector of quantiles used in estimation.}
#'  \item{call }{ the function call}
#' @export
#' @examples
#' ages = runif(20, 10000, 25000) #simulating some random data
#' sd = runif(20, 50, 100)
#'
#' # for a point estimate plus 95% CI
#' minmi(ages=ages, sd=sd, K=22000, alpha=0.05) 
#' 
#' # finding just a point estimate, but using 200 Monte Carlo samples
#' minmi(ages=ages, sd=sd, K=22000, q=0.5, B = 200) 
#' 
#' # now using large enough B to keep Monte Carlo standard error less than 1 
#' minmi(ages=ages, sd=sd, K=22000, alpha=0.05, A = 1) 
#' 
#' @importFrom stats dnorm pnorm runif var
minmi <- function (ages, sd, K, alpha = 0.05, q = c(lower = alpha/2, point = 0.5, upper = 1-alpha/2), B = NULL, A = 0.1 * min(sd[sd>0]), .B_init = 500) {
  
  # set up result list.  Note for future work - this could include Monte Carlo standard errors too
  tmp = rep(NA,length(q))
  if(is.null(names(q)))
    names(tmp)=paste0("q=",q)
  else
    names(tmp)=names(q)
  if(is.null(B))
    Btmp = tmp
  else
  {
    Btmp = rep(B,length(q))
    names(Btmp) = names(tmp)
  }
  result <- list(theta=tmp,B=Btmp)

  n <- length(ages)
  m <- min(ages)

  flag.delta_model <- (all(sd == 0))

  if (flag.delta_model) {
    for (i in 1:length(q)) {
      result$theta[[i]] <- estimate_quantile.minmi(K = K, W = ages, u=NULL, eps.sigma = sd, q = q[i])
    }
    return(result)
  }
  else {
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
      result$theta[[i]] <- estimate_quantile.minmi(K = K, W = ages, u = mc.samples[, 1:result$B[[i]]], eps.sigma = sd, q = q[i])
    }
    result$q=q
    result$call <- match.call()
    class(result)="minmi"
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
