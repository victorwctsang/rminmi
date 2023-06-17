#' Estimating the Required Number of Monte Carlo Samples for MINMI
#'
#' @param A Numeric greater than 0 specifying the maximum Monte Carlo error standard error we aim to have associated with our MINMI estimates.
#' @param K Numeric upper bound for fossil ages.
#' @param m Numeric. Sample minimum from fossil record.
#' @param n Number of fossils.
#' @param u Matrix of numbers sampled from a uniform distribution from 0 to 1. Rows correspond to fossils and Columns correspond to Monte Carlo samples
#' @param eps.sigma Vector of standard errors for each fossil. Length n.
#' @param q Numeric between 0 and 1 specifying the quantile we want to estimate.
#' @param theta Numeric estimate of extinction time to use in finding \code{B}. Defaults to the estimate of extinction time if there were no measurement error.
#'
#' @details This function is used internally by \link{minmi} to choose the appropriate number of Monte Carlo samples \code{B} to use, and should rarely be called directly.
#' Given a desired Monte Carlo standard error \code{A}, this function approximates the number of Monte Carlo samples \code{B} that would be needed to achieve this level
#' of precision. It assumes the true extinction time is \code{theta} in these calculations, if \code{theta} is not provided, this defaults to the estimate if we ignore measurement error.
#' If the required value of \code{B} is less than 100, it is reset to 100 to be on the safe side.
#' @return An integer representing the number of Monte Carlo samples to use.
#'
#' @importFrom stats dnorm pnorm var
choose_B <- function (A, K, m, n, u, eps.sigma, q, theta = K - q ^ (-1 / n) * (K - m)) {
  # Initial estimate of theta_q using no measurement error case

  # PDF and CDF evaluations (for convenience)
  f_eps.m <- dnorm(m - theta, mean=0, sd=eps.sigma)
  F_eps.m <- pnorm(m - theta, mean=0, sd=eps.sigma)

  f_eps.K <- dnorm(K - theta, mean=0, sd=eps.sigma)
  F_eps.K <- pnorm(K - theta, mean=0, sd=eps.sigma)

  # Monte Carlo Samples
  e <- uniform_to_tnorm(u, eps.sigma, a = -Inf, b = m - theta)

  # Estimate \hat\psi and \hat\psi prime
  psi_hat <- apply((m - e - theta) / (K - e - theta), MARGIN=1, FUN=mean)
  psi_hat.prime <- - apply((K - m) / (K - e - theta) ^ 2, MARGIN=1, FUN=mean)

  # calculate B
  du.dtheta <- sum( (F_eps.m/F_eps.K * (f_eps.K/F_eps.K * psi_hat + psi_hat.prime))/(1 - F_eps.m/F_eps.K * psi_hat) )
  sigma.psi <- apply((m - e - theta) / (K - e - theta), MARGIN=1, FUN=var)
  du.dpsi <- F_eps.m / (F_eps.m * psi_hat - F_eps.K)

  B <- ceiling(1/A^2 * du.dtheta ^ (-2) * sum( du.dpsi^2 * sigma.psi))

  if (B < 100) {
    warning(sprintf('Estimated number of Monte Carlo samples for q = %.3f is a bit small (B = %i), using 100 instead. Consider using a smaller target MC standard error (currently using A = %.2f).\n', q, B, A))
    B = 100
  }

  return(B)
}
