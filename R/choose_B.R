#' Function for Finding the Number of Monte Carlo Samples for MINMI
#'
#' @param A Numeric greater than 0 specifying the maximum Monte Carlo error variance we aim to have associated with our MINMI estimates.
#' @param K Numeric upper bound for fossil ages.
#' @param m Numeric. Sample minimum from fossil record.
#' @param n Number of fossils.
#' @param u Matrix of numbers sampled from a uniform distribution from 0 to 1. Rows correspond to fossils and Columns correspond to Monte Carlo samples
#' @param eps.sigma Vector of standard errors for each fossil. Length n.
#' @param q Numeric between 0 and 1 specifying the quantile we want to estimate.
#'
#' @return An integer representing the number of Monte Carlo samples to use.
#'
#' @importFrom stats dnorm pnorm var
choose_B <- function (A, K, m, n, u, eps.sigma, q) {
  # Initial estimate of theta_q using no measurement error case
  theta_q.hat.init <- K - q ^ (-1 / n) * (K - m)

  # PDF and CDF evaluations (for convenience)
  f_eps.m <- dnorm(m - theta_q.hat.init, mean=0, sd=eps.sigma)
  F_eps.m <- pnorm(m - theta_q.hat.init, mean=0, sd=eps.sigma)

  f_eps.K <- dnorm(K - theta_q.hat.init, mean=0, sd=eps.sigma)
  F_eps.K <- pnorm(K - theta_q.hat.init, mean=0, sd=eps.sigma)

  # Monte Carlo Samples
  e <- uniform_to_tnorm(u, eps.sigma, a = -Inf, b = m - theta_q.hat.init)

  # Estimate \hat\psi and \hat\psi prime
  psi_hat <- apply((m - e - theta_q.hat.init) / (K - e - theta_q.hat.init), MARGIN=1, FUN=mean)
  psi_hat.prime <- - apply((K - m) / (K - e - theta_q.hat.init) ^ 2, MARGIN=1, FUN=mean)

  # calculate B
  du.dtheta <- sum( (F_eps.m/F_eps.K * (f_eps.K/F_eps.K * psi_hat + psi_hat.prime))/(1 - F_eps.m/F_eps.K * psi_hat) )
  sigma.psi <- apply((m - e - theta_q.hat.init) / (K - e - theta_q.hat.init), MARGIN=1, FUN=var)
  du.dpsi <- F_eps.m / (F_eps.m * psi_hat - F_eps.K)

  B <- ceiling(1/A * du.dtheta ^ (-2) * sum( du.dpsi^2 * sigma.psi^2))

  if (B < 100) {
    warning(sprintf('Estimated number of Monte Carlo samples for q = %.3f is too small (B = %i), using 100 instead. Consider using a smaller target MCE variance (currently using A = %.3f).\n\n', q, B, A))
    B = 100
  }

  return(B)
}
