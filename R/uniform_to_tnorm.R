#' Helper function to transform samples from a uniform dist to a truncated normal
#'
#' @param u Matrix of numbers. Rows = fossils; columns = MC samples
#' @param eps.sigma Vector of numbers. Standard deviations for each fossil.
#' @param a Lower bound for truncated normal
#' @param b Upper bound for truncated normal
#'
#' @return Matrix of numbers. Same dimensions as u
uniform_to_tnorm <- function (u, eps.sigma, a, b) {
  mc.samples <- matrix(extraDistr::qtnorm(p = c(u), mean = 0, sd = eps.sigma, a = a, b = b), ncol=ncol(u)) # TODO fix this janky matrix unpacking
  return(mc.samples)
}
