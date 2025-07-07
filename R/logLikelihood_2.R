#' Compute Log-Likelihood Penalty for Random Effects
#'
#' This function computes the second component of the log-likelihood (penalty term)
#' associated with the random effects (frailty terms) in a mixed-effects model.
#'
#' The penalty is derived from the Gaussian assumption on the random effects and
#' involves both the dimension of the random effects vector and the variance parameter \code{theta}.
#'
#' @param u A numeric vector of length \code{K} containing the random effects (frailties).
#' @param theta A positive numeric scalar, representing the variance component of the random effects.
#' @param K An integer representing the number of clusters (i.e., the length of \code{u}).
#'
#' @return A numeric scalar representing the log-likelihood penalty contribution from the random effects.
#'
#' @examples
#' u <- rnorm(5, 0, 1)
#' theta <- 0.5
#' K <- length(u)
#' logLikelihood_2(u, theta, K)
#'
#' @export
logLikelihood_2 = function(u, theta, K) {
  penalty <- -(1/2) * (K * log(2 * pi * theta) + (1 / theta) * (t(u) %*% u))
  return(penalty)
}
