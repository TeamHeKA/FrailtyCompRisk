#' Log-Likelihood Contribution for Cause 1
#'
#' Computes the contribution to the log-likelihood from individuals who experienced event type 1,
#' using a weighting matrix `W` and an IPCW matrix `M`.
#'
#' @param status A numeric vector indicating the event type:
#' `0` for censored, `1` for cause 1, `2` for other causes.
#' @param M A numeric matrix (usually computed with \code{\link{compute_M_optimized}}), representing
#' IPCW weights between individuals.
#' @param W A sparse matrix (e.g., of class \code{dgCMatrix} from the \pkg{Matrix} package)
#' representing estimated weights or hazards.
#'
#' @return A numeric scalar: the log-likelihood contribution from all individuals who failed from cause 1.
#'
#' @details
#' For each individual `i` such that `status[i] == 1`, the function computes:
#' \deqn{
#' \log W_{ii} - \log\left(\sum_j M_{ij} W_{ij}\right)
#' }
#' and sums this across all such individuals.
#'
#' @seealso [compute_M_optimized()], [Matrix::Matrix-class]


logLikelihood_1 = function(status, M, W)
{
  loglik <- 0
  for (i in which(status == 1)) {
    v <- log(sum(M[i, ] * W@x))
    loglik = loglik + log(W[i, i]) - v
  }
  return(loglik)
}
