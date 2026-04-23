#' Maximum Likelihood Estimation for Cause 1 in Competing Risks
#'
#' Performs maximum likelihood estimation (MLE) for the effect of covariates on the cause-specific hazard function of cause 1 in a competing risks framework.
#'
#' @param data A data frame containing:
#' \describe{
#'   \item{times}{Observed survival or censoring times.}
#'   \item{status}{Event indicator: 0 = censored, 1 = cause 1 (event of interest), 2+ = other causes.}
#'   \item{Covariates}{All columns from the 4th onward are considered as covariates.}
#' }
#' @param max_iter Maximum number of Newton-Raphson iterations (default = 100).
#' @param tol Convergence tolerance for the change in parameter estimates (default = 1e-6).
#'
#' @return A named list with one element:
#' \describe{
#'   \item{beta}{A numeric vector of estimated regression coefficients for the
#'   subdistribution hazard of cause 1. Each coefficient represents the log
#'   effect of the corresponding covariate on the subdistribution hazard.}
#' }
#'
#' @details
#' The function fits a Cox-type model for the subdistribution hazard of cause 1, treating all other causes (including censoring) as censored (i.e., status != 1 becomes 0).
#' It uses a Newton-Raphson procedure to maximize the partial likelihood.
#'
#' The design matrix is reordered internally by event time to facilitate risk set computation.
#'
#' @seealso \code{\link{Reml_Cox_frailty}}, \code{\link{Ml_Cox}}, \code{\link{logLikelihood_1}}
#'
#' @examples
#' \donttest{
#' data <- data.frame(
#'  times = c(0.099, 0.006, 0.260, 0.151, 0.262,
#'            0.019, 0.026, 0.241, 0.017, 0.195,
#'            0.007, 0.057, 0.241, 0.132, 0.044,
#'            0.242, 0.416, 0.096, 0.172, 0.015),
#'  status = c(1, 0, 1, 1, 0,
#'             0, 0, 0, 0, 0,
#'             0, 0, 0, 0, 0,
#'             0, 0, 0, 0, 0),
#'  clusters = c(1,1,1,1,1,
#'               1,1,1,1,1,
#'              2,2,2,2,2,
#'               2,2,2,2,2),
#'  V1 = c(0.720, 0.289, 0.382, 0.995, 0.045,
#'         0.048, 0.042, 0.810, 0.964, 0.781,
#'         0.443, 0.303, 0.902, 0.596, 0.462,
#'         0.521, 0.656, 0.164, 0.261, 0.635)
#')
#'   beta_hat <- Ml_CompRisk(data)
#' }
#'
#' @export
Ml_CompRisk <- function(data, max_iter = 100, tol = 1e-6)
{
  times <- data$times
  status <- data$status
  if (("clusters" %in% names(data)) && ncol(data) == 3){stop("With no covariables, there is nothing to estimate")}
  if (!("clusters" %in% names(data)) && ncol(data) == 2){stop("With no covariables, there is nothing to estimate")}
  if (("clusters" %in% names(data)) && ncol(data) > 3){X <- as.matrix(data[, 4:ncol(data)])}
  if (!("clusters" %in% names(data)) && ncol(data) > 2){X <- as.matrix(data[, 3:ncol(data)])}
  event_of_interest <- (status == 1)
  status_mod <- ifelse(event_of_interest, 1, 0)
  N <- length(times)
  p <- ncol(X)
  beta <- rep(0, p)
  ord <- order(times)
  times_ord <- times[ord]
  status_ord <- status_mod[ord]
  X_ord <- X[ord, , drop = FALSE]
  for (iter in 1:max_iter) {
    eta <- as.vector(X_ord %*% beta)
    exp_eta <- exp(eta)

    U <- rep(0, p)
    I <- matrix(0, p, p)

    for (i in which(status_ord == 1)) {
      R <- which(times_ord >= times_ord[i])

      sum_exp_eta <- sum(exp_eta[R])
      sum_exp_eta_X <- colSums(exp_eta[R] * X_ord[R, , drop = FALSE])
      sum_exp_eta_XX <- t(X_ord[R, , drop = FALSE]) %*% (exp_eta[R] * X_ord[R, , drop = FALSE])

      U <- U + (X_ord[i, ] - sum_exp_eta_X / sum_exp_eta)
      I <- I + (sum_exp_eta_XX / sum_exp_eta) - (sum_exp_eta_X %*% t(sum_exp_eta_X)) / (sum_exp_eta^2)
    }

    step <- solve(I, U)
    beta_new <- beta + step

    if (sqrt(sum(step^2)) < tol) {
      break
    }
    beta <- beta_new
  }

  return(list(beta = beta))
}
