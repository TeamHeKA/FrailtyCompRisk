#' Maximum Likelihood Estimation for the Cox Model
#'
#' This function computes the maximum likelihood estimates of the regression coefficients in the Cox proportional hazards model using Newton-Raphson iterations.
#'
#' @param data A data frame containing the columns:
#'   - `times`: observed survival or censoring times,
#'   - `status`: event indicators (1 if event occurred, 0 otherwise),
#'   - Covariates in columns 4 to the end.
#' @param max_iter Maximum number of iterations for the Newton-Raphson algorithm (default = 100).
#' @param tol Convergence tolerance based on the Euclidean norm of coefficient changes (default = 1e-6).
#'
#' @return A numeric vector of estimated regression coefficients \code{beta}.
#'
#' @examples
#' \dontrun{
#' n <- 200
#' G <- rep(1:10, each = 20)
#' Z <- matrix(rnorm(n), ncol = 1)
#' df <- simulate_data(G, Z, prop = 0.6, beta = c(0.5), theta = 0.01, cens = TRUE, pcens = 0.25)
#' check_data_format(df)
#' Ml_Cox(df)
#' }
#'
#' @seealso \code{\link{simulate_data}}, \code{\link{check_data_format}}
#' @export
Ml_Cox <- function(data, max_iter=100, tol=1e-6)
{
  times <- data$times
  status <- data$status
  X <- as.matrix(data[,4:(length(data[1,]))])

  N <- length(times)
  p <- ncol(X)
  beta <- rep(0, p)

  ord <- order(times)
  ord_times <- times[ord]
  ord_by_time_status <- status[ord]
  X_bis <- X[ord, , drop=FALSE]

  for (iter in 1:max_iter) {
    eta <- as.vector(X_bis %*% beta)
    exp_eta <- exp(eta)

    loglik <- 0
    U <- rep(0, p)
    I <- matrix(0, p, p)

    for (i in which(ord_by_time_status == 1)) {
      R <- which(ord_times >= ord_times[i])

      sum_exp_eta <- sum(exp_eta[R])
      sum_exp_eta_X <- colSums(exp_eta[R] * X_bis[R, , drop=FALSE])
      sum_exp_eta_XX <- t(X_bis[R, , drop=FALSE]) %*% (exp_eta[R] * X_bis[R, , drop=FALSE])

      loglik <- loglik + (eta[i] - log(sum_exp_eta))
      U <- U + (X_bis[i, ] - sum_exp_eta_X / sum_exp_eta)
      I <- I + (sum_exp_eta_XX / sum_exp_eta) - (sum_exp_eta_X %*% t(sum_exp_eta_X)) / (sum_exp_eta^2)
    }

    step <- solve(I, U)
    beta_new <- beta + step

    if (sqrt(sum(step^2)) < tol) {
      break
    }
    beta <- beta_new
  }

  return(beta = beta)
}
