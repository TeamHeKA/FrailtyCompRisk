#' Maximum Likelihood Estimation for the Cox Model (Optimized)
#'
#' This function computes the maximum likelihood estimates of the regression coefficients in the Cox proportional hazards model using Newton-Raphson iterations.
#'
#' @param data A data frame with columns: `times`, `status`, and covariates starting from the 4th column.
#' @param max_iter Maximum number of iterations (default = 100).
#' @param tol Convergence tolerance for the Euclidean norm (default = 1e-6).
#' @return A numeric vector of estimated regression coefficients.
#' @export
Ml_Cox <- function(data, max_iter = 100, tol = 1e-6) {
  times <- data$times
  status <- data$status
  if (("clusters" %in% names(data)) && ncol(data) == 3){stop("With no covariables, there is nothing to estimate")}
  if (!("clusters" %in% names(data)) && ncol(data) == 2){stop("With no covariables, there is nothing to estimate")}
  if (("clusters" %in% names(data)) && ncol(data) > 3){X <- as.matrix(data[, 4:ncol(data)])}
  if (!("clusters" %in% names(data)) && ncol(data) > 2){X <- as.matrix(data[, 3:ncol(data)])}
  N <- length(times)
  p <- ncol(X)

  beta <- rep(0, p)

  # Tri par temps croissants
  ord <- order(times)
  times <- times[ord]
  status <- status[ord]
  X <- X[ord, , drop = FALSE]

  for (iter in seq_len(max_iter)) {
    eta <- drop(X %*% beta)
    exp_eta <- exp(eta)

    rev_exp_eta <- rev(cumsum(rev(exp_eta)))
    rev_exp_eta_X <- apply(X * exp_eta, 2, function(col) rev(cumsum(rev(col))))

    # Calcule vectorisé de S2
    X_eta <- X * sqrt(exp_eta)
    S2_array <- array(0, dim = c(p, p, N))
    for (i in seq_len(p)) {
      for (j in seq_len(p)) {
        S2_array[i, j, ] <- rev(cumsum(rev(X[, i] * X[, j] * exp_eta)))
      }
    }

    I <- matrix(0, p, p)
    U <- numeric(p)

    for (i in which(status == 1)) {
      S0 <- rev_exp_eta[i]
      S1 <- rev_exp_eta_X[i, ]
      S2 <- S2_array[, , i]

      U <- U + (X[i, ] - S1 / S0)
      I <- I + (S2 / S0) - (tcrossprod(S1) / S0^2)
    }

    step <- tryCatch(solve(I, U), error = function(e) rep(0, p))
    beta_new <- beta + step

    if (sqrt(sum(step^2)) < tol) {
      break
    }
    beta <- beta_new
  }

  return(list(beta = beta))
}

