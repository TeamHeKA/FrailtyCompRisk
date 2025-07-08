#' Cox Model with Random Frailty (REML Estimation)
#'
#' Estimates the regression coefficients, random effects, and frailty variance \eqn{\theta} in a Cox proportional hazards model with shared frailty using Restricted Maximum Likelihood (REML).
#'
#' @param data A data frame containing:
#'   - `times`: observed survival or censoring times,
#'   - `status`: event indicator (1 = event, 0 = censored),
#'   - `clusters`: cluster identifiers (e.g., center, group),
#'   - covariates from column 4 onward.
#' @param max_iter Maximum number of Newton-Raphson iterations (default = 300).
#' @param tol Tolerance for convergence (default = 1e-6).
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{beta}{Estimated fixed effects (regression coefficients).}
#'   \item{u}{Estimated random effects (frailties) for each cluster.}
#'   \item{theta}{Estimated variance of the random effects.}
#'   \item{p_value}{P-value for testing the null hypothesis \eqn{H_0: \theta = 0} (no frailty).}
#' }
#'
#' @details
#' The function uses a REML-based Newton-Raphson approach to estimate the parameters.
#' A penalized likelihood including a log-likelihood contribution for the frailty variance \eqn{\theta} is maximized.
#'
#' The hypothesis test of \eqn{H_0: \theta = 0} is based on a Wald-type statistic. A high p-value (e.g. > 0.05) may indicate that the random effects (frailties) are negligible.
#'
#' @importFrom Matrix Matrix Diagonal
#' @importFrom stats pchisq
#'
#' @seealso \code{\link{Ml_Cox}}, \code{\link{logLikelihood_1}}, \code{\link{logLikelihood_2}}
#'
#' @export

Reml_Cox_frailty <- function(data,max_iter=300, tol = 1e-6)
{
  times <- data$times
  status <- data$status
  clusters <- data$clusters
  p <- (length(data[1,])-3)
  N <- length(times)
  K <- length(unique(clusters))

  if (p > 0) {
    X <- Matrix(as.matrix(data[, 4:(p + 3)]), sparse = TRUE)
  } else {
    X <- Matrix(0, N, 0, sparse = TRUE)
  }

  if (p>0){
    data_ini <- data.frame(times = times,status = status,clusters = clusters, as.matrix(X))
    gamma_Cox <- Ml_Cox(data_ini)
    gamma_0 <- gamma_Cox
  }else{
    gamma_0 <- c()
  }

  u_0 <- rep(0, K)
  theta_0 <- 1
  theta_1 <- 2

  Z <- X
  Q <- Matrix(0, nrow = N, ncol = K, sparse = TRUE)
  for (i in 1:N) {
    Q[i, clusters[i]] <- 1
  }

  if (p>0){
    Y <- cbind(Z, Q)
    Y_t <- rbind(Matrix::t(Z), Matrix::t(Q))
  }else{
    Y <- Q
    Y_t <- Matrix::t(Q)
  }

  D <- Matrix(as.numeric(status == 1), ncol = 1, sparse = TRUE)

  eta <- if (p > 0) as.vector(Z %*% gamma_0 + Q %*% u_0) else as.vector(Q %*% u_0)
  W_diag <- exp(eta)
  W <- Diagonal(x = W_diag)

  M <- Matrix(0,N,N)
  for (i in 1:N){
    for (j in 1:N){
      if (times[i] <= times[j])
        M[j,i] <- 1
    }
  }

  Mt <- Matrix::t(M)

  v <- as.vector(colSums(as.matrix(M * W_diag)))
  a <- numeric(N)
  a[status == 1] <- 1 / v[status == 1]
  A <- Diagonal(x = a)

  b_vec <- as.vector(M %*% a)
  B <- Diagonal(x = b_vec)

  lik_0 <- logLikelihood_1(status, M, W) + logLikelihood_2(u_0,theta_0,K)
  lik_1 <- lik_0 + 1
  iter <- 1

  not_converged <- ((abs(lik_1 - lik_0) >= tol) && (abs(theta_1 - theta_0) >= tol))

  while ((not_converged) && (iter <= max_iter)) {
    grad <- as.vector(D) - W %*% (M %*% a)

    a_sq <- a^2
    A2 <- Diagonal(x = a_sq)
    negHess <- W %*% (B - M %*% A2 %*% Mt %*% W)

    d <- c(rep(1, p), rep(1 / theta_0, K))
    M_pen <- Diagonal(x = d)

    V <- Y_t %*% negHess %*% Y + M_pen
    inv_V <- solve(V)

    update_direction <- inv_V %*% (Y_t %*% grad) - inv_V %*% c(rep(0, p), (1 / theta_0) * u_0)

    new_par <- c(gamma_0, u_0) + update_direction

    if (p>0){
      gamma_0 <- new_par[1:p]
    }
    u_0 <- new_par[(p + 1):(p + K)]

    theta_1 <- theta_0
    trace_term <- sum(diag(inv_V[(p + 1):(p + K), (p + 1):(p + K)]))
    theta_0 <- as.numeric((t(u_0) %*% u_0) / (K - (trace_term / theta_0)))

    eta <- if (p > 0) as.vector(Z %*% gamma_0 + Q %*% u_0) else as.vector(Q %*% u_0)
    W_diag <- exp(eta)

    W <- Diagonal(x = W_diag)

    v <- as.vector(colSums(as.matrix(M * W_diag)))
    a <- numeric(N)
    a[status == 1] <- 1 / v[status == 1]
    A <- Diagonal(x = a)

    b_vec <- as.vector(M %*% a)
    B <- Diagonal(x = b_vec)

    lik_0 <- lik_1
    lik_1 <- logLikelihood_1(status, M, W) + logLikelihood_2(u_0, theta_0,K)

    not_converged <- ((abs(lik_1 - lik_0) >= tol) && (abs(theta_1 - theta_0) >= tol))
    iter <- iter + 1
  }

  A22 <- inv_V[(p + 1):(p + K), (p + 1):(p + K)]
  wald_stat <- as.numeric(t(u_0) %*% solve(A22) %*% u_0)

  df_eff <- K - sum(diag(A22)) / theta_0

  p_value <- pchisq(wald_stat, df = df_eff, lower.tail = FALSE)
  if (p < 0.05){
    message("A p-value greater than 0.05 in the test of theta=0 suggests that the cluster effect may be negligible.\n","p-value = ",p_value,"\n")
  }
  return(list(beta = gamma_0,u = u_0, theta = theta_0,p_value = p_value))
}
