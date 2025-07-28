#' REML Estimation for Cause 1 in a Competing Risks Model with Shared Frailty
#'
#' Fits a cause-specific hazard model for cause 1 in a competing risks framework using restricted maximum likelihood (REML) with a shared frailty term for clusters.
#'
#' @param data A data frame with at least the following columns:
#' \describe{
#'   \item{times}{Observed times to event or censoring.}
#'   \item{status}{Event type: 0 = censored, 1 = cause of interest, >1 = other competing risks.}
#'   \item{clusters}{Cluster membership (e.g., centers or subjects).}
#'   \item{covariates}{(Optional) Columns from the 4th onward are used as covariates.}
#' }
#' @param cluster_censoring Logical. If \code{TRUE}, accounts for center-specific censoring using a frailty model on censoring times. Default is \code{FALSE}.
#' @param max_iter Maximum number of iterations for the Newton-Raphson algorithm (default = 300).
#' @param tol Convergence tolerance for the likelihood and frailty variance (default = 1e-6).
#' @param threshold Lower bound for the frailty variance parameter \eqn{\theta}. If the estimated value falls below this threshold, frailty is considered negligible (default = 1e-5).
#'
#' @return A list with:
#' \describe{
#'   \item{beta}{Estimated regression coefficients for the cause 1 hazard.}
#'   \item{u}{Estimated cluster-specific random effects (frailties).}
#'   \item{theta}{Estimated frailty variance.}
#'   \item{p_value}{P-value testing the null hypothesis \eqn{H_0: \theta = 0}.}
#' }
#'
#' @details
#' The function fits a proportional hazards model for cause 1, accounting for unobserved heterogeneity via a shared frailty term on clusters. When \code{cluster_censoring = TRUE}, a frailty-adjusted survival curve is used to correct for informative censoring using either a Kaplan-Meier or Nelson-Aalen-based estimator.
#'
#' If the estimated frailty variance \eqn{\theta} becomes smaller than the provided \code{threshold}, the model reverts to a standard Cox model for cause 1 using \code{\link{Ml_CompRisk}}.
#'
#' @seealso \code{\link{Ml_CompRisk}}, \code{\link{Reml_Cox_frailty}}, \code{\link{logLikelihood_1}}, \code{\link{compute_M_optimized}}
#'
#' @examples
#' \dontrun{
#' data <- data.frame(
#'   times = c(4, 6, 10, 12, 3),
#'   status = c(1, 0, 2, 1, 0),
#'   clusters = c(1, 1, 2, 2, 3),
#'   x1 = rnorm(5),
#'   x2 = sample(0:1, 5, replace = TRUE)
#' )
#' result <- Reml_CompRisk_frailty(data)
#' result$beta
#' result$p_value
#' }
#'
#' @importFrom stats pchisq
#' @importFrom Matrix Matrix Diagonal
#' @importFrom Matrix sparseMatrix
#' @export

Reml_CompRisk_frailty <- function(data, cluster_censoring = FALSE, max_iter = 300, tol = 1e-6, threshold = 1e-5) {
  p <- ncol(data) - 3
  times <- data$times
  status <- data$status
  N <- length(times)

  if (p > 0) {
    X <- Matrix(as.matrix(data[, 4:(p + 3)]), sparse = TRUE)
  } else {
    X <- Matrix(0, N, 0, sparse = TRUE)
  }

  clusters <- data$clusters
  K <- length(unique(clusters))
  if (K==1)
    {stop("There is only one cluster, therefore there are no cluster effect, you might use: \n - method = 'CompRisk' if you are using 'Parameters_estimation()'\n - 'Ml_CompRisk()' if you are using 'Reml_CompRisk_frailty()'")}

  if (p > 0) {
    data_ini <- data.frame(times = times, status = status, clusters = clusters, as.matrix(X))
    gamma_0 <- Ml_Cox(data_ini)$beta
  } else {
    gamma_0 <- numeric(0)
  }

  u_0 <- rep(0, K)
  theta_0 <- 1
  theta_1 <- 2

  Z <- X
  Q <- sparseMatrix(i = 1:N, j = clusters, x = 1, dims = c(N, K))

  if (p > 0) {
    Y <- cbind(Z, Q)
    Y_t <- rbind(Matrix::t(Z), Matrix::t(Q))
  } else {
    Y <- Q
    Y_t <- Matrix::t(Q)
  }

  D <- Matrix(as.numeric(status == 1), ncol = 1, sparse = TRUE)
  eta <- if (p > 0) as.vector(Z %*% gamma_0 + Q %*% u_0) else as.vector(Q %*% u_0)
  W_diag <- exp(eta)
  W <- Diagonal(x = W_diag)

  if (cluster_censoring && !(0 %in% status)) {
    cluster_censoring <- FALSE
  }

  if (!cluster_censoring) {
    M <- compute_M_optimized(times, status, "KaplanMeier_Censoring_vectorized", u_bis = rep(0, 1), clusters)
  } else {
    d_cens <- as.integer(status == 0)
    data_bis <- data.frame(times = times, status = d_cens, clusters = clusters, matrix(0, N, 0))
    l <- Reml_Cox_frailty(data_bis)
    u_bis <- l$u
    M <- compute_M_optimized(times, status, "Nelson_Censoring_vectorized", u_bis, clusters)
  }

  Mt <- Matrix::t(M)
  v <- as.vector(Matrix::colSums(M * W_diag))
  a <- numeric(N)
  a[status == 1] <- 1 / v[status == 1]
  A <- Diagonal(x = a)
  b_vec <- as.vector(M %*% a)
  B <- Diagonal(x = b_vec)

  lik_0 <- logLikelihood_1(status, M, W) + logLikelihood_2(u_0, theta_0, K)
  lik_1 <- lik_0 + 1
  iter <- 1

  while ((abs(lik_1 - lik_0) >= tol) && (abs(theta_1 - theta_0) >= tol) && iter <= max_iter) {
    grad <- as.vector(D) - W %*% (M %*% a)
    A2 <- Diagonal(x = a^2)
    negHess <- W %*% (B - M %*% A2 %*% Mt %*% W)

    d <- c(rep(1, p), rep(1 / theta_0, K))
    M_pen <- Diagonal(x = d)

    V <- Y_t %*% negHess %*% Y + M_pen
    rhs_vec <- as.numeric(Y_t %*% grad) - c(rep(0, p), (1 / theta_0) * u_0)
    update_direction <- solve(V, rhs_vec)


    new_par <- c(gamma_0, u_0) + update_direction
    if (p > 0) gamma_0 <- new_par[1:p]
    u_0 <- new_par[(p + 1):(p + K)]

    theta_1 <- theta_0
    invV_diag <- diag(solve(V)[(p + 1):(p + K), (p + 1):(p + K)])
    trace_term <- sum(invV_diag)
    theta_0 <- as.numeric((sum(u_0^2)) / (K - (trace_term / theta_0)))

    if (theta_0 < threshold) {
      stop("The variance of the cluster effect is inferior to the threshold ('threshold' =",threshold,"), we suggest you use the method without frailty:\n - method = 'CompRisk' if you are using 'Parameters_estimation()'\n - 'Ml_CompRisk()' if you are using 'Reml_CompRisk_frailty()'")
    }

    eta <- if (p > 0) as.vector(Z %*% gamma_0 + Q %*% u_0) else as.vector(Q %*% u_0)
    W_diag <- exp(eta)
    W <- Diagonal(x = W_diag)

    v <- as.vector(Matrix::colSums(M * W_diag))
    a <- numeric(N)
    a[status == 1] <- 1 / v[status == 1]
    A <- Diagonal(x = a)
    b_vec <- as.vector(M %*% a)
    B <- Diagonal(x = b_vec)

    lik_0 <- lik_1
    lik_1 <- logLikelihood_1(status, M, W) + logLikelihood_2(u_0, theta_0, K)
    iter <- iter + 1
  }

  V_inv <- solve(V)
  A22 <- V_inv[(p + 1):(p + K), (p + 1):(p + K), drop = FALSE]
  wald_stat <- as.numeric(t(u_0) %*% solve(A22, u_0))
  trace_A22 <- sum(diag(as.matrix(A22)))
  df_eff <- K - trace_A22 / theta_0


  suppressWarnings({
      p_value <- pchisq(wald_stat, df = df_eff, lower.tail = FALSE)
  })
  if (is.nan(p_value)){
    message("The estimated value of theta is too small to allow a reliable statistical test of whether theta = 0, hence the p-value is 'NaN'. \nThis suggests that the cluster effect may be negligible.")
  }

  if (!is.na(p_value) && p_value >= 0.05) {
    message("A p-value greater than 0.05 in the test of theta=0 suggests that the cluster effect may be negligible.\n p-value = ", p_value)
  }

  return(list(beta = gamma_0, u = u_0, theta = theta_0, p_value = p_value))
}
