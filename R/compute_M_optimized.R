#' Compute the Weight Matrix M (Optimized)
#'
#' This function computes the matrix of inverse probability of censoring weights (IPCW),
#' used in multicenter competing risks models with frailty. It supports both Kaplan-Meier and Nelson-Aalen estimators
#' for the censoring distribution.
#'
#' @param times A numeric vector of observed times (events or censoring).
#' @param status A numeric vector of status indicators (0 = censored, 1 or 2 = cause of failure).
#' @param f A string, either `"KaplanMeier_Censoring_vectorized"` or `"Nelson_Censoring_vectorized"`,
#' specifying the method to estimate the censoring distribution.
#' @param u_bis A numeric vector of cluster-level random effects (used only if `f == "Nelson_Censoring_vectorized"`).
#' @param clusters A vector indicating cluster membership for each subject (same length as `times`).
#'
#' @return A square matrix `M` of size N x N, where `M[i, j]` is the IPCW weight between subject i and subject j.
#'
#' @details
#' - If `f == "KaplanMeier_Censoring_vectorized"`, then `u_bis` is ignored.
#' - The weight matrix accounts for pairwise contributions, depending on event times and causes.
#'
#' @seealso [KaplanMeier_Censoring_vectorized()], [Nelson_Censoring_vectorized()], [create_G_function()]
#'
#' @examples
#' n <- 100
#' times <- rexp(n, 0.2)
#' status <- sample(0:2, n, replace = TRUE)
#' clusters <- sample(1:5, n, replace = TRUE)
#' u_bis <- rnorm(5)
#' M1 <- compute_M_optimized(times, status, "KaplanMeier_Censoring_vectorized", NULL, clusters)
#' M2 <- compute_M_optimized(times, status, "Nelson_Censoring_vectorized", u_bis, clusters)
#'
#' @export
compute_M_optimized <- function(times, status, f, u_bis, clusters)
{
  N <- length(times)
  unique_times <- sort(unique(times))

  # Censoring survival estimation
  if (f == "KaplanMeier_Censoring_vectorized") {
    G_vals <- KaplanMeier_Censoring_vectorized(times, status, unique_times)
  }
  if (f == "Nelson_Censoring_vectorized") {
    G_vals <- Nelson_Censoring_vectorized(times, unique_times, status, clusters, u_bis)
  }
  G_func <- create_G_function(unique_times, G_vals)

  # Matrix form for efficient vectorized comparison
  T_mat_i <- matrix(rep(times, each = N), nrow = N)
  T_mat_j <- matrix(rep(times, times = N), nrow = N)
  fail_mat_j <- matrix(rep(status, times = N), nrow = N)

  cond1 <- (T_mat_j <= T_mat_i) & (fail_mat_j > 1)
  cond2 <- (T_mat_j >= T_mat_i)
  cond3 <- (T_mat_j < T_mat_i) & (fail_mat_j <= 1)

  G_T_i <- G_func(times)
  G_T_j <- G_func(times)
  G_T_i_mat <- matrix(rep(G_T_i, each = N), nrow = N)
  G_T_j_mat <- matrix(rep(G_T_j, times = N), nrow = N)

  M <- Matrix(0, nrow = N, ncol = N)
  M[cond1] <- G_T_i_mat[cond1] / G_T_j_mat[cond1]
  M[cond2] <- 1
  M[cond3] <- 0

  return(M)
}
