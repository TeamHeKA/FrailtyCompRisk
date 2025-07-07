#' Nelson-Aalen Estimator of Censoring Distribution (Vectorized)
#'
#' Computes the estimated censoring survival probabilities using a Nelson-Aalen-type estimator,
#' accounting for shared frailty through a cluster-specific random effect.
#'
#' @param times A numeric vector of observed times.
#' @param unique_times A numeric vector of distinct, sorted time points.
#' @param status A numeric vector where 0 indicates censoring.
#' @param clusters A vector indicating the cluster/group for each observation.
#' @param u_bis A numeric vector of cluster-specific random effects (e.g., frailty values).
#'
#' @return A numeric vector of estimated censoring survival probabilities for each individual.
#'
#' @examples
#' \dontrun{
#' n_cov <- 1
#' n_cluster <- 5
#' n_per_cluster <- 100
#' n <- n_cluster * n_per_cluster
#' G <- rep(1:n_cluster, each = n_per_cluster)
#' Z <- matrix(runif(n * n_cov), ncol = n_cov)
#' df <- simulate_data(G, Z, prop = 0.6, beta = c(0), theta = 0.01, cens = TRUE, pcens = 0.25, tau = 0)
#' check_data_format(df)
#' unique_times <- sort(unique(df$times))
#' u_bis <- rnorm(n_cluster, 0, 0.1)
#' S_C <- Nelson_Censoring_vectorized(df$times, unique_times, df$status, df$clusters, u_bis)
#' }
#'
#' @export
Nelson_Censoring_vectorized <- function(times, unique_times, status, clusters, u_bis) {
  check_data_format(data.frame(times = times, status = status, clusters = clusters))
  
  n_times <- length(unique_times)
  N <- length(times)
  
  # Step 1: Estimate cumulative hazard function H(t) at each unique time
  d_vec <- sapply(unique_times, function(t) sum(times == t & status == 0))  # censored at t
  R_vec <- sapply(unique_times, function(t) {
    ind <- which(times >= t)
    sum(exp(u_bis[clusters[ind]]))
  })
  
  lambda_hat <- ifelse(R_vec > 0, d_vec / R_vec, 0)
  H_vec <- cumsum(lambda_hat)  # Estimated Nelson-Aalen cumulative hazard
  
  # Step 2: Map each individual's time to H(t_i)
  H_ti <- H_vec[findInterval(times, unique_times)]
  
  # Step 3: Compute censoring survival probability for each individual
  S_C <- exp(-H_ti * exp(u_bis[clusters]))
  
  return(S_C)
}
