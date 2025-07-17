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
