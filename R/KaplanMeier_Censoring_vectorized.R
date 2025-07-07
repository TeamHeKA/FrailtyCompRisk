#' Kaplan-Meier Estimator for Censoring Distribution (Vectorized)
#'
#' Computes the Kaplan-Meier estimator of the censoring distribution at given time points.
#'
#' @param times A numeric vector of observed times.
#' @param status A numeric vector of status indicators (0 = censored, >0 = event).
#' @param unique_times A numeric vector of time points at which to evaluate the estimator.
#'
#' @return A numeric vector of estimated survival probabilities at each unique time.
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
#' KM <- KaplanMeier_Censoring_vectorized(df$times, df$status, unique_times)
#' }
#'
#' @export
KaplanMeier_Censoring_vectorized <- function(times, status, unique_times) {

  n <- length(unique_times)
  surv <- numeric(n)
  surv[1] <- 1

  for (i in seq_len(n)) {
    t <- unique_times[i]
    d <- sum(times == t & status == 0)  # number of censored at time t
    R <- sum(times >= t)                # number at risk at time t
    if (R > 0) {
      surv[i] <- ifelse(i == 1, 1, surv[i - 1] * (1 - d / R))
    } else {
      surv[i] <- ifelse(i == 1, 1, surv[i - 1])
    }
  }

  return(surv)
}
