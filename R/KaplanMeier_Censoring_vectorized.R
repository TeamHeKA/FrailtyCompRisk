#' Kaplan-Meier Estimator for Censoring Distribution (Vectorized)
#'
#' Computes the Kaplan-Meier estimator of the censoring distribution at given time points.
#'
#' @param times A numeric vector of observed times.
#' @param status A numeric vector of status indicators (0 = censored, >0 = event).
#' @param unique_times A numeric vector of time points at which to evaluate the estimator.
#'
#' @return A numeric vector of estimated survival probabilities at each unique time.

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
