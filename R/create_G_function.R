#' Create Step Function for Estimated Survival (G)
#'
#' This function creates a right-continuous step function (vectorized) based on
#' estimated values of the censoring survival function at given unique times.
#'
#' The resulting function can be used to evaluate \eqn{G(t)} for any vector of time points.
#'
#' @param unique_times A numeric vector of increasing unique times (e.g., censoring times).
#' @param G_vals A numeric vector of the same length as `unique_times`, corresponding to estimated values of the survival function at each time.
#'
#' @return A function `G(t)` that can be evaluated on a numeric vector of time points.

create_G_function <- function(unique_times, G_vals)
{
  function(t) {
    idx <- findInterval(t, unique_times)
    idx[idx == 0] <- 1  # If t < min(unique_times), return first survival value (usually 1)
    idx[idx > length(G_vals)] <- length(G_vals)  # If t > max(unique_times), return last known value
    return(G_vals[idx])
  }
}
