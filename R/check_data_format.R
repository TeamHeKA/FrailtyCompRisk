#' Check the Format of a Data Frame for Multicentre Competing Risks with Frailty
#'
#' Validates whether the submitted data frame is correctly structured for the analysis
#' of multicentre competing risks with frailty.
#'
#' @param df A data frame to be checked.
#'
#' @return Returns \code{TRUE} if the data frame is correctly structured; otherwise, returns
#' a character string describing the cause of failure.
#'
#' @examples
#' \dontrun{
#' n_cov = 2
#' n_cluster = 5
#' n_per_cluster = 100
#' n = n_cluster * n_per_cluster
#' a1 = 0
#' b1 = 1
#' G <- rep(1:n_cluster, each = n_per_cluster)
#' Z <- matrix(runif(n * n_cov, a1, b1), ncol = n_cov)
#' df = simulate_data(
#'   G,
#'   Z,
#'   prop = 0.6,
#'   beta = c(0),
#'   theta = 0.01,
#'   cens = TRUE,
#'   pcens = 0.25,
#'   tau = 0
#' )
#' check_data_format(df)
#' }
#'
#' @export


check_data_format <- function(df)
{
  if (!is.data.frame(df)) {
    stop("The argument must be a data frame.")
  }

  n_cols <- ncol(df)
  if (n_cols < 3) {
    stop(
      "The data frame must contain at least 3 columns:\n - 'times': the observed time,\n - 'status': 0 for right-censoring, i in [1,n] otherwise, where i is the cause of failure,\n - 'clusters': cluster/group indicator."
    )
  }

  col_names <- colnames(df)
  expected_names <- c("times", "status", "clusters")
  for (i in 1:3) {
    if (col_names[i] != expected_names[i]) {
      stop(paste0("Column ", i, " must be named '", expected_names[i], "'."))
    }
  }

  if (anyNA(df$times))
    stop("Some values in 'times' are missing (NA).")
  if (anyNA(df$status))
    stop("Some values in 'status' are missing (NA).")
  if (anyNA(df$clusters))
    stop("Some values in 'clusters' are missing (NA).")
  if (n_cols > 3) {
    covariate_cols <- df[, 4:n_cols, drop = FALSE]
    if (!all(sapply(covariate_cols, is.numeric)) || any(is.na(covariate_cols))) {
      stop("All covariates (columns 4 to end) must be numeric and contain no missing values (NA).")
    }
  }
  return(TRUE)
}
