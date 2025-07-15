#' Unified Interface for Parameter Estimation in (Competing Risks) Frailty Models
#'
#' Provides a unified wrapper to estimate model parameters under four different modeling strategies:
#' \itemize{
#'   \item Competing risks with shared frailty for cause 1 (`"CompRisk_frailty"`)
#'   \item Competing risks without frailty for cause 1 (`"CompRisk"`)
#'   \item Standard Cox model with shared frailty (`"Cox_frailty"`)
#'   \item Standard Cox model without frailty (`"Cox"`)
#' }
#'
#' @param data A data frame containing at least the following columns:
#' \describe{
#'   \item{times}{Event or censoring times.}
#'   \item{status}{Event indicator: 0 = censored, 1 = event of interest, >1 = other types (for competing risks).}
#'   \item{clusters}{Cluster/group variable (e.g., center ID).}
#'   \item{covariates}{(Optional) Additional covariates from the 4th column onward.}
#' }
#' @param method Character string indicating the estimation method. One of:
#' \itemize{
#'   \item \code{"CompRisk_frailty"} (default): Competing risks model for cause 1 with shared frailty.
#'   \item \code{"CompRisk"}: Competing risks model for cause 1 without frailty.
#'   \item \code{"Cox_frailty"}: Cox proportional hazards model with shared frailty.
#'   \item \code{"Cox"}: Standard Cox proportional hazards model.
#' }
#' @param cluster_censoring Logical. If \code{TRUE}, adjusts for cluster-specific censoring. Only applicable when \code{method = "CompRisk_frailty"} (default = \code{FALSE}).
#' @param max_iter Maximum number of iterations for the selected model estimation (default = 300).
#' @param tol Convergence tolerance (default = 1e-6).
#' @param threshold Lower bound for the frailty variance parameter \eqn{\theta}. If the estimated value falls below this threshold, frailty is considered negligible (default = 1e-5).

#'
#' @return A list of results from the selected estimation method, typically including estimated regression coefficients, frailty variance (if applicable), random effects (if applicable), and a p-value testing the frailty variance.
#' Returns \code{NULL} in case of error.
#'
#' @details
#' This wrapper allows seamless switching between various modeling frameworks. It automatically calls:
#' \itemize{
#'   \item \code{\link{Reml_CompRisk_frailty}} for `"CompRisk_frailty"`
#'   \item \code{\link{Ml_CompRisk}} for `"CompRisk"`
#'   \item \code{\link{Reml_Cox_frailty}} for `"Cox_frailty"`
#'   \item \code{\link{Ml_Cox}} for `"Cox"`
#' }
#'
#' Data format is validated via the helper function \code{\link{check_data_format}}. An error is raised if input format is invalid or if `cluster_censoring` is used with an incompatible method.
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
#' result <- Parameters_estimation(data, method = "CompRisk_frailty")
#' result$beta
#' }
#'
#' @seealso \code{\link{Reml_CompRisk_frailty}}, \code{\link{Ml_CompRisk}}, \code{\link{Reml_Cox_frailty}}, \code{\link{Ml_Cox}}, \code{\link{check_data_format}}
#'
#' @export

Parameters_estimation <- function(data,method = "CompRisk_frailty",cluster_censoring = F,max_iter = 300, tol = 1e-6,threshold = 1e-6)
{
  if (!(method %in% c("CompRisk_frailty","CompRisk","Cox_frailty","Cox")))
  {
    stop("The method chosen must be : 'CompRisk_frailty','CompRisk','Cox_frailty' or 'Cox'.")
  }
  if (!(check_data_format(data)) == TRUE)
  {
    stop()
  }
  if (cluster_censoring && (method != "CompRisk_frailty"))
  {stop("'cluster_censoring' is only for 'CompRisk_frailty' method.")}

  if (method == "CompRisk_frailty")
  {
    res = tryCatch(Reml_CompRisk_frailty(data,cluster_censoring,max_iter = 300, tol = 1e-6,threshold = 1e-6),
                   error = function(e) {
                     message(e$message)
                     return(NULL)})
  }
  if (method == "CompRisk")
  {
    res = tryCatch(Ml_CompRisk(data,max_iter = 300, tol = 1e-6),
                   error = function(e) {
                     message(e$message)
                     return(NULL)})
  }
  if (method == "Cox_frailty")
  {
    res = tryCatch(Reml_Cox_frailty(data,max_iter = 300, tol = 1e-6),
                   error = function(e) {
                     message(e$message)
                     return(NULL)})
  }
  if (method == "Cox")
  {
    res = tryCatch(Ml_Cox(data,max_iter = 300, tol = 1e-6),
                   error = function(e) {
                     message(e$message)
                     return(NULL)})
  }
  return(res)
}
