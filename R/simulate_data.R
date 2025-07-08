#' Simulate clustered competing risks data with shared frailty
#'
#' This function simulates clustered time-to-event data under a competing risks framework
#' with shared frailty and the possibility of random censoring. Cause 1 is modeled through
#' an inversion method of the subdistribution function, and cause 2 is introduced when inversion fails.
#'
#' @param G A vector of group or cluster identifiers (length \code{N}). Each value indicates which cluster the individual belongs to.
#' @param Z A matrix of covariates (dimensions \code{N x p}). Can be set to \code{NULL} if no covariates are used.
#' @param prop Proportion of individuals susceptible to cause 1 (default: \code{0.6}).
#' @param beta A numeric vector of regression coefficients, one per covariate (length \code{p}).
#' @param theta Variance of the shared frailty term for event times (cause 1).
#' @param cens Logical, indicating whether censoring should be simulated (default: \code{FALSE}).
#' @param pcens Target censoring proportion (default: \code{0.25}).
#' @param tau Variance of the shared frailty term for censoring times (default: \code{0}).
#'
#' @return A \code{data.frame} with columns:
#' \describe{
#'   \item{times}{Observed time (either event time or censoring time).}
#'   \item{status}{Event indicator: \code{0} for censored, \code{1} for cause 1, \code{2} for cause 2.}
#'   \item{clusters}{Cluster ID for each individual.}
#'   \item{Cov}{Optional covariate matrix (if any).}
#' }
#'
#' @examples
#' \dontrun{
#' n_cov <- 1
#' n_cluster <- 5
#' n_per_cluster <- 100
#' n <- n_cluster * n_per_cluster
#'
#' G <- rep(1:n_cluster, each = n_per_cluster)
#' Z <- matrix(runif(n * n_cov), ncol = n_cov)
#'
#' df <- simulate_data(
#'   G = G,
#'   Z = Z,
#'   prop = 0.6,
#'   beta = c(0.5),
#'   theta = 0.01,
#'   cens = TRUE,
#'   pcens = 0.25,
#'   tau = 0.01
#' )
#'
#' head(df)
#' table(df$status)
#' }
#'
#' @importFrom stats rexp rnorm runif
#' @export


simulate_data <- function(G, Z = NULL, prop = 0.6, beta = NULL, theta = 0.5, cens = FALSE, pcens = 0.25, tau = 0.5) {
  N <- length(G)
  k <- length(unique(G))
  p <- if (!is.null(Z)) if (is.null(dim(Z))) 1 else ncol(Z) else 0

  if (p > 0 && is.null(beta)) stop("You must provide beta when Z is given.")
  if (p > 0 && length(beta) != p) stop("Length of beta must match the number of columns in Z.")

  # Frailty
  cluster_id <- as.numeric(factor(G))
  frailty <- rnorm(k, mean = 0, sd = sqrt(theta))[cluster_id]

  # Predictors
  linpred <- if (p > 0) as.vector(Z %*% beta) + frailty else frailty

  # Simulate uniform random variable
  u <- runif(N)

  # Inverse CDF of cause 1 failure time
  tf <- -log(1 - u) / (prop * exp(linpred))  # Exponential approximation

  # Censure
  if (cens) {
    frailty_cens <- rnorm(k, mean = 0, sd = sqrt(tau))[cluster_id]

    # Calibrate censoring distribution
    # Assume base rate such that proportion censored ~ pcens
    base_rate <- uniroot(function(r) {
      mean(rexp(N, rate = r * exp(frailty_cens)) > tf) - pcens
    }, interval = c(0.001, 100))$root

    censure_times <- rexp(N, rate = base_rate * exp(frailty_cens))
  } else {
    censure_times <- rep(max(tf) + 1, N)
  }

  # Observed data
  X <- pmin(tf, censure_times)
  status <- as.integer(tf <= censure_times)

  # Build output
  df <- data.frame(
    times = X,
    status = status,
    clusters = G
  )
  if (p > 0) df <- cbind(df, as.data.frame(Z))

  return(df)
}

