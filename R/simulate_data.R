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


simulate_data <- function(G, Z, prop = 0.6, beta, theta, cens = FALSE, pcens = 0.25, tau = 0)
{

  k <- length(unique(G))
  N <- length(G)
  p <- ncol(Z)

  if (is.null(p)){p <- 0}
  if ((length(beta)!=0) && (length(beta) != p)) {
    stop("The length of beta must match the number of columns in Z.")
  }

  # Initialization: primary failure cause (1)
  cf <- rep(1, N)
  ux <- runif(N)

  # Simulation of random effects (shared frailty)
  G.seed <- rnorm(k, mean = 0, sd = sqrt(theta))
  fragile <- factor(G)
  levels(fragile) <- G.seed
  fragile <- as.numeric(as.character(fragile))  # vector of size N

  # Inversion of the cumulative incidence function for cause 1
  suppressWarnings({
    if (p>0){
      linpred <- as.vector(Z %*% beta) + fragile   # linear predictor + frailty
    }else{
      linpred <- fragile
    }
    tf <- log(prop) - log(prop - 1 + exp(log(1 - ux) / exp(linpred)))
  })

  # Case where tf is NaN (e.g. log of negative value)
  cf[is.na(tf)] <- 2
  l.na <- sum(is.na(tf))

  # Type II durations (for NAs in tf)
  if (p > 0) {
    taux <- exp(Z %*% beta)      # here we assume the same beta is used
    taux <- as.vector(taux)      # vector of size N
    tf[is.na(tf)] <- rexp(l.na, rate = taux[is.na(tf)])
  }else{
    tf[is.na(tf)] <- rexp(l.na,1)
  }

  # Censoring
  if (cens) {
    G.seed2 <- rnorm(k, mean = 0, sd = sqrt(tau))
    fragile2 <- factor(G)
    levels(fragile2) <- G.seed2
    fragile2 <- as.numeric(as.character(fragile2))

    # Lambda computation to reach pcens
    if (p == 0){
      ex <- 1
    }else{
      ex <- exp(sum(beta)) / p # rough approximation of E[exp(Xβ)] if unknown
    }

    lamt <- 1
    b <- 0.5 * (1 - ex) - pcens * (1 + ex) + ex
    racstatus <- lamt * sqrt(b^2 + 4 * (1 - pcens) * pcens * ex)
    lamc <- (-b * lamt + racstatus) / (2 * (1 - pcens))

    censure1 <- rexp(N, rate = lamc * exp(fragile2))
  } else {
    censure1 <- rep(max(tf) + 1e-6, N)
  }

  # Construction of observed time variables
  X <- pmin(tf, censure1)
  epsilon <- cf
  epsilon[tf > X] <- 0
  status <- ifelse(epsilon == 1, 1, 0)
  Y <- censure1 * (status == 0) + X * (status == 1)

  # Output data.frame construction
  data <- data.frame(tf = tf,
                     censure1 = censure1,
                     X = X,
                     Y = Y,
                     status = status,
                     epsilon = epsilon,
                     cf = cf,
                     G = G,
                     Z = Z)

  N <- length(data[,1])
  n <- length(data[1,])
  times <- rep(0,N)
  status <- rep(0,N)
  clusters <- data$G

  for (i in 1:N)
  {
    if (data$status[i] == 0)
    {
      status[i] = 0
      times[i] = data$censure1[i]
    }
    if (data$status[i] != 0)
    {
      status[i] = data$cf[i]
      times[i] = data$tf[i]
    }
  }
  if (p >= 1){
    Cov = as.matrix(data[,(n+1-p):n])
  }else{
    Cov = matrix(0,N,0)
  }

  df = data.frame(times = times,status = status,clusters = clusters,Cov)
  return(df)
}
