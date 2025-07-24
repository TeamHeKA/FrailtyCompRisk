test_that("Ml_Cox works", {
  set.seed(123)

  n_cov = 0
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(),theta = 0.3,cens = TRUE)
  expect_error(Ml_Cox(df),"With no covariables, there is nothing to estimate")

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)
  res = Ml_Cox(df)
  expect_equal(length(res$beta),n_cov)
  expect_lte(abs(res$beta[1]-1),1)
  expect_lte(abs(res$beta[2]-1.2),1)
})
