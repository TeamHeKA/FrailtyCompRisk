test_that("Reml_CompRisk_frailty works without covariables", {
  n_cov = 0
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(),theta = 0.3,cens = TRUE)
  res = Reml_CompRisk_frailty(df)
  expect_equal(length(res$beta),n_cov)
  expect_equal(length(res$theta),1)
  expect_equal(length(res$u),n_cluster)
  expect_lte(abs(res$theta - 0.3),0.5)
})
