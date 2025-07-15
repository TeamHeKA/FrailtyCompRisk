test_that("Ml_CompRisk works", {
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)
  res = Ml_CompRisk(df)
  expect_equal(length(res$beta),n_cov)
  expect_lte(abs(res$beta[1]-1),0.5)
  expect_lte(abs(res$beta[2]-1.2),0.5)
})
