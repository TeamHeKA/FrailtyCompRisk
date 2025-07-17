test_that("Reml_Cox_frailty works", {
  suppressMessages({
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.6,cens = TRUE)
  res = Reml_Cox_frailty(df)
  expect_equal(length(res$beta),n_cov)
  expect_equal(length(res$theta),1)
  expect_equal(length(res$u),n_cluster)
  expect_lte(abs(res$beta[1]-1),1)
  expect_lte(abs(res$beta[2]-1.2),1)
  expect_lte(abs(res$theta - 0.3),1)

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0,cens = TRUE)
  res = Reml_Cox_frailty(df)
  expect_equal(length(res$beta),n_cov)
  expect_equal(length(res$theta),1)
  expect_equal(length(res$u),n_cluster)
  expect_lte(abs(res$beta[1]-1),1)
  expect_lte(abs(res$beta[2]-1.2),1)
  expect_lte(abs(res$theta),1)

  n_cov = 4
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2,0,-0.2),theta = 0.6,cens = TRUE)
  res = Reml_Cox_frailty(df)
  expect_equal(length(res$beta),n_cov)
  expect_equal(length(res$theta),1)
  expect_equal(length(res$u),n_cluster)
  expect_lte(abs(res$beta[1]-1),1)
  expect_lte(abs(res$beta[2]-1.2),1)
  expect_lte(abs(res$beta[3]),1)
  expect_lte(abs(res$beta[4]+0.2),1)
  expect_lte(abs(res$theta - 0.3),1)

  n_cov = 0
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(),theta = 0.6,cens = TRUE)
  res = Reml_Cox_frailty(df)
  expect_equal(length(res$beta),n_cov)
  expect_equal(length(res$theta),1)
  expect_equal(length(res$u),n_cluster)
  expect_lte(abs(res$theta - 0.6),1)
  })
})
