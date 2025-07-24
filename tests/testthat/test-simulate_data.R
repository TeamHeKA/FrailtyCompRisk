test_that("simulated data with simulate_data works", {
  n_cov = 0
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(),theta = 0.3,cens = TRUE)
  expect_equal(check_data_format(df), T)

  n_cov = 0
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(),theta = 0.3,cens = FALSE)
  expect_equal(check_data_format(df), T)

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)
  expect_equal(check_data_format(df), T)

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = FALSE)
  expect_equal(check_data_format(df), T)

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  expect_error(simulate_data(G,Z,prop = 0.6,beta = NULL,theta = 0.3,cens = FALSE), "You must provide beta when Z is given.")

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  expect_error(simulate_data(G,Z,prop = 0.6,beta = c(1),theta = 0.3,cens = FALSE), "Length of beta must match the number of columns in Z.")
})
