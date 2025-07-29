test_that("Reml_Cox_frailty works in a typical case.", {
  set.seed(123)
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
  })
})

test_that("Reml_Cox_frailty works when the variance of the frailty effect is close to 0.",{
  set.seed(123)
  suppressMessages({

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
  })
})

test_that("Reml_Cox_frailty works without covariables.",{
  set.seed(123)
  suppressMessages({

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

test_that("The p_value for the test theta = 0 might be computationnaly unavailable and therefore equals 'NaN'.",{
  set.seed(1)

  n_cov = 5
  n_per_cluster = 20
  n_cluster = 15
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2,0,0,1),theta = 0,cens = TRUE)
  expect_message(Reml_Cox_frailty(df),"The estimated value of theta is too small to allow a reliable statistical test of whether theta = 0, hence the p-value is 'NaN'. \\nThis suggests that the cluster effect may be negligible.")
})

test_that("When the estimated variance for the frailty effect is too small, it suggests to use a Cox model.",{
  set.seed(2)

  n_cov = 5
  n_per_cluster = 20
  n_cluster = 15
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2,0,0,1),theta = 0,cens = TRUE)
  expect_error(Reml_Cox_frailty(df, threshold = 1e-5),"The estimated variance of the cluster effect is inferior to the threshold \\('threshold' =1e-05\\), we suggest you use the method without frailty:\\n - method = 'Cox' if you are using 'Parameters_estimation\\(\\)'\\n - 'Ml_Cox\\(\\)' if you are using 'Reml_Cox_frailty\\(\\)'")
})

test_that("When there is only one cluster, it suggests to use a Cox model.",{
  set.seed(123)

  n_cov = 2
  n_per_cluster = 20
  n_cluster = 1
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(0,-0.7),theta = 0,cens = TRUE)
  expect_error(Reml_Cox_frailty(df),"There is only one cluster, therefore there are no cluster effect, you might use: \\n - method = 'Cox' if you are using 'Parameters_estimation\\(\\)'\\n - 'Ml_Cox\\(\\)' if you are using 'Reml_Cox_frailty\\(\\)'")
})
