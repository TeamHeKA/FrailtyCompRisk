test_that("Parameters_estimation needs a data.frame", {
  expect_error(Parameters_estimation(1,"Clearly a wrong method"),"The method chosen must be : 'CompRisk_frailty','CompRisk','Cox_frailty' or 'Cox'.")
})

test_that("Parameters_estimation needs a well structured data.frame.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)
  df_bis = data.frame(status = df$status,clusters = df$clusters,as.matrix(df[,4:length(df[1,])]))
  expect_error(Parameters_estimation(df_bis),"")
})

test_that("cluster_censoring is specific to the 'CompRisk_frailty' method.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)
  expect_error(Parameters_estimation(df,"CompRisk",cluster_censoring = TRUE),"'cluster_censoring' is only for 'CompRisk_frailty' method.")
})

test_that("The method 'CompRisk_frailty' works.",{
  set.seed(123)
  suppressMessages({

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)
  res = Parameters_estimation(df,"CompRisk_frailty")
  expect_equal(length(res$beta),n_cov)
  expect_equal(length(res$theta),1)
  expect_equal(length(res$u),n_cluster)
  expect_lte(abs(res$beta[1]-1),1)
  expect_lte(abs(res$beta[2]-1.2),1)
  expect_lte(abs(res$theta - 0.3),1)
  })
})

test_that("The method 'CompRisk' works.",{
  set.seed(123)
  suppressMessages({

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)
  res = Parameters_estimation(df,"CompRisk")
  expect_equal(length(res$beta), n_cov)
  expect_lte(abs(res$beta[1] - 1), 1)
  expect_lte(abs(res$beta[2] - 1.2), 1)
  })
})

test_that("The method 'Cox_frailty' works.",{
  set.seed(123)
  suppressMessages({

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n * n_cov, 0, 1), ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1, 1.2),theta = 0.6,cens = TRUE)
  res = Parameters_estimation(df, "Cox_frailty")
  expect_equal(length(res$beta), n_cov)
  expect_equal(length(res$theta), 1)
  expect_equal(length(res$u), n_cluster)
  expect_lte(abs(res$beta[1] - 1), 1)
  expect_lte(abs(res$beta[2] - 1.2), 1)
  expect_lte(abs(res$theta - 0.3), 1)
  })
})

test_that("The method 'Cox' works.",{
  set.seed(123)
  suppressMessages({

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)
  res = Parameters_estimation(df,"Cox")
  expect_equal(length(res$beta),n_cov)
  expect_lte(abs(res$beta[1]-1),1)
  expect_lte(abs(res$beta[2]-1.2),1)
  })
})
