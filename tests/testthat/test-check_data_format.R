test_that("check_data_format works", {

  expect_error(check_data_format(1),"The argument must be a data frame.")
})

test_that("the data_frame cannot be empty.",{
  m = matrix(0, nrow = 0, ncol = 4)
  expect_error(check_data_format(data.frame(m)),"The dataframe must contain at least one row, otherwise, these is nothing to analyse.")
})

test_that("The data.frame require at least 3 columns.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df_1 = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)

  df_2 = data.frame(times = df_1$times,status = df_1$status)
  expect_error(check_data_format(df_2),"The data frame must contain at least 3 columns:\\n - 'times': the observed time,\\n - 'status': 0 for right-censoring, i in \\[1,n\\] otherwise, where i is the cause of failure,\\n - 'clusters': cluster/group indicator.")
})

test_that("The first column has to be named 'times'.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df_1 = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)

  df_3 = data.frame(status = df_1$status,clusters = df_1$clusters,as.matrix(df_1[,4:length(df_1[1,])]))
  expect_error(check_data_format(df_3),"Column 1 must be named 'times'.")
})

test_that("The second column has to be named 'status'.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df_1 = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)

  df_4 = data.frame(times = df_1$times,clusters = df_1$clusters,as.matrix(df_1[,4:length(df_1[1,])]))
  expect_error(check_data_format(df_4),"Column 2 must be named 'status'.")
})

test_that("The third column has to be named 'clusters'.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df_1 = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)

  df_5 = data.frame(times = df_1$times,status = df_1$status,as.matrix(df_1[,4:length(df_1[1,])]))
  expect_error(check_data_format(df_5),"Column 3 must be named 'clusters'.")
})

test_that("'times' is a vector of positive numeric values.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df_1 = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)

  df_6 = df_1
  df_6$times[1] <- NA
  expect_error(check_data_format(df_6), "Some values in 'times' are missing \\(NA\\) or non-positive numeric.")
})

test_that("'status' is a vector of positive integer values.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df_1 = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)

  df_7 = df_1
  df_7$status[1] <- NA
  expect_error(check_data_format(df_7), "Some values in 'status' are missing \\(NA\\) or non-positive integer.")
})

test_that("'clusters' is a vector of positive integer values.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df_1 = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)

  df_8 = df_1
  df_8$clusters[1] <- NA
  expect_error(check_data_format(df_8), "Some values in 'clusters' are missing \\(NA\\) or non-positive integer.")
})

test_that("Covariables are numeric.",{
  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df_1 = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)

  df_9 = df_1
  df_9[1,length(df_9)] <- NA
  expect_error(check_data_format(df_9), "All covariates \\(columns 4 to end\\) must be numeric and contain no missing values \\(NA\\).")
})
