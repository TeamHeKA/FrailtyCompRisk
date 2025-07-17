test_that("check_data_format works", {
  expect_error(check_data_format(1),"The argument must be a data frame.")

  n_cov = 2
  n_per_cluster = 15
  n_cluster = 20
  n = n_cluster * n_per_cluster
  G = rep(1:n_cluster, each = n_per_cluster)
  Z = matrix(rnorm(n*n_cov,0,1),ncol = n_cov)
  df_1 = simulate_data(G,Z,prop = 0.6,beta = c(1,1.2),theta = 0.3,cens = TRUE)

  df_2 = data.frame(times = df_1$times,status = df_1$status,as.matrix(df_1[,4:length(df_1[1,])]))
  expect_error(check_data_format(df_2),"Column 3 must be named 'clusters'.")

  df_3 = df_1
  df_3$times[1] <- NA
  expect_error(check_data_format(df_3), "Some values in 'times' are missing \\(NA\\).")

  df_4 = df_1
  df_4[1,length(df_4)] <- NA
  expect_error(check_data_format(df_4), "All covariates \\(columns 4 to end\\) must be numeric and contain no missing values \\(NA\\).")
})
