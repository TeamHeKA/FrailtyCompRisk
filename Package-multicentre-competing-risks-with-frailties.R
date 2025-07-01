library(Matrix)

check_data_format <- function(df)
{
  if (!is.data.frame(df)) {
    stop("The argument must be a data frame.")
}
  
  n_cols <- ncol(df)
  if (n_cols < 3) {
    stop(
      "The data frame must contain at least 3 columns:\n - 'times': the observed time,\n - 'status': 0 for right-censoring, i in [1,n] otherwise, where i is the cause of failure,\n - 'clusters': cluster/group indicator."
    )
  }
  
  col_names <- colnames(df)
  expected_names <- c("times", "status", "clusters")
  for (i in 1:3) {
    if (col_names[i] != expected_names[i]) {
      stop(paste0("Column ", i, " must be named '", expected_names[i], "'."))
    }
  }
  
  if (anyNA(df$times))
    stop("Some values in 'times' are missing (NA).")
  if (anyNA(df$status))
    stop("Some values in 'status' are missing (NA).")
  if (anyNA(df$clusters))
    stop("Some values in 'clusters' are missing (NA).")
  if (n_cols > 3) {
    covariate_cols <- df[, 4:n_cols, drop = FALSE]
    if (!all(sapply(covariate_cols, is.numeric))) {
      stop("All covariates (columns 4 to end) must be numeric.")
    }
  }
  
  message("The data frame is well structured.")}

#check_data_format() est validé


simulate_data <- function(G, Z, prop = 0.6, beta, theta, cens = FALSE, pcens = 0.25, tau = 0) 
{
  
  k <- length(unique(G))      # nombre de groupes
  N <- length(G)              # nombre d'observations
  p <- ncol(Z)                # nombre de covariables
  
  if (is.null(p)){p <- 0}
  if ((length(beta)!=0) && (length(beta) != p)) {
    stop("La longueur de beta doit correspondre au nombre de colonnes de Z.")
  }
  
  # Initialisation : cause de défaillance primaire (1)
  cf <- rep(1, N)
  ux <- runif(N)
  
  # Simulation des effets aléatoires (fragilité partagée)
  G.seed <- rnorm(k, mean = 0, sd = sqrt(theta))
  fragile <- factor(G)
  levels(fragile) <- G.seed
  fragile <- as.numeric(as.character(fragile))  # vecteur de taille N
  
  # Inversion de la fonction de répartition de la subdistribution pour la cause 1
  suppressWarnings({
    if (p>0){
      linpred <- as.vector(Z %*% beta) + fragile   # prédicteur linéaire + fragilité
    }else{
      linpred <- fragile
    }
    tf <- log(prop) - log(prop - 1 + exp(log(1 - ux) / exp(linpred)))
  })
  
  # Cas où tf est NaN (log de valeur négative par exemple)
  cf[is.na(tf)] <- 2
  l.na <- sum(is.na(tf))
  
  # Durées de type II (pour les NA dans tf)
  if (p > 0) {
    taux <- exp(Z %*% beta)      # ici on suppose que même beta est utilisé
    taux <- as.vector(taux)      # vecteur de taille N
    tf[is.na(tf)] <- rexp(l.na, rate = taux[is.na(tf)])
  }else{
    tf[is.na(tf)] <- rexp(l.na,1)
  }
  
  # Censure
  if (cens) {
    G.seed2 <- rnorm(k, mean = 0, sd = sqrt(tau))
    fragile2 <- factor(G)
    levels(fragile2) <- G.seed2
    fragile2 <- as.numeric(as.character(fragile2))
    
    # Calcul de lambda pour atteindre pcens
    if (p == 0){
      ex <- 1
    }else{
      ex <- exp(sum(beta)) / p # approximation grossière du E[exp(Xβ)] si inconnu
    }
    
    lamt <- 1
    b <- 0.5 * (1 - ex) - pcens * (1 + ex) + ex
    racstatus <- lamt * sqrt(b^2 + 4 * (1 - pcens) * pcens * ex)
    lamc <- (-b * lamt + racstatus) / (2 * (1 - pcens))
    
    censure1 <- rexp(N, rate = lamc * exp(fragile2))
  } else {
    censure1 <- rep(max(tf) + 1e-6, N)
  }
  
  # Construction des variables de temps observé
  X <- pmin(tf, censure1)
  epsilon <- cf
  epsilon[tf > X] <- 0
  status <- ifelse(epsilon == 1, 1, 0)
  Y <- censure1 * (status == 0) + X * (status == 1)
  
  # Construction du data.frame de sortie
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
#simulate_data() est validé

KaplanMeier_Censoring_vectorized <- function(times, status, unique_times)
{
  n <- length(unique_times)
  surv <- numeric(n)
  surv[1] <- 1
  
  for (i in 1:n) {
    t <- unique_times[i]
    d <- sum(times == t & status == 0)  # nb censures à t
    R <- sum(times >= t)                 # nb à risque à t
    if (R > 0) {
      surv[i] <- ifelse(i == 1, 1, surv[i-1] * (1 - d/R))
    } else {
      surv[i] <- ifelse(i == 1, 1, surv[i-1])
    }
  }
  
  return(surv)
}
#KaplanMeier_Censoring_vectorized() est validé

Nelson_Censoring_vectorized <- function(times, unique_times, status, clusters, u_bis)
{
  n_times <- length(unique_times)
  N <- length(times)
  
  # Étape 1 : Calcul de H(t) à chaque temps unique
  d_vec <- sapply(unique_times, function(t) sum(times == t & status == 0))
  R_vec <- sapply(unique_times, function(t) {
    ind <- which(times >= t)
    sum(exp(u_bis[clusters[ind]]))
  })
  
  lambda_hat <- ifelse(R_vec > 0, d_vec / R_vec, 0)
  H_vec <- cumsum(lambda_hat)  # Fonction de Nelson-Aalen estimée
  
  # Étape 2 : Associer à chaque times[i] la bonne valeur H(times[i])
  # Pour cela, on crée un vecteur indexé par times
  H_ti <- H_vec[findInterval(times, unique_times)]
  
  # Étape 3 : Calcul de la survie de la censure pour chaque individu
  S_C <- exp(-H_ti * exp(u_bis[clusters]))
  
  return(S_C)
}
#Nelson_Censoring_vectorized() est validé

create_G_function <- function(unique_times, G_vals)
{
  function(t) {
    idx <- findInterval(t, unique_times)
    idx[idx == 0] <- 1  # Si t < min(unique_times), surv = 1
    idx[idx > length(G_vals)] <- length(G_vals) # Si t > max, surv = 0 ou dernière valeur
    return(G_vals[idx])
  }
}
#create_G_function() est validé

compute_M_optimized <- function(times, status,f,u_bis,clusters) #Calcul optimisé de la matrice des poids
{
  N <- length(times)
  unique_times <- sort(unique(times))
  
  # Calcul vectorisé de G(t) pour les temps uniques
  if (f == "KaplanMeier_Censoring_vectorized"){
  G_vals <- KaplanMeier_Censoring_vectorized(times, status, unique_times)}
  if (f == "Nelson_Censoring_vectorized"){
  G_vals <- Nelson_Censoring_vectorized(times,unique_times,status,clusters,u_bis)}
  G_func <- create_G_function(unique_times, G_vals)
  
  # Matrices répétées pour comparaison vectorisée
  T_mat_i <- matrix(rep(times, each = N), nrow = N)
  T_mat_j <- matrix(rep(times, times = N), nrow = N)
  fail_mat_j <- matrix(rep(status, times = N), nrow = N)
  
  # Évaluation vectorisée des conditions
  cond1 <- (T_mat_j <= T_mat_i) & (fail_mat_j > 1)
  cond2 <- (T_mat_j >= T_mat_i)
  cond3 <- (T_mat_j < T_mat_i) & (fail_mat_j <= 1)
  
  # Vecteur G(t) évalué en chaque T_i et T_j
  G_T_i <- G_func(times)
  G_T_j <- G_func(times)
  G_T_i_mat <- matrix(rep(G_T_i, each = N), nrow = N)
  G_T_j_mat <- matrix(rep(G_T_j, times = N), nrow = N)
  
  # Construction de M
  M <- matrix(0, nrow = N, ncol = N)
  M[cond1] <- G_T_i_mat[cond1] / G_T_j_mat[cond1]
  M[cond2] <- 1
  M[cond3] <- 0
  
  return(M)
}
#compute_M_optimized() est validé

logLikelihood_1 = function(status,M,W)
{
  loglik <- 0
  for (i in which(status == 1)) #formule de l1
  {
    v <- log(sum(M[i,] * W@x))
    loglik = loglik + log(W[i,i]) - v
  }
  return(loglik)
}
#logLikelihood_1() est validé
  
logLikelihood_2 = function(u,theta,K) 
{
  penalty <- -(1/2)*(K*log(2*pi*theta) + (1/theta)*(t(u)%*%u))
  return(penalty)
}
#logLikelihood_2() est validé

Ml_Cox <- function(data, max_iter=100, tol=1e-6)
{
  times <- data$times
  status <- data$status
  X <- as.matrix(data[,4:(length(data[1,]))])
  
  N <- length(times)
  p <- ncol(X)
  beta <- rep(0, p)
  
  ord <- order(times)
  ord_times <- times[ord]
  ord_by_time_status <- status[ord]
  X_bis <- X[ord, , drop=FALSE]
  
  for (iter in 1:max_iter) {
    eta <- as.vector(X_bis %*% beta)
    exp_eta <- exp(eta)
    
    loglik <- 0
    U <- rep(0, p)
    I <- matrix(0, p, p)
    
    for (i in which(ord_by_time_status == 1)) {
      R <- which(ord_times >= ord_times[i])
      
      sum_exp_eta <- sum(exp_eta[R])
      sum_exp_eta_X <- colSums(exp_eta[R] * X_bis[R, , drop=FALSE])
      sum_exp_eta_XX <- t(X_bis[R, , drop=FALSE]) %*% (exp_eta[R] * X_bis[R, , drop=FALSE])
      
      loglik <- loglik + (eta[i] - log(sum_exp_eta))
      
      U <- U + (X_bis[i, ] - sum_exp_eta_X / sum_exp_eta)
      
      I <- I + (sum_exp_eta_XX / sum_exp_eta) - (sum_exp_eta_X %*% t(sum_exp_eta_X)) / (sum_exp_eta^2)
    }
    
    status <- solve(I, U)
    beta_new <- beta + status
    
    if (sqrt(sum(status^2)) < tol) {
      break
    }
    beta <- beta_new
  }
  return(beta = beta)
}
#Reml() est validé

Reml_Cox_frailty <- function(data,max_iter=300, tol = 1e-6) 
{
  times <- data$times
  status <- data$status
  clusters <- data$clusters
  X <- as.matrix(data[,4:length(data[1,])])
  
  p <- ncol(X)
  N <- length(times)
  K <- length(unique(clusters))
  
  if (is.null(p)){p <- 0}
  if (p>0){
    data_ini <- data.frame(times = times,status = status,clusters = clusters, X)
    gamma_Cox <- Ml_Cox(data_ini)
    gamma_0 <- gamma_Cox
  }else{
    gamma_0 <- c()
  }
  
  u_0 <- rep(0, K)
  theta_0 <- 1
  theta_1 <- 2
  
  Z <- X
  Q <- Matrix(0, nrow = N, ncol = K, sparse = TRUE)
  for (i in 1:N) {
    Q[i, clusters[i]] <- 1
  }
  
  if (p>0){
    Y <- cbind(Z, Q)
    Y_t <- rbind(t(Z), t(Q))
  }else{
    Y <-Q
    Y_t <- t(Q)
  }
  
  D <- Matrix(as.numeric(status == 1), ncol = 1, sparse = TRUE)
  
  if (p>0){
    eta <- as.vector(Z %*% gamma_0 + Q %*% u_0)
    W_diag <- exp(eta)
  }else{
    eta <- as.vector(Q %*% u_0)
    W_diag <- exp(eta)
  }
  W <- Diagonal(x = W_diag)

  M <- matrix(0,N,N)
  for (i in 1:N){
    for (j in 1:N){
      if (times[i] <= times[j])
      M[j,i] <- 1
    }
  }
  
  Mt <- t(M)
  
  v <- as.vector(colSums(M * W_diag))
  a <- numeric(N)
  a[status == 1] <- 1 / v[status == 1]
  A <- Diagonal(x = a)
  
  b_vec <- as.vector(M %*% a)
  B <- Diagonal(x = b_vec)
  
  lik_0 <- logLikelihood_1(status, M, W) + logLikelihood_2(u_0,theta_0,K)
  lik_1 <- lik_0 + 1
  iter <- 1
  
  not_converged <- ((abs(lik_1 - lik_0) >= tol) && (abs(theta_1 - theta_0) >= tol))
  
  while ((not_converged) && (iter <= max_iter)) {
    grad <- as.vector(D) - W %*% (M %*% a)
    
    a_sq <- a^2
    A2 <- Diagonal(x = a_sq)
    negHess <- W %*% (B - M %*% A2 %*% Mt %*% W)
    
    M_pen <- Diagonal(p + K)
    diag(M_pen)[(p + 1):(p + K)] <- 1 / theta_0
    
    V <- Y_t %*% negHess %*% Y + M_pen
    inv_V <- solve(V)
    
    update_direction <- inv_V %*% (Y_t %*% grad) - inv_V %*% c(rep(0, p), (1 / theta_0) * u_0)
    
    new_par <- c(gamma_0, u_0) + update_direction
    
    if (p>0){
      gamma_0 <- new_par[1:p]
    }
    u_0 <- new_par[(p + 1):(p + K)]
    
    theta_1 <- theta_0
    trace_term <- sum(diag(inv_V[(p + 1):(p + K), (p + 1):(p + K)]))
    theta_0 <- as.numeric((t(u_0) %*% u_0) / (K - (trace_term / theta_0)))
    
    if (p>0){
      eta <- as.vector(Z %*% gamma_0 + Q %*% u_0)
      W_diag <- exp(eta)
    }else{
      eta <- as.vector(Q %*% u_0)
      W_diag <- exp(eta)
    }
    
    W <- Diagonal(x = W_diag)
    
    v <- as.vector(colSums(M * W_diag))
    a <- numeric(N)
    a[status == 1] <- 1 / v[status == 1]
    A <- Diagonal(x = a)
    
    b_vec <- as.vector(M %*% a)
    B <- Diagonal(x = b_vec)
    
    lik_0 <- lik_1
    lik_1 <- logLikelihood_1(status, M, W) + logLikelihood_2(u_0, theta_0,K)
    
    not_converged <- ((abs(lik_1 - lik_0) >= tol) && (abs(theta_1 - theta_0) >= tol))
    iter <- iter + 1
  }
  
  A22 <- inv_V[(p + 1):(p + K), (p + 1):(p + K)]
  wald_stat <- as.numeric(t(u_0) %*% solve(A22) %*% u_0)
  
  df_eff <- K - sum(diag(A22)) / theta_0
  
  p_value <- pchisq(wald_stat, df = df_eff, lower.tail = FALSE)
  if (p < 0.05){
    message("A p-value greater than 0.05 in the test of θ=0 suggests that the cluster effect may be negligible.\n p-value = ",p_value)
  }
  return(list(beta = gamma_0,u = u_0, theta = theta_0,p_value = p_value))
}
#Reml_Cox_frailty() est validé

Ml_CompRisk <- function(data, max_iter = 100, tol = 1e-6) 
{
  times <- data$times
  status <- data$status
  X <- as.matrix(data[,4:length(data[1,])])
  event_of_interest <- (status == 1)
  status_mod <- ifelse(event_of_interest, 1, 0)
  N <- length(times)
  p <- ncol(X)
  beta <- rep(0, p)
  ord <- order(times)
  times_ord <- times[ord]
  status_ord <- status_mod[ord]
  X_ord <- X[ord, , drop = FALSE]
  for (iter in 1:max_iter) {
    eta <- as.vector(X_ord %*% beta)
    exp_eta <- exp(eta)
    
    U <- rep(0, p)
    I <- matrix(0, p, p)
    
    for (i in which(status_ord == 1)) {
      R <- which(times_ord >= times_ord[i])
      
      sum_exp_eta <- sum(exp_eta[R])
      sum_exp_eta_X <- colSums(exp_eta[R] * X_ord[R, , drop = FALSE])
      sum_exp_eta_XX <- t(X_ord[R, , drop = FALSE]) %*% (exp_eta[R] * X_ord[R, , drop = FALSE])
      
      U <- U + (X_ord[i, ] - sum_exp_eta_X / sum_exp_eta)
      I <- I + (sum_exp_eta_XX / sum_exp_eta) - (sum_exp_eta_X %*% t(sum_exp_eta_X)) / (sum_exp_eta^2)
    }
    
    step <- solve(I, U)
    beta_new <- beta + step
    
    if (sqrt(sum(step^2)) < tol) {
      break
    }
    beta <- beta_new
  }
  
  return(beta = beta)
}
#Reml_CompRisk() est validé

Reml_CompRisk_frailty <- function(data,cluster_censoring = F, max_iter=300, tol = 1e-6, threshold = 1e-5) 
{
  p <- (ncol(data) - 3)
  times <- data$times
  status <- data$status
  N <- length(times)
  
  if (p>0){
  X <- as.matrix(data[, 4:(p + 3)])
  }else{
    X <- matrix(0,N,0)
  }
  
  clusters <- data$clusters
  
  K <- length(unique(clusters))
  
  if (p>0){
    data_ini <- data.frame(times = data$times,status = data$status,clusters = clusters,X)
    gamma_Cox <- Ml_Cox(data_ini)
    gamma_0 <- gamma_Cox
  }else{
    gamma_0 <- c()
  }
  u_0 <- rep(0, K)
  theta_0 <- 1
  theta_1 <- 2
  
  Z <- X
  Q <- Matrix(0, nrow = N, ncol = K, sparse = TRUE)
  for (i in 1:N) {
    Q[i, clusters[i]] <- 1
  }
  
  if (p>0){
    Y <- cbind(Z, Q)
    Y_t <- rbind(t(Z), t(Q))
  }else{
    Y <-Q
    Y_t <- t(Q)
  }
  
  D <- Matrix(as.numeric(status == 1), ncol = 1, sparse = TRUE)
  
  if (p>0){
    eta <- as.vector(Z %*% gamma_0 + Q %*% u_0)
    W_diag <- exp(eta)
  }else{
    eta <- as.vector(Q %*% u_0)
    W_diag <- exp(eta)
  }
  W <- Diagonal(x = W_diag)
  
  if (cluster_censoring && !(0%in%status)){cluster_censoring <- F}
  
  if (!cluster_censoring){
    M <- compute_M_optimized(times, status,"KaplanMeier_Censoring_vectorized",u_bis=rep(0,1),clusters)}
  else{
    d_cens = as.integer(status == 0)
    data_bis = data.frame(times = times,status = d_cens,clusters = clusters,matrix(0,N,0))
    l = Reml_Cox_frailty(data_bis)
    u_bis <- l$u
    M <- compute_M_optimized(times, status,"Nelson_Censoring_vectorized",u_bis,clusters)}
  
  Mt <- t(M)
  
  v <- as.vector(colSums(M * W_diag))
  a <- numeric(N)
  a[status == 1] <- 1 / v[status == 1]
  A <- Diagonal(x = a)
  
  b_vec <- as.vector(M %*% a)
  B <- Diagonal(x = b_vec)
  
  lik_0 <- logLikelihood_1(status, M, W) + logLikelihood_2(u_0,theta_0,K)
  lik_1 <- lik_0 + 1
  iter <- 1
  
  not_converged <- ((abs(lik_1 - lik_0) >= tol) && (abs(theta_1 - theta_0) >= tol))
  
  while ((not_converged) && (iter <= max_iter)) {
    grad <- as.vector(D) - W %*% (M %*% a)
    a_sq <- a^2
    A2 <- Diagonal(x = a_sq)
    negHess <- W %*% (B - M %*% A2 %*% Mt %*% W)
    
    M_pen <- Diagonal(p + K)
    diag(M_pen)[(p + 1):(p + K)] <- 1 / theta_0
    
    V <- Y_t %*% negHess %*% Y + M_pen
    inv_V <- solve(V)
    
    update_direction <- inv_V %*% (Y_t %*% grad) - inv_V %*% c(rep(0, p), (1 / theta_0) * u_0)

    new_par <- c(gamma_0, u_0) + update_direction
    
    if (p>0){
    gamma_0 <- new_par[1:p]
    }
    u_0 <- new_par[(p + 1):(p + K)]
    
    theta_1 <- theta_0
    trace_term <- sum(diag(inv_V[(p + 1):(p + K), (p + 1):(p + K)]))
    theta_0 <- as.numeric((t(u_0) %*% u_0) / (K - (trace_term / theta_0)))
    
    if (theta_0 < threshold)
    {
      theta_0 <- threshold
      gamma_0 <- Ml_CompRisk(data)
      u_0 <- rep(0,K)
      message("treshold \n")
      break
    }
    
    if (p>0){
      eta <- as.vector(Z %*% gamma_0 + Q %*% u_0)
      W_diag <- exp(eta)
    }else{
      eta <- as.vector(Q %*% u_0)
      W_diag <- exp(eta)
    }
    
    W <- Diagonal(x = W_diag)
    
    v <- as.vector(colSums(M * W_diag))
    a <- numeric(N)
    a[status == 1] <- 1 / v[status == 1]
    A <- Diagonal(x = a)
    
    b_vec <- as.vector(M %*% a)
    B <- Diagonal(x = b_vec)
    
    lik_0 <- lik_1
    lik_1 <- logLikelihood_1(status, M, W) + logLikelihood_2(u_0, theta_0,K)
    
    not_converged <- ((abs(lik_1 - lik_0) >= tol) && (abs(theta_1 - theta_0) >= tol))
    
    iter <- iter + 1
  }
  
  A22 <- inv_V[(p + 1):(p + K), (p + 1):(p + K)]
  wald_stat <- as.numeric(t(u_0) %*% solve(A22) %*% u_0)
  
  df_eff <- K - sum(diag(A22)) / theta_0
  
  if (!(threshold == theta_0)){
    p_value <- pchisq(wald_stat, df = df_eff, lower.tail = FALSE)
  }else{
    message("theta might be inferior to the threshold")
    p_value <- NA
  }
  if (!(is.na(p_value)) && (p_value >= 0.05))
  {
    message("A p-value greater than 0.05 in the test of θ=0 suggests that the cluster effect may be negligible.\n p-value = ",p_value)
  }
  return(list(beta = gamma_0,u = u_0, theta = theta_0,p_value = p_value))
}
#Reml_CompRisk_frailty() est validé

Parameters_estimation <- function(data,method = "CompRisk_frailty",cluster_censoring = F,max_iter = 300, tol = 1e-6)
{
  if (!(method %in% c("Comprisk_frailty","Comprisk","Cox_Frailty","Cox")))
  {
    stop("The method chosen must be : 'Comprisk_frailty','Comprisk','Cox_Frailty' or 'Cox'.")
  }
  check_data_format(data)
  if (cluster_censoring && (method != "CompRisk_frailty"))
    {stop("'cluster_censoring' is only for 'Comprisk_frailty' method.")}
  
  if (method == "CompRisk_frailty")
  {
    res = Reml_CompRisk_frailty(data,cluster_censoring)
  }
  if (method == "Comprisk")
  {
    res = Ml_CompRisk(data)
  }
  if (method == "Cox_frailty")
  {
    res = Reml_Cox_frailty(data)
  }
  if (method == "Cox")
  {
    res = Ml_Cox(data)
  }
  return(res)
}