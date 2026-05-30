sandwich_cov = function(Sigma, y, X, beta, Diagonal = TRUE){
  e = y - X%*%beta
  Sigma_inv = solve(Sigma)
  if (Diagonal == TRUE){
    middle = diag(as.vector(e)^2)
  } else {
    middle = e%*%t(e)
  }
  part1 = solve(t(X)%*%Sigma_inv%*%X)
  part2 = t(X)%*%Sigma_inv%*%(middle)%*%Sigma_inv%*%X
  VarB = part1%*%part2%*%part1
  list(VarB = VarB)
}

mat3Dplot = function(matrix){
  n_row <- nrow(matrix)
  n_col <- ncol(matrix)
  x <- 1:n_row
  y <- 1:n_col
  z <- matrix
  grid <- expand.grid(x = x, y = y)
  z_values <- as.vector(z)
  library(plotly)
  plot_ly(x = ~grid$x, y = ~grid$y, z = ~z_values, type = "mesh3d")
}

estimate_params_cov <- function(D, Sigma_full, Sigma_UB, x0, 
                                constrained = TRUE, constraint_type = "leq") {
  sum_M <- sum(Sigma_full) - sum(Sigma_UB)
  
  objective_fn <- function(params) {
    phi <- params[1]
    sigma2 <- params[2]
    C_est <- sigma2 * exp(-D / phi)
    (sum(C_est) - sum_M)^2
  }
  
  # General constraint function
  constraint_fn <- function(params) {
    phi <- params[1]
    sigma2 <- params[2]
    C_est <- sigma2 * exp(-D / phi)
    val <- sum(C_est) - sum_M
    if (constraint_type == "leq") {
      return(val)           # Enforce sum(C_est) - sum_M ≤ 0
    } else if (constraint_type == "geq") {
      return(-val)          # Enforce sum(C_est) - sum_M ≥ 0 → -val ≤ 0
    } else {
      stop("Invalid constraint_type. Use 'leq' or 'geq'.")
    }
  }
  
  if (constrained) {
    result <- nloptr(
      x0 = x0,
      eval_f = objective_fn,
      eval_g_ineq = constraint_fn,
      lb = c(1e-6, 1e-6),
      opts = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6)
    )
  } else {
    result <- nloptr(
      x0 = x0,
      eval_f = objective_fn,
      lb = c(1e-6, 1e-6),
      opts = list("algorithm" = "NLOPT_LN_BOBYQA", "xtol_rel" = 1e-6)
    )
  }
  
  return(result$solution)
}


Q_mat <- function(Sigma, y, X, beta, Diagonal = TRUE) {
  e <- y - X %*% beta
  
  # Define E matrix (empirical covariance of residuals)
  if (Diagonal) {
    E <- diag(as.vector(e)^2)
  } else {
    E <- e %*% t(e)
  }
  
  V_inv <- solve(Sigma)
  XtVinv <- t(X) %*% V_inv
  A_inv <- solve(XtVinv %*% X)
  B <- XtVinv %*% E %*% V_inv %*% X
  
  Q <- X %*% A_inv %*% B %*% A_inv %*% t(X)
  
  return(list(Q = Q))
}

Simulation_res = function(Sigma, y, X, beta, B){
  n <- length(y)
  p <- ncol(X)
  Sigma_inv <- solve(Sigma)
  residual_list <- list()
  for (b in 1:B) {
    y_sim <- as.vector(X %*% beta + MASS::mvrnorm(1, mu = rep(0, n), Sigma = Sigma))
    beta_sim <- solve(t(X) %*% Sigma_inv %*% X) %*% t(X) %*% Sigma_inv %*% y_sim
    e_sim <- y_sim - X %*% beta_sim
    
    residual_list[[b]] <- e_sim
  }
  V_hat <- matrix(0, n, n)
  for (e_sim in residual_list) {
    V_hat <- V_hat + e_sim %*% t(e_sim)
  }
  V_hat <- V_hat/B
  Xt_Sinv <- t(X) %*% Sigma_inv
  middle <- Xt_Sinv %*% V_hat %*% Sigma_inv %*% X
  Var_beta_hat_sandwich <- solve(Xt_Sinv %*% X) %*% middle %*% solve(Xt_Sinv %*% X)
  
  
  list(VarB = Var_beta_hat_sandwich)
}

Find_expand_var = function(Sigma1, Sigma2, X){
  Cov = solve(t(X) %*% solve(Sigma1) %*% X) %*%
    t(X) %*% solve(Sigma1) %*%
    Sigma2 %*%
    solve(Sigma1) %*% X %*%
    solve(t(X) %*% solve(Sigma1) %*% X)
  list(cov = Cov)
}

# Delta method first order
DT_Delta_1st_jacobianfun = function(w, Sigma_w, g){
  J <- jacobian(g, w)
  mean_Z <- g(w)
  cov_Z <- J %*% Sigma_w %*% t(J)
  list(z = mean_Z, cov_z = cov_Z)
}

# Delta method first order
DT_Delta_1st = function(w, Sigma_w, g, dg){
  diag_vals = dg(w)
  J = diag(c(diag_vals))
  mean_Z <- g(w)
  cov_Z <- J %*% Sigma_w %*% t(J)
  list(z = mean_Z, cov_z = cov_Z)
}


# MCT, Monte Carlo Transformation
DT_MC = function(w, Sigma_w, g, iterations){
  Z_set = NULL
  W_set = NULL
  for (i in 1:iterations){
    w_i = t(w) + rmvnorm(1, sigma = Sigma_w)
    z_i = g(w_i)
    Z_set = rbind(Z_set, z_i)
    W_set = rbind(W_set, w_i)
  }
  Cov_Z = matrix(0, length(idx), length(idx))
  for (j in 1:iterations){
    Cov_Z = Cov_Z + (Z_set[j,] - colMeans(Z_set))%*%t(Z_set[j,] - colMeans(Z_set))
  }
  CovU_mc = 1/(iterations-1)*Cov_Z
  z_mc = colMeans(Z_set)
  list(cov_z = CovU_mc, z = z_mc)
}

# Find Bures distance between matrix A and matrix B:
bures_distance <- function(A, B) {
  A_half <- sqrtm(A)
  inner <- A_half %*% B %*% A_half
  inner_sqrt <- sqrtm(inner)
  trace_A <- sum(diag(A))
  trace_B <- sum(diag(B))
  trace_middle <- sum(diag(inner_sqrt))
  val <- trace_A + trace_B - 2 * trace_middle
  # For Numerical stable
  val <- Re(val)
  val <- max(val, 0)
  return(sqrt(val))
}

# Find the transformation matrix based on Bures distance:
find_T_matrix <- function(Sigma_x, Sigma_z) {
  S_z_sqrt <- sqrtm(Sigma_z)
  middle <- S_z_sqrt %*% Sigma_x %*% S_z_sqrt
  middle_sqrt_inv <- solve(sqrtm(middle))
  Tm <- S_z_sqrt %*% middle_sqrt_inv %*% S_z_sqrt
  return(Tm)
}

# Use to check the positive definite of the matrix:
is_positive_definite <- function(A) {
  if (!isSymmetric(A)) return(FALSE)
  result <- tryCatch({
    chol(A)
    TRUE
  }, error = function(e) FALSE)
  return(result)
}




# The REML for RBEGLS method modified version that can estimate the nugget effect now:
reml_loglik_modifed <- function(theta, alpha_est, Dall, idx, y, W, cov_model = "exponential", nug = TRUE) {
  n <- length(y)
  p <- ncol(W)
  
  # Covariance parameters
  range <- theta[1]
  sill <- theta[2]
  nugget <- if (nug) theta[3] else 0  
  
  xi0 = theta[3]
  xi1 = theta[4]
  
  # covariance 
  if (cov_model == "exponential") {
    Sigma <- sill * exp(-Dall[idx,idx] / range)
  } else {
    stop("This is a Unsupported covariance model")
  }
  
  # Add nugget effect if needed
  diag(Sigma) <- diag(Sigma) + nugget
  
  # kSigma matrix
  cov_xiall = xi1*(1 + 3^.5*Dall/xi0)*exp(-3^0.5*Dall/xi0)
  
  w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 
  W = cbind(1, w)
  X = W
  fols = lm(y~W-1)
  beta_ols = fols$coef
  
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  kSigma = CovU*beta_ols^2 
  
  
  # Total covariance matrix: spatial + known measurement error component
  covmat <- Sigma + kSigma
  
  # Cholesky decomposition and terms
  L <- t(chol(covmat))
  log_det_covmat <- 2 * sum(log(diag(L)))
  cov_inv <- solve(covmat)
  
  # Term 2: log determinant 
  XtSiX <- t(X) %*% cov_inv %*% X
  L1 <- t(chol(XtSiX))
  log_det_X_covmat_inv_X <- 2 * sum(log(diag(L1)))
  
  # Term 3: projection part
  beta_hat <- solve(XtSiX, t(X) %*% cov_inv %*% y)
  r <- y - X %*% beta_hat
  part3 <- t(r) %*% cov_inv %*% r
  
  # REML negative log-likelihood
  nll <- (log_det_covmat + log_det_X_covmat_inv_X + part3 + (n - p) * log(2 * pi)) / 2
  return(as.numeric(nll))
}


# # Objective to maximize
# neg_loglik <- function(theta) reml_loglik_modified(theta, alpha_est, Dall, idx, y, W, cov_model = "exponential", nug = FALSE)
# 
# # Optimization
# fit <- optim(
#   par = c(3.5, 1.2, 4, 4),           # initial values: range, sill, nugget
#   fn = neg_loglik,
#   method = "L-BFGS-B",
#   lower = c(1e-3, 1e-3, 1e-3, 1e-3),    # avoid singularities
#   upper = c(20, 20, 20, 20)
# )
# 
# # Estimated parameters
# fit$par

check_Q_properties <- function(Q) {
  # 1. Rank
  matrix_rank <- Matrix::rankMatrix(Q)[1]  # rankMatrix returns a matrix, extract the value
  
  # 2. Symmetry
  is_symmetric <- all.equal(Q, t(Q), tolerance = 1e-8)
  
  # 3. Positive semi-definiteness
  eigenvalues <- eigen(Q, symmetric = TRUE, only.values = TRUE)$values
  is_psd <- all(eigenvalues >= -1e-8)  # allow small numerical error
  
  # Print results
  cat("Rank of Q:", matrix_rank, "\n")
  cat("Is symmetric?:", isTRUE(is_symmetric), "\n")
  cat("Is positive semi-definite?:", is_psd, "\n")
  
  invisible(list(rank = matrix_rank, symmetric = is_symmetric, psd = is_psd, eigenvalues = eigenvalues))
}