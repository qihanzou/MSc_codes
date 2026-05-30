reml_loglik <- function(theta, y, X, dist_mat, cov_model = "exponential", kSigma, nug = FALSE) {
  n <- length(y)
  p <- ncol(X)
  
  # Covariance parameters
  range <- theta[1]
  sill <- theta[2]
  nugget <- if (nug) theta[3] else 0  # Set to 0 if nug == FALSE
  
  # Build covariance matrix
  if (cov_model == "exponential") {
    Sigma <- sill * exp(-dist_mat / range)
    
  } else if (cov_model == "gaussian") {
    Sigma <- sill * exp(-(dist_mat / range)^2)
    
  } else {
    stop("Unsupported covariance model")
  }
  
  # Add nugget effect if enabled
  diag(Sigma) <- diag(Sigma) + nugget
  
  # Total covariance matrix: spatial + known measurement error component
  covmat <- Sigma + kSigma
  
  # Cholesky decomposition and terms
  L <- t(chol(covmat))
  log_det_covmat <- 2 * sum(log(diag(L)))
  cov_inv <- solve(covmat)
  
  # Term 2: log determinant of Xt Σ-1 X
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

reml_loglik_modifed <- function(theta, xi, beta, Dall, idx, y, X, dist_mat, cov_model = "exponential", kSigma, nug = TRUE) {
  n <- length(y)
  p <- ncol(X)
  
  # Covariance parameters
  range <- theta[1]
  sill <- theta[2]
  nugget <- if (nug) theta[3] else 0  # Set to 0 if nug == FALSE
  
  xi0 = xi[1]
  xi1 = xi[2]
  
  # Build covariance matrix
  if (cov_model == "exponential") {
    Sigma <- sill * exp(-dist_mat / range)
  } else {
    stop("Unsupported covariance model")
  }
  
  # Add nugget effect if enabled
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
  
  # Term 2: log determinant of Xt Σ-1 X
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

# Objective to maximize
neg_loglik <- function(theta) reml_loglik(theta, y, W, Dall[idx,idx], cov_model = "exponential", CovUB, nug = FALSE)

# Optimization
fit <- optim(
  par = c(3.5, 1.2),           # initial values: range, sill, nugget
  fn = neg_loglik,
  method = "L-BFGS-B",
  lower = c(1e-3, 1e-3),    # avoid singularities
  upper = c(20, 20)
)

# Estimated parameters
fit$par