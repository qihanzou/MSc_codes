MLE.fit_GLS = function(y, X, D, kSigma, cov_model, theta.ini, nug, opt, lo_bound, up_bound){
  n = length(y)
  q = length(theta.ini) 

  if (opt == "LB"){
    prof.min = optim(theta.ini, profile.nll_GLS, method = "L-BFGS-B", lower = lo_bound, upper = up_bound, 
                     y = y, X = X, kSigma = kSigma, cov_model = cov_model, D = D, nug = nug)
  } else if (opt == "NM"){
    prof.min = optim(theta.ini, profile.nll_GLS, method = "Nelder-Mead",
                     y = y, X = X, kSigma = kSigma, cov_model = cov_model, D = D, nug = nug)
  } else if (opt == "CG"){
    prof.min = optim(theta.ini, profile.nll_GLS, method = "CG",
                     y = y, X = X, kSigma = kSigma, cov_model = cov_model, D = D, nug = nug)
  }

  # Return theta
  theta = prof.min$par
  suc = prof.min$convergence + 1
  nll = prof.min$value
  
  if (nug == TRUE){
    cratio = theta[q]/(theta[q-1]+theta[q])
    eta = c(theta[1:(q-2)], cratio)
  } else if (nug == FALSE){
    eta = theta[1:(q-2)]
  }
  
  rlist = list(theta = theta, suc, nll = nll, eta = eta)
  
  invisible(rlist)
}


 
profile.nll_GLS = function(theta, y, X, D, kSigma, cov_model, nug){
  n = length(y)
  p = dim(X)[2]
  q = length(theta) 
  
  if (nug == TRUE){
    sill = theta[q-1]+theta[q]
    cratio = theta[q]/sill
    eta = c(theta[1:(q-2)], cratio)
  } else if (nug == FALSE){
    eta = theta[1:(q-1)]
    sill = theta[q]
  }
  
  
  # Full covariance matrix:
  cormat_e = cor.mat(D, eta, cov_model, nug) 
  covmat = sill*cormat_e + kSigma
  
  # term 1
  L = t(chol(covmat))
  log_det_covmat = 2*sum(log(diag(L))) 
  cov_inv = solve(covmat)
  
  # term 2
  L1 = t(chol(t(X)%*%cov_inv%*%X))
  log_det_X_covmat_inv_X = 2*sum(log(diag(L1)))
  
  # term 3
  r = y - X%*%solve(t(X)%*%cov_inv%*%X)%*%t(X)%*%cov_inv%*%y
  part3 = t(r)%*%cov_inv%*%r
  
  # term 4
  log_det_XX = 2*sum(log(diag(t(chol(t(X)%*%X))))) 

  
  nll2 = (log_det_covmat + log_det_X_covmat_inv_X + part3 - log_det_XX + (n-p)*log(2*pi))
  return(nll2/2)
}






MLE.fit_IREML = function(y, X, D, beta1, kSigma, theta.ini, nug, r1, r1s, idx, opt, lo_bound, up_bound){
  n = length(y)
  q = length(theta.ini) 
  
  if (opt == "LB"){
    prof.min = optim(theta.ini, profile.nll_IREML, method = "L-BFGS-B", lower = lo_bound, upper = up_bound, 
                     y = y, X = X, beta1 = beta1, kSigma = kSigma, D = D, nug = nug, r1 = r1, r1s = r1s, idx = idx)
  } 
  
  # Return theta
  theta = prof.min$par
  suc = prof.min$convergence + 1
  nll = prof.min$value
  
  rlist = list(theta = theta, suc, nll = nll)
  
  invisible(rlist)
}



profile.nll_IREML = function(theta, beta1, y, X, D, kSigma, nug, r1, r1s, idx){
  n = length(y)
  p = dim(X)[2]
  q = length(theta) 
  
  phi = theta[1]
  sigma2 = theta[2]
  
  # Full covariance matrix:
  cov_xiall = sigma2*(1 + 3^.5*Dall/phi)*exp(-3^0.5*Dall/phi)
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  
  
  covmat = beta1^2*CovU + kSigma
  
  # term 1
  L = t(chol(covmat))
  log_det_covmat = 2*sum(log(diag(L))) 
  cov_inv = solve(covmat)
  
  # term 2
  L1 = t(chol(t(X)%*%cov_inv%*%X))
  log_det_X_covmat_inv_X = 2*sum(log(diag(L1)))
  
  # term 3
  r = y - X%*%solve(t(X)%*%cov_inv%*%X)%*%t(X)%*%cov_inv%*%y
  part3 = t(r)%*%cov_inv%*%r
  
  # term 4
  log_det_XX = 2*sum(log(diag(t(chol(t(X)%*%X))))) 
  
  
  nll2 = (log_det_covmat + log_det_X_covmat_inv_X + part3 - log_det_XX + (n-p)*log(2*pi))
  return(nll2/2)
}



LOUP = function(D, paraest, cov_model, nug){
  sdout = SDCal(D, paraest, cov_model, nug)
  mu_xi <- paraest
  Sigma_xi <- sdout$covmat
  eig <- eigen(Sigma_xi)
  eig_values <- eig$values
  eig_vectors <- eig$vectors
  conf_level <- 0.95
  c2 <- qchisq(conf_level, df = 2)
  c <- sqrt(c2)
  vertex1 <- mu_xi + c * sqrt(eig_values[1]) * eig_vectors[,1]
  vertex2 <- mu_xi - c * sqrt(eig_values[1]) * eig_vectors[,1]
  upb = c(vertex1[1], vertex1[2])
  lob = c(vertex2[1], vertex2[2])
  
  if (lob[1]<=0){
    lob[1] = 0.001
  } 
  if (lob[2]<=0){
    lob[2] = 0.001
  }
  list(lob = lob, upb = upb)
}


IREML_loop_final2 = function(iter, btor, ttor, ator, xtor, 
                            betaest, thetaest, alphaest, xiest, CovE_GLS, 
                            y, W, r1, r1s, x1s, nug, Dall, idx,
                            cov_ep, cov_xi){
  
  iter0 = 1
  W_set = list()
  beta_set = NULL
  theta_set  = NULL
  alpha_set = NULL
  xi_set = NULL
  
  W_set[[iter0]] = W
  beta_set = cbind(beta_set, betaest)
  theta_set = cbind(theta_set, thetaest)
  alpha_set = cbind(alpha_set, alphaest)
  xi_set = cbind(xi_set, xiest)
  
  b_diff = btor + 0.1
  t_diff = ttor + 0.1
  a_diff = ator + 0.1
  x_diff = xtor + 0.1
  
  
  while ((b_diff > btor | t_diff > ttor | a_diff > ator | x_diff > xtor)&iter0 < iter){
    
    LOPU_xi = LOUP(Dall[-idx,-idx], xiest, cov_xi, nug)
    LOPU_theta = LOUP(Dall[idx,idx], thetaest, cov_ep, nug)
    
    out1_j = MLE.fit_IREML(y, W, Dall, betaest[2], CovE_GLS, xiest, nug, r1, r1s, idx, "LB", LOPU_xi$lob, LOPU_xi$upb)
    cov_xiall_j = out1_j$theta[2]*cor.mat(Dall, out1_j$theta[1], cov_xi, nug = 0)
    xi_j = out1_j$theta
    alpha_j =  solve(t(r1s)%*%solve(cov_xiall_j[-idx,-idx])%*%r1s)%*%t(r1s)%*%solve(cov_xiall_j[-idx,-idx])%*%x1s
    w_j = r1%*%alpha_j + cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx])%*%(x1s - r1s%*%alpha_j) 
    W_j = cbind(1, w_j)
    pc_j = r1 - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% r1s
    CovU_j = cov_xiall_j[idx, idx] - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% t(cov_xiall_j[idx, -idx]) + (pc_j)%*%solve(t(r1s)%*%solve(cov_xiall_j[-idx,-idx])%*%r1s)%*%t(pc_j)
    CovUB_j = CovU_j*betaest[2]^2 
    
    out2_j = MLE.fit_GLS(y, W_j, Dall[idx,idx], CovUB_j, cov_ep, thetaest, nug, "LB", LOPU_theta$lob, LOPU_theta$upb)
    CovE_GLS_j = out2_j$theta[2]*cor.mat(Dall[idx,idx], out2_j$eta, cov_ep, nug = 0)
    CovUBE_GLS_j = CovUB_j + CovE_GLS_j
    beta_RBEGLS_j = solve(t(W_j)%*%solve(CovUBE_GLS_j)%*%W_j)%*%t(W_j)%*%solve(CovUBE_GLS_j)%*%y
    theta_RBEGLS_j = out2_j$theta
    
    b_diff = max(abs(beta_RBEGLS_j - betaest))
    t_diff = max(abs(theta_RBEGLS_j - thetaest))
    a_diff = max(abs(alpha_j - alphaest))
    x_diff = max(abs(xi_j - xiest))
    
    betaest = beta_RBEGLS_j
    thetaest = theta_RBEGLS_j
    alphaest = alpha_j
    xiest = xi_j
    iter0 = iter0 + 1
    
    beta_set = cbind(beta_set, betaest)
    theta_set = cbind(theta_set, thetaest)
    alpha_set = cbind(alpha_set, alphaest)
    xi_set = cbind(xi_set, xiest)
    W_set[[iter0]] = W_j
    
  }
  
  rlist = list(beta_set = beta_set, theta_set = theta_set, CovUBE_IREML = CovUBE_GLS_j, 
               alpha_set = alpha_set, xi_set = xi_set, W_set = W_set,
               iter = iter0, b_diff = b_diff, t_diff = t_diff, a_diff = a_diff, x_diff = x_diff)
  invisible(rlist)
}


IREML_loop_final = function(iter, btor, ttor, ator, xtor, 
                             betaest, thetaest, alphaest, xiest, CovE_GLS, 
                             y, W, r1, r1s, x1s, nug, Dall, idx,
                             cov_ep, cov_xi){
  
  iter0 = 1
  W_set = list()
  beta_set = NULL
  theta_set  = NULL
  alpha_set = NULL
  xi_set = NULL
  
  W_set[[iter0]] = W
  beta_set = cbind(beta_set, betaest)
  theta_set = cbind(theta_set, thetaest)
  alpha_set = cbind(alpha_set, alphaest)
  xi_set = cbind(xi_set, xiest)
  
  b_diff = btor + 0.1
  t_diff = ttor + 0.1
  a_diff = ator + 0.1
  x_diff = xtor + 0.1
  
  LOPU_xi = LOUP(Dall[-idx,-idx], xiest, cov_xi, nug)
  LOPU_theta = LOUP(Dall[idx,idx], thetaest, cov_ep, nug)
  
  
  while ((b_diff > btor | t_diff > ttor | a_diff > ator | x_diff > xtor)&iter0 < iter){
    
    
    out1_j = MLE.fit_IREML(y, W, Dall, betaest[2], CovE_GLS, xiest, nug, r1, r1s, idx, "LB", LOPU_xi$lob, LOPU_xi$upb)
    cov_xiall_j = out1_j$theta[2]*cor.mat(Dall, out1_j$theta[1], cov_xi, nug = 0)
    xi_j = out1_j$theta
    alpha_j =  solve(t(r1s)%*%solve(cov_xiall_j[-idx,-idx])%*%r1s)%*%t(r1s)%*%solve(cov_xiall_j[-idx,-idx])%*%x1s
    w_j = r1%*%alpha_j + cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx])%*%(x1s - r1s%*%alpha_j) 
    W_j = cbind(1, w_j)
    pc_j = r1 - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% r1s
    CovU_j = cov_xiall_j[idx, idx] - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% t(cov_xiall_j[idx, -idx]) + (pc_j)%*%solve(t(r1s)%*%solve(cov_xiall_j[-idx,-idx])%*%r1s)%*%t(pc_j)
    CovUB_j = CovU_j*betaest[2]^2 
    
    out2_j = MLE.fit_GLS(y, W_j, Dall[idx,idx], CovUB_j, cov_ep, thetaest, nug, "LB", LOPU_theta$lob, LOPU_theta$upb)
    CovE_GLS_j = out2_j$theta[2]*cor.mat(Dall[idx,idx], out2_j$eta, cov_ep, nug = 0)
    CovUBE_GLS_j = CovUB_j + CovE_GLS_j
    beta_RBEGLS_j = solve(t(W_j)%*%solve(CovUBE_GLS_j)%*%W_j)%*%t(W_j)%*%solve(CovUBE_GLS_j)%*%y
    theta_RBEGLS_j = out2_j$theta
    
    b_diff = max(abs(beta_RBEGLS_j - betaest))
    t_diff = max(abs(theta_RBEGLS_j - thetaest))
    a_diff = max(abs(alpha_j - alphaest))
    x_diff = max(abs(xi_j - xiest))
    
    betaest = beta_RBEGLS_j
    thetaest = theta_RBEGLS_j
    alphaest = alpha_j
    xiest = xi_j
    iter0 = iter0 + 1
    
    beta_set = cbind(beta_set, betaest)
    theta_set = cbind(theta_set, thetaest)
    alpha_set = cbind(alpha_set, alphaest)
    xi_set = cbind(xi_set, xiest)
    W_set[[iter0]] = W_j
    
  }
  
  rlist = list(beta_set = beta_set, theta_set = theta_set, CovUBE_IREML = CovUBE_GLS_j, 
               alpha_set = alpha_set, xi_set = xi_set, W_set = W_set,
               iter = iter0, b_diff = b_diff, t_diff = t_diff, a_diff = a_diff, x_diff = x_diff)
  invisible(rlist)
}






IREML = function(iter, btor, ttor, ator, xtor, 
                            betaest, thetaest, alphaest, xiest, CovE_GLS, 
                            y, W, r1, r1s, x1s, nug, Dall, idx,
                            cov_ep, cov_xi){
  
  iter0 = 1
  W_set = list()
  beta_set = NULL
  theta_set  = NULL
  alpha_set = NULL
  xi_set = NULL
  
  W_set[[iter0]] = W
  beta_set = cbind(beta_set, betaest)
  theta_set = cbind(theta_set, thetaest)
  alpha_set = cbind(alpha_set, alphaest)
  xi_set = cbind(xi_set, xiest)
  
  b_diff = btor + 0.1
  t_diff = ttor + 0.1
  a_diff = ator + 0.1
  x_diff = xtor + 0.1
  
  LOPU_xi = LOUP(Dall[-idx,-idx], xiest, cov_xi, nug)
  LOPU_theta = LOUP(Dall[idx,idx], thetaest, cov_ep, nug)
  
  
  while ((b_diff > btor | t_diff > ttor | a_diff > ator | x_diff > xtor)&iter0 < iter){
    
    
    out1_j = MLE.fit_IREML(y, W, Dall, betaest[2], CovE_GLS, xiest, nug, r1, r1s, idx, "LB", LOPU_xi$lob, LOPU_xi$upb)
    cov_xiall_j = out1_j$theta[2]*cor.mat(Dall, out1_j$theta[1], cov_xi, nug = 0)
    xi_j = out1_j$theta
    
    alpha_j =  alpha_est
    W_j = W
    pc_j = r1 - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% r1s
    CovU_j = cov_xiall_j[idx, idx] - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% t(cov_xiall_j[idx, -idx]) + (pc_j)%*%solve(t(r1s)%*%solve(cov_xiall_j[-idx,-idx])%*%r1s)%*%t(pc_j)
    CovUB_j = CovU_j*betaest[2]^2 
    
    out2_j = MLE.fit_GLS(y, W_j, Dall[idx,idx], CovUB_j, cov_ep, thetaest, nug, "LB", LOPU_theta$lob, LOPU_theta$upb)
    CovE_GLS_j = out2_j$theta[2]*cor.mat(Dall[idx,idx], out2_j$eta, cov_ep, nug = 0)
    CovUBE_GLS_j = CovUB_j + CovE_GLS_j
    beta_RBEGLS_j = solve(t(W_j)%*%solve(CovUBE_GLS_j)%*%W_j)%*%t(W_j)%*%solve(CovUBE_GLS_j)%*%y
    theta_RBEGLS_j = out2_j$theta
    
    b_diff = max(abs(beta_RBEGLS_j - betaest))
    t_diff = max(abs(theta_RBEGLS_j - thetaest))
    a_diff = max(abs(alpha_j - alphaest))
    x_diff = max(abs(xi_j - xiest))
    
    betaest = beta_RBEGLS_j
    thetaest = theta_RBEGLS_j
    alphaest = alpha_j
    xiest = xi_j
    iter0 = iter0 + 1
    
    beta_set = cbind(beta_set, betaest)
    theta_set = cbind(theta_set, thetaest)
    alpha_set = cbind(alpha_set, alphaest)
    xi_set = cbind(xi_set, xiest)
    W_set[[iter0]] = W_j
    
  }
  
  rlist = list(beta_set = beta_set, theta_set = theta_set, CovUBE_IREML = CovUBE_GLS_j, 
               alpha_set = alpha_set, xi_set = xi_set, W_set = W_set,
               iter = iter0, b_diff = b_diff, t_diff = t_diff, a_diff = a_diff, x_diff = x_diff)
  invisible(rlist)
}