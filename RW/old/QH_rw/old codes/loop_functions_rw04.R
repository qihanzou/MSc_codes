loop_funs = function(w, alpha, theta, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xiall){
  alpha = alpha
  theta = theta
  
  W = cbind(1, w)
  fols = lm(y~W-1)
  beta_ols = fols$coef
  varbeta_ols = diag(vcov(fols))
  
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  CovUB = CovU*beta_ols[2]^2 
  
  out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS = out2$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y
  
  rlist = list(theta = out2$theta, nll = out2$nll, beta = beta_RBEGLS, alpha = alpha, xi = theta)
  invisible(rlist)
}

loop_funs2 = function(w, alpha, theta, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xi){
  alpha = alpha
  theta = theta
  
  cov_xiall = theta[2]*cor.mat(Dall, theta[1], cov_xi, nug = 0)
  
  W = cbind(1, w)
  fols = lm(y~W-1)
  beta_ols = fols$coef
  varbeta_ols = diag(vcov(fols))
  
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  CovUB = CovU*beta_ols[2]^2 
  
  out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS = out2$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y
  
  rlist = list(theta = out2$theta, nll = out2$nll, beta = beta_RBEGLS, alpha = alpha, xi = theta)
  invisible(rlist)
}

loop_funs3 = function(w, alpha, theta, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xi){
  alpha = alpha
  theta = theta
  
  W = cbind(1, w)
  fols = lm(y~W-1)
  beta_ols = fols$coef
  varbeta_ols = diag(vcov(fols))
  
  cov_xiall = theta[2]*cor.mat(Dall, theta[1], cov_xi, nug = 0)
  
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  CovUB = CovU*beta_ols[2]^2 
  
  out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS = out2$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y
  
  CovUBE_GLS_rw_inv = solve(CovUBE_GLS)
  A1_rw = solve(t(W)%*%CovUBE_GLS_rw_inv%*%W)
  varbeta_rw = diag(A1_rw)
  
  rlist = list(theta = out2$theta, nll = out2$nll, beta = beta_RBEGLS, alpha = alpha, xi = theta, varbeta_rw = varbeta_rw)
  invisible(rlist)
}



loop_funs4 = function(w, alpha, theta, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xi, cov_inv){
  alpha = alpha
  theta = theta
  
  W = cbind(1, w)
  fols = lm(y~W-1)
  beta_ols = fols$coef
  varbeta_ols = diag(vcov(fols))
  
  cov_xiall = theta[2]*cor.mat(Dall, theta[1], cov_xi, nug = 0)
  
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  CovUB = CovU*beta_ols[2]^2 
  
  out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS = out2$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y
  
  CovUBE_GLS_rw_inv = solve(CovUBE_GLS)
  A1_rw = solve(t(W)%*%CovUBE_GLS_rw_inv%*%W)
  varbeta_rw = diag(A1_rw)
  
  
  varbeta_rw_v2m = solve(t(W)%*%cov_inv%*%W)                       
  B2_e = t(W)%*%cov_inv%*%CovUBE_GLS%*%cov_inv%*%W     

  varbeta_rw_v3m = varbeta_rw_v2m%*%B2_e%*%varbeta_rw_v2m                    
  
  varbeta_rw_v2 = diag(varbeta_rw_v2m)
  varbeta_rw_v3 = diag(varbeta_rw_v3m)
  
  
  rlist = list(theta = out2$theta, nll = out2$nll, beta = beta_RBEGLS, alpha = alpha, xi = theta, varbeta_rw = varbeta_rw, varbeta_rw_v2 = varbeta_rw_v2, varbeta_rw_v3 = varbeta_rw_v3)
  invisible(rlist)
}



Comp_W = function(alpha_i, theta_i, r1, r1s, x1s, cov_xi, Dall){
  cov_xiall = theta_i[2]*cor.mat(Dall, theta_i[1], cov_xi, nug = 0)
  w = r1%*%t(alpha_i) + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%t(alpha_i))
}

Expand_Var = function(alpha0_est_i, alpha1_est_i, theta0_est_i, theta1_est_i, dalpha0, dalpha1, dtheta0, dtheta1){
  
  alpha0_unc = alpha0_est_i
  alpha0_up = alpha0_est_i + 1*dalpha0
  alpha0_down = alpha0_est_i - 1*dalpha0
  
  alpha1_unc = alpha1_est_i
  alpha1_up = alpha1_est_i + 1*dalpha1
  alpha1_down = alpha1_est_i - 1*dalpha1
  
  theta0_unc = theta0_est_i
  theta0_up = theta0_est_i + 1*dtheta0
  theta0_down = theta0_est_i - 1*dtheta0
  if (theta0_down <= 0){
    theta0_down = 0.01
  }
  
  theta1_unc = theta1_est_i
  theta1_up = theta1_est_i + 1*dtheta1
  theta1_down = theta1_est_i - 1*dtheta1
  if (theta1_down <= 0){
    theta1_down = 0.01
  }
  
  Var_alpha0 = cbind(alpha0_unc, alpha0_up, alpha0_down)
  Var_alpha1 = cbind(alpha1_unc, alpha1_up, alpha1_down)
  Var_theta0 = cbind(theta0_unc, theta0_up, theta0_down)
  Var_theta1 = cbind(theta1_unc, theta1_up, theta1_down)
  
  Var_data = expand.grid(Var_alpha0, Var_alpha1, Var_theta0, Var_theta1)
}

