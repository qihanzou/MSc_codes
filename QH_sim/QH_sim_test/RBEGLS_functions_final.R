
RBEGLS_loop_final = function(iter, btor, ttor, betaest, thetaest, theta_ep.ini, CovU1, y, W, Dall, cov_ep, lo.bounde,up.bounde){
  
  iter0 = 0
  thetaest = theta_ep.ini
  b_diff = btor + 0.1
  t_diff = ttor + 0.1
  while ((b_diff > btor | t_diff > ttor)&iter0 < iter){
    CovUB1 = CovU1*betaest[2]^2
    CovUB = CovUB1
    out = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, thetaest, nug, "LB", lo.bounde,up.bounde)
    CorE = cor.mat(Dall[idx,idx], out$eta, cov_ep, nug = 0)
    CovE = out$theta[2]*CorE
    CovUBE = CovUB + CovE
    beta_update = solve(t(W)%*%solve(CovUBE)%*%W)%*%t(W)%*%solve(CovUBE)%*%y
    theta_update = out$theta
    
    b_diff = max(abs(beta_update - betaest))
    t_diff = max(abs(theta_update - thetaest))
    
    betaest = beta_update
    thetaest = theta_update
    
    iter0 = iter0 + 1
    
  }
  
  rlist = list(beta_est = betaest, theta_est = thetaest, iter = iter0, b_diff = b_diff, t_diff = t_diff, nll = out$nll)
  invisible(rlist)
}


Comp_W = function(alpha_i, theta_i, r1, r1s, x1s, cov_xi, Dall){
  cov_xiall = theta_i[2]*cor.mat(Dall, theta_i[1], cov_xi, nug = 0)
  w = r1%*%t(alpha_i) + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%t(alpha_i))
}

loop_funs_final = function(w, xi_j, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xi, cov_inv){
  
  W_j = cbind(1, w)
  fols_j = lm(y~W_j-1)
  beta_ols_j = fols_j$coef
  cov_xiall_j = xi_j[2]*cor.mat(Dall, xi_j[1], cov_xi, nug = 0)
  
  pc_j = r1 - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% r1s
  CovU_j = cov_xiall_j[idx, idx] - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% t(cov_xiall_j[idx, -idx]) + (pc_j)%*%solve(t(r1s)%*%solve(cov_xiall_j[-idx,-idx])%*%r1s)%*%t(pc_j)
  CovUB_j = CovU_j*beta_ols_j[2]^2 
  
  out_gls = MLE.fit_GLS(y, W_j, Dall[idx,idx], CovUB_j, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CovE_GLS_j = out_gls$theta[2]*cor.mat(Dall[idx,idx], out_gls$eta, cov_ep, nug = 0)
  CovUBE_GLS_j = CovUB_j + CovE_GLS_j
  
  
  varbeta_rw_v2m = solve(t(W_j)%*%cov_inv%*%W_j)                       
  B2_e = t(W_j)%*%cov_inv%*%CovUBE_GLS_j%*%cov_inv%*%W_j  
  varbeta_rw_v3 = diag(varbeta_rw_v2m%*%B2_e%*%varbeta_rw_v2m)               
  
  
  rlist = list(varbeta_rw_v3 = varbeta_rw_v3, beta_ols_j = beta_ols_j)
  invisible(rlist)
}


Expand_theta_final = function(theta0_m, theta1_m, dtheta0, dtheta1){
  
  theta0_up = theta0_m + 1*dtheta0
  theta0_down = theta0_m - 1*dtheta0
  if (theta0_down <= 0){
    theta0_down = 0.01
  }
  
  theta1_up = theta1_m + 1*dtheta1
  theta1_down = theta1_m - 1*dtheta1
  if (theta1_down <= 0){
    theta1_down = 0.01
  }
  
  set_theta0 = cbind(theta0_up, theta0_down)
  set_theta1 = cbind(theta1_up, theta1_down)
  
  set_data = expand.grid(set_theta0, set_theta1)
}