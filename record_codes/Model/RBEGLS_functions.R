
RBEGLS_loop = function(iter, btor, ttor, betaest, thetaest, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde){
  
  iter0 = 0
  thetaest = theta_ep.ini
  b_diff = btor + 0.1
  t_diff = ttor + 0.1
  while ((b_diff > btor | t_diff > ttor)&iter0 < iter){
    CovUB = CovU*betaest[2]^2
    out = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, thetaest, nug, "LB", lo.bounde,up.bounde)
    CorE = cor.mat(Dall[idx,idx], out$eta, cov_ep, nug = 0)
    CovE = out$theta[2]*CorE
    CovUBE = CovUB + CovE
    beta_update = solve(t(W)%*%solve(CovUBE)%*%W)%*%t(W)%*%solve(CovUBE)%*%y
    theta_update = out$theta
    
    b_diff = max(abs(beta_update - betaest))
    t_diff = max(abs(theta_update - thetaest))
    #all_diff = max(b_diff, t_diff)
    
    betaest = beta_update
    thetaest = theta_update
    
    iter0 = iter0 + 1
    
  }
  
  rlist = list(beta_est = betaest, theta_est = thetaest, iter = iter0, b_diff = b_diff, t_diff = t_diff)
  invisible(rlist)
}


RBEGLS_loop2 = function(iter, btor, ttor, betaest, thetaest, CovU, y, W, D, cov_ep, lo.bounde,up.bounde){
  
  iter0 = 0
  #thetaest = theta_ep.ini
  b_diff = btor + 0.1
  t_diff = ttor + 0.1
  while ((b_diff > btor | t_diff > ttor)&iter0 < iter){
    CovUB = CovU*betaest[2]^2
    out = MLE.fit_GLS(y, W, D, CovUB, cov_ep, thetaest, nug, "LB", lo.bounde,up.bounde)
    CorE = cor.mat(D, out$eta, cov_ep, nug = 0)
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
  
  rlist = list(beta_est = betaest, theta_est = thetaest, iter = iter0, b_diff = b_diff, t_diff = t_diff)
  invisible(rlist)
}


RBEGLS_loop3 = function(iter, btor, ttor, betaest, thetaest, CovU, y, W, D, cov_ep, lo.bounde,up.bounde){
  
  iter0 = 0
  b_diff = btor + 0.1
  t_diff = ttor + 0.1
  while ((b_diff > btor | t_diff > ttor)&iter0 < iter){
    CovUB = CovU*betaest[2]^2
    out = MLE_fit_GLS(y = y, X = W, D = D, kSigma = CovUB, cov_model = cov_ep)
    CovE = out[[1]][2]*exp(-D/out[[1]][1]) 
    CovUBE = CovUB + CovE
    beta_update = solve(t(W)%*%solve(CovUBE)%*%W)%*%t(W)%*%solve(CovUBE)%*%y
    theta_update = out[[1]]
    
    b_diff = max(abs(beta_update - betaest))
    t_diff = max(abs(theta_update - thetaest))
    
    betaest = beta_update
    thetaest = theta_update
    
    iter0 = iter0 + 1
    
  }
  
  rlist = list(beta_est = betaest, theta_est = thetaest, iter = iter0, b_diff = b_diff, t_diff = t_diff)
  invisible(rlist)
}
