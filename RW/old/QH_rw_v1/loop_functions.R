loop_funs = function(w, alpha, theta){
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


Comp_W = function(alpha_i, theta_i){
  cov_xiall = theta_i[2]*cor.mat(Dall, theta_i[2], cov_xi, nug = 0)
  w = r1%*%t(alpha_i) + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%t(alpha_i))
}



