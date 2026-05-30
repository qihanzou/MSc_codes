

fit_function <- function(y, W, Dall, cov_xiall, designx1, designx1s, beta_est, ff, out1){
  # -------------------------------------------------- #
  # 3. Proposed Method 1 and variance
  # -------------------------------------------------- #
  pc = designx1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% designx1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx])
  + (pc)%*%solve(t(designx1s)%*%solve(cov_xiall[-idx,-idx])%*%designx1s)%*%t(pc)
  CovUB = CovU*beta_est[2]^2 
  H = W%*%solve(t(W)%*%W)%*%t(W) 
  tmp_e = ff$resid/(1-diag(H)) 
  eet = as.matrix(tmp_e%*%t(tmp_e) - CovUB) 
  out2 = MLE.fit_tr(eet, Dall[idx,idx], cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE = out2$theta[2]*CorE
  CovUBE = CovUB + CovE
  beta_update = solve(t(W)%*%solve(CovUBE)%*%W)%*%t(W)%*%solve(CovUBE)%*%y
  ######################################################################################
  CovUB2 = CovU*beta_update[2]^2
  CovUBE2 = CovUB2 + CovE 
  tmp = t(W)%*%solve(CovUBE2)%*%W
  varboot2 = solve(tmp)
  varbeta2  = diag(varboot2)   # proposed 1
  
  alist = list(beta_update = beta_update,  
       thetae = out2$theta, 
       varbeta2 = varbeta2)
  return(alist)
}



