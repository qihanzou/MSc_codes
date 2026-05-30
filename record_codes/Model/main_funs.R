

main_cal = function(Local_gas, Local_data){
  plot(Local_gas$X, Local_gas$Y, cex = 0.3, pch = 16)
  points(Local_data$X, Local_data$Y, cex = 0.3, col = "red", pch = 16)

  D_gas = deprojected_distances(Local_gas$RA, Local_gas$DEC, Local_gas$RA, Local_gas$DEC)
  D_data = deprojected_distances(Local_data$RA, Local_data$DEC, Local_data$RA, Local_data$DEC)
  D_data_to_gas = deprojected_distances(Local_data$RA, Local_data$DEC, Local_gas$RA, Local_gas$DEC)
  x1s = as.matrix(select(Local_gas, c("CO21_mgd")))
  r1s = as.matrix(cbind(1, select(Local_gas, c("proj_dist"))))
  
  
  N = dim(D_data)[1]
  M = dim(D_gas)[1]
  idx = 1:N
  Dall = matrix(0, N + M, N + M)
  Dall[1:N, 1:N] = D_data
  Dall[1:N, (1:M) + N] = D_data_to_gas
  Dall[(1:M) + N, 1:N] = t(D_data_to_gas)
  Dall[(1:M) + N, (1:M) + N] = D_gas
  
  
  xi.ini = c(0.5)
  lo.bound = c(1e-5)
  up.bound = c(4)
  nug = 0
  cov_xi = "Exp"
  out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
  
  
  xi_est = out1$theta     #. res1[[1]] 
  alpha_est  = out1$beta     #res1[[4]] 
  Sigma1 = xi_est[2]*exp(-D_gas/xi_est[1]) 
  diag(Sigma1) = xi_est[2] 
  r1 = cbind(1, Local_data$proj_dist)
  cmat1 = xi_est[2]*exp(-D_data_to_gas/xi_est[1]) 
  pl1 = as.matrix(r1)%*%alpha_est 
  pe1 = cmat1%*% solve(Sigma1)%*%(x1s - r1s%*%alpha_est)
  w = pl1 + pe1
  
  
  # ---------------------------------------------------------------------------- #
  # Training set and testing set
  # ---------------------------------------------------------------------------- #
  set.seed(2024)
  idx = sample.int(n = nrow(new_data), size = floor(0.2*nrow(new_data)), replace = F)
  test= new_data[idx, ]
  train = new_data[-idx, ] 
  Wtrain = as.matrix(cbind(1, select(train, c("w", "proj_dist"))))
  Ytrain = as.matrix(select(train, c("Z_N2S2Ha")))
  
  
  
  
  
  
  
  
  
  inv_Sigma = solve(Sigma1)
  Cz0 = xi_est[2]*exp(-D_data/xi_est[1]) 
  pc = r1 - cmat1%*%inv_Sigma%*%r1s
  CovU = Cz0 - cmat1%*%inv_Sigma%*%t(cmat1) + (pc)%*%solve(t(r1s)%*%inv_Sigma%*%r1s)%*%t(pc)
  
  
  new_data = Local_data
  new_data$w = w
  train = new_data 
  W = Wtrain = as.matrix(cbind(1, select(train, c("w", "proj_dist"))))
  Ytrain = y = as.matrix(select(train, c("Z_N2S2Ha")))
  
  beta_ols = solve(t(Wtrain)%*%Wtrain)%*%t(Wtrain)%*%Ytrain
  
  
  theta.ini = c(0.1)
  lo.bound = c(1e-5)
  up.bound = c(2)
  cov_ep = "Exp"
  kr = MLE.fit(y, W, Dall[idx, idx], cov_ep, theta.ini, nug = 0, "LB", lo.bound, up.bound)
  
  # ---------------------------------------------------------------------------- #
  # KR RMSE and MAD
  # ---------------------------------------------------------------------------- #
  thetahat_kr = kr[[1]] 
  betahat_kr  = kr[[4]] 
  Sigma_kr = thetahat_kr[2]*exp(-D_data[-idx,-idx]/thetahat2[1]) 
  d1 = cbind(1, test$w, test$proj_dist)
  cmat2 = thetahat2[2]*exp(-D_data[idx,-idx]/thetahat2[1]) 
  pl2 = as.matrix(d1)%*%betahat2       
  pe2 = cmat2%*% solve(Sigma2)%*%(Ytrain - Wtrain%*%betahat2) 
  pred_Y = pl2 + pe2
  (RMSE_geo_dist_kr <- sqrt(mean((pred_Y - test$Z_N2S2Ha)^2)))
  (MAD_geo_dist_kr <- mean(abs(pred_Y - test$Z_N2S2Ha)))
  
  
  
  
  
  kr2 = MLE.fit(y, W[,-2], Dall[idx, idx], cov_ep, theta.ini, nug = 0, "LB", lo.bound, up.bound)
  
  
  
  theta_ep.ini = c(0.5, 0.01)
  lo.bounde = c(1e-5, 1e-10)
  up.bounde = c(2,2)
  nug = 0
  RBE_res = RBEGLS_loop(50, 0.001, 0.001, beta_ols, theta_ep.ini, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
  
  beta_RBEGLS = RBE_res$beta_est
  theta_RBEGLS  = RBE_res$theta_est
  CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
  CovE_GLS = theta_RBEGLS[2]*CorE_GLS
  CovUB = CovU*beta_RBEGLS[2]^2 
  CovUBE_GLS = CovUB + CovE_GLS
  CovUBE_GLS_inv = solve(CovUBE_GLS)
  
  
  sdout = SDCal(Dall[-idx,-idx], out1$theta, 'Exp', nug = 0)
  
  boot_iter = 100
  ETAallboot_all = SimuData_D(out1$theta, Dall, "Exp", nug, boot_iter)
  EVar = EVarKR = EVarKRT = matrix(0,3, 3)
  r1all = rbind(r1, r1s)
  x2 = W[,3]
  
  
  for (v in 1:boot_iter){
    
    
    etaallboot = ETAallboot_all$Ymat[v,]
    x1allboot = r1all%*%alpha_est + etaallboot
    x1boot = x1allboot[idx]
    x1sboot = x1allboot[-idx]
    Xboot = cbind(1, x1boot)
    
    
    xie = out1$theta + rmvnorm(1,sigma=sdout$covmat)
    xie = out1$theta
    cov_xialle = cor.mat(Dall, xie[1], cov_xi, nug = 0)
    
    tmp_inv = solve(cov_xialle[-idx,-idx])
    alphaboot = solve(t(r1s) %*% tmp_inv %*% r1s)%*%t(r1s) %*% tmp_inv %*% x1s
    wboot_e = r1%*%alphaboot + cov_xialle[idx, -idx]%*%tmp_inv%*%(x1sboot - r1s%*%alphaboot)
    
    Wboot_e = cbind(1, wboot_e, x2)
    A1_e = solve(t(Wboot_e)%*%CovUBE_GLS_inv%*%Wboot_e)
    EVar = EVar + A1_e
  }
  
  
  beta_rb = beta_RBEGLS
  A1 = solve(t(W)%*%CovUBE_GLS_inv%*%W)
  varbeta_rb = diag(A1)
  varbeta_rbm = diag(EVar/boot_iter)
  beta_rb/varbeta_rb^.5
  beta_rb/varbeta_rbm^.5
  
  beta_kr = kr$beta
  varbeta_kr = diag(kr$beta_var)
  beta_kr/varbeta_kr^.5
  
  
  beta_18 = kr2$beta
  varbeta_18 = diag(kr2$beta_var)
  
  line1 = c(beta_18[-1], beta_kr[-1], beta_rb[-1])
  line2 = c(varbeta_18[-1], varbeta_kr[-1], varbeta_rb[-1], varbeta_rbm[-1])^.5
  
  round(line1, 4)
  round(line2, 4)
  
  return()
}




