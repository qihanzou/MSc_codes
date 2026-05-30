rm(list = ls(all = TRUE))

resTrf = function(res){
  BetaOLS  = BetaKR = BetaRB = XiE = ThetaE = NULL
  VarOLS = VarKR = VarKRm = VarRB = VarRBm = VarKRT = VarKRTm = NULL
  for(i in 1:nreps){
    BetaOLS = rbind(BetaOLS, res[[i]]$beta_ols)
    BetaKR = cbind(BetaKR, res[[i]]$beta_kr)
    BetaRB = cbind(BetaRB, res[[i]]$beta_RBEGLS)
    
    VarOLS = rbind(VarOLS, res[[i]]$varbeta_ols)
    VarKR = rbind(VarKR, res[[i]]$varbeta_kr)
    VarKRm = rbind(VarKRm, res[[i]]$varbeta_krm)
    VarRB = rbind(VarRB, res[[i]]$varbeta_rb)
    VarRBm = rbind(VarRBm, res[[i]]$varbeta_rbm)
    VarKRT = rbind(VarKRT, res[[i]]$varbeta_krt)
    VarKRTm = rbind(VarKRTm, res[[i]]$varbeta_krtm)
    #Var3 = rbind(Var3, res[[i]]$varbeta_cl2)
    #ThetaX = rbind(ThetaX, res[[i]]$thetax)
    XiE = rbind(XiE, res[[i]]$xi_est)
    ThetaE = rbind(ThetaE, res[[i]]$theta_est)
    
  }
  return(list(BetaOLS = BetaOLS,BetaKR = t(BetaKR), BetaRB = t(BetaRB), 
              VarOLS = VarOLS, VarKR = VarKR, VarKRm = VarKRm, VarRB = VarRB, VarRBm = VarRBm,
              VarKRT = VarKRT, VarKRTm = VarKRTm, XiE = XiE, ThetaE = ThetaE))
}

outcover = function(TP, TSD, betas0){
  nsim = dim(TP)[1]
  bmean = apply(TP,2,mean)
  rMSE = sqrt(apply(TP,2,sd)^2 + (apply(TP,2,mean) - betas0)^2)
  SEm = apply(TSD, 2, mean)
  
  Betas0 = matrix(betas0, nsim, dim(TP)[2], byrow = T)
  UP = TP + qnorm(0.95)*TSD
  LO = TP - qnorm(0.95)*TSD
  CPT =  (UP > Betas0 & LO<Betas0)
  CP90 = apply(CPT, 2, mean)
  
  Betas0 = matrix(betas0, nsim, dim(TP)[2], byrow = T)
  UP = TP + qnorm(0.975)*TSD
  LO = TP - qnorm(0.975)*TSD
  CPT =  (UP > Betas0 & LO<Betas0)
  CP95 = apply(CPT, 2, mean)
  
  Tab = cbind(betas0, bmean, rMSE, SEm, CP90, CP95)
  return(Tab)
}




setwd("C:/Users/qihan/Desktop/QH_new")

  n = 1600
  nameRdata=strwrap(paste("new3_X1000n_", n, "_", n/2, ".Rdata", sep=""))
  load(nameRdata)
  betas0 = c(2,beta1, 1, 3, 5)
  out = resTrf(res)
  T0 = outcover(out$BetaOLS, out$VarOLS^.5, betas0)
  TK = outcover(out$BetaKR, out$VarKR^.5, betas0)
  TK2 = outcover(out$BetaKR, out$VarKRm^.5, betas0)
  TR1 = outcover(out$BetaRB, out$VarRB^.5, betas0)
  TR2 = outcover(out$BetaRB, out$VarRBm^.5, betas0)
  #TKT = outcover(out$BetaKR, out$VarKRT^.5, betas0)
  #TKT2 = outcover(out$BetaKR, out$VarKRTm^.5, betas0)
  
  
  
  print(T0)
  print(TK)
  #print(TK2)
  print(TR1)
  print(TR2)