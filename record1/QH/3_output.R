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




setwd("C:/Users/qihan/Desktop/QH")
#setwd("C:/Users/qihan/Desktop/QH/B4 folder")

#load("X1000n_200_100.Rdata")

Tab0 = TabK = TabKm = TabR = TabRm = TabKT = TabKTm  = matrix(0,8,6)
TabXi= TabThe = matrix(0,8,3)
for (fi in 1:4){
  nsam = c(200,400, 800,1600)
  n = nsam[fi]
  #n = nsam[1]
  # seed2024_test_X1000n_200_100
  nameRdata=strwrap(paste("seed2024_beta1_X1000n_", n, "_", n/2, ".Rdata", sep=""))
  #load("seed2024_test_X1000n_200_100.Rdata")
  load(nameRdata)
  
  
  betas0 = c(2,beta1, 1)
  out = resTrf(res)
  T0 = outcover(out$BetaOLS, out$VarOLS^.5, betas0)
  TK = outcover(out$BetaKR, out$VarKR^.5, betas0)
  TK2 = outcover(out$BetaKR, out$VarKRm^.5, betas0)
  TR1 = outcover(out$BetaRB, out$VarRB^.5, betas0)
  TR2 = outcover(out$BetaRB, out$VarRBm^.5, betas0)
  TKT = outcover(out$BetaKR, out$VarKRT^.5, betas0)
  TKT2 = outcover(out$BetaKR, out$VarKRTm^.5, betas0)
  
  print(TR2)
  
  Tab0[1:2+(fi-1)*2, ] = T0[-1,]
  TabK[1:2+(fi-1)*2, ] = TK[-1,]
  TabKm[1:2+(fi-1)*2, ] = TK2[-1,]
  TabR[1:2+(fi-1)*2, ] = TR1[-1,]
  TabRm[1:2+(fi-1)*2, ] = TR2[-1,]
  TabKT[1:2+(fi-1)*2, ] = TKT[-1,]
  TabKTm[1:2+(fi-1)*2, ] = TKT2[-1,]
  
  # A indicate MLE in one step
  x1 = apply(out$XiE,2,mean)
  x2 = apply(out$XiE,2,sd)
  TabXi[1:2+(fi-1)*2, ] = cbind(xi0, x1, x2)
  e1 = apply(out$ThetaE,2,mean)
  e2 = apply(out$ThetaE,2,sd)
  TabThe[1:2+(fi-1)*2, ] = cbind(ep_theta0, e1, e2)
}

library("xtable")  
T1 = cbind(TabK, TabR[,-1], TabRm[,-(1:3)])

print(xtable(t(T1), digits = 3))



