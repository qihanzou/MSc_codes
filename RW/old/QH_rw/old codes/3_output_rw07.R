rm(list = ls(all = TRUE))

resTrf = function(res){
  BetaOLS = BetaRW = BetaRB =  BetaKR = NULL
  VarOLS = VarRB = VarRW1 = VarRW2 = VarRW3 = VarRW4 = VarRW5 = VarRW_v2 = VarRW_v3 = VarMC = VarKR = VarSI = NULL
  XiRB  = ALPHARB = NULL
  ThetaRB = ThetaKR = ThetaRW = NULL
  
  for(z in 1:1000){
    BetaOLS  = cbind(BetaOLS, res[[z]]$beta_ols)
    BetaRB   = cbind(BetaRB, res[[z]]$beta_RBEGLS)
    BetaKR   = cbind(BetaKR, res[[z]]$beta_kr)
    BetaRW   = cbind(BetaRW, res[[z]]$beta_rw)
    
    VarOLS   = cbind(VarOLS, res[[z]]$varbeta_ols)
    VarKR    = cbind(VarKR, res[[z]]$varbeta_kr)
    VarRB    = cbind(VarRB, res[[z]]$varbeta_rb)
    VarMC    = cbind(VarMC, res[[z]]$varbeta_mc)
    VarRW1   = cbind(VarRW1, res[[z]]$varbeta_rw_mean)
    VarRW_v2 = cbind(VarRW_v2, res[[z]]$varbeta_rw_v2_mean)
    VarRW_v3 = cbind(VarRW_v3, res[[z]]$varbeta_rw_v3_mean)
    
    XiRB     = cbind(XiRB, res[[z]]$xi_est)
    ALPHARB  = cbind(ALPHARB, res[[z]]$alpha_est)
    ThetaRB  = cbind(ThetaRB, res[[z]]$theta_est)
    ThetaKR  = cbind(ThetaKR, res[[z]]$theta_kr)
    ThetaRW  = cbind(ThetaRW, res[[z]]$theta_rw)
    
  }
  return(list(BetaOLS = BetaOLS, BetaRB = BetaRB, BetaKR = BetaKR, BetaRW = BetaRW,
              VarOLS = VarOLS, VarKR = VarKR, VarRB = VarRB, VarRW1 = VarRW1, VarMC = VarMC,
              VarRW_v2 = VarRW_v2, VarRW_v3 = VarRW_v3,
              XiRB = XiRB,
              ALPHARB = ALPHARB, 
              ThetaRB = ThetaRB, ThetaKR = ThetaKR, ThetaRW = ThetaRW))
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




setwd("C:/Users/qihan/Desktop/QH_rw")
n = 800
nameRdata=strwrap(paste("Qihan_rw07_1000n_", n, "_", n/2, "_mc100_ci975.Rdata", sep=""))
load(nameRdata)
betas0 = c(2,beta1)
out = resTrf(res)
TO = outcover(t(out$BetaOLS), t(out$VarOLS)^.5, betas0)
TKR = outcover(t(out$BetaKR), t(out$VarKR)^.5, betas0)
TRB = outcover(t(out$BetaRB), t(out$VarRB)^.5, betas0)
TMC = outcover(t(out$BetaRB), t(out$VarMC)^.5, betas0)


TRW_rb_v3 = outcover(t(out$BetaRB), t(out$VarRW_v3)^.5, betas0)

TRW_rw_v3 = outcover(t(out$BetaRW), t(out$VarRW_v3)^.5, betas0)


rowMeans(out$XiRB)

rowMeans(out$ALPHARB)

rowMeans(out$ThetaRB)



print(TO)
print(TKR)
print(TRB)
print(TMC)
print(TRW_rb_v3)
print(TRW_rw_v3)

library("xtable")  
T1 = cbind(t(TO), t(TKR), t(TRB), t(TMC), t(TRW_rb_v3))

print(xtable(t(T1), digits = 3))



rowMeans(out$BetaRB)
rowMeans(out$BetaRW)
