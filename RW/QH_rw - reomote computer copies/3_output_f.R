rm(list = ls(all = TRUE))

resTrf = function(res){
  BetaOLS  = BetaRB =  BetaKR = NULL
  VarOLS = VarRB = VarRW1  = Var_v3_85 = Var_v3_90 = Var_v3_95 = Var_v3_m = VarMC = VarKR  = NULL
  XiRB  = ALPHARB = Var_v3_m859095 = Var_v3_m9095 = Var_v3_m8590 = NULL
  ThetaRB = ThetaKR  = NULL
  
  for(z in 1:1000){
    BetaOLS  = cbind(BetaOLS, res[[z]]$beta_ols)
    BetaRB   = cbind(BetaRB, res[[z]]$beta_RBEGLS)
    BetaKR   = cbind(BetaKR, res[[z]]$beta_kr)
    
    VarOLS   = cbind(VarOLS, res[[z]]$varbeta_ols)
    VarKR    = cbind(VarKR, res[[z]]$varbeta_kr)
    VarRB    = cbind(VarRB, res[[z]]$varbeta_rb)
    VarMC    = cbind(VarMC, res[[z]]$varbeta_mc)
    Var_v3_85 = cbind(Var_v3_85, res[[z]]$varbeta_v3_85)
    Var_v3_90 = cbind(Var_v3_90, res[[z]]$varbeta_v3_90)
    Var_v3_95 = cbind(Var_v3_95, res[[z]]$varbeta_v3_95)
    Var_v3_m859095 = cbind(Var_v3_m859095, res[[z]]$varbeta_v3_zmean859095)
    Var_v3_m9095 = cbind(Var_v3_m9095, res[[z]]$varbeta_v3_zmean9095)
    Var_v3_m8590 = cbind(Var_v3_m8590, res[[z]]$varbeta_v3_zmean8590)
    
    XiRB     = cbind(XiRB, res[[z]]$xi_est)
    ALPHARB  = cbind(ALPHARB, res[[z]]$alpha_est)
    ThetaRB  = cbind(ThetaRB, res[[z]]$theta_est)
    ThetaKR  = cbind(ThetaKR, res[[z]]$theta_kr)
    
  }
  return(list(BetaOLS = BetaOLS, BetaRB = BetaRB, BetaKR = BetaKR,
              VarOLS = VarOLS, VarKR = VarKR, VarRB = VarRB, VarMC = VarMC,
              Var_v3_85 = Var_v3_85, Var_v3_90 = Var_v3_90, Var_v3_95 = Var_v3_95,
              XiRB = XiRB, Var_v3_m859095 = Var_v3_m859095, Var_v3_m9095 = Var_v3_m9095, Var_v3_m8590 = Var_v3_m8590,
              ALPHARB = ALPHARB, 
              ThetaRB = ThetaRB, ThetaKR = ThetaKR))
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




setwd("C:/Users/qihan/Desktop/QH_rw - reomote computer copies")
n = 1600
nameRdata=strwrap(paste("Qihan_f_1000n_", n, "_", n/2, "_mc100.Rdata", sep=""))
load(nameRdata)
betas0 = c(2,beta1)
out = resTrf(res)
TO = outcover(t(out$BetaOLS), t(out$VarOLS)^.5, betas0)
TKR = outcover(t(out$BetaKR), t(out$VarKR)^.5, betas0)
TRB = outcover(t(out$BetaRB), t(out$VarRB)^.5, betas0)
TMC = outcover(t(out$BetaRB), t(out$VarMC)^.5, betas0)
TRB_v3_85 = outcover(t(out$BetaRB), t(out$Var_v3_85)^.5, betas0)
TRB_v3_90 = outcover(t(out$BetaRB), t(out$Var_v3_90)^.5, betas0)
TRB_v3_95 = outcover(t(out$BetaRB), t(out$Var_v3_95)^.5, betas0)
TRB_v3_m859095 = outcover(t(out$BetaRB), t(out$Var_v3_m859095)^.5, betas0)
TRB_v3_m9095 = outcover(t(out$BetaRB), t(out$Var_v3_m9095)^.5, betas0)
TRB_v3_m8590 = outcover(t(out$BetaRB), t(out$Var_v3_m8590)^.5, betas0)


# ---------------------------------------------------------------------------- #
print(TO)
print(TKR)
print(TRB)
print(TMC)

print(TRB_v3_85)
print(TRB_v3_90)
print(TRB_v3_95)

print(TRB_v3_m859095)
print(TRB_v3_m9095)
print(TRB_v3_m8590)
# ---------------------------------------------------------------------------- #

 
library("xtable")  
T1 = cbind(t(TO), t(TKR), t(TRB), t(TMC), t(TRB_v3_85), t(TRB_v3_90), t(TRB_v3_95), t(TRB_v3_m8590), t(TRB_v3_m9095), t(TRB_v3_m859095))
print(xtable(t(T1), digits = 3))



