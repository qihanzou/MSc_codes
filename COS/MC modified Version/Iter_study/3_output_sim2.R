rm(list = ls(all = TRUE))

resTrf = function(res){
  BetaOLS  = BetaRB =  BetaKR = NULL
  VarOLS = VarRB = VarRW1  = Var_v3_85 = Var_v3_90 = Var_v3_95 = Var_v3_m = VarMC = VarKR  = NULL
  XiRB  = ALPHARB = VarSI = Var_v3_m9095 = Var_v3_m8590 = NULL
  ThetaRB = ThetaKR  = Nout = XiSI = ThetaSI = NULL
  
  for(z in 1:1000){
    BetaOLS  = cbind(BetaOLS, res[[z]]$beta_ols)
    BetaRB   = cbind(BetaRB, res[[z]]$beta_RBEGLS)
    BetaKR   = cbind(BetaKR, res[[z]]$beta_kr)
    
    VarOLS   = cbind(VarOLS, res[[z]]$varbeta_ols)
    VarKR    = cbind(VarKR, res[[z]]$varbeta_kr)
    VarRB    = cbind(VarRB, res[[z]]$varbeta_rb)
    VarMC    = cbind(VarMC, res[[z]]$varbeta_mc2)
    VarSI    = cbind(VarSI, res[[z]]$varbeta_A1 ) 
    
    XiRB     = cbind(XiRB, res[[z]]$xi_est)
    XiSI     = cbind(XiSI, res[[z]]$samples_outliers_mean)
    ALPHARB  = cbind(ALPHARB, res[[z]]$alpha_est)
    ThetaRB  = cbind(ThetaRB, res[[z]]$theta_est)
    ThetaKR  = cbind(ThetaKR, res[[z]]$theta_kr)
    ThetaSI  = cbind(ThetaSI, res[[z]]$theta_gls)
    Nout = cbind(Nout, res[[z]]$num_outliers)
    
  }
  return(list(BetaOLS = BetaOLS, BetaRB = BetaRB, BetaKR = BetaKR,
              VarOLS = VarOLS, VarKR = VarKR, VarRB = VarRB, VarMC = VarMC,
              ALPHARB = ALPHARB, VarSI = VarSI, XiRB = XiRB, XiSI = XiSI,
              ThetaRB = ThetaRB, ThetaKR = ThetaKR, Nout = Nout, ThetaSI = ThetaSI))
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




setwd("C:/Users/qihan/Desktop/QH_sim2_ab4_iter_study")

nameRdata=strwrap(paste("SIM_ab4_MCiter_100_reps_1000n_200_100.Rdata"))
load(nameRdata)
betas0 = c(2,beta1)
out = resTrf(res)
TO = outcover(t(out$BetaOLS), t(out$VarOLS)^.5, betas0)
TKR = outcover(t(out$BetaKR), t(out$VarKR)^.5, betas0)
TRB = outcover(t(out$BetaRB), t(out$VarRB)^.5, betas0)
TMC = outcover(t(out$BetaRB), t(out$VarMC)^.5, betas0)

print(TO)
print(TKR)
print(TRB)
print(TMC)
#print(TSI)


rowMeans(out$XiRB)

rowMeans(out$ThetaKR)
rowMeans(out$ThetaRB)

rowMeans(out$ALPHARB)


library("xtable")  
T1 = cbind(t(TO), t(TKR), t(TRB), t(TMC))
print(xtable(t(T1), digits = 3))



