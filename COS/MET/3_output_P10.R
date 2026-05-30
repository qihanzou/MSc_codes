rm(list = ls(all = TRUE))

resTrf = function(res){
  BetaOLS  = BetaRB_dm = BetaRB_mc = BetaRB_mean = BetaKR = NULL
  VarOLS = VarRB_mean  = VarRB_dm = VarKR  = VarRB_mc = VarRB_dm_mc = NULL
  XiRB  = ALPHARB  = NULL
  ThetaRB_dm = ThetaRB_mc = ThetaRB_mean = ThetaKR  = Bdist  = NULL
  VarRB_dm_mc_mean = VarRB_mc_mc_mean = VarRB_dm_mc_median = VarRB_mc_mc_median = NULL
  VarRB_dm_mc_mean2 = VarRB_mc_mc_mean2 = VarRB_dm_mc_median2 = VarRB_mc_mc_median2 = NULL
  
  for(z in 1:1000){
    XiRB     = cbind(XiRB, res[[z]]$xi_est)
    ALPHARB  = cbind(ALPHARB, res[[z]]$alpha_est)
    
    BetaOLS  = cbind(BetaOLS, res[[z]]$beta_ols)
    BetaKR   = cbind(BetaKR, res[[z]]$beta_kr)
    BetaRB_dm   = cbind(BetaRB_dm, res[[z]]$beta_RBEGLS)
    BetaRB_mc   = cbind(BetaRB_mc, res[[z]]$beta_RBEGLS_mc)
    
    ThetaKR  = cbind(ThetaKR, res[[z]]$theta_kr)
    ThetaRB_dm  = cbind(ThetaRB_dm, res[[z]]$theta_est)
    ThetaRB_mc  = cbind(ThetaRB_mc, res[[z]]$theta_est_mc)
    
    VarOLS   = cbind(VarOLS, res[[z]]$varbeta_ols)
    VarKR    = cbind(VarKR, res[[z]]$varbeta_kr)
    
    VarRB_mc    = cbind(VarRB_mc, res[[z]]$varbeta_rb_mct)
    VarRB_dm    = cbind(VarRB_dm, res[[z]]$varbeta_rb_dm)
    
    VarRB_dm_mc_mean = cbind(VarRB_dm_mc_mean, res[[z]]$varbeta_mc_mean_dm_M)
    VarRB_mc_mc_mean = cbind(VarRB_mc_mc_mean, res[[z]]$varbeta_mc_mean_mc_M) 
    VarRB_dm_mc_median = cbind(VarRB_dm_mc_median, res[[z]]$varbeta_mc_median_dm_M)
    VarRB_mc_mc_median = cbind(VarRB_mc_mc_median, res[[z]]$varbeta_mc_median_mc_M) 
    
    VarRB_dm_mc_mean2 = cbind(VarRB_dm_mc_mean2, res[[z]]$varbeta_mc_mean_dm2_M)
    VarRB_mc_mc_mean2 = cbind(VarRB_mc_mc_mean2, res[[z]]$varbeta_mc_mean_mc2_M) 
    VarRB_dm_mc_median2 = cbind(VarRB_dm_mc_median2, res[[z]]$varbeta_mc_median_dm2_M)
    VarRB_mc_mc_median2 = cbind(VarRB_mc_mc_median2, res[[z]]$varbeta_mc_median_mc2_M) 

  }
  return(list(BetaOLS = BetaOLS, BetaKR = BetaKR,
              BetaRB_dm = BetaRB_dm, BetaRB_mc = BetaRB_mc,
              VarOLS = VarOLS, VarKR = VarKR, 
              VarRB_dm = VarRB_dm, VarRB_mc = VarRB_mc,
              VarRB_dm_mc_mean = VarRB_dm_mc_mean, VarRB_mc_mc_mean = VarRB_mc_mc_mean, 
              VarRB_dm_mc_median = VarRB_dm_mc_median, VarRB_mc_mc_median = VarRB_mc_mc_median,
              VarRB_dm_mc_mean2 = VarRB_dm_mc_mean2, VarRB_mc_mc_mean2 = VarRB_mc_mc_mean2, 
              VarRB_dm_mc_median2 = VarRB_dm_mc_median2, VarRB_mc_mc_median2 = VarRB_mc_mc_median2,
              ALPHARB = ALPHARB, XiRB = XiRB,
              ThetaRB_dm = ThetaRB_dm, ThetaRB_mc = ThetaRB_mc, ThetaKR = ThetaKR))
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




setwd("C:/Users/qihan/Desktop/MET")

nameRdata=strwrap(paste("MET_cos_1000n_800_400.Rdata"))
load(nameRdata)
betas0 = c(2,beta1)
out = resTrf(res)
TO = outcover(t(out$BetaOLS), t(out$VarOLS)^.5, betas0)
TKR = outcover(t(out$BetaKR), t(out$VarKR)^.5, betas0)


TRB_dm = outcover(t(out$BetaRB_dm), t(out$VarRB_dm)^.5, betas0)
TRB_mc = outcover(t(out$BetaRB_mc), t(out$VarRB_mc)^.5, betas0)

TRB_dm_mc_mean = outcover(t(out$BetaRB_dm), t(out$VarRB_dm_mc_mean)^.5, betas0)
TRB_dm_mc_median = outcover(t(out$BetaRB_dm), t(out$VarRB_dm_mc_median)^.5, betas0)
TRB_mc_mc_mean = outcover(t(out$BetaRB_mc), t(out$VarRB_mc_mc_mean)^.5, betas0)
TRB_mc_mc_median = outcover(t(out$BetaRB_mc), t(out$VarRB_mc_mc_median)^.5, betas0)

TRB_dm_mc_mean2 = outcover(t(out$BetaRB_dm), t(out$VarRB_dm_mc_mean2)^.5, betas0)
TRB_dm_mc_median2 = outcover(t(out$BetaRB_dm), t(out$VarRB_dm_mc_median2)^.5, betas0)
TRB_mc_mc_mean2 = outcover(t(out$BetaRB_mc), t(out$VarRB_mc_mc_mean2)^.5, betas0)
TRB_mc_mc_median2 = outcover(t(out$BetaRB_mc), t(out$VarRB_mc_mc_median2)^.5, betas0)


# rowMeans(out$ALPHARB)
# rowMeans(out$XiRB)
# rowMeans(out$ThetaKR)
# rowMeans(out$ThetaRB_dm)
# rowMeans(out$ThetaRB_mc)

print(TO)
print(TKR)
print(TRB_dm)
print(TRB_mc)
print(TRB_dm_mc_mean)
print(TRB_dm_mc_median)
print(TRB_mc_mc_mean)
print(TRB_mc_mc_median)
print(TRB_dm_mc_mean2)
print(TRB_dm_mc_median2)
print(TRB_mc_mc_mean2)
print(TRB_mc_mc_median2)


library("xtable")  
T1 = cbind(t(TKR), t(TRB_dm), t(TRB_mc), 
           t(TRB_dm_mc_mean), t(TRB_dm_mc_median), t(TRB_mc_mc_mean), t(TRB_mc_mc_median),
           t(TRB_dm_mc_mean2), t(TRB_dm_mc_median2), t(TRB_mc_mc_mean2), t(TRB_mc_mc_median2))

print(xtable(t(T1), digits = 3))

