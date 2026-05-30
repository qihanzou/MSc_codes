rm(list = ls(all = TRUE))

resTrf = function(res){
  BetaOLS  = BetaRB_dm = BetaRB_mc = BetaRB_mean = BetaKR = NULL
  VarOLS = VarRB_mean  = VarRB_dm = VarKR  = VarRB_mc = VarRB_dm_mc = NULL
  XiRB  = ALPHARB  = NULL
  ThetaRB_dm = ThetaRB_mc = ThetaRB_mean = ThetaKR  = Bdist  = NULL
  
  for(z in 1:1000){
    XiRB     = cbind(XiRB, res[[z]]$xi_est)
    ALPHARB  = cbind(ALPHARB, res[[z]]$alpha_est)
    
    BetaOLS  = cbind(BetaOLS, res[[z]]$beta_ols)
    BetaKR   = cbind(BetaKR, res[[z]]$beta_kr)
    BetaRB_dm   = cbind(BetaRB_dm, res[[z]]$beta_RBEGLS)
    BetaRB_mc   = cbind(BetaRB_mc, res[[z]]$beta_RBEGLS_mc)
    BetaRB_mean = cbind(BetaRB_mean, res[[z]]$beta_mean)
    
    ThetaKR  = cbind(ThetaKR, res[[z]]$theta_kr)
    ThetaRB_dm  = cbind(ThetaRB_dm, res[[z]]$theta_est)
    ThetaRB_mc  = cbind(ThetaRB_mc, res[[z]]$theta_est_mc)
    ThetaRB_mean  = cbind(ThetaRB_mean, res[[z]]$theta_mean)
    
    VarOLS   = cbind(VarOLS, res[[z]]$varbeta_ols)
    VarKR    = cbind(VarKR, res[[z]]$varbeta_kr)
    VarRB_mean    = cbind(VarRB_mean, res[[z]]$varbeta_rb_mean)
    VarRB_dm    = cbind(VarRB_dm, res[[z]]$varbeta_rb_dm)
    VarRB_dm_mc    = cbind(VarRB_dm_mc, res[[z]]$varbeta_rb_dm_m)
    VarRB_mc    = cbind(VarRB_mc, res[[z]]$varbeta_rb_mc) 
    
    Bdist = cbind(Bdist, res[[z]]$bdist)
    
  }
  return(list(BetaOLS = BetaOLS, BetaKR = BetaKR,
              BetaRB_dm = BetaRB_dm, BetaRB_mc = BetaRB_mc, BetaRB_mean = BetaRB_mean,
              VarOLS = VarOLS, VarKR = VarKR, 
              VarRB_mean = VarRB_mean, VarRB_dm = VarRB_dm, VarRB_mc = VarRB_mc,VarRB_dm_mc = VarRB_dm_mc,
              ALPHARB = ALPHARB, XiRB = XiRB, Bdist = Bdist,
              ThetaRB_dm = ThetaRB_dm, ThetaRB_mc = ThetaRB_mc, ThetaRB_mean = ThetaRB_mean, ThetaKR = ThetaKR))
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




setwd("C:/Users/qihan/Desktop/Transformation_and_ME")

nameRdata=strwrap(paste("nl2_exp_1000n_200_100.Rdata"))
load(nameRdata)
betas0 = c(2,beta1)
out = resTrf(res)
TO = outcover(t(out$BetaOLS), t(out$VarOLS)^.5, betas0)
TKR = outcover(t(out$BetaKR), t(out$VarKR)^.5, betas0)


TRB_dm = outcover(t(out$BetaRB_dm), t(out$VarRB_dm)^.5, betas0)
TRB_mc = outcover(t(out$BetaRB_mc), t(out$VarRB_mc)^.5, betas0)
TRB_mean = outcover(t(out$BetaRB_mean), t(out$VarRB_mean)^.5, betas0)
TRB_dm_mc = outcover(t(out$BetaRB_dm), t(out$VarRB_dm_mc)^.5, betas0)

rowMeans(out$ALPHARB)
rowMeans(out$XiRB)
rowMeans(out$ThetaKR)
rowMeans(out$ThetaRB_dm)
rowMeans(out$ThetaRB_mc)
rowMeans(out$ThetaRB_mean)
rowMeans(out$Bdist)

print(TO)
print(TKR)
print(TRB_dm)
print(TRB_mc)
print(TRB_mean)
print(TRB_dm_mc)



