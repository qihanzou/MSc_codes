rm(list = ls(all = TRUE))

resTrf = function(res){
  BetaIR  = BetaRB  = NULL
  VarIR = VarRB = VarW  = NULL
  XiRB  = XiIR = AlphaRB = AlphaIR = NULL
  ThetaRB = ThetaIR = RMSE_IR = RMSE_RB = NULL
  
  for(z in 1:1000){
    BetaRB   = cbind(BetaRB, res[[z]]$beta_RBEGLS)
    BetaIR   = cbind(BetaIR, res[[z]]$beta_IREML)
    
    VarRB    = cbind(VarRB, res[[z]]$varbeta_RBEGLS)
    VarIR    = cbind(VarIR, res[[z]]$varbeta_IREML)
    VarW    = cbind(VarW, res[[z]]$varbeta_w ) 
    
    XiRB     = cbind(XiRB, res[[z]]$xi_est)
    XiIR     = cbind(XiIR, res[[z]]$xi_IREML)
    
    AlphaRB  = cbind(AlphaRB, res[[z]]$alpha_est)
    AlphaIR  = cbind(AlphaIR, res[[z]]$alpha_IREML)
    
    ThetaRB  = cbind(ThetaRB, res[[z]]$theta_RBEGLS)
    ThetaIR  = cbind(ThetaIR, res[[z]]$theta_IREML)
    
    RMSE_IR  = cbind(RMSE_IR, res[[z]]$RMSE_IREML)
    RMSE_RB  = cbind(RMSE_RB, res[[z]]$RMSE_RBEGLS)
  }
  return(list(BetaRB = BetaRB, BetaIR = BetaIR, VarRB = VarRB, VarIR = VarIR, VarW = VarW,
              XiRB = XiRB, XiIR = XiIR, AlphaRB = AlphaRB, AlphaIR = AlphaIR,
              ThetaRB = ThetaRB, ThetaIR = ThetaIR, 
              RMSE_IR = RMSE_IR, RMSE_RB = RMSE_RB))
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




setwd("C:/Users/qihan/Desktop/rebound")

nameRdata=strwrap(paste("rebound_1000n_400_200.Rdata"))
load(nameRdata)
betas0 = c(2,beta1)
out = resTrf(res)

TRB = outcover(t(out$BetaRB), t(out$VarRB)^.5, betas0)
TIR = outcover(t(out$BetaIR), t(out$VarIR)^.5, betas0)
TW  = outcover(t(out$BetaIR), t(out$VarIR)^.5, betas0)


print(TRB)
print(TIR)
print(TW)


rowMeans(out$XiRB)
rowMeans(out$XiIR)

rowMeans(out$ThetaRB)
rowMeans(out$ThetaIR)

rowMeans(out$AlphaRB)
rowMeans(out$AlphaIR)

mean(out$RMSE_RB)
mean(out$RMSE_IR)

