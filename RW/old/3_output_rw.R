rm(list = ls(all = TRUE))

resTrf = function(res){
  BetaOLS = BetaRW = BetaRB =  NULL
  VarOLS = VarRB = VarRW = NULL
  XiRB = XiRW = ALPHARB = ALPHARW = NULL
  ThetaRB = ThetaRW = NumEpo = NllRB = NllRW = SucRW = NULL
  
  for(i in 1:nreps){
    BetaOLS = cbind(BetaOLS, res[[i]]$beta_ols)
    BetaRB = cbind(BetaRB, res[[i]]$beta_RBEGLS)
    BetaRW = cbind(BetaRW, res[[i]]$beta_rw)
    
    VarOLS = cbind(VarOLS, res[[i]]$varbeta_ols)
    VarRB = cbind(VarRB, res[[i]]$varbeta_rb)
    VarRW = cbind(VarRW, res[[i]]$varbeta_rw)
    
    
    XiRB = cbind(XiRB, res[[i]]$xi_est)
    XiRW = cbind(XiRW, res[[i]]$xi_rw)
    
    ALPHARB = cbind(ALPHARB, res[[i]]$alpha_est)
    ALPHARW = cbind(ALPHARW, res[[i]]$alpha_rw)
    
    ThetaRB = cbind(ThetaRB, res[[i]]$theta_est)
    ThetaRW = cbind(ThetaRW, res[[i]]$theta_rw)
    
    NumEpo = rbind(NumEpo, res[[i]]$num_epoches)
    NllRB = rbind(NllRB, res[[i]]$nll_rb_full)
    NllRW = rbind(NllRW, res[[i]]$nll_rw)
    SucRW = rbind(SucRW, res[[i]]$rw_update_exist)
    

  }
  return(list(BetaOLS = BetaOLS, BetaRB = BetaRB, BetaRW = BetaRW, 
              VarOLS = VarOLS, VarRB = VarRB, VarRW = VarRW, 
              XiRB = XiRB, XiRW = XiRW,
              ALPHARB = ALPHARB, ALPHARW = ALPHARW,
              ThetaRB = ThetaRB, ThetaRW = ThetaRW,
              NumEpo = NumEpo, NllRB = NllRB, NllRW = NllRW, SucRW = SucRW))
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
 #rw_test1000n_200_100
  n = 400
  nameRdata=strwrap(paste("rw_test_2_1000n_", n, "_", n/2, ".Rdata", sep=""))
  load(nameRdata)
  betas0 = c(2,beta1)
  out = resTrf(res)
  TO = outcover(t(out$BetaOLS), t(out$VarOLS)^.5, betas0)
  TRB = outcover(t(out$BetaRB), t(out$VarRB)^.5, betas0)
  TRW = outcover(t(out$BetaRW), t(out$VarRW)^.5, betas0)

  

  
  rowMeans(out$XiRB)
  rowMeans(out$XiRW)
  
  rowMeans(out$ALPHARB)
  rowMeans(out$ALPHARW)
  
  rowMeans(out$ThetaRB)
  rowMeans(out$ThetaRW)

  mean(out$NumEpo)
  mean(out$NllRB)
  mean(out$NllRW)
  mean(out$SucRW)
  
  
  print(TO)
  print(TRB)
  print(TRW)
  