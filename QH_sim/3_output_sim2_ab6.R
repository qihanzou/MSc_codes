rm(list = ls(all = TRUE))

resTrf = function(res){
  BetaRB = VarRB = XiRB = Xisim = Xiout = Varxi = Varxisim = Varxiout = Nout = NULL
  
  for(z in 1:1000){
    BetaRB   = cbind(BetaRB, res[[z]]$beta_RBEGLS)
    VarRB    = cbind(VarRB, res[[z]]$varbeta_rb)
    
    XiRB     = cbind(XiRB, res[[z]]$xi_est)
    Xisim    = cbind(Xisim, res[[z]]$xi_sim)
    Xiout    = cbind(Xiout, res[[z]]$xi_out)
    
    Varxi  = cbind(Varxi, res[[z]]$varxi)
    Varxisim  = cbind(Varxisim, res[[z]]$varxi_sim)
    Varxiout  = cbind(Varxiout, res[[z]]$varxi_out)
    Nout = cbind(Nout, res[[z]]$n_out)
    
  }
  return(list(BetaRB = BetaRB, VarRB = VarRB, XiRB=XiRB, Xisim = Xisim, 
              Xiout = Xiout, Varxi = Varxi, Varxisim = Varxisim,
              Varxiout = Varxiout, Nout = Nout))
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




setwd("C:/Users/qihan/Desktop/QH_sim2_ab6")

nameRdata=strwrap(paste("SIM_ab6_1000n_200_100.Rdata"))
load(nameRdata)
betas0 = c(3,4)
out = resTrf(res)
XIRB =  outcover(t(out$XiRB),  t(out$Varxi)^.5, betas0)
XISIM = outcover(t(out$Xisim), t(out$Varxi)^.5, betas0)
XIOUT = outcover(t(out$Xiout), t(out$Varxi)^.5, betas0)

print(XIRB)
print(XISIM)
print(XIOUT)


rowMeans(out$XiRB)
rowMeans(out$Xisim)
rowMeans(out$Xiout)





