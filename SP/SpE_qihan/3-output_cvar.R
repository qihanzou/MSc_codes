resTrf = function(res){
  Beta1  = Beta2 = Beta3 = Var1 = Var2 = Var3 =Var4 = cvar = ThetaE0 = ThetaCL = ThetaKR = NULL
  Alpha = Xi = Bols = NULL
  for(i in 1:nreps){
    Beta1 = rbind(Beta1, res[[i]]$beta_update)
    
    Var1 = rbind(Var1, res[[i]]$varbeta_cl0)
    Var2 = rbind(Var2, res[[i]]$varbeta_cl1)
    
    cvar = rbind(cvar, res[[i]]$cvar)
  
  }
  return(list(Beta1=Beta1, 
              Var1 =Var1, Var2=Var2,
              cvar = cvar))
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


w1_address = "C:/Users/qihan/Desktop/SpE_ncl"

setwd(w1_address)
load("cvar_nreps1000_n_1000_500_ver_1.Rdata")
out = resTrf(res)

betas0 = c(2,4)
(Tab1 = outcover(out$Beta1, out$Var1^.5, betas0))
(Tab2 = outcover(out$Beta1, out$Var2^.5, betas0))


colMeans(out$cvar)








