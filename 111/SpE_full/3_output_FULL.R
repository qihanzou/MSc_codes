library(xtable)

resTrf = function(res){
  Beta1  = Beta2 = Beta3 = Beta4 = Beta5 = Beta6 = Var1 = Var2 = Var3 = Var4 = Var5 = Var6= ThetaX = ThetaE = ThetaE1 = ThetaE2 = ThetaW2 = ThetaW3 = NULL
  for(i in 1:nreps){
    Beta1 = cbind(Beta1, res[[i]]$beta_est)
    Beta2 = cbind(Beta2, res[[i]]$beta_update)
    Beta3 = cbind(Beta3, t(res[[i]]$beta_RBEGLS))
    Beta4 = cbind(Beta4, res[[i]]$beta_uk)
    Beta5 = cbind(Beta5, res[[i]]$beta_w2)
    Beta6 = cbind(Beta6, res[[i]]$beta_w3)
    
    Var1 = rbind(Var1, res[[i]]$varbeta1)
    Var2 = rbind(Var2, res[[i]]$varbeta2)
    Var3 = rbind(Var3, res[[i]]$varbeta3)
    Var4 = rbind(Var4, res[[i]]$varbeta_uk)
    Var5 = rbind(Var5, t(res[[i]]$verbeta_w2))
    Var6 = rbind(Var6, t(res[[i]]$verbeta_w2))
    
    ThetaX = rbind(ThetaX, res[[i]]$thetax)
    ThetaE = rbind(ThetaE, res[[i]]$thetae)
    ThetaE1 = rbind(ThetaE1, res[[i]]$theta_RBEGLS)
    ThetaE2 = rbind(ThetaE2, res[[i]]$theta_uk)
    ThetaW2 = rbind(ThetaW2, res[[i]]$theta_w2)
    ThetaW3 = rbind(ThetaW3, res[[i]]$theta_w3)
    
  }
  return(list(Beta1=t(Beta1), Beta2=t(Beta2), Beta3 = t(Beta3), Beta4 = t(Beta4), Beta5 = t(Beta5), Beta6 = t(Beta6),
              Var1 = Var1, Var2 =Var2, Var3 = Var3, Var4 = Var4, Var5 = Var5, Var6 = Var6,
              ThetaX = ThetaX, ThetaE = ThetaE, ThetaE2 = ThetaE2, ThetaE1 = ThetaE1, ThetaW2 = ThetaW2, ThetaW3 = ThetaW3))
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

outcover_bind8 = function(TP1, TP2, TP3, TSD1, TSD2, TSD3, betas0){
  TP = TP3
  bmean = apply(TP,2,mean)
  rMSE = sqrt(apply(TP,2,sd)^2 + (apply(TP,2,mean) - betas0)^2)
  nsim = dim(TP)[1]
  Betas0 = matrix(betas0, nsim, dim(TP)[2], byrow = T)
  bias = bmean - betas0
  
  b0UP = NULL
  b0LO = NULL
  b1UP = NULL
  b1LO = NULL
  
  b0UP1 = TP1[,1] + qnorm(0.95)*TSD1[,1]
  b0UP2 = TP2[,1] + qnorm(0.95)*TSD2[,1]
  b0UP3 = TP3[,1] + qnorm(0.95)*TSD3[,1]
  
  b1UP1 = TP1[,2] + qnorm(0.95)*TSD1[,2]
  b1UP2 = TP2[,2] + qnorm(0.95)*TSD2[,2]
  b1UP3 = TP3[,2] + qnorm(0.95)*TSD3[,2]
  
  for (i in 1:length(b0UP1)){
    largest <- mean(c(b0UP1[i], b0UP2[i], b0UP3[i]))
    b0UP = cbind(b0UP, largest)
  }
  
  for (i in 1:length(b1UP1)){
    largest <- mean(c(b1UP1[i], b1UP2[i], b1UP3[i]))
    b1UP = cbind(b1UP, largest)
  }
  
  UP = cbind(t(b0UP),t(b1UP))
  
  
  b0LO1 = TP1[,1] - qnorm(0.95)*TSD1[,1]
  b0LO2 = TP2[,1] - qnorm(0.95)*TSD2[,1]
  b0LO3 = TP3[,1] - qnorm(0.95)*TSD3[,1]
  
  b1LO1 = TP1[,2] - qnorm(0.95)*TSD1[,2]
  b1LO2 = TP2[,2] - qnorm(0.95)*TSD2[,2]
  b1LO3 = TP3[,2] - qnorm(0.95)*TSD3[,2]
  
  
  for (i in 1:length(b0LO1)){
    smallest <- mean(c(b0LO1[i], b0LO2[i], b0LO3[i]))
    b0LO = cbind(b0LO, smallest)
  }
  
  for (i in 1:length(b1LO1)){
    smallest <- mean(c(b1LO1[i], b1LO2[i], b1LO3[i]))
    b1LO = cbind(b1LO, smallest)
  }
  
  LO = cbind(t(b0LO),t(b1LO))
  
  b0TP90 = (b0LO+b0UP)/2
  b1TP90 = (b1LO+b1UP)/2
  bmean_90 = t(cbind(mean(b0TP90),mean(b1TP90)))
  
  b0TSD90 = as.vector(as.numeric((b0UP - mean(b0TP90))/qnorm(0.95)))
  b1TSD90 = as.vector(as.numeric((b1UP - mean(b1TP90))/qnorm(0.95)))
  
  SEm_90 = t(cbind(mean(b0TSD90),mean(b1TSD90)))
  
  rMSE_90 = t(cbind(sqrt((sd(b0TP90))^2 + (mean(b0TP90) - betas0[1])^2), 
                    sqrt((sd(b1TP90))^2 + (mean(b1TP90) - betas0[2])^2)))
  
  bias_90 = bmean_90 - betas0
  
  
  CPT =  (UP > Betas0 & LO<Betas0)
  CP90 = apply(CPT, 2, mean)
  
  
  Betas0 = matrix(betas0, nsim, dim(TP)[2], byrow = T)
  b0UP = NULL
  b0LO = NULL
  b1UP = NULL
  b1LO = NULL
  
  b0UP1 = TP1[,1] + qnorm(0.975)*TSD1[,1]
  b0UP2 = TP2[,1] + qnorm(0.975)*TSD2[,1]
  b0UP3 = TP3[,1] + qnorm(0.975)*TSD3[,1]
  
  b1UP1 = TP1[,2] + qnorm(0.975)*TSD1[,2]
  b1UP2 = TP2[,2] + qnorm(0.975)*TSD2[,2]
  b1UP3 = TP3[,2] + qnorm(0.975)*TSD3[,2]
  
  for (i in 1:length(b0UP1)){
    largest <- mean(c(b0UP1[i], b0UP2[i], b0UP3[i]))
    b0UP = cbind(b0UP, largest)
  }
  
  for (i in 1:length(b1UP1)){
    largest <- mean(c(b1UP1[i], b1UP2[i], b1UP3[i]))
    b1UP = cbind(b1UP, largest)
  }
  
  UP = cbind(t(b0UP),t(b1UP))
  
  
  b0LO1 = TP1[,1] - qnorm(0.975)*TSD1[,1]
  b0LO2 = TP2[,1] - qnorm(0.975)*TSD2[,1]
  b0LO3 = TP3[,1] - qnorm(0.975)*TSD3[,1]
  
  b1LO1 = TP1[,2] - qnorm(0.975)*TSD1[,2]
  b1LO2 = TP2[,2] - qnorm(0.975)*TSD2[,2]
  b1LO3 = TP3[,2] - qnorm(0.975)*TSD3[,2]
  
  
  for (i in 1:length(b0LO1)){
    smallest <- mean(c(b0LO1[i], b0LO2[i], b0LO3[i]))
    b0LO = cbind(b0LO, smallest)
  }
  
  for (i in 1:length(b1LO1)){
    smallest <- mean(c(b1LO1[i], b1LO2[i], b1LO3[i]))
    b1LO = cbind(b1LO, smallest)
  }
  
  LO = cbind(t(b0LO),t(b1LO))
  
  
  b0TP95 = (b0LO+b0UP)/2
  b1TP95 = (b1LO+b1UP)/2
  bmean_95 = t(cbind(mean(b0TP95),mean(b1TP95)))
  
  b0TSD95 = as.vector(as.numeric((b0UP - mean(b0TP95))/qnorm(0.975)))
  b1TSD95 = as.vector(as.numeric((b1UP - mean(b1TP95))/qnorm(0.975)))
  
  SEm_95 = t(cbind(mean(b0TSD95),mean(b1TSD95)))
  
  
  rMSE_95 = t(cbind(sqrt((sd(b0TP95))^2 + (mean(b0TP95) - betas0[1])^2), 
                    sqrt((sd(b1TP95))^2 + (mean(b1TP95) - betas0[2])^2)))
  
  bias_95 = bmean_95 - betas0
  
  CPT =  (UP > Betas0 & LO<Betas0)
  CP95 = apply(CPT, 2, mean)
  
  
  Tab = cbind(betas0 = c(betas0), bmean = c(bmean_95), rMSE = c(rMSE_95), SEm = c(SEm_95), CP90 = c(CP90), CP95 = c(CP95))
  return(Tab)
}


w1_address = "C:/Users/qihan/Desktop/SpE_full"
setwd(w1_address)


load("n_500_100_Ver1_FULL.Rdata")
out = resTrf(res)


betas0 = c(2,4) # ver2
Tab1 = outcover(out$Beta1, out$Var1^.5, betas0)
Tab2 = outcover(out$Beta2, out$Var2^.5, betas0)
Tab3 = outcover(out$Beta3, out$Var3^.5, betas0)
Tab4 = outcover(out$Beta4, out$Var4^.5, betas0)
Tab5 = outcover_bind8(out$Beta5, out$Beta6, out$Beta2, out$Var5^.5, out$Var6^.5, out$Var2^.5, betas0)



print(Tab1)
print(Tab2)
print(Tab3)
print(Tab4)
print(Tab5)


Taball = rbind(Tab1, Tab2, Tab3, Tab4)
xtable(Taball, digits = 3)


Taball = data.frame(rbind(Tab1, Tab2, Tab3, Tab4, Tab5)) #, row.names = c("OLS1","OLS2","Updated1","Updated2","UK1","UK2","Corr1","Corr2"))
xt = xtable(Taball, digits = 2)
