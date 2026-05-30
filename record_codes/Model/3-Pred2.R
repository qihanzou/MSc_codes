rm(list = ls(all = TRUE)) 
library(pracma)
library(mvtnorm) 
library(geoR)


setwd("C:/Users/qihan/Desktop/Model")
source("BFuns.R")
source("SDCal_modified.R")
source("SDCal_RMLE.R")
source("BFuns_TR.R")
source("BFuns_RBEGLS.R")
source("RBEGLS_functions.R")
source("SimuData.R")
source("R_gas_functions.R")


load("1-Data.Rdata")

plot(Local_gas$X, Local_gas$Y, cex = 0.3, pch = 16)
points(Local_data$X, Local_data$Y, cex = 0.3, col = "red", pch = 16)

dim(Local_gas)
dim(Local_data)




D_gas = deprojected_distances(Local_gas$RA, Local_gas$DEC, Local_gas$RA, Local_gas$DEC)
D_data = deprojected_distances(Local_data$RA, Local_data$DEC, Local_data$RA, Local_data$DEC)
D_data_to_gas = deprojected_distances(Local_data$RA, Local_data$DEC, Local_gas$RA, Local_gas$DEC)
x1s = as.matrix(select(Local_gas, c("CO21_mgd")))
r1s = as.matrix(cbind(1, select(Local_gas, c("proj_dist"))))


# ---------------------------------------------------------------------------- #
# Applying Universal Kriging
# ---------------------------------------------------------------------------- #
xi.ini = c(0.5)
lo.bound = c(1e-5)
up.bound = c(4)
nug = 0
cov_xi = "Exp"
out = MLE.fit(x1s, r1s, D_gas, cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
xi_est =     out$theta    
alpha_est  = out$theta
Sigma1 = xi_est[2]*exp(-D_gas/xi_est[1]) 
r1 = cbind(1, Local_data$proj_dist)
cmat1 = xi_est[2]*exp(-D_data_to_gas/xi_est[1]) 
pl1 = as.matrix(r1)%*%alpha_est 
pe1 = cmat1%*% solve(Sigma1)%*%(x1s - r1s%*%alpha_est)
w = pl1 + pe1
# ---------------------------------------------------------------------------- #
# Create data with the predicted w
# ---------------------------------------------------------------------------- #
new_data = Local_data
new_data$w = w

inv_Sigma = solve(Sigma1)
Cz0 = xi_est[2]*exp(-D_data/xi_est[1]) 
pc = r1 - cmat1%*%inv_Sigma%*%r1s
CovU = Cz0 - cmat1%*%inv_Sigma%*%t(cmat1) + (pc)%*%solve(t(r1s)%*%inv_Sigma%*%r1s)%*%t(pc)

# ---------------------------------------------------------------------------- #
# k-folds Cross-Validation
# ---------------------------------------------------------------------------- #
x.comb = as.matrix(cbind(1, new_data$w, new_data$proj_dist))
y.comb = as.matrix(new_data$Z_N2S2Ha)

cv_k = 5
n.comb = length(y.comb)
if (n.comb%%cv_k == 0){
  split.ind = rep(1:cv_k, each = floor(n.comb/cv.k))
}else{
  split.ind = c(rep(1:cv_k, each = floor(n.comb/cv_k)),
                1:(n.comb%%cv_k))
}
set.seed(1174116)
split.ind = sample(split.ind)



RMSE_OLS_data = c()
MAD_OLS_data = c()
RMSE_OLS_data2 = c()
MAD_OLS_data2 = c()
RMSE_KR_data = c()
MAD_KR_data = c()
RMSE_KR2_data = c()
MAD_KR2_data = c()
RMSE_RBE_data = c()
MAD_RBE_data = c()
RMSE_RBE_data2 = c()
MAD_RBE_data2 = c()
for(j in 1:cv_k){
  Wtest = x.comb[split.ind == j,]
  Ytest = y.comb[split.ind == j]
  W = Wtrain = x.comb[split.ind != j,]
  y = Ytrain = y.comb[split.ind != j]
  
  train_TF = split.ind != j
  test_TF = split.ind == j
  
  idx_train = which(train_TF == TRUE)
  idx_test = which(test_TF == TRUE)
  # ---------------------------------------------------------------------------- #
  # OLS:
  beta_ols = solve(t(Wtrain[,1:2])%*%Wtrain[,1:2])%*%t(Wtrain[,1:2])%*%Ytrain
  pred_ols = Wtest[,1:2]%*%beta_ols
  RMSE_ols <- sqrt(mean((pred_ols - Ytest)^2))
  MAD_ols <- mean(abs(pred_ols - Ytest))
  RMSE_OLS_data[j] = RMSE_ols
  MAD_OLS_data[j] = MAD_ols
  
  beta_ols2 = solve(t(Wtrain[,-2])%*%Wtrain[,-2])%*%t(Wtrain[,-2])%*%Ytrain
  pred_ols2 = Wtest[,-2]%*%beta_ols2
  RMSE_ols2 <- sqrt(mean((pred_ols2 - Ytest)^2))
  MAD_ols2 <- mean(abs(pred_ols2 - Ytest))
  RMSE_OLS_data2[j] = RMSE_ols2
  MAD_OLS_data2[j] = MAD_ols2
  
  # KR: 
  theta.ini = c(0.1)
  lo.bound = c(1e-10)
  up.bound = c(2)
  cov_ep = "Exp"
  
  kr = MLE.fit(Ytrain, Wtrain[,1:2], D_data[idx_train,idx_train], cov_ep, theta.ini, nug = 0, "LB", lo.bound, up.bound)
  theta_kr = kr$theta
  beta_kr  = kr$beta
  Sigma_kr = theta_kr[2]*exp(-D_data[idx_train,idx_train]/theta_kr[1]) 
  cmat_kr = theta_kr[2]*exp(-D_data[idx_test,idx_train]/theta_kr[1]) 
  pl_kr = as.matrix(Wtest[,1:2])%*%beta_kr
  pe_kr = cmat_kr%*% solve(Sigma_kr)%*%(Ytrain - Wtrain[,1:2]%*%beta_kr) 
  pred_kr = pl_kr + pe_kr
  RMSE_kr <- sqrt(mean((pred_kr - Ytest)^2))
  MAD_kr <- mean(abs(pred_kr - Ytest))
  RMSE_KR_data[j] = RMSE_kr
  MAD_KR_data[j] = MAD_kr
  
  
  kr2 = MLE.fit(Ytrain, Wtrain[,-2], D_data[idx_train,idx_train], cov_ep, theta.ini, nug = 0, "LB", lo.bound, up.bound)
  theta_kr2 = kr2$theta 
  beta_kr2  = kr2$beta
  Sigma_kr2 = theta_kr2[2]*exp(-D_data[idx_train,idx_train]/theta_kr2[1]) 
  cmat_kr2 = theta_kr2[2]*exp(-D_data[idx_test,idx_train]/theta_kr2[1]) 
  pl_kr2 = as.matrix(Wtest[,-2])%*%beta_kr2    
  pe_kr2 = cmat_kr2%*% solve(Sigma_kr2)%*%(Ytrain - W[,-2]%*%beta_kr2) 
  pred_kr2 = pl_kr2 + pe_kr2
  RMSE_kr2 <- sqrt(mean((pred_kr2 - Ytest)^2))
  MAD_kr2 <- mean(abs(pred_kr2 - Ytest))
  RMSE_KR2_data[j] = RMSE_kr2
  MAD_KR2_data[j] = MAD_kr2
  
  
  # RBEGLS
  theta_ep.ini = c(0.5, 0.001)
  lo.bounde = c(1e-5, 1e-10)
  up.bounde = c(2,2)
  nug = 0
  cov_ep1 = "Exp"
  
  RBE_res = RBEGLS_loop2(50, 0.001, 0.001, beta_ols, theta_ep.ini, 
                         CovU[idx_train, idx_train], Ytrain, Wtrain[,1:2], D_data[idx_train,idx_train], 
                         cov_ep1, lo.bounde, up.bounde)
  
  beta_RBEGLS = RBE_res$beta_est
  theta_RBEGLS  = RBE_res$theta_est
  CorE_GLS = cor.mat(D_data[idx_train,idx_train], theta_RBEGLS[1], cov_ep1, nug = 0)
  CovE_GLS = theta_RBEGLS[2]*CorE_GLS
  CovUB = CovU[idx_train, idx_train]*beta_RBEGLS[2]^2 
  CovUBE_GLS = CovUB + CovE_GLS
  
  CovE_GLS2 = theta_RBEGLS[2]*exp(-D_data[idx_test,idx_train]/theta_RBEGLS[1]) 
  CovUB2 = CovU[idx_test, idx_train]*beta_RBEGLS[2]^2 
  cc1 = CovUB2 + CovE_GLS2
  pl_RBE = as.matrix(Wtest[,1:2])%*%beta_RBEGLS
  pe_RBE = cc1%*% solve(CovUBE_GLS)%*%(Ytrain - Wtrain[,1:2]%*%beta_RBEGLS) 
  pred_RBE = pl_RBE + pe_RBE
  
  RMSE_RBE <- sqrt(mean((pred_RBE - Ytest)^2))
  MAD_RBE <- mean(abs(pred_RBE - Ytest))
  RMSE_RBE_data[j] = RMSE_RBE
  MAD_RBE_data[j] = MAD_RBE
  
  
  
  RBE_res2 = RBEGLS_loop2(50, 0.001, 0.001, beta_ols, theta_ep.ini, 
                         CovU[idx_train, idx_train], Ytrain, Wtrain[,-2], D_data[idx_train,idx_train], 
                         cov_ep1, lo.bounde, up.bounde)
  
  beta_RBEGLS2 = RBE_res2$beta_est
  theta_RBEGLS2  = RBE_res2$theta_est
  CorE_GLS2 = cor.mat(D_data[idx_train,idx_train], theta_RBEGLS2[1], cov_ep1, nug = 0)
  CovE_GLS2 = theta_RBEGLS2[2]*CorE_GLS2
  CovUB2 = CovU[idx_train, idx_train]*beta_RBEGLS2[2]^2 
  CovUBE_GLS2 = CovUB2 + CovE_GLS2
  
  CovE_GLS22 = theta_RBEGLS2[2]*exp(-D_data[idx_test,idx_train]/theta_RBEGLS2[1]) 
  CovUB22 = CovU[idx_test, idx_train]*beta_RBEGLS[2]^2 
  cc12 = CovUB22 + CovE_GLS22
  pl_RBE2 = as.matrix(Wtest[,-2])%*%beta_RBEGLS2
  pe_RBE2 = cc12%*% solve(CovUBE_GLS2)%*%(Ytrain - Wtrain[,-2]%*%beta_RBEGLS2) 
  pred_RBE2 = pl_RBE2 + pe_RBE2
  
  RMSE_RBE2 <- sqrt(mean((pred_RBE2 - Ytest)^2))
  MAD_RBE2 <- mean(abs(pred_RBE2 - Ytest))
  RMSE_RBE_data2[j] = RMSE_RBE2
  MAD_RBE_data2[j] = MAD_RBE2
  
  print(RMSE_ols)
  print(RMSE_ols2)
  print(RMSE_kr)
  print(RMSE_kr2)
  print(RMSE_RBE)
  print(RMSE_RBE2)
}

