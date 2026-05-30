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
out1 = MLE.fit(x1s, r1s, D_gas, cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)

xi_est =     out1$theta     #. res1[[1]] 
alpha_est  = out1$beta     #res1[[4]] 
Sigma1 = xi_est[2]*exp(-D_gas/xi_est[1]) 
diag(Sigma1) = xi_est[2] 
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
k = 10
set.seed(2024)
new_data = new_data[sample(1:nrow(new_data)),]
index_folds <- cut(seq(1, nrow(new_data)), breaks = k, labels = FALSE)
RMSE_OLS_data = c()
MAD_OLS_data = c()
RMSE_KR_data = c()
MAD_KR_data = c()
RMSE_KR2_data = c()
MAD_KR2_data = c()
RMSE_RBE_data = c()
MAD_RBE_data = c()
for(i in 1:k){
  idx <- which(index_folds==i, arr.ind=TRUE)
  test <- new_data[idx, ]
  train <- new_data[-idx, ]
  
  W = Wtrain = as.matrix(cbind(1, select(train, c("w", "proj_dist"))))
  Ytrain = y = as.matrix(select(train, c("Z_N2S2Ha")))
  # ---------------------------------------------------------------------------- #
  # OLS:
  beta_ols = solve(t(Wtrain)%*%Wtrain)%*%t(Wtrain)%*%Ytrain
  d_ols = cbind(1, test$w, test$proj_dist)
  pred_ols = as.matrix(d_ols)%*%beta_ols
  RMSE_ols <- sqrt(mean((pred_ols - test$Z_N2S2Ha)^2))
  MAD_ols <- mean(abs(pred_ols - test$Z_N2S2Ha))
  
  RMSE_OLS_data[i] = RMSE_ols
  MAD_OLS_data[i] = MAD_ols
  
  # KR: 
  theta.ini = c(0.1)
  lo.bound = c(1e-5)
  up.bound = c(2)
  cov_ep = "Exp"
  
  kr = MLE.fit(y, W, D_data[-idx,-idx], cov_ep, theta.ini, nug = 0, "LB", lo.bound, up.bound)
  theta_kr = kr[[1]] 
  beta_kr  = kr[[4]] 
  Sigma_kr = theta_kr[2]*exp(-D_data[-idx,-idx]/theta_kr[1]) 
  d_kr = cbind(1, test$w, test$proj_dist)
  cmat_kr = theta_kr[2]*exp(-D_data[idx,-idx]/theta_kr[1]) 
  pl_kr = as.matrix(d_kr)%*%beta_kr      
  pe_kr = cmat_kr%*% solve(Sigma_kr)%*%(Ytrain - Wtrain%*%beta_kr) 
  pred_kr = pl_kr + pe_kr
  RMSE_kr <- sqrt(mean((pred_kr - test$Z_N2S2Ha)^2))
  MAD_kr <- mean(abs(pred_kr - test$Z_N2S2Ha))
  RMSE_KR_data[i] = RMSE_kr
  MAD_KR_data[i] = MAD_kr
  
  
  kr2 = MLE.fit(y, W[,-2], D_data[-idx,-idx], cov_ep, theta.ini, nug = 0, "LB", lo.bound, up.bound)
  theta_kr2 = kr2[[1]] 
  beta_kr2  = kr2[[4]] 
  Sigma_kr2 = theta_kr2[2]*exp(-D_data[-idx,-idx]/theta_kr2[1]) 
  d_kr2 = cbind(1, test$proj_dist)
  cmat_kr2 = theta_kr2[2]*exp(-D_data[idx,-idx]/theta_kr2[1]) 
  pl_kr2 = as.matrix(d_kr2)%*%beta_kr2    
  pe_kr2 = cmat_kr2%*% solve(Sigma_kr2)%*%(Ytrain - W[,-2]%*%beta_kr2) 
  pred_kr2 = pl_kr2 + pe_kr2
  RMSE_kr2 <- sqrt(mean((pred_kr2 - test$Z_N2S2Ha)^2))
  MAD_kr2 <- mean(abs(pred_kr2 - test$Z_N2S2Ha))
  RMSE_KR2_data[i] = RMSE_kr2
  MAD_KR2_data[i] = MAD_kr2
  
  
  # RBEGLS
  theta_ep.ini = c(0.5, 0.01)
  lo.bounde = c(1e-5, 1e-10)
  up.bounde = c(2,2)
  nug = 0

  RBE_res = RBEGLS_loop2(50, 0.001, 0.001, beta_ols, theta_ep.ini, 
                         CovU[-idx, -idx], Ytrain, Wtrain, D_data[-idx,-idx], 
                         cov_ep, lo.bounde, up.bounde)
  
  beta_RBEGLS = RBE_res$beta_est
  theta_RBEGLS  = RBE_res$theta_est
  CorE_GLS = cor.mat(D_data[-idx,-idx], theta_RBEGLS[1], cov_ep, nug = 0)
  CovE_GLS = theta_RBEGLS[2]*CorE_GLS
  CovUB = CovU[-idx, -idx]*beta_RBEGLS[2]^2 
  CovUBE_GLS = CovUB + CovE_GLS
  
  CorE_GLS2 = cor.mat(D_data[idx,-idx], theta_RBEGLS[1], cov_ep, nug = 0)
  CovE_GLS2 = theta_RBEGLS[2]*CorE_GLS2
  CovUB2 = CovU[idx, -idx]*beta_RBEGLS[2]^2 
  CovUBE_GLS2 = CovUB2 + CovE_GLS2
  
  d_RBE = cbind(1, test$w, test$proj_dist)
  pl_RBE = as.matrix(d_RBE)%*%beta_RBEGLS
  pe_RBE = CovUBE_GLS2%*% solve(CovUBE_GLS)%*%(Ytrain - Wtrain%*%beta_RBEGLS) 
  pred_RBE = pl_RBE + pe_RBE
  
  RMSE_RBE <- sqrt(mean((pred_RBE - test$Z_N2S2Ha)^2))
  MAD_RBE <- mean(abs(pred_RBE - test$Z_N2S2Ha))
  RMSE_RBE_data[i] = RMSE_RBE
  MAD_RBE_data[i] = MAD_RBE
}


