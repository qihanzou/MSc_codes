# # # # test use:
rm(list = ls(all = TRUE)) 
nall = 200
idx = 1:100
beta1 = 4
R_address = "C:/Users/qihan/Desktop/Rebound"


setwd(R_address)
library(mvtnorm)  
library(moments)
library(revss)
library(asbio)
library(MASS)
library(r2spss)
source("SimuData.R")
source("BFuns.R")
source("BFuns_IREML.R")
source("SDCal_modified.R")
source("RBEGLS_functions_final.R")


# -------------------------------------------------------------------------- #
locRdata=strwrap(paste("locs/loc", nall, ".Rdata",sep=""))
load(locRdata)

Dall = matrix(0, nrow = nall, ncol = nall)            
for(i in 1:nall) {
  Dall[i, ] <- sqrt((locall$x[i] - locall$x)^2 + (locall$y[i] - locall$y)^2)
}
#------------------------------------------------------------------------------#

alpha0 = c(1,0.5)
xi0 = c(3, 4)  # initial values of (phi, sigma^2)
cov_xi = "Mat32"
xi.ini = c(4)
up.bound = c(20)
lo.bound = c(0.01)

cov_ep = "Exp"
ep_theta0 = c(2,1) 
theta_ep.ini = c(3.5,1.2)
up.bounde = c(20, 20)
lo.bounde = c(0.01, 0.01)


beta0 = 2

nug = 0


#run.sim <- function(){
  
  
  etaall = SimuData(xi0, rate, locall, cov_xi, nug, 1)
  epall = SimuData(ep_theta0, rate, locall, cov_ep, nug, 1)
  
  k=1
  
  distall = sqrt((locall$x - 0)^2 + (locall$y - 0)^2)
  
  # -------------------------------------------------------------------------- #
  r1all = as.matrix(cbind(1, distall)) 
  r1s = as.matrix(cbind(1, distall[-idx]))
  r1 = as.matrix(cbind(1,distall[idx]))
  
  x1all = r1all%*%alpha0 + etaall$Ymat[k,]  
  x1 = x1all[idx] 
  x1s = x1all[-idx] 
  N = length(x1) 
  M = length(x1s)
  
  
  #  simple model
  ep_error = epall$Ymat[k,idx]        
  y = beta0 + beta1*x1 + ep_error     
  
  #--------------------------------------------------------------------------------------------------------------
  # Initial iteration:
  out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
  alpha_est = out1$beta
  cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
  w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 
  W = cbind(1, w)
  fols = lm(y~W-1)
  beta_ols = fols$coef
  
  # REML
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  CovUB = CovU*beta_ols[2]^2 
  out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS1 = out2$theta[2]*CorE_GLS
  CovUBE_GLS1 = CovUB + CovE_GLS1
  CovUBE_GLS1_inv = solve(CovUBE_GLS1)
  beta_est1 = solve(t(W)%*%CovUBE_GLS1_inv%*%W)%*%t(W)%*%CovUBE_GLS1_inv%*%y
  
  res1 = RBEGLS_loop_final(50, 0.001, 0.001, beta_est1, out2$theta, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
  beta_RBEGLS = res1$beta_est
  theta_RBEGLS  = res1$theta_est
  CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
  CovE_GLS = theta_RBEGLS[2]*CorE_GLS
  CovUB = CovU*beta_RBEGLS[2]^2
  CovUBE_GLS = CovUB + CovE_GLS
  CovUBE_GLS_inv = solve(CovUBE_GLS)
  A1 = solve(t(W)%*%CovUBE_GLS_inv%*%W)
  varbeta_rb = diag(A1)
  
  
  LOPU_xi = LOUP(Dall[-idx,-idx], out1$theta, cov_xi, nug)
  LOPU_theta = LOUP(Dall[idx,idx], theta_RBEGLS, cov_ep, nug)
  
  out1_j = MLE.fit_IREML(y, W, Dall, beta_RBEGLS[2], CovE_GLS, out1$theta, nug, r1, r1s, idx, "LB", LOPU_xi$lob, LOPU_xi$upb)
  out1_j
  out1_j$theta
  out1$theta
  theta_RBEGLS
  
  (out1_j$theta + out1$theta)/theta_RBEGLS
  
  
  cov_xiall_j = out1_j$theta[2]*cor.mat(Dall, out1_j$theta[1], cov_xi, nug = 0)
  pc_j = r1 - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% r1s
  CovU_j = cov_xiall_j[idx, idx] - cov_xiall_j[idx, -idx]%*%solve(cov_xiall_j[-idx,-idx]) %*% t(cov_xiall_j[idx, -idx]) + (pc_j)%*%solve(t(r1s)%*%solve(cov_xiall_j[-idx,-idx])%*%r1s)%*%t(pc_j)
  CovUB_j = CovU_j*beta_RBEGLS[2]^2 
  
  res2 = RBEGLS_loop_final(50, 0.001, 0.001, beta_RBEGLS, out1_j$theta, out1_j$theta, CovU_j, y, W, Dall, cov_ep, lo.bounde,up.bounde)
  
  res1$nll
  out1_j$nll
  res2$nll
  
  
  res2 = IREML_loop_final(50, 0.001, 0.001, 0.001, 0.001, beta_est1, out2$theta, alpha_est, out1$theta, CovE_GLS1, y, W, r1, r1s, x1s, nug, Dall, idx, cov_ep, cov_xi)
  beta_IREML = res2$beta_set[,res2$iter]
  theta_IREML = res2$theta_set[,res2$iter]
  alpha_IREML = res2$alpha_set[,res2$iter]
  xi_IREML = res2$xi_set[,res2$iter]
  CovUBE_IREML = res2$CovUBE_IREML
  CovUBE_IREML_inv = solve(CovUBE_IREML)
  W_IREML = res2$W_set[[res2$iter]]
  A2 = solve(t(W_IREML)%*%CovUBE_IREML_inv%*%W_IREML)
  varbeta_IREML = diag(A2)
  
  A3 = solve(t(W)%*%CovUBE_IREML_inv%*%W)
  varbeta_IREML2 = diag(A3)
  
  
  EVar  = matrix(0, 2, 2)
  for (k in 1:res2$iter){
     W_IREML_k = res2$W_set[[k]]
     A1_e = solve(t(W_IREML_k)%*%CovUBE_IREML_inv%*%W_IREML_k)
     EVar = EVar + A1_e   
   }
   varbeta_w = diag(EVar/res2$iter)
  
  
  RMSE_IREML = sqrt(mean((x1 - W_IREML)^2))
  RMSE_RBEGLS = sqrt(mean((x1 - W)^2))
  RMSE_IREML_RBEGLS = sqrt(mean((W_IREML - W)^2))
  
  MAD_IREML = mean(abs(x1 - W_IREML))
  MAD_RBEGLS = mean(abs(x1 - W))
  MAD_IREML_RBEGLS = mean(abs(W_IREML - W))


  #--------------------------------------------------------------------------------------------------------------
  list(alpha_est = out1$beta,
       alpha_IREML = as.matrix(alpha_IREML),
       xi_est = as.matrix(out1$theta),
       xi_IREML = as.matrix(xi_IREML),
       
       beta_RBEGLS = beta_RBEGLS,  
       beta_IREML = as.matrix(beta_IREML),
       
       theta_RBEGLS = as.matrix(theta_RBEGLS),
       theta_IREML = as.matrix(theta_IREML),
       
       
       RMSE_IREML = RMSE_IREML,
       RMSE_RBEGLS = RMSE_RBEGLS,
       RMSE_IREML_RBEGLS = RMSE_IREML_RBEGLS,
       MAD_IREML = MAD_IREML,
       MAD_RBEGLS = MAD_RBEGLS,
       MAD_IREML_RBEGLS = MAD_IREML_RBEGLS,
       
       varbeta_RBEGLS = as.matrix(varbeta_rb),
       varbeta_IREML = as.matrix(varbeta_IREML),
       varbeta_IREML2 = as.matrix(varbeta_IREML2),
       varbeta_w = as.matrix(varbeta_w)
  )
  
  
#}


