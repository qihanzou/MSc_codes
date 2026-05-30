# # # # test use:
# rm(list = ls(all = TRUE)) 
# nall = 800
# idx = 1:400
# beta1 = 4
# MC_iter = 100
# R_address = "C:/Users/qihan/Desktop/QH_nonlinear"


setwd(R_address)
library(mvtnorm)  
library(moments)
library(revss)
library(asbio)
library(MASS)
library(r2spss)
source("SimuData.R")
source("BFuns.R")
source("BFuns_RBEGLS.R")
source("SDCal_modified.R")
source("RBEGLS_functions_final.R")
library(sigmoid)


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


N_sim = 1000
# define the non-linear transformation here:
nl_trans = function(input){
  output = sigmoid(input) # change this line
}



run.sim <- function(){
  
  
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
  #y = beta0 + beta1*x1 + ep_error     
  z1 = nl_trans(x1) 
  y = beta0 + beta1*z1 + ep_error 
  #--------------------------------------------------------------------------------------------------------------
  # w
  out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
  alpha_est = out1$beta
  cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
  
  w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  
  
  # MCT, Monte Carlo Transformation
  Z_set = NULL
  for (i in 1:N_sim){
    w_i = t(w) + rmvnorm(1, sigma = CovU)
    z_i = nl_trans(w_i)
    Z_set = rbind(Z_set, z_i)
  }
  Cov_Z = matrix(0, length(idx), length(idx))
  for (j in 1:N_sim){
    Cov_Z = Cov_Z + (Z_set[j,] - colMeans(Z_set))%*%t(Z_set[j,] - colMeans(Z_set))
  }
  CovU_nl = 1/(N_sim-1)*Cov_Z
  z = colMeans(Z_set)

  
  #--------------------------------------------------------------------------------------------------------------
  # OLS
  Z = cbind(1, z)
  fols = lm(y~Z-1)
  beta_ols = fols$coef
  varbeta_ols = diag(vcov(fols))
  
  #--------------------------------------------------------------------------------------------------------------
  # KR
  kr = MLE.fit(y, Z, Dall[idx,idx], cov_ep, c(2.5), nug, "LB", lo.bound,up.bound)
  theta_kr = kr[[1]] 
  beta_kr  = kr[[4]] 
  varbeta_kr = diag(kr[[5]])
  
  #--------------------------------------------------------------------------------------------------------------
  # RBEGLS
  CovUB = CovU_nl*beta_ols[2]^2 
  
  out2 = MLE.fit_GLS(y, Z, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS = out2$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(Z)%*%solve(CovUBE_GLS)%*%Z)%*%t(Z)%*%solve(CovUBE_GLS)%*%y
  
  
  # loop:
  res1 = RBEGLS_loop_final(50, 0.001, 0.001, beta_RBEGLS, out2$theta, theta_ep.ini, CovU_nl, y, Z, Dall, cov_ep, lo.bounde,up.bounde)
  beta_RBEGLS = res1$beta_est
  theta_RBEGLS  = res1$theta_est
  CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
  CovE_GLS = theta_RBEGLS[2]*CorE_GLS
  CovUB = CovU_nl*beta_RBEGLS[2]^2
  CovUBE_GLS = CovUB + CovE_GLS
  CovUBE_GLS_inv = solve(CovUBE_GLS)
  A1 = solve(t(Z)%*%CovUBE_GLS_inv%*%Z)
  varbeta_rb = diag(A1)
  #--------------------------------------------------------------------------------------------------------------
  # MC under construction
  
  sdout = SDCal(Dall[-idx,-idx], out1$theta, cov_xi, nug = 0)
  EVar = EVarKR = EVarKRT = matrix(0, 2, 2)
  for (v in 1:MC_iter){
    
    xie = out1$theta + rmvnorm(1, sigma=sdout$covmat)
    while ((xie[1] < 0)|(xie[2] < 0)) {
      xie = out1$theta + rmvnorm(1, sigma=sdout$covmat)
    }
    cov_xialle = xie[2]*cor.mat(Dall, xie[1], cov_xi, nug = 0)

    tmp_inv = solve(cov_xialle[-idx,-idx])
    alphaboot = solve(t(r1s) %*% tmp_inv %*% r1s)%*%t(r1s) %*% tmp_inv %*% x1s
    wboot_e = r1%*%alphaboot + cov_xialle[idx, -idx]%*%tmp_inv%*%(x1s - r1s%*%alphaboot)
    wboot_e1 = nl_trans(wboot_e)
    Wboot_e = cbind(1, wboot_e1)
    
    A1_e = solve(t(Wboot_e)%*%CovUBE_GLS_inv%*%Wboot_e)
    EVar = EVar + A1_e                      
  }
  varbeta_mc = diag(EVar/MC_iter)
  

  #--------------------------------------------------------------------------------------------------------------
  list(alpha_est = as.matrix(out1$beta),
       xi_est = as.matrix(out1$theta),
       
       beta_ols = as.matrix(as.numeric(beta_ols)), 
       beta_RBEGLS = as.matrix(as.numeric(beta_RBEGLS)),  
       beta_kr = as.matrix(as.numeric(beta_kr)),
       
       theta_est = as.matrix(as.numeric(theta_RBEGLS)),
       theta_kr = as.matrix(as.numeric(theta_kr)),
       
       varbeta_ols = as.matrix(as.numeric(varbeta_ols)), 
       varbeta_kr = as.matrix(as.numeric(varbeta_kr)),
       varbeta_rb = as.matrix(as.numeric(varbeta_rb)),
       varbeta_mc = as.matrix(as.numeric(varbeta_mc))
  )
  
  
}


