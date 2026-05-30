# test use:
#boot_iter = 100
#nall = 400
#idx = 1:200
#beta1 = 4
# 
# 
#R_address = "C:/Users/qihan/Desktop/QH_w"

setwd(R_address)
library(mvtnorm)               
source("SimuData.R")
source("BFuns.R")
source("SDCal_modified.R")
source("SDCal_RMLE.R")
source("BFuns_TR.R")
source("BFuns_RBEGLS.R")
source("RBEGLS_functions.R")
source("rho_est.R")




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
beta2 = 1


nug = 0



run.sim <- function(){


  etaall = SimuData(xi0, rate, locall, cov_xi, nug, 1)
  epall = SimuData(ep_theta0, rate, locall, cov_ep, nug, 1)
  
  k=1
  
  distall = sqrt((locall$x - 0)^2 + (locall$y - 0)^2)
  
  r1all = as.matrix(cbind(1, distall)) 
  r1s = as.matrix(cbind(1, distall[-idx]))
  r1 = as.matrix(cbind(1,distall[idx]))
  
  x1all = r1all%*%alpha0 + etaall$Ymat[k,]  
  x1 = x1all[idx] 
  x1s = x1all[-idx] 
  N = length(x1) 
  M = length(x1s)
  
  #  simple model
  ep_error = epall$Ymat[k,idx]        # error of s1,...,sN
  x2 = distall[idx]
  y = beta0 + beta1*x1 + beta2*x2 + ep_error     # create y values for s1,...,sN
  
  
  # -------------------------------------------------- #
  # 1. OLS beta estimates
  # -------------------------------------------------- #
  out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
  alpha_est = out1$beta
  cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
  w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 
  W = cbind(1, w, x2)
  fols = lm(y~W-1)
  beta_ols = fols$coef
  varbeta_ols = diag(vcov(fols))

  # -------------------------------------------------- #
  # 2. kr:
  # -------------------------------------------------- #
  kr = MLE.fit(y, W, Dall[idx,idx], cov_ep, c(2.5), nug, "LB", lo.bound,up.bound)
  theta_kr = kr[[1]] 
  beta_kr  = kr[[4]] 
  beta_var_kr = diag(kr[[5]])
  covKR = theta_kr[2]*cor.mat(Dall[idx,idx], kr$eta, cov_ep, nug = 0)
  covKR.inv = solve(covKR)

  # -------------------------------------------------- #
  # 3. RBEGLS method
  # -------------------------------------------------- #
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  CovUB = CovU*beta_ols[2]^2 
  
  
  out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS = out2$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y
  
  # loop start here:
  res1 = RBEGLS_loop(50, 0.001, 0.001, beta_RBEGLS, out2$theta, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
  
  beta_RBEGLS = res1$beta_est
  theta_RBEGLS  = res1$theta_est
  CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
  CovE_GLS = theta_RBEGLS[2]*CorE_GLS
  CovUB = CovU*beta_RBEGLS[2]^2 
  CovUBE_GLS = CovUB + CovE_GLS
  CovUBE_GLS_inv = solve(CovUBE_GLS)
  
  
  
  

 boot_iter = 100
 Bbar = matrix(0, 3, boot_iter)
 
  for (v in 1:boot_iter){
    wboot_e = w + t(rmvnorm(1, sigma = CovU))
    Wboot_e = cbind(1, wboot_e, x2)
    beta_j = solve(t(Wboot_e)%*%Wboot_e)%*%t(Wboot_e)%*%y
    
    Bbar[, v] = beta_j                    
  }

  varbet = apply(Bbar, 1, var)
  b0 = Bbar[1,]
  b1 = Bbar[2,]
  b2 = Bbar[3,]
    
  nB0 = boot_iter*beta_RBEGLS[1]
  D0 = nB0 - sum(b0)
  D0n = D0/boot_iter
  nB1 = boot_iter*beta_RBEGLS[2]
  D1 = nB1 - sum(b1)
  D1n = D1/boot_iter
  nB2 = boot_iter*beta_RBEGLS[3]
  D2 = nB2 - sum(b2)
  D2n = D2/boot_iter
  
  mb0 = b0 + D0n
  mb1 = b1 + D1n
  mb2 = b2 + D2n
  
  rho_b0 = est_rho(5000, beta_RBEGLS[1], as.matrix(mb0), varbet[1], boot_iter)
  rho_b1 = est_rho(5000, beta_RBEGLS[2], as.matrix(mb1), varbet[2], boot_iter)
  rho_b2 = est_rho(5000, beta_RBEGLS[3], as.matrix(mb2), varbet[3], boot_iter)

  A1 = solve(t(W)%*%CovUBE_GLS_inv%*%W)
  varbeta_rb = diag(A1)
  varbeta_rbc = varbeta_rb + varbet
  rho_varbet = c(rho_b0*varbet[1], rho_b1*varbet[2], rho_b2*varbet[3])
  varbeta_rbc_rho = varbeta_rb + rho_varbet
  varbeta_kr = beta_var_kr 
  
  list(xi_est = out1$theta, alpha_est = out1$beta,
       beta_ols = beta_ols, beta_kr = beta_kr, beta_RBEGLS = beta_RBEGLS,  
       theta_est = theta_RBEGLS,
       varbeta_ols = varbeta_ols, varbeta_kr = varbeta_kr, 
       varbeta_rb = varbeta_rb,
       varbeta_rbc = varbeta_rbc, varbeta_rbc_rho = varbeta_rbc_rho)
       

}

