#test use:
boot_iter = 100

R_address = "C:/Users/qihan/Desktop/SpE_qihan"


setwd(R_address)
library(mvtnorm)               
source("SimuData.R")
source("BFuns.R")
source("SDCal_modified.R")
source("SDCal_RMLE.R")
source("BFuns_RBEGLS.R")
source("RBEGLS_functions.R")
source("rho_est.R")

nall = 200
idx = 1:100

locRdata=strwrap(paste("locs/loc", nall, ".Rdata",sep=""))
load(locRdata)

Dall = matrix(0, nrow = nall, ncol = nall)            
for(i in 1:nall) {
  Dall[i, ] <- sqrt((locall$x[i] - locall$x)^2 + (locall$y[i] - locall$y)^2)
}

alpha0 = c(1,0.5)

xi0 = c(3, 4) 
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
beta1 = 4

rate = 0.5
nug = 0

#run.sim <- function(){
  xall = SimuData(xi0, rate, locall, cov_xi, nug, 1)
  epall = SimuData(ep_theta0, rate, locall, cov_ep, nug, 1)
  
  k=1
  
  distall = matrix(0, nrow = nall, ncol = 1)            
  for(j in 1:nall) {
    distall[j] <- sqrt((locall$x[j] - 0)^2 + (locall$y[j] - 0)^2)
  }
  
  designx1all = as.matrix(cbind(1, distall)) 
  designx1s = as.matrix(cbind(1, distall[-idx]))
  designx1 = as.matrix(cbind(1,distall[idx]))
  
  x1all = designx1all%*%alpha0 + xall$Ymat[k,]  
  x1 = x1all[idx] 
  x1s = x1all[-idx] 
  N = length(x1) 
  M = length(x1s)
  
  #  simple model
  ep_error = epall$Ymat[k,idx] # error of s1,...,sN
  y = beta0 + beta1*x1 + ep_error # create y values for s1,...,sN
  #--------------------------------------------------------------------------------------------------------------#
  
  ## OLS beta estimates
  
  # step1: estimate alpha and theta_n: from x* to w (estimator of x)
  out1 = MLE.fit(x1s, designx1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
  betaout1 = out1$beta
  cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
  
  # step2: Derive W from Universal Kriging and use it in place of X:
  w = designx1%*%betaout1 + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - designx1s%*%betaout1) 
  W = cbind(1, w)
  beta_est = solve(t(W)%*%W)%*%t(W)%*%y 
  
  pc = designx1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% designx1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(designx1s)%*%solve(cov_xiall[-idx,-idx])%*%designx1s)%*%t(pc)
  CovUB = CovU*beta_est[2]^2 
  
  out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS = out2$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y
  res1 = RBEGLS_loop(50, 0.001, 0.001, 0.001, 0.001, beta_RBEGLS, out2$theta, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
  
  beta_RBEGLS = cbind(res1[1], res1[2])
  theta_RBEGLS  = cbind(res1[3], res1[4])
  CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
  CovE_GLS = theta_RBEGLS[2]*CorE_GLS
  CovUB = CovU*beta_RBEGLS[2]^2 
  CovUBE_GLS = CovUB + CovE_GLS
  
  
  b0 = matrix(0, boot_iter)
  b1 = matrix(0, boot_iter)
  var_boot = c()
  inv_Cov_UBE_GLS = solve(CovUBE_GLS)
  for (v in 1:boot_iter){
    w_i = w + t(rmvnorm(1, sigma = CovU))
    W_i = cbind(1, w_i)
    beta_i = solve(t(W_i)%*%inv_Cov_UBE_GLS%*%W_i)%*%t(W_i)%*%inv_Cov_UBE_GLS%*%y
    b0[v] = beta_i[1]
    b1[v] = beta_i[2]
    
    A_i = solve(t(W_i)%*%inv_Cov_UBE_GLS%*%W_i)
    varbeta_i = diag(A_i)
    var_boot = rbind(var_boot, varbeta_i)
  }
  varb0 = sd(b0)^2
  varb1 = sd(b1)^2
  cvar = cbind(varb0, varb1)
  SEm = apply(var_boot^.5, 2, mean)
  A1 = solve(t(W)%*%inv_Cov_UBE_GLS%*%W)
  varbeta_cl0 = diag(A1)
  varbeta_cl1 = SEm^2
  
  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  nB0 = boot_iter*beta_RBEGLS[1]
  D0 = nB0 - sum(b0)
  D0n = D0/boot_iter
  nB1 = boot_iter*beta_RBEGLS[2]
  D1 = nB1 - sum(b1)
  D1n = D1/boot_iter
  mb0 = b0 + D0n
  mb1 = b1 + D1n
  TP0 = mb0
  TP1 = mb1
  
  nV0 = boot_iter*varbeta_cl0[1]
  DV0 = nV0 - sum(var_boot[,1])
  DV0n = DV0/boot_iter
  nV1 = boot_iter*varbeta_cl0[2]
  DV1 = nV1 - sum(var_boot[,2])
  DV1n = DV1/boot_iter
  
  mv0 = as.matrix(as.numeric(var_boot[,1] + DV0n))
  mv1 = as.matrix(as.numeric(var_boot[,2] + DV1n))
  
  SEm2s = apply(cbind(mv0,mv1)^.5, 2, mean)^2
  

  rho_b0 = est_rho(5000, beta_RBEGLS[1], mb0, varb0, boot_iter)
  rho_b1 = est_rho(5000, beta_RBEGLS[2], mb1, varb1, boot_iter)
  
  rho_v0 = est_rho(5000, varbeta_cl0[1], mv0, sd(mv0)^2, boot_iter)
  rho_v1 = est_rho(5000, varbeta_cl0[2], mv1, sd(mv1)^2, boot_iter)
  
  cvar_rho = cbind(rho_b0*varb0 + (varb0/boot_iter)*(1 - rho_b0), rho_b1*varb1 + (varb1/boot_iter)*(1 - rho_b1))
  
  
  vcl0_b0 = varbeta_cl0[1]
  vcl0_b1 = varbeta_cl0[2]
  varbeta_cl0_mod = cbind(rho_v0*vcl0_b0 + (vcl0_b0/boot_iter)*(1 - rho_v0), rho_v1*vcl0_b1 + (vcl0_b1/boot_iter)*(1 - rho_v1))
  
  
  vcl1_b0 = varbeta_cl1[1]
  vcl1_b1 = varbeta_cl1[2]
  varbeta_cl1_mod = cbind(rho_v0*vcl1_b0 + (vcl1_b0/boot_iter)*(1 - rho_v0), rho_v1*vcl1_b1 + (vcl1_b1/boot_iter)*(1 - rho_v1))
  

  # -------------------------------------------------------------------------- #
  # -------------------------------------------------------------------------- #
  varbeta_cl2 = varbeta_cl0 + cvar
  varbeta_cl3 = varbeta_cl0 + cvar_rho
  varbeta_cl4 = varbeta_cl1 + cvar
  varbeta_cl5 = varbeta_cl1 + cvar_rho
  
  varbeta_cl6 = varbeta_cl0_mod + cvar
  varbeta_cl7 = varbeta_cl0_mod + cvar_rho
  varbeta_cl8 = varbeta_cl1_mod + cvar
  varbeta_cl9 = varbeta_cl1_mod + cvar_rho
  
  
  list(theta_ols = out1$theta, 
       alpha_est = out1$beta,
       beta_ols = beta_est,
       beta_update = beta_RBEGLS, 
       theta_update = out2$theta,
       cvar = cvar,
       rho_b0 = rho_b0,
       rho_b1 = rho_b1,
       rho_v0 = rho_v0,
       rho_v1 = rho_v1,
       varbeta_cl0 = varbeta_cl0, 
       varbeta_cl1 = varbeta_cl1,
       varbeta_cl2 = varbeta_cl2,
       varbeta_cl3 = varbeta_cl3,
       varbeta_cl4 = varbeta_cl4,
       varbeta_cl5 = varbeta_cl5,
       varbeta_cl6 = varbeta_cl6,
       varbeta_cl7 = varbeta_cl7,
       varbeta_cl8 = varbeta_cl8,
       varbeta_cl9 = varbeta_cl9)

#}

  varbeta_cl0
  varbeta_cl1
  varbeta_cl0_mod
  
  
  
  