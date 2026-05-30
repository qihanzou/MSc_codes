library(mvtnorm) 
library(geoR)
source("SimuData.R")
source("BFuns.R")
source("BFuns_TR.R")
source("BFuns_RBEGLS.R")
source("RBEGLS_functions.R")
source("fit_function_FULL.R")
setwd(R_address)
setwd(w1_address)


# -------------------------------------------------- #
#            Locations and Distance matrix
# -------------------------------------------------- #
locRdata=strwrap(paste("locs/loc", nall, ".Rdata",sep=""))
load(locRdata) 
plot(locall$x, locall$y)
Dall = matrix(0, nrow = nall, ncol = nall)            
for(i in 1:nall) {
  Dall[i, ] <- sqrt((locall$x[i] - locall$x)^2 + (locall$y[i] - locall$y)^2)
}

# -------------------------------------------------- #
#                  Initialization
# -------------------------------------------------- #
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
beta1 = 4

rate = 0.5
nug = 0


run.sim <- function(){
  # -------------------------------------------------- #
  # Simulation for two data sets
  # -------------------------------------------------- #
  xall = SimuData(xi0, rate, locall, cov_xi, nug, 1)
  epall = SimuData(ep_theta0, rate, locall, cov_ep, nug, 1)
  
  k = 1
  
  # -------------------------------------------------- #
  # pre-processing:
  # -------------------------------------------------- #
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
  
  
 
  
  
  # -------------------------------------------------- #
  # 1. OLS beta estimates and variance
  # -------------------------------------------------- #
  out1 = MLE.fit(x1s, designx1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
  betaout1 = out1$beta
  cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
  w = designx1%*%betaout1 + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - designx1s%*%betaout1) 
  W = cbind(1, w)
  beta_est = solve(t(W)%*%W)%*%t(W)%*%y 
  ff = lm(y~W-1) 
  
  # -------------------------------------------------- #
  # 2. original kriging method:
  # -------------------------------------------------- #
  out3 = MLE.fit(y, W, Dall[idx,idx], cov_ep, xi.ini, nug, "LB", lo.bound,up.bound)
  thetahat_UK = out3[[1]] 
  betahat_UK  = out3[[4]] 
  betahat_var_UK = diag(out3[[5]])

  
  # -------------------------------------------------- #
  # 3. Proposed Method and variance
  # -------------------------------------------------- #
  pc = designx1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% designx1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(designx1s)%*%solve(cov_xiall[-idx,-idx])%*%designx1s)%*%t(pc)
  CovUB = CovU*beta_est[2]^2 
  
  H = W%*%solve(t(W)%*%W)%*%t(W) 
  tmp_e = ff$resid/(1-diag(H)) 
  eet = as.matrix(tmp_e%*%t(tmp_e) - CovUB) #  ee^T - beta1^2 * Sigma(n)

  out2 = MLE.fit_tr(eet, Dall[idx,idx], cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  
  CorE = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE = out2$theta[2]*CorE
  CovUBE = CovUB + CovE
  beta_update = solve(t(W)%*%solve(CovUBE)%*%W)%*%t(W)%*%solve(CovUBE)%*%y
  
  # --------------------------------------------------#
  # 4. RBEGLS method
  # --------------------------------------------------#
  out4 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out4$eta, cov_ep, nug = 0)
  CovE_GLS = out4$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y
  
  # loop start here:
  res1 = RBEGLS_loop(30, 0.001, 0.001,0.001,0.001, beta_RBEGLS, out4$theta, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
  print(res1)
  
  beta_RBEGLS = cbind(res1[1], res1[2])
  theta_RBEGLS  = cbind(res1[3], res1[4])
  
  # --------------------------------------------------#
  # 5. MCI method
  # --------------------------------------------------#
  varW = CovU
  w2 = w + diag(varW)
  W2 = cbind(1, w2)
  beta_est_up = solve(t(W2)%*%W2)%*%t(W2)%*%y 
  
  w3 = w - diag(varW)
  W3 = cbind(1, w3)
  beta_est_down = solve(t(W3)%*%W3)%*%t(W3)%*%y 
  
  
  list2 = fit_function(y, W2, Dall, cov_xiall, designx1, designx1s, beta_est_up, ff, out1)
  list3 = fit_function(y, W3, Dall, cov_xiall, designx1, designx1s, beta_est_down, ff, out1)
  
  
  
  # CovUB_w2 = CovU*beta_est_up[2]^2 
  # H2 = W2%*%solve(t(W2)%*%W2)%*%t(W2) 
  # tmp_e2 = ff$resid/(1-diag(H2))
  # eet2 = as.matrix(tmp_e2%*%t(tmp_e2) - CovUB_w2) 
  # out_w2 = MLE.fit_tr(eet2, Dall[idx,idx], cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  # CorE_w2 = cor.mat(Dall[idx,idx], out_w2$eta, cov_ep, nug = 0)
  # CovE_w2 = out_w2$theta[2]*CorE_w2
  # CovUBE_w2 = CovUB_w2 + CovE_w2
  # beta_update_w2 = solve(t(W2)%*%solve(CovUBE_w2)%*%W2)%*%t(W2)%*%solve(CovUBE_w2)%*%y
  # 
  # CovUB_w3 = CovU*beta_est_down[2]^2 
  # H3 = W3%*%solve(t(W3)%*%W3)%*%t(W3) 
  # tmp_e3 = ff$resid/(1-diag(H3))
  # eet3 = as.matrix(tmp_e3%*%t(tmp_e3) - CovUB_w3) 
  # out_w3 = MLE.fit_tr(eet3, Dall[idx,idx], cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  # CorE_w3 = cor.mat(Dall[idx,idx], out_w3$eta, cov_ep, nug = 0)
  # CovE_w3 = out_w3$theta[2]*CorE_w3
  # CovUBE_w3 = CovUB_w3 + CovE_w3
  # beta_update_w3 = solve(t(W3)%*%solve(CovUBE_w3)%*%W3)%*%t(W3)%*%solve(CovUBE_w3)%*%y
  # 
  # 
  # CovUB_w2 = CovU*beta_update_w2[2]^2
  # CovUBE_w2 = CovUB_w2 + CovE_w2
  # tmp_w2 = t(W2)%*%solve(CovUBE_w2)%*%W2
  # varboot_w2 = solve(tmp_w2)
  # varbeta_w2  = diag(varboot_w2)   
  # 
  # CovUB_w3 = CovU*beta_update_w3[2]^2
  # CovUBE_w3 = CovUB_w3 + CovE_w3
  # tmp_w3 = t(W3)%*%solve(CovUBE_w3)%*%W3
  # varboot_w3 = solve(tmp_w3)
  # varbeta_w3  = diag(varboot_w3)
  
  
  
  
  # --------------------------------------------------#
  # 6. Variances:
  # --------------------------------------------------#
  Q2 = t(W)%*%W
  Q1 = t(W)%*%CovUBE%*%W
  varboot = solve(Q2)%*%Q1%*%solve(Q2)
  varbeta1 = diag(varboot)     # least square
  
  CovUB2 = CovU*beta_update[2]^2
  CovUBE2 = CovUB2 + CovE 
  tmp = t(W)%*%solve(CovUBE2)%*%W
  varboot2 = solve(tmp)
  varbeta2  = diag(varboot2)   # proposed
  
  CovUB3 = CovU*beta_RBEGLS[2]^2
  CovUBE3 = CovUB3 + CovE_GLS
  tmp3 = t(W)%*%solve(CovUBE3)%*%W
  varboot3 = solve(tmp3)
  varbeta3  = diag(varboot3)   # RBEGLS
  
  #####################################################
  
  
  
  list(beta_est = beta_est, beta_update = beta_update, beta_uk = betahat_UK, beta_RBEGLS = beta_RBEGLS, 
       beta_w2 = list2$beta_update, beta_w3 = list3$beta_update,
       thetax = out1$theta, thetae = out2$theta, theta_uk = out3$theta, theta_RBEGLS = theta_RBEGLS,
       theta_w2 = list2$thetae, theta_w3 = list3$thetae,
       varbeta1 = varbeta1, varbeta2 = varbeta2, varbeta_uk = betahat_var_UK, varbeta3 = varbeta3, 
       verbeta_w2 = list2$varbeta2, verbeta_w3 = list3$varbeta2,
       Ep = ep_error[1:5])
  
}  
