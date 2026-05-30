# test use:
#boot_iter = 100
#nall = 400
#idx = 1:200
#beta1 = 4
#R_address = "C:/Users/qihan/Desktop/QH_new"

setwd(R_address)
library(mvtnorm)               
source("SimuData.R")
source("BFuns.R")
source("SDCal_modified.R")
source("SDCal_RMLE.R")
source("BFuns_TR.R")
source("BFuns_RBEGLS.R")
source("RBEGLS_functions_New3.R")


locRdata=strwrap(paste("locs/loc", nall, ".Rdata",sep=""))
load(locRdata)


Dall = matrix(0, nrow = nall, ncol = nall)            
for(i in 1:nall) {
  Dall[i, ] <- sqrt((locall$x[i] - locall$x)^2 + (locall$y[i] - locall$y)^2)
}

#------------------------------------------------------------------------------#

alpha0 = c(1,0.5)
xi0 = c(3, 4)  # (phi, sigma^2)
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


gamma0 = c(2,1)    # new
kappa0 = c(2,3)    # new
cov_kappa = "Gau"  # new
kappa_ini = c(2.5) # new
beta3 = 3          # new


zeta0 = c(2,2,2)   # new
iota0 = c(2.5,2.5) # new
cov_iota = "Sph"   # new
iota_ini = c(2)    # new
beta4 = 5          # new

  
  
nug = 0



run.sim <- function(){


  etaall = SimuData(xi0, rate, locall, cov_xi, nug, 1)
  epall = SimuData(ep_theta0, rate, locall, cov_ep, nug, 1)
  etaall3 = SimuData(kappa0, rate, locall, cov_kappa, nug, 1) # new!
  etaall4 = SimuData(iota0, rate, locall, cov_iota, nug, 1) # new!
  
  k=1
  
  distall = sqrt((locall$x - 0)^2 + (locall$y - 0)^2)
  
  # ------------------------------------------------------------------#
  r1all = as.matrix(cbind(1, distall)) 
  r1s = as.matrix(cbind(1, distall[-idx]))
  r1 = as.matrix(cbind(1,distall[idx]))
  
  x1all = r1all%*%alpha0 + etaall$Ymat[k,]  
  x1 = x1all[idx] # complete misaligned
  x1s = x1all[-idx] # that I have
  N = length(x1) 
  M = length(x1s)
  # -------------------------------------------------- #
  xy = locall$x*locall$y
  p3all = as.matrix(cbind(1, xy))
  x3all = p3all%*%gamma0 + etaall3$Ymat[k,]
  x3 = as.matrix(x3all[idx])
  x3s = as.matrix(x3all[-idx]) 
  
  idx_new = 1:(length(idx) + length(-idx)/2)
  x3tar = as.matrix(x3all[idx_new])
  x3obs = as.matrix(x3all[-idx_new])
  p3tar = as.matrix(cbind(1, xy[idx_new]))
  p3obs = as.matrix(cbind(1, xy[-idx_new]))
  # -------------------------------------------------- #
  xdy = locall$x/locall$y
  ydx = locall$y/locall$x
  p4all = as.matrix(cbind(1, xdy, ydx))
  x4all = p4all%*%zeta0 + etaall4$Ymat[k,]
  x4 = as.matrix(x4all[idx])
  x4s = as.matrix(x4all[-idx]) 

  idx_new4 = 1:(length(idx) - length(-idx)/2)
  x4tar = as.matrix(x4all[idx_new4])
  x4obs = as.matrix(x4all[-idx_new4])
  p4tar = as.matrix(cbind(1, xdy[idx_new4], ydx[idx_new4]))
  p4obs = as.matrix(cbind(1, xdy[-idx_new4], ydx[-idx_new4]))
  # -------------------------------------------------- #

  ep_error = epall$Ymat[k,idx]
  x2 = distall[idx]
  y = beta0 + beta1*x1 + beta2*x2 + beta3*x3 + beta4*x4 + ep_error    
  # -------------------------------------------------- #
  # 1. compute w1, w3, w4 and OLS beta estimates
  # -------------------------------------------------- #
  # for w1.
  # -------------------------------------------------- #
  out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
  alpha_est = out1$beta
  cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
  w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 
  # -------------------------------------------------- #
  # for w3.
  # -------------------------------------------------- #
  out_w3 = MLE.fit(x3obs, p3obs, Dall[-idx_new,-idx_new], cov_kappa, kappa_ini, nug, "LB", lo.bound,up.bound)
  gamma_est = out_w3$beta
  cov_kappaall = out_w3$theta[2]*cor.mat(Dall, out_w3$eta, cov_kappa, nug = 0)
  w3_all = p3tar%*%gamma_est + cov_kappaall[idx_new, -idx_new]%*%solve(cov_kappaall[-idx_new,-idx_new])%*%(x3obs - p3obs%*%gamma_est)
  w3 = w3_all[idx]
  # -------------------------------------------------- #
  # for w4.
  # -------------------------------------------------- #
  out_w4 = MLE.fit(x4obs, p4obs, Dall[-idx_new4,-idx_new4], cov_iota, iota_ini, nug, "LB", lo.bound,up.bound)
  zeta_est = out_w4$beta
  cov_iotaall = out_w4$theta[2]*cor.mat(Dall, out_w4$eta, cov_iota, nug = 0)
  w4_partial = p4tar%*%zeta_est + cov_iotaall[idx_new4, -idx_new4]%*%solve(cov_iotaall[-idx_new4,-idx_new4])%*%(x4obs - p4obs%*%zeta_est)
  w4 = as.matrix(c(w4_partial, x4obs[1:(length(idx)/2)]))
  # -------------------------------------------------- #
  # OLS
  # -------------------------------------------------- #
  W = cbind(1, w, x2, w3, w4)
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
  # pc w1
  # -------------------------------------------------- #
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  CovUB_w1 = CovU*beta_ols[2]^2 
  # -------------------------------------------------- #
  # pc w3
  # -------------------------------------------------- #
  pc_w3 = p3tar - cov_kappaall[idx_new, -idx_new]%*%solve(cov_kappaall[-idx_new,-idx_new]) %*% p3obs
  CovU_w3_all = cov_kappaall[idx_new, idx_new] - cov_kappaall[idx_new, -idx_new]%*%solve(cov_kappaall[-idx_new,-idx_new]) %*% t(cov_kappaall[idx_new, -idx_new]) + (pc_w3)%*%solve(t(p3obs)%*%solve(cov_kappaall[-idx_new,-idx_new])%*%p3obs)%*%t(pc_w3)
  CovU_w3 = CovU_w3_all[idx,idx]
  CovUB_w3 = CovU_w3*beta_ols[4]^2
  # -------------------------------------------------- #
  # pc w4
  # -------------------------------------------------- #
  pc_w4 = p4tar - cov_iotaall[idx_new4, -idx_new4]%*%solve(cov_iotaall[-idx_new4,-idx_new4]) %*% p4obs
  CovU_w4_partial = cov_iotaall[idx_new4, idx_new4] - cov_iotaall[idx_new4, -idx_new4]%*%solve(cov_iotaall[-idx_new4,-idx_new4]) %*% t(cov_iotaall[idx_new4, -idx_new4]) + (pc_w4)%*%solve(t(p4obs)%*%solve(cov_iotaall[-idx_new4,-idx_new4])%*%p4obs)%*%t(pc_w4)
  zero_matrix = CovU_w4_partial*0
  mm1 = cbind(CovU_w4_partial, zero_matrix)
  mm2 = cbind(zero_matrix, zero_matrix)
  CovU_w4 = rbind(mm1,mm2)
  CovUB_w4 = CovU_w4*beta_ols[5]^2
  # -------------------------------------------------- #
  # compute CovUB
  # -------------------------------------------------- #
  CovUB = CovUB_w1 + CovUB_w3 + CovUB_w4
  # -------------------------------------------------- #
  out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS = out2$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y
  
  source("RBEGLS_functions_New3.R")
  # -------------------------------------------------- #
  # loop start here:
  # -------------------------------------------------- #
  res1 = RBEGLS_loop_New3(50, 0.001, 0.001, beta_RBEGLS, out2$theta, theta_ep.ini, CovU, CovU_w3, CovU_w4, y, W, Dall, cov_ep, lo.bounde,up.bounde)
  
  beta_RBEGLS = res1$beta_est
  theta_RBEGLS  = res1$theta_est
  CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
  CovE_GLS = theta_RBEGLS[2]*CorE_GLS
  CovUB = CovU*beta_RBEGLS[2]^2 + CovU_w3*beta_RBEGLS[4]^2 + CovU_w4*beta_RBEGLS[5]^2
  CovUBE_GLS = CovUB + CovE_GLS
  CovUBE_GLS_inv = solve(CovUBE_GLS)
  
  
  
  # -------------------------------------------------- #
  ## Bootstrap initialization:
  sdout = SDCal(Dall[-idx,-idx], out1$theta, cov_xi, nug = 0)
  sdout_w3 = SDCal(Dall[-idx_new,-idx_new], out_w3$theta, cov_kappa, nug = 0)
  sdout_w4 = SDCal(Dall[-idx_new4,-idx_new4], out_w4$theta, cov_iota, nug = 0)

  ETAallboot_all    = SimuData(out1$theta, rate, locall, cov_xi, nug, boot_iter)
  ETAallboot_all_w3 = SimuData(out_w3$theta, rate, locall, cov_kappa, nug, boot_iter)
  ETAallboot_all_w4 = SimuData(out_w4$theta, rate, locall, cov_iota,  nug, boot_iter)
  
  EVar = EVarKR = EVarKRT = matrix(0, 5, 5)
  # -------------------------------------------------- #
  # step3: parameter bootstrap start here:
  for (v in 1:boot_iter){
    # -------------------------------------------------- #
    # for x1, x3, x4
    # -------------------------------------------------- #
    x1allboot = r1all%*%alpha_est + ETAallboot_all$Ymat[v,]
    x1sboot = x1allboot[-idx]
    
    x3allboot = p3all%*%gamma_est + ETAallboot_all_w3$Ymat[v,]
    x3obsboot = as.matrix(x3allboot[-idx_new])
    
    x4allboot = p4all%*%zeta_est + ETAallboot_all_w4$Ymat[v,] # here wrong. should be v
    x4obsboot = as.matrix(x4allboot[-idx_new4])
    # -------------------------------------------------- #

    #(b) simulate xi, kappa, iota
    xie = out1$theta + rmvnorm(1,sigma=sdout$covmat)
    cov_xialle = cor.mat(Dall, xie[1], cov_xi, nug = 0)
    
    kappae = out_w3$theta + rmvnorm(1, sigma = sdout_w3$covmat)
    cov_kappaalle = cor.mat(Dall, kappae[1], cov_kappa, nug = 0)
    
    iotae = out_w4$theta + rmvnorm(1, sigma = sdout_w4$covmat)
    cov_iotaalle = cor.mat(Dall, iotae[1], cov_iota, nug = 0)
    
    
    
    # (c): plug parameters into kriging equation to derive Wj:
    tmp_inv = solve(cov_xialle[-idx,-idx])
    alphaboot = solve(t(r1s) %*% tmp_inv %*% r1s)%*%t(r1s) %*% tmp_inv %*% x1s # why we use x1s rather than x1sboot here?
    wboot_e = r1%*%alphaboot + cov_xialle[idx, -idx]%*%tmp_inv%*%(x1sboot - r1s%*%alphaboot)
    
    tmp_inv_w3 = solve(cov_kappaalle[-idx_new,-idx_new])
    gammaboot = solve(t(p3obs) %*% tmp_inv_w3 %*% p3obs)%*%t(p3obs) %*% tmp_inv_w3 %*% x3obs
    wboot3_all = p3tar%*%gammaboot + cov_kappaalle[idx_new, -idx_new]%*%tmp_inv_w3%*%(x3obsboot - p3obs%*%gammaboot)
    w3boote = wboot3_all[idx]
    
    tmp_inv_w4 = solve(cov_iotaalle[-idx_new4,-idx_new4])
    zetaboot = solve(t(p4obs) %*% tmp_inv_w4 %*% p4obs)%*%t(p4obs) %*% tmp_inv_w4 %*% x4obs
    wboot4_partial = p4tar%*%zetaboot + cov_iotaalle[idx_new4, -idx_new4]%*%tmp_inv_w4%*%(x4obsboot - p4obs%*%zetaboot)
    w4boote = as.matrix(c(wboot4_partial, x4obs[1:(length(idx)/2)]))
    
    
    Wboot_e = cbind(1, wboot_e, x2, w3boote, w4boote)
    # ------------------------------------------------------------------------- #
    A1_e = solve(t(Wboot_e)%*%CovUBE_GLS_inv%*%Wboot_e)
    EVar = EVar + A1_e
    
    B1_e = solve(t(Wboot_e)%*%covKR.inv%*%Wboot_e)                       #
    B2_e = t(Wboot_e)%*%covKR.inv%*%CovUBE_GLS%*%covKR.inv%*%Wboot_e     #
    EVarKR = EVarKR + solve(t(Wboot_e)%*%covKR.inv%*%Wboot_e)   #
    EVarKRT = EVarKRT + B1_e%*%B2_e%*%B1_e                      #   
  }
  
  # step4: Calculate the parametric bootstrap SE as the empirical standard deviation of beta:
  # -------------------------------------------------- #
  # -------------------------------------------------- #
  A1 = solve(t(W)%*%CovUBE_GLS_inv%*%W)
  varbeta_rb = diag(A1)
  varbeta_rbm = diag(EVar/boot_iter)
  
  varbeta_kr = beta_var_kr                             #
  varbeta_krm = diag(EVarKR/boot_iter)                 #
  
  B1 = solve(t(W)%*%covKR.inv%*%W)                     #
  B2 = t(W)%*%covKR.inv%*%CovUBE_GLS%*%covKR.inv%*%W   #
  varbeta_krt = diag(B1%*%B2%*%B1)                     #
  varbeta_krtm = diag(EVarKRT/boot_iter)               #
   
  # -------------------------------------------------- #

  list(xi_est = out1$theta, alpha_est = out1$beta,
       kappa_est = out_w3$theta, gamma_est = out_w3$beta,
       iota_est = out_w4$theta, zeta_est = out_w4$beta,
       
       beta_ols = beta_ols, 
       beta_kr = beta_kr, theta_kr = theta_kr,
       beta_RBEGLS = beta_RBEGLS,  theta_est = theta_RBEGLS,
       
       varbeta_ols = varbeta_ols, varbeta_kr = varbeta_kr, varbeta_krm = varbeta_krm, 
       varbeta_krt = varbeta_krt, varbeta_krtm = varbeta_krtm, 
       varbeta_rb = varbeta_rb, varbeta_rbm = varbeta_rbm)
       

}

