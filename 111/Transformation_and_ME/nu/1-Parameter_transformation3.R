# # # # test use:
rm(list = ls(all = TRUE)) 
nall = 200
idx = 1:100
beta1 = 4
MC_iter = 100
R_address = "C:/Users/qihan/Desktop/Transformation_and_ME"


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
source("Matrix_Funs.R")
library(sigmoid)
library(numDeriv)


# Define the transformation here: 
g = function(w) 1/w


# -------------------------------------------------------------------------- #
locRdata=strwrap(paste("locs/loc", nall, ".Rdata",sep=""))
load(locRdata)

Dall = matrix(0, nrow = nall, ncol = nall)            
for(i in 1:nall) {
  Dall[i, ] <- sqrt((locall$x[i] - locall$x)^2 + (locall$y[i] - locall$y)^2)
}
#plot(locall$x, locall$y)
#------------------------------------------------------------------------------#

alpha0 = c(1,0.5)
xi0 = c(3, 4)  # initial values of (phi, sigma^2) (3, 4)
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
y = beta0 + beta1*g(x1) + ep_error 

#--------------------------------------------------------------------------------------------------------------
# w
out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
alpha_est = out1$beta
cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)

w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 
pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)

#gdist_mc = DT_MC(w, CovU, g, 2000) # if use MCT method
gdist = DT_Delta_1st(w, CovU, g) # 1st order Delta method
CovU_nl = gdist$cov_z
z = gdist$z

#round(gdist$z - gdist_mc$z)
# round(gdist$cov_z - CovU_mc)
#bures_distance(gdist$cov_z, gdist_mc$cov_z)
# bures_distance(CovU, CovU_mc)
# bures_distance(CovU, gdist$cov_z)

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

theta_RBEGLS  =  res1$theta_est 
CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
CovE_GLS = theta_RBEGLS[2]*CorE_GLS
CovUB = CovU_nl*beta_RBEGLS[2]^2
CovUBE_GLS = CovUB + CovE_GLS
CovUBE_GLS_inv = solve(CovUBE_GLS)
A1 = solve(t(Z)%*%CovUBE_GLS_inv%*%Z)
varbeta_rb = diag(A1)


#--------------------------------------------------------------------------------------------------------------
# MC 
sdout = SDCal(Dall[-idx,-idx], out1$theta, cov_xi, nug = 0)
ETAallboot_all = SimuData(out1$theta, rate, locall, cov_xi, nug, MC_iter)
EVar = EVar_2 = EVar_3 = matrix(0, 2, 2)
for (v in 1:MC_iter){
  
  x1allboot = r1all%*%alpha_est + ETAallboot_all$Ymat[v,]
  x1sboot = x1allboot[-idx]
  
  xie = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  while ((xie[1] < 0)|(xie[2] < 0)) {
    xie = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  }
  cov_xialle = xie[2]*cor.mat(Dall, xie[1], cov_xi, nug = 0)
  tmp_inv = solve(cov_xialle[-idx,-idx])
  alphaboot = solve(t(r1s) %*% tmp_inv %*% r1s)%*%t(r1s) %*% tmp_inv %*% x1sboot 
  wboot_e = r1%*%alphaboot + cov_xialle[idx, -idx]%*%tmp_inv%*%(x1sboot - r1s%*%alphaboot)
  Zboot_e = cbind(1, g(wboot_e))
  
  A1_e = solve(t(Zboot_e)%*%CovUBE_GLS_inv%*%Zboot_e)
  EVar = EVar + A1_e                     
  
  
  pc_e = r1 - cov_xialle[idx, -idx]%*%solve(cov_xialle[-idx,-idx]) %*% r1s
  CovU_e = cov_xialle[idx, idx] - cov_xialle[idx, -idx]%*%solve(cov_xialle[-idx,-idx]) %*% t(cov_xialle[idx, -idx]) + (pc_e)%*%solve(t(r1s)%*%solve(cov_xialle[-idx,-idx])%*%r1s)%*%t(pc_e)
  gdist_e = DT_Delta_1st(wboot_e, CovU_e, g) # 1st order Delta method
  CovU_nl_e = gdist_e$cov_z
  z_e = gdist_e$z
  Z_e = cbind(1, z_e)
  
  CovUB_e = CovU_nl_e*beta_RBEGLS[2]^2
  CovUBE_GLS_e = CovUB_e + CovE_GLS
  
  #CovUBE_GLS_e = CovUBE_GLS_e + diag(1e-3, nrow(CovUBE_GLS_e))
  
  CovUBE_GLS_inv_e = solve(CovUBE_GLS_e)
  
  A1_e_2 = solve(t(Z_e)%*%CovUBE_GLS_inv_e%*%Z_e)
  EVar_2 = EVar_2 + A1_e_2
  
  EVar_3 = EVar_3 + Find_expand_var(CovUBE_GLS, CovUBE_GLS_e, Z_e)$cov
}
varbeta_mc = diag(EVar/MC_iter)
varbeta_mc_2 = diag(EVar_2/MC_iter)
varbeta_mc_3 = diag(EVar_3/MC_iter)

# ---------------------------------------------------------------------------- #
# MC modified
EVar2  = EVar2_2 = EVar2_3 = matrix(0, 2, 2)

for (v2 in 1:MC_iter){
  
  xie2 = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  while ((xie2[1] < 0)|(xie2[2] < 0)) {
    xie2 = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  }
  
  ETAallboot_all2 = SimuData(xie2, rate, locall, cov_xi, nug, 1)
  x1allboot2 = r1all%*%alpha_est + ETAallboot_all2$Ymat[1,]
  x1sboot2 = x1allboot2[-idx]
  
  cov_xialle2 = xie2[2]*cor.mat(Dall, xie2[1], cov_xi, nug = 0)
  tmp_inv2 = solve(cov_xialle2[-idx,-idx])
  alphaboot2 = solve(t(r1s) %*% tmp_inv2 %*% r1s)%*%t(r1s) %*% tmp_inv2 %*% x1sboot2 
  wboot_e2 = r1%*%alphaboot2 + cov_xialle2[idx, -idx]%*%tmp_inv2%*%(x1sboot2 - r1s%*%alphaboot2)
  Zboot_e2 = cbind(1, g(wboot_e2))
  
  A1_e2 = solve(t(Zboot_e2)%*%CovUBE_GLS_inv%*%Zboot_e2)
  EVar2 = EVar2 + A1_e2    
  
  pc_e2 = r1 - cov_xialle2[idx, -idx]%*%solve(cov_xialle2[-idx,-idx]) %*% r1s
  CovU_e2 = cov_xialle2[idx, idx] - cov_xialle2[idx, -idx]%*%solve(cov_xialle2[-idx,-idx]) %*% t(cov_xialle2[idx, -idx]) + (pc_e2)%*%solve(t(r1s)%*%solve(cov_xialle2[-idx,-idx])%*%r1s)%*%t(pc_e2)
  gdist_e2 = DT_Delta_1st(wboot_e2, CovU_e2, g) # 1st order Delta method
  CovU_nl_e2 = gdist_e2$cov_z
  z_e2 = gdist_e2$z
  Z_e2 = cbind(1, z_e2)
  
  CovUB_e2 = CovU_nl_e2*beta_RBEGLS[2]^2
  CovUBE_GLS_e2 = CovUB_e2 + CovE_GLS
  
  #CovUBE_GLS_e2 = CovUBE_GLS_e2 + diag(1e-3, nrow(CovUBE_GLS_e2))
  
  CovUBE_GLS_inv_e2 = solve(CovUBE_GLS_e2)
  
  A1_e2_2 = solve(t(Z_e2)%*%CovUBE_GLS_inv_e2%*%Z_e2)
  EVar2_2 = EVar2_2 + A1_e2_2
  
  EVar2_3 = EVar2_3 + Find_expand_var(CovUBE_GLS, CovUBE_GLS_e2, Z_e2)$cov

  
}
varbeta_mc2 = diag(EVar2/MC_iter)
varbeta_mc2_2 = diag(EVar2_2/MC_iter)
varbeta_mc2_3 = diag(EVar2_3/MC_iter)

# ---------------------------------------------------------------------------- #
# MC modified 2


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
     varbeta_mc = as.matrix(as.numeric(varbeta_mc)),
     varbeta_mc_2 = as.matrix(as.numeric(varbeta_mc_2)),
     varbeta_mc_3 = as.matrix(as.numeric(varbeta_mc_3)),
     varbeta_mc2 = as.matrix(as.numeric(varbeta_mc2)),
     varbeta_mc2_2 = as.matrix(as.numeric(varbeta_mc2_2)),
     varbeta_mc2_3 = as.matrix(as.numeric(varbeta_mc2_3))
)


#}


