# # # # test use:
rm(list = ls(all = TRUE)) 
nall = 200
idx = 1:100
beta1 = 4
MC_iter = 100
N_sim = 1000
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
library(sigmoid)
library(expm)

find_T_matrix <- function(Sigma_x, Sigma_z) {
  S_z_sqrt <- sqrtm(Sigma_z)
  middle <- S_z_sqrt %*% Sigma_x %*% S_z_sqrt
  middle_sqrt_inv <- solve(sqrtm(middle))
  Tm <- S_z_sqrt %*% middle_sqrt_inv %*% S_z_sqrt
  return(Tm)
}

# define the non-linear transformation here:
nl_trans = function(input){
  output = input^2
  #output = input
  #output = sigmoid(input, method = "logistic")  # logistic
  #output = sigmoid(input, method = "Gompertz")  # Gompertz
  #output = sigmoid(input, method = "tanh")      # tanh
  #output = sigmoid(input, method = "ReLU")      # ReLU linear in our case
  #output = sigmoid(input, method = "leakyReLU") # leakyReLU
  #output = sin(input)                           # Sinusoid
}


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
y = beta0 + beta1*nl_trans(x1) + ep_error 
y2 = beta0 + beta1*x1 + ep_error
hist(y, breaks = 100)
hist(y2, breaks = 100)










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
W_set = NULL
for (i in 1:N_sim){
  w_i = t(w) + rmvnorm(1, sigma = CovU)
  z_i = nl_trans(w_i)
  Z_set = rbind(Z_set, z_i)
  W_set = rbind(W_set, w_i)
}
Cov_Z = matrix(0, length(idx), length(idx))
for (j in 1:N_sim){
  Cov_Z = Cov_Z + (Z_set[j,] - colMeans(Z_set))%*%t(Z_set[j,] - colMeans(Z_set))
}
CovU_nl = 1/(N_sim-1)*Cov_Z
z = colMeans(Z_set)
wbar = colMeans(W_set)


Act1 = find_T_matrix(CovU, CovU_nl)
Act2 = solve(CovU)%*%CovU_nl
#bures_distance(CovU_nl, CovU)
#bures_distance(CovU, CovU)
#AA = Tm %*% CovU %*% t(Tm)
#bures_distance(AA, CovU_nl)
#all.equal(Tm %*% CovU %*% t(Tm), CovU_nl)

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

#--------------------------------------------------------------------------------------------------------------
# RBEGLS
up.bounde = c(100,100)

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
# MC 
sdout = SDCal(Dall[-idx,-idx], out1$theta, cov_xi, nug = 0)
ETAallboot_all = SimuData(out1$theta, rate, locall, cov_xi, nug, MC_iter)
EVar = EVar_act1 = EVar_act2 = matrix(0, 2, 2)
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
  Zboot_e = cbind(1, nl_trans(wboot_e))
  A1_e = solve(t(Zboot_e)%*%CovUBE_GLS_inv%*%Zboot_e)
  EVar = EVar + A1_e         
  
  pce1 = r1 - cov_xialle[idx, -idx]%*%solve(cov_xialle[-idx,-idx]) %*% r1s
  CovUe1 = cov_xialle[idx, idx] - cov_xialle[idx, -idx]%*%solve(cov_xialle[-idx,-idx]) %*% t(cov_xialle[idx, -idx]) + (pce1)%*%solve(t(r1s)%*%solve(cov_xialle[-idx,-idx])%*%r1s)%*%t(pce1)
  CovU_nl_e1_act1 = Act1 %*% CovUe1 %*% t(Act1)
  CovU_nl_e1_act2 = CovUe1%*%Act2
  
  CovUB_e1_act1 = CovU_nl_e1_act1*beta_RBEGLS[2]^2
  CovUBE_GLS_e1_act1 = CovUB_e1_act1 + CovE_GLS
  A1_e1_act1 = solve(t(Zboot_e)%*%solve(CovUBE_GLS_e1_act1)%*%Zboot_e)
  EVar_act1 = EVar_act1 + A1_e1_act1
  
  CovUB_e1_act2 = CovU_nl_e1_act2*beta_RBEGLS[2]^2
  CovUBE_GLS_e1_act2 = CovUB_e1_act2 + CovE_GLS
  A1_e1_act2 = solve(t(Zboot_e)%*%solve(CovUBE_GLS_e1_act2)%*%Zboot_e)
  EVar_act2 = EVar_act2 + A1_e1_act2
}
varbeta_mc = diag(EVar/MC_iter)
varbeta_mc_act1 = diag(EVar_act1/MC_iter)
varbeta_mc_act2 = diag(EVar_act2/MC_iter)

# ---------------------------------------------------------------------------- #
# MC modified
EVar2  = EVar2_act1 = EVar2_act2 = matrix(0, 2, 2)

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
  Zboot_e2 = cbind(1, nl_trans(wboot_e2))
  A1_e2 = solve(t(Zboot_e2)%*%CovUBE_GLS_inv%*%Zboot_e2)
  EVar2 = EVar2 + A1_e2           
  
  pce2 = r1 - cov_xialle2[idx, -idx]%*%solve(cov_xialle2[-idx,-idx]) %*% r1s
  CovUe2 = cov_xialle2[idx, idx] - cov_xialle2[idx, -idx]%*%solve(cov_xialle2[-idx,-idx]) %*% t(cov_xialle2[idx, -idx]) + (pce2)%*%solve(t(r1s)%*%solve(cov_xialle2[-idx,-idx])%*%r1s)%*%t(pce2)
  CovU_nl_e2_act1 = Act1 %*% CovUe2 %*% t(Act1)
  CovU_nl_e2_act2 = CovUe2%*%Act2
  
  CovUB_e2_act1 = CovU_nl_e2_act1*beta_RBEGLS[2]^2
  CovUBE_GLS_e2_act1 = CovUB_e2_act1 + CovE_GLS
  A1_e2_act1 = solve(t(Zboot_e2)%*%solve(CovUBE_GLS_e2_act1)%*%Zboot_e2)
  EVar2_act1 = EVar2_act1 + A1_e2_act1
  
  CovUB_e2_act2 = CovU_nl_e2_act2*beta_RBEGLS[2]^2
  CovUBE_GLS_e2_act2 = CovUB_e2_act2 + CovE_GLS
  A1_e2_act2 = solve(t(Zboot_e2)%*%solve(CovUBE_GLS_e2_act2)%*%Zboot_e2)
  EVar2_act2 = EVar2_act2 + A1_e2_act2
}
varbeta_mc2 = diag(EVar2/MC_iter)
varbeta_mc2_act1 = diag(EVar2_act1/MC_iter)
varbeta_mc2_act2 = diag(EVar2_act2/MC_iter)



#--------------------------------------------------------------------------------------------------------------
# MC under construction
ETAallboot_all = SimuData(out1$theta, rate, locall, cov_xi, nug, MC_iter)
EVar_krig = EVarKR = EVarKRT = matrix(0, 2, 2)
for (v in 1:MC_iter){
  x1allboot = r1all%*%alpha_est + ETAallboot_all$Ymat[v,]
  x1sboot = x1allboot[-idx]
  out_e = MLE.fit(x1sboot, r1s, Dall[-idx,-idx], cov_xi, out1$theta, nug, "LB", lo.bound,up.bound)
  alpha_est_e = out_e$beta
  cov_xialle = out_e$theta[2]*cor.mat(Dall, out_e$eta, cov_xi, nug = 0)
  
  tmp_inv = solve(cov_xialle[-idx,-idx])
  alphaboot = solve(t(r1s) %*% tmp_inv %*% r1s)%*%t(r1s) %*% tmp_inv %*% x1s
  wboot_e = r1%*%alphaboot + cov_xialle[idx, -idx]%*%tmp_inv%*%(x1s - r1s%*%alphaboot)
  Zboot_e = cbind(1, nl_trans(wboot_e))
  
  pc_e = r1 - cov_xialle[idx, -idx]%*%solve(cov_xialle[-idx,-idx]) %*% r1s
  CovU_e = cov_xialle[idx, idx] - cov_xialle[idx, -idx]%*%solve(cov_xialle[-idx,-idx]) %*% t(cov_xialle[idx, -idx]) + (pc_e)%*%solve(t(r1s)%*%solve(cov_xialle[-idx,-idx])%*%r1s)%*%t(pc_e)
  CovU_nl_e = Act1 %*% CovU_e %*% t(Act1)
  
  CovUB_e = CovU_nl_e*beta_RBEGLS[2]^2
  CovUBE_GLS_e = CovUB_e + CovE_GLS
  CovUBE_GLS_inv_e = solve(CovUBE_GLS_e)
  
  A1_e = solve(t(Zboot_e)%*%CovUBE_GLS_inv_e%*%Zboot_e)
  EVar_krig = EVar_krig + A1_e                      
}
varbeta_mc_krig = diag(EVar_krig/MC_iter)
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
     varbeta_mc_act1 = as.matrix(as.numeric(varbeta_mc_act1)),
     varbeta_mc_act2 = as.matrix(as.numeric(varbeta_mc_act2)),
     varbeta_mc2 = as.matrix(as.numeric(varbeta_mc2)),
     varbeta_mc2_act1 = as.matrix(as.numeric(varbeta_mc2_act1)),
     varbeta_mc2_act2 = as.matrix(as.numeric(varbeta_mc2_act2)),
     varbeta_mc_krig = as.matrix(as.numeric(varbeta_mc_krig))
)


#}


