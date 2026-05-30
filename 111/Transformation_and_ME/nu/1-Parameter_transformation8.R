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
library(nloptr)

# Define the transformation here: 
# logistic transformation:
#g = function(w) 1/(1+exp(-w))
#dg = function(w) (1/(1+exp(-w)))*(1-1/(1+exp(-w)))

# identity transformation g(x) = x:
g = function(w) w
#dg = function(w) w-w+1

# exp transformation g(x) = exp(x):
#g = function(w) exp(w)
#dg = function(w) exp(w)

# power of 2 transformation g(x) = x^2:
# g = function(w) w^2
# dg = function(w) 2*w

# power of 3 transformation g(x) = x^3:
# g = function(w) w^3
# dg = function(w) 3*w^2

# tanh
#g = function(w) tanh(w)
#dg = function(w) 1 - tanh(w)^2


# -------------------------------------------------------------------------- #
locRdata=strwrap(paste("locs/loc", nall, ".Rdata",sep=""))
load(locRdata)

Dall = matrix(0, nrow = nall, ncol = nall)            
for(i in 1:nall) {
  Dall[i, ] <- sqrt((locall$x[i] - locall$x)^2 + (locall$y[i] - locall$y)^2)
}
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

#  simple model be applied g(x)
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

gdist_mc = DT_MC(w, CovU, g, 1000) # if use MCT method
CovU_nl_mc = gdist_mc$cov_z
z_mc = gdist_mc$z

gdist = DT_Delta_1st_jacobianfun(w, CovU, g) # 1st order Delta method
CovU_nl = gdist$cov_z
z = gdist$z

bdist = bures_distance(gdist$cov_z, gdist_mc$cov_z)
#--------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------
# KR
Z = cbind(1, z)
kr = MLE.fit(y, Z, Dall[idx,idx], cov_ep, c(2.5), nug, "LB", lo.bound,up.bound)
theta_kr = kr[[1]] 
beta_kr  = kr[[4]] 
varbeta_kr = diag(kr[[5]])
#--------------------------------------------------------------------------------------------------------------
# OLS dm
fols = lm(y~Z-1)
beta_ols = fols$coef
varbeta_ols = diag(vcov(fols))

# RBEGLS dm
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
varbeta_rb_dm = diag(A1)

#--------------------------------------------------------------------------------------------------------------
# OLS mc
Z_mc = cbind(1, z_mc)
fols_mc = lm(y~Z_mc-1)
beta_ols_mc = fols_mc$coef

# RBEGLS mc
CovUB_mc = CovU_nl_mc*beta_ols_mc[2]^2 
out2_mc = MLE.fit_GLS(y, Z_mc, Dall[idx,idx], CovUB_mc, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
CorE_GLS_mc = cor.mat(Dall[idx,idx], out2_mc$eta, cov_ep, nug = 0)
CovE_GLS_mc = out2_mc$theta[2]*CorE_GLS_mc
CovUBE_GLS_mc = CovUB_mc + CovE_GLS_mc
beta_RBEGLS_mc = solve(t(Z_mc)%*%solve(CovUBE_GLS_mc)%*%Z_mc)%*%t(Z_mc)%*%solve(CovUBE_GLS_mc)%*%y
# loop:
res1_mc = RBEGLS_loop_final(50, 0.001, 0.001, beta_RBEGLS_mc, out2_mc$theta, theta_ep.ini, CovU_nl_mc, y, Z_mc, Dall, cov_ep, lo.bounde,up.bounde)
beta_RBEGLS_mc = res1_mc$beta_est
theta_RBEGLS_mc  =  res1_mc$theta_est 
CorE_GLS_mc = cor.mat(Dall[idx,idx], theta_RBEGLS_mc[1], cov_ep, nug = 0)
CovE_GLS_mc = theta_RBEGLS_mc[2]*CorE_GLS_mc
CovUB_mc = CovU_nl_mc*beta_RBEGLS_mc[2]^2
CovUBE_GLS_mc = CovUB_mc + CovE_GLS_mc
CovUBE_GLS_inv_mc = solve(CovUBE_GLS_mc)
A1_mc = solve(t(Z_mc)%*%CovUBE_GLS_inv_mc%*%Z_mc)
varbeta_rb_mc = diag(A1_mc)

# Mean of two models:
varbeta_rb_mean = (varbeta_rb_dm + varbeta_rb_mc)/2
beta_mean = (beta_RBEGLS + beta_RBEGLS_mc)/2
theta_mean = (theta_RBEGLS + theta_RBEGLS_mc)/2
CovE_GLS_mean = theta_mean[2]*cor.mat(Dall[idx,idx], theta_mean[1], cov_ep, nug = 0)
CovUBE_GLS_mean = (CovUB_mc + CovUB)/2 + CovE_GLS_mean

# Robust sandwich variance for beta:
varbeta_sdm = diag(sandwich_cov(CovUBE_GLS, y, Z, beta_RBEGLS, Diagonal = TRUE)$VarB)
varbeta_smc = diag(sandwich_cov(CovUBE_GLS_mc, y, Z_mc, beta_RBEGLS_mc, Diagonal = TRUE)$VarB)
varbeta_smean = diag(sandwich_cov(CovUBE_GLS_mean, y, (Z_mc+Z)/2, beta_mean, Diagonal = TRUE)$VarB)

# Base estimates
theta_base_constrained_leq      <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS, CovUB, x0 = c(3.5, 1.2), constrained = TRUE, constraint_type = "leq")
theta_base_constrained_geq      <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS, CovUB, x0 = c(3.5, 1.2), constrained = TRUE, constraint_type = "geq")
theta_base_unconstrained        <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS, CovUB, x0 = c(3.5, 1.2), constrained = FALSE)

# Base refit estimates
theta_base_refit_constrained_leq   <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS, CovUB, x0 = c(theta_RBEGLS), constrained = TRUE, constraint_type = "leq")
theta_base_refit_constrained_geq   <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS, CovUB, x0 = c(theta_RBEGLS), constrained = TRUE, constraint_type = "geq")
theta_base_refit_unconstrained     <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS, CovUB, x0 = c(theta_RBEGLS), constrained = FALSE)

# MC estimates
theta_mc_constrained_leq        <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mc, CovUB_mc, x0 = c(3.5, 1.2), constrained = TRUE, constraint_type = "leq")
theta_mc_constrained_geq        <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mc, CovUB_mc, x0 = c(3.5, 1.2), constrained = TRUE, constraint_type = "geq")
theta_mc_unconstrained          <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mc, CovUB_mc, x0 = c(3.5, 1.2), constrained = FALSE)

# MC refit estimates
theta_mc_refit_constrained_leq   <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mc, CovUB_mc, x0 = c(theta_RBEGLS_mc), constrained = TRUE, constraint_type = "leq")
theta_mc_refit_constrained_geq   <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mc, CovUB_mc, x0 = c(theta_RBEGLS_mc), constrained = TRUE, constraint_type = "geq")
theta_mc_refit_unconstrained     <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mc, CovUB_mc, x0 = c(theta_RBEGLS_mc), constrained = FALSE)

# Mean estimates
theta_mean_constrained_leq      <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mean, (CovUB_mc + CovUB)/2, x0 = c(3.5, 1.2), constrained = TRUE, constraint_type = "leq")
theta_mean_constrained_geq      <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mean, (CovUB_mc + CovUB)/2, x0 = c(3.5, 1.2), constrained = TRUE, constraint_type = "geq")
theta_mean_unconstrained        <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mean, (CovUB_mc + CovUB)/2, x0 = c(3.5, 1.2), constrained = FALSE)

# Mean refit estimates
theta_mean_refit_constrained_leq   <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mean, (CovUB_mc + CovUB)/2, x0 = c(theta_mean), constrained = TRUE, constraint_type = "leq")
theta_mean_refit_constrained_geq   <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mean, (CovUB_mc + CovUB)/2, x0 = c(theta_mean), constrained = TRUE, constraint_type = "geq")
theta_mean_refit_unconstrained     <- estimate_params_cov(Dall[idx,idx], CovUBE_GLS_mean, (CovUB_mc + CovUB)/2, x0 = c(theta_mean), constrained = FALSE)

#--------------------------------------------------------------------------------------------------------------
list(
  alpha_est = as.matrix(out1$beta),
  xi_est = as.matrix(out1$theta),
  
  beta_ols = as.matrix(as.numeric(beta_ols)), 
  beta_RBEGLS = as.matrix(as.numeric(beta_RBEGLS)),  
  beta_RBEGLS_mc = as.matrix(as.numeric(beta_RBEGLS_mc)), 
  beta_mean = as.matrix(as.numeric(beta_mean)),
  beta_kr = as.matrix(as.numeric(beta_kr)),
  
  theta_est = as.matrix(as.numeric(theta_RBEGLS)),
  theta_est_mc = as.matrix(as.numeric(theta_RBEGLS_mc)),
  theta_mean = as.matrix(as.numeric(theta_mean)),
  theta_kr = as.matrix(as.numeric(theta_kr)),
  
  # --- Expanded theta_* entries ---
  theta_base_constrained_leq = as.matrix(as.numeric(theta_base_constrained_leq)),
  theta_base_constrained_geq = as.matrix(as.numeric(theta_base_constrained_geq)),
  theta_base_unconstrained = as.matrix(as.numeric(theta_base_unconstrained)),
  theta_base_refit_constrained_leq = as.matrix(as.numeric(theta_base_refit_constrained_leq)),
  theta_base_refit_constrained_geq = as.matrix(as.numeric(theta_base_refit_constrained_geq)),
  theta_base_refit_unconstrained = as.matrix(as.numeric(theta_base_refit_unconstrained)),
  
  theta_mc_constrained_leq = as.matrix(as.numeric(theta_mc_constrained_leq)),
  theta_mc_constrained_geq = as.matrix(as.numeric(theta_mc_constrained_geq)),
  theta_mc_unconstrained = as.matrix(as.numeric(theta_mc_unconstrained)),
  theta_mc_refit_constrained_leq = as.matrix(as.numeric(theta_mc_refit_constrained_leq)),
  theta_mc_refit_constrained_geq = as.matrix(as.numeric(theta_mc_refit_constrained_geq)),
  theta_mc_refit_unconstrained = as.matrix(as.numeric(theta_mc_refit_unconstrained)),
  
  theta_mean_constrained_leq = as.matrix(as.numeric(theta_mean_constrained_leq)),
  theta_mean_constrained_geq = as.matrix(as.numeric(theta_mean_constrained_geq)),
  theta_mean_unconstrained = as.matrix(as.numeric(theta_mean_unconstrained)),
  theta_mean_refit_constrained_leq = as.matrix(as.numeric(theta_mean_refit_constrained_leq)),
  theta_mean_refit_constrained_geq = as.matrix(as.numeric(theta_mean_refit_constrained_geq)),
  theta_mean_refit_unconstrained = as.matrix(as.numeric(theta_mean_refit_unconstrained)),
  
  # --- Variances ---
  varbeta_ols = as.matrix(as.numeric(varbeta_ols)), 
  varbeta_kr = as.matrix(as.numeric(varbeta_kr)),
  varbeta_rb_dm = as.matrix(as.numeric(varbeta_rb_dm)),
  varbeta_rb_mc = as.matrix(as.numeric(varbeta_rb_mc)),
  varbeta_rb_mean = as.matrix(as.numeric(varbeta_rb_mean)),
  varbeta_sdm = as.matrix(as.numeric(varbeta_sdm)),
  varbeta_smc = as.matrix(as.numeric(varbeta_smc)),
  varbeta_smean = as.matrix(as.numeric(varbeta_smean)),
  
  # --- Distance metric ---
  bdist = bdist
)



#}

















