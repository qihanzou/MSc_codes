# test use:
rm(list = ls(all = TRUE)) 
nall = 200
idx = 1:100
beta1 = 4
ci_level = 0.95
MC_iter = 100


R_address = "C:/Users/qihan/Desktop/QH_rw"

setwd(R_address)
library(mvtnorm)               
source("SimuData.R")
source("BFuns.R")
source("SDCal_modified.R")
source("SDCal_RMLE.R")
source("BFuns_RBEGLS.R")
source("RBEGLS_functions_rw.R")
source("loop_functions_rw06.R")


if (as.numeric(exists("update_est")) == 1  ){
  remove(update_est)
} 
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
# w
out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
alpha_est = out1$beta
cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 

#--------------------------------------------------------------------------------------------------------------
# OLS
W = cbind(1, w)
fols = lm(y~W-1)
beta_ols = fols$coef
varbeta_ols = diag(vcov(fols))

#--------------------------------------------------------------------------------------------------------------
# KR
kr = MLE.fit(y, W, Dall[idx,idx], cov_ep, c(2.5), nug, "LB", lo.bound,up.bound)
theta_kr = kr[[1]] 
beta_kr  = kr[[4]] 
varbeta_kr = diag(kr[[5]])

#--------------------------------------------------------------------------------------------------------------
# RBEGLS
pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
CovUB = CovU*beta_ols[2]^2 

out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
CovE_GLS = out2$theta[2]*CorE_GLS
CovUBE_GLS = CovUB + CovE_GLS
beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y


# loop:
res1 = RBEGLS_loop_mod2(50, 0.001, 0.001, beta_RBEGLS, out2$theta, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
beta_RBEGLS = res1$beta_est
theta_RBEGLS  = res1$theta_est
CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
CovE_GLS = theta_RBEGLS[2]*CorE_GLS
CovUB = CovU*beta_RBEGLS[2]^2
CovUBE_GLS = CovUB + CovE_GLS
CovUBE_GLS_inv = solve(CovUBE_GLS)
A1 = solve(t(W)%*%CovUBE_GLS_inv%*%W)
varbeta_rb = diag(A1)
#--------------------------------------------------------------------------------------------------------------
# MC 

sdout = SDCal(Dall[-idx,-idx], out1$theta, cov_xi, nug = 0)
ETAallboot_all = SimuData(out1$theta, rate, locall, cov_xi, nug, MC_iter)
EVar = EVarKR = EVarKRT = matrix(0, 2, 2)

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
  Wboot_e = cbind(1, wboot_e)
  
  A1_e = solve(t(Wboot_e)%*%CovUBE_GLS_inv%*%Wboot_e)
  EVar = EVar + A1_e                      
}
varbeta_mc = diag(EVar/MC_iter)

#--------------------------------------------------------------------------------------------------------------
# RW

dtheta0 = qnorm(ci_level)*sdout$se[1]
dtheta1 = qnorm(ci_level)*sdout$se[2]
dtheta = cbind(dtheta0, dtheta1)

Var_data = Expand_theta(out1$theta[1], out1$theta[2], dtheta0, dtheta1)

cov_inv = solve(CovUBE_GLS)

RES_set = list()
Var_simple_set = list()
for (j in 1:dim(Var_data)[1]){
  aaj = as.matrix(as.numeric(cbind(alpha_est[1], alpha_est[2])))
  ttj = as.matrix(as.numeric(cbind(Var_data[j,][1], Var_data[j,][2])))
  w_j = Comp_W(t(aaj), ttj, r1, r1s, x1s, cov_xi, Dall)
  
  res_j = loop_funs4(w_j, aaj, ttj, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xi, cov_inv)
  RES_set[[j]] = res_j
  
  W_j = cbind(1, w_j)
  Var_simple_set[[j]] = diag(solve(t(W_j)%*%CovUBE_GLS_inv%*%W_j))
  }

Var_Res = NULL
Var_Res_v2 = NULL
Var_Res_v3 = NULL
beta_Res = NULL
theta_Res = NULL
alpha_Res = NULL
xi_Res = NULL
Var_simple = NULL
for (k in 1:dim(Var_data)[1]){
  Var_Res = cbind(Var_Res, RES_set[[k]]$varbeta_rw)
  Var_Res_v2 = cbind(Var_Res_v2, RES_set[[k]]$varbeta_rw_v2)
  Var_Res_v3 = cbind(Var_Res_v3, RES_set[[k]]$varbeta_rw_v3)
  Var_simple = cbind(Var_simple, Var_simple_set[[k]])
  
  beta_Res = cbind(beta_Res, RES_set[[k]]$beta)
  theta_Res = cbind(theta_Res, RES_set[[k]]$theta)
  alpha_Res = cbind(alpha_Res, RES_set[[k]]$alpha)
  xi_Res = cbind(xi_Res, RES_set[[k]]$xi)
}

beta_rw = rowMeans(beta_Res)
theta_rw = rowMeans(theta_Res)
alpha_rw = rowMeans(alpha_Res)
xi_rw = rowMeans(xi_Res)

varbeta_rw_simple = rowMeans(Var_simple)



varbeta_rw_mean = rowMeans(Var_Res)
varbeta_rw_v2_mean = rowMeans(Var_Res_v2)
varbeta_rw_v3_mean = rowMeans(Var_Res_v3)

#--------------------------------------------------------------------------------------------------------------
list(alpha_est = out1$beta,
     alpha_rw = as.matrix(alpha_rw),
     xi_est = as.matrix(out1$theta),
     xi_rw = as.matrix(xi_rw),
     
     beta_ols = as.matrix(as.numeric(beta_ols)), 
     beta_RBEGLS = beta_RBEGLS,  
     beta_kr = beta_kr,
     beta_rw = as.matrix(beta_rw),
     
     theta_est = as.matrix(theta_RBEGLS),
     theta_kr = as.matrix(theta_kr),
     theta_rw = as.matrix(theta_rw),
     
     varbeta_ols = as.matrix(as.numeric(varbeta_ols)), 
     varbeta_kr = as.matrix(varbeta_kr),
     varbeta_rb = as.matrix(varbeta_rb),
     varbeta_mc = as.matrix(varbeta_mc),
     
     varbeta_rw_mean = as.matrix(varbeta_rw_mean),
     
     varbeta_rw_v2_mean = as.matrix(varbeta_rw_v2_mean),
     
     varbeta_rw_v3_mean = as.matrix(varbeta_rw_v3_mean),
     
     varbeta_rw_simple = as.matrix(varbeta_rw_simple)
)


#}





