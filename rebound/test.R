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
pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)

res1 = RBEGLS_loop_final(50, 0.001, 0.001, beta_ols, theta_ep.ini, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
beta_RBEGLS = res1$beta_est
theta_RBEGLS  = res1$theta_est
CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
CovE_GLS = theta_RBEGLS[2]*CorE_GLS
CovUB = CovU*beta_RBEGLS[2]^2
CovUBE_GLS = CovUB + CovE_GLS
CovUBE_GLS_inv = solve(CovUBE_GLS)
A1 = solve(t(W)%*%CovUBE_GLS_inv%*%W)
varbeta_rb = diag(A1)


cov_xiall = xi0[2]*cor.mat(Dall, xi0[1], cov_xi, nug = 0)
w = r1%*%alpha0 + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha0) 
W = cbind(1, w)

res1 = RBEGLS_loop_final(50, 0.001, 0.001, beta_ols, theta_ep.ini, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
beta_RBEGLS = res1$beta_est
theta_RBEGLS  = res1$theta_est
CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
CovE_GLS = theta_RBEGLS[2]*CorE_GLS
CovUB = CovU*beta_RBEGLS[2]^2
CovUBE_GLS = CovUB + CovE_GLS
CovUBE_GLS_inv = solve(CovUBE_GLS)
A1 = solve(t(W)%*%CovUBE_GLS_inv%*%W)
varbeta_rb2 = diag(A1)


varbeta_rb
varbeta_rb2

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


