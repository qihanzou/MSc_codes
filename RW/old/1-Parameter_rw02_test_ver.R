# test use:
# nall = 200
# idx = 1:100
# beta1 = 4
# ci_level = 0.95
# 
#
#
# R_address = "C:/Users/qihan/Desktop/QH_rw"

setwd(R_address)
library(mvtnorm)               
source("SimuData.R")
source("BFuns.R")
source("SDCal_modified.R")
source("SDCal_RMLE.R")
source("BFuns_RBEGLS.R")
source("RBEGLS_functions_rw.R")
source("loop_functions02.R")


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
ep_error = epall$Ymat[k,idx]        # error of s1,...,sN
y = beta0 + beta1*x1 + ep_error     # create y values for s1,...,sN



# -------------------------------------------------------------------------- #
out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
alpha_est = out1$beta
cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 

W = cbind(1, w)
fols = lm(y~W-1)
beta_ols = fols$coef
varbeta_ols = diag(vcov(fols))

pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
CovUB = CovU*beta_ols[2]^2 

out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
CovE_GLS = out2$theta[2]*CorE_GLS
CovUBE_GLS = CovUB + CovE_GLS
beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y


# loop start here:
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


# initialization
Var_aplha_est = out1$beta_var
dalpha0 = qnorm(ci_level)*sqrt(diag(Var_aplha_est)[1])
dalpha1 = qnorm(ci_level)*sqrt(diag(Var_aplha_est)[2])
dalpha = cbind(dalpha0, dalpha1)

sdout = SDCal(Dall[-idx,-idx], out1$theta, cov_xi, nug = 0)
dtheta0 = qnorm(ci_level)*sdout$se[1]
dtheta1 = qnorm(ci_level)*sdout$se[2]
dtheta = cbind(dtheta0, dtheta1)

Var_data = Expand_Var(alpha_est[1], alpha_est[2], out1$theta[1], out1$theta[2], dalpha0, dalpha1, dtheta0, dtheta1)
  

RES_set = list()
for (j in 1:dim(Var_data)[1]){
  aaj = as.matrix(as.numeric(cbind(Var_data[j,][1], Var_data[j,][2])))
  ttj = as.matrix(as.numeric(cbind(Var_data[j,][3], Var_data[j,][4])))
  w_j = Comp_W(t(aaj), ttj, r1, r1s, x1s, cov_xi, Dall)
  res_j = loop_funs3(w_j, aaj, ttj, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xi)
  RES_set[[j]] = res_j
  }

Var_Res = NULL
for (k in 1:dim(Var_data)[1]){
  Var_Res = cbind(Var_Res, RES_set[[k]]$varbeta_rw)
}

varbeta_rw_mean = rowMeans(Var_Res)
varbeta_rw_max = c(max(Var_Res[1,]), max(Var_Res[2,]))

varbeta_rw_mean_middle95 = c(mean(sort(Var_Res[1,])[3:79]), mean(sort(Var_Res[2,])[3:79]))
varbeta_rw_mean_middle90 = c(mean(sort(Var_Res[1,])[5:77]), mean(sort(Var_Res[2,])[5:77]))
varbeta_rw_mean_cut_minmax = c(mean(sort(Var_Res[1,])[2:80]), mean(sort(Var_Res[2,])[2:80]))



list(alpha_est = out1$beta,
     xi_est = as.matrix(out1$theta),
     
     beta_ols = as.matrix(as.numeric(beta_ols)), 
     beta_RBEGLS = beta_RBEGLS,  
     theta_est = as.matrix(theta_RBEGLS),
     
     varbeta_ols = as.matrix(as.numeric(varbeta_ols)), 
     varbeta_rb = as.matrix(varbeta_rb),
     varbeta_rw_mean = as.matrix(varbeta_rw_mean),
     varbeta_rw_max = as.matrix(varbeta_rw_max),
     
     varbeta_rw_mean_cut_minmax = as.matrix(varbeta_rw_mean_cut_minmax),
     varbeta_rw_mean_middle95 = as.matrix(varbeta_rw_mean_middle95),
     varbeta_rw_mean_middle90 = as.matrix(varbeta_rw_mean_middle90)
     )


}





