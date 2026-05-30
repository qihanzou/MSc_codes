# test use:
rm(list = ls(all = TRUE)) 
nall = 200
idx = 1:100
beta1 = 4
# 
MC_iter = 100
N_sim_out = 10000
N_sim_in = 10000
# 
R_address = "C:/Users/qihan/Desktop/QH_rw"


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
res1 = RBEGLS_loop_final(50, 0.001, 0.001, beta_RBEGLS, out2$theta, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
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
XI1_list = c()
XI2_list = c()
for (v in 1:MC_iter){
  
  x1allboot = r1all%*%alpha_est + ETAallboot_all$Ymat[v,]
  x1sboot = x1allboot[-idx]
  
  xie = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  while ((xie[1] < 0)|(xie[2] < 0)) {
    xie = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  }
  cov_xialle = xie[2]*cor.mat(Dall, xie[1], cov_xi, nug = 0)
  XI1_list[v] = xie[1]
  XI2_list[v] = xie[2]
  tmp_inv = solve(cov_xialle[-idx,-idx])
  alphaboot = solve(t(r1s) %*% tmp_inv %*% r1s)%*%t(r1s) %*% tmp_inv %*% x1sboot 
  wboot_e = r1%*%alphaboot + cov_xialle[idx, -idx]%*%tmp_inv%*%(x1sboot - r1s%*%alphaboot)
  Wboot_e = cbind(1, wboot_e)
  
  A1_e = solve(t(Wboot_e)%*%CovUBE_GLS_inv%*%Wboot_e)
  EVar = EVar + A1_e                      
}
varbeta_mc = diag(EVar/MC_iter)

# --------------------------------------------------------------------------- #
# 4m method (min_max_max_min):

max0 = matrix(NA, nrow=N_sim_out, ncol=1)
max1 = matrix(NA, nrow=N_sim_out, ncol=1)
min0 = matrix(NA, nrow=N_sim_out, ncol=1)
min1 = matrix(NA, nrow=N_sim_out, ncol=1)

for (t in 1:N_sim_out) {
  valid_samples <- matrix(NA, nrow=N_sim_in, ncol=2)
  count <- 0
  
  while (count < N_sim_in) {
    remaining <- N_sim_in - count
    samples <- rmvnorm(remaining, mean = out1$theta, sigma = sdout$covmat)  
    valid <- rowSums(samples >= 0) == 2  # Check non-negative
    
    num_valid <- sum(valid)
    if (num_valid > 0) {
      valid_samples[(count + 1):(count + num_valid), ] <- samples[valid, ]
      count <- count + num_valid
    }
  }
  
  theta_sim <- valid_samples
  theta_sim0 = theta_sim[,1]
  theta_sim1 = theta_sim[,2]
  
  max0[t] = max(theta_sim0)
  max1[t] = max(theta_sim1)
  min0[t] = min(theta_sim0)
  min1[t] = min(theta_sim1)
  
}


min_max0 = min(max0)
min_max1 = min(max1)
max_min0 = max(min0)
max_min1 = max(min1)

mm_sim_t0 = cbind(min_max0, max_min0, out1$theta[1])
mm_sim_t1 = cbind(min_max1, max_min1, out1$theta[2])
mm_sim_data = expand.grid(mm_sim_t0, mm_sim_t1)


mm_sim_RES_set = list()
for (j3 in 1:dim(mm_sim_data)[1]){
  mm_sim_ttj = as.matrix(as.numeric(cbind(mm_sim_data[j3,][1], mm_sim_data[j3,][2])))
  mm_sim_w_j = Comp_W(t(alpha_est), mm_sim_ttj, r1, r1s, x1s, cov_xi, Dall)
  mm_sim_RES_set[[j3]] = loop_funs_final(mm_sim_w_j, mm_sim_ttj, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xi, CovUBE_GLS_inv)
}


mm_Var_Res_sim = NULL
for (k3 in 1:dim(mm_sim_data)[1]){
  mm_Var_Res_sim = cbind(mm_Var_Res_sim, mm_sim_RES_set[[k3]]$varbeta_rw_v3)
}

varbeta_4m_mean = rowMeans(mm_Var_Res_sim)
varbeta_4m_median = apply(mm_Var_Res_sim, 1, FUN = median)
varbeta_4m_robloc = apply(mm_Var_Res_sim, 1, FUN = robLoc)
varbeta_4m_skew = apply(mm_Var_Res_sim, 1, FUN = skewness)


varbeta_4m_huber.mu = apply(mm_Var_Res_sim, 1, FUN = huber.mu)
varbeta_4m_huber.one.step = apply(mm_Var_Res_sim, 1, FUN = huber.one.step)
varbeta_4m_huber = as.numeric(cbind(apply(mm_Var_Res_sim, 1, FUN = huber)[[1]]$mu, apply(mm_Var_Res_sim, 1, FUN = huber)[[2]]$mu))


AA = trimmed_mean(mm_Var_Res_sim[1,], trim = 0.5 - 0.5*exp(-abs(varbeta_4m_skew[1])))
BB = trimmed_mean(mm_Var_Res_sim[2,], trim = 0.5 - 0.5*exp(-abs(varbeta_4m_skew[2])))
varbeta_4m_tri = as.numeric(cbind(AA, BB))



#--------------------------------------------------------------------------------------------------------------
list(alpha_est = out1$beta,
     xi_est = as.matrix(out1$theta),
     
     beta_ols = as.matrix(as.numeric(beta_ols)), 
     beta_RBEGLS = beta_RBEGLS,  
     beta_kr = beta_kr,
     
     theta_est = as.matrix(theta_RBEGLS),
     theta_kr = as.matrix(theta_kr),
     
     varbeta_ols = as.matrix(as.numeric(varbeta_ols)), 
     varbeta_kr = as.matrix(varbeta_kr),
     varbeta_rb = as.matrix(varbeta_rb),
     varbeta_mc = as.matrix(varbeta_mc),
     varbeta_4m_mean = as.matrix(varbeta_4m_mean),
     varbeta_4m_median = as.matrix(varbeta_4m_median),
     varbeta_4m_robloc = as.matrix(varbeta_4m_robloc),
     varbeta_4m_huber.mu = as.matrix(varbeta_4m_huber.mu),
     varbeta_4m_huber.one.step = as.matrix(varbeta_4m_huber.one.step),
     varbeta_4m_huber = as.matrix(varbeta_4m_huber),
     varbeta_4m_tri = as.matrix(varbeta_4m_tri),
     varbeta_4m_skew = as.matrix(varbeta_4m_skew),
     mm_sim_data = mm_sim_data
)


#}


