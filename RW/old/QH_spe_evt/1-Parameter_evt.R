# test use:
rm(list = ls(all = TRUE)) 
nall = 200
idx = 1:100
beta1 = 4

MC_iter = 100

N_sim_out = 10000
N_sim_in = 10000


R_address = "C:/Users/qihan/Desktop/QH_rw"

setwd(R_address)
library(mvtnorm)               
source("SimuData.R")
source("BFuns.R")
source("BFuns_RBEGLS.R")
source("SDCal_modified.R")
source("RBEGLS_functions_final.R")
library(extRemes)

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



# v3
ci_level = c(0.925, 0.95, 0.975)
V_list = matrix(0, 3, 2)
for (l in 1:3){
  dtheta0 = qnorm(ci_level[l])*sdout$se[1]
  dtheta1 = qnorm(ci_level[l])*sdout$se[2]
  Var_data = Expand_theta_final(out1$theta[1], out1$theta[2], dtheta0, dtheta1)
  
  RES_set = list()
  for (j in 1:dim(Var_data)[1]){
    ttj = as.matrix(as.numeric(cbind(Var_data[j,][1], Var_data[j,][2])))
    w_j = Comp_W(t(alpha_est), ttj, r1, r1s, x1s, cov_xi, Dall)
    RES_set[[j]] = loop_funs_final(w_j, ttj, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xi, CovUBE_GLS_inv)
  }
  
  
  Var_Res_v3 = NULL
  for (k in 1:dim(Var_data)[1]){
    Var_Res_v3 = cbind(Var_Res_v3, RES_set[[k]]$varbeta_rw_v3)
  }
  
  varbeta_v3 = rowMeans(Var_Res_v3)
  V_list[l,] = varbeta_v3
}
varbeta_v3_zmean859095 = colMeans(V_list)
varbeta_v3_zmean9095 = colMeans(V_list[2:3,])
varbeta_v3_zmean8590 = colMeans(V_list[1:2,])
# ---------------------------------------------------------------------------- #


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
  
  # --------------------------------------------------------------------------- #
  
  # min_max_max_min method:
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
  
  mm_varbeta_evt = rowMeans(mm_Var_Res_sim)
  # --------------------------------------------------------------------------- #
  
  # library(climextRemes)
  # fit_gev(max0)
  # fit_gev(max1)

  # fit_max0 = fevd(max0,type="GEV", method = "Lmoments")
  # fit_max1 = fevd(max1,type="GEV", method = "Lmoments")
  # fit_min0 = fevd(-min0,type="GEV", method = "Lmoments")
  # fit_min1 = fevd(-min1,type="GEV", method = "Lmoments")
  # 
  # 
  # # # Mean:
  # # mean_up_theta0 = fit_max0$results[1] + fit_max0$results[2]*(gamma(1 - 1*fit_max0$results[3]) - 1)/fit_max0$results[3]
  # # mean_down_theta0 = -fit_min0$results[1] + fit_min0$results[2]*(gamma(1 - 1*fit_min0$results[3]) - 1)/fit_min0$results[3]
  # # mean_up_theta1 = fit_max1$results[1] + fit_max1$results[2]*(gamma(1 - 1*fit_max1$results[3]) - 1)/fit_max1$results[3]
  # # mean_down_theta1 = -fit_min1$results[1] + fit_min1$results[2]*(gamma(1 - 1*fit_min1$results[3]) - 1)/fit_min1$results[3]
  # # 
  # # mean_sim_t0 = cbind(mean_up_theta0, mean_down_theta0)
  # # mean_sim_t1 = cbind(mean_up_theta1, mean_down_theta1)
  # # mean_sim_data = expand.grid(mean_sim_t0, mean_sim_t1)
  # 
  # # Median:
  # median_up_theta0 =  fit_max0$results[1] + fit_max0$results[2]*((log(2))^(-fit_max0$results[3]) - 1)/fit_max0$results[3]
  # median_down_theta0 =  -fit_min0$results[1] + fit_min0$results[2]*((log(2))^(-fit_min0$results[3]) - 1)/fit_min0$results[3]
  # median_up_theta1 =  fit_max1$results[1] + fit_max1$results[2]*((log(2))^(-fit_max1$results[3]) - 1)/fit_max1$results[3]
  # median_down_theta1 =  -fit_min1$results[1] + fit_min1$results[2]*((log(2))^(-fit_min1$results[3]) - 1)/fit_min1$results[3]
  # 
  # med_sim_t0 = cbind(median_up_theta0, median_down_theta0)
  # med_sim_t1 = cbind(median_up_theta1, median_down_theta1)
  # med_sim_data = expand.grid(med_sim_t0, med_sim_t1)
  # 
  # # # Location:
  # # sim_set_theta0 = cbind(fit_max0$results[1], -fit_min0$results[1])
  # # sim_set_theta1 = cbind(fit_max1$results[1], -fit_min1$results[1])
  # # sim_set_data = expand.grid(sim_set_theta0, sim_set_theta1)
  # 
  # sim_set_data = med_sim_data
  # # -----------------------------------------------------------------------------------#
  # sim_RES_set = list()
  # for (j2 in 1:dim(sim_set_data)[1]){
  #   sim_ttj = as.matrix(as.numeric(cbind(sim_set_data[j2,][1], sim_set_data[j2,][2])))
  #   sim_w_j = Comp_W(t(alpha_est), sim_ttj, r1, r1s, x1s, cov_xi, Dall)
  #   sim_RES_set[[j2]] = loop_funs_final(sim_w_j, sim_ttj, y, r1, r1s, Dall, cov_ep, theta_ep.ini, nug, lo.bounde, up.bounde, cov_xi, CovUBE_GLS_inv)
  # }
  # 
  # Var_Res_sim = NULL
  # for (k2 in 1:dim(sim_set_data)[1]){
  #   Var_Res_sim = cbind(Var_Res_sim, sim_RES_set[[k2]]$varbeta_rw_v3)
  # }
  # 
  # varbeta_evt = rowMeans(Var_Res_sim)




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
     varbeta_v3_85 = as.matrix(V_list[1,]),
     varbeta_v3_90 = as.matrix(V_list[2,]),
     varbeta_v3_95 = as.matrix(V_list[3,]),
     varbeta_v3_zmean9095 = as.matrix(varbeta_v3_zmean9095),
     varbeta_v3_zmean8590 = as.matrix(varbeta_v3_zmean8590),
     varbeta_v3_zmean859095 = as.matrix(varbeta_v3_zmean859095),
     mm_varbeta_evt = as.matrix(mm_varbeta_evt),
     mm_sim_data = mm_sim_data
)


#}

print(varbeta_kr)
print(varbeta_mc)
print(t(as.matrix(V_list[2,])))
print(mm_varbeta_evt)

Var_data
mm_sim_data
