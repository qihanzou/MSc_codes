# test use:
rm(list = ls(all = TRUE)) 
nall = 200
idx = 1:100
beta1 = 4
# 
MC_iter = 100
num_sim_samples = 100
# 
R_address = "C:/Users/qihan/Desktop/QH_sim"


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
# 
samples <- rmvnorm(num_sim_samples, mean = out1$theta, sigma = sdout$covmat)  
for (q in 1:dim(samples)[1]){
  while ((samples[q,1] < 0)|(samples[q,2] < 0)) {
    samples[q,] = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  }
}

# ---------------------------------------------------------------------------- #
# https://stackoverflow.com/questions/45289225/removing-multivariate-outliers-with-mvoutlier
# ---------------------------------------------------------------------------- #

n_samples <- nrow(samples)
p_samples <- ncol(samples)
Beta <- n_samples/(n_samples-1)^2*mahalanobis(samples, center = colMeans(samples), cov = cov(samples))
F_stat <- ((n_samples-p_samples-1)/p_samples)*(Beta/(1-Beta))
outliers <- which(F_stat > qf(0.05, df1 = p_samples, df2 = n_samples-p_samples-1, lower.tail=FALSE))
n_out = length(outliers)

sample_lib = cbind(samples, F_stat, mahalanobis(samples, center = colMeans(samples), cov = cov(samples)))
sample_lib_out = sample_lib[outliers,]

# ---------------------------------------------------------------------------- #


plot(samples[,1], samples[,2], xlab = "Phi", ylab  ="Sigma2")
points(x = samples[outliers,1], y = samples[outliers,2],col = "red")
points(x = out1$theta[1], y = out1$theta[2], col = "blue", pch = 19)
points(x = xi0[1], y = xi0[2], col = "green", pch = 19)
points(x = colMeans(samples[outliers,])[1], 
       y = colMeans(samples[outliers,])[2], col = "red", pch = 19)



EVar_sim  = EVarGLS_sim = EVar_sim2 = NULL

for (c in 1:dim(samples[outliers,])[1]){

  xi_sim = samples[outliers,][c,]
  cov_xiall_sim = xi_sim[2]*cor.mat(Dall, xi_sim[1], cov_xi, nug = 0)
  tmp_inv_sim = solve(cov_xiall_sim[-idx,-idx])
  
  alpha_j = solve(t(r1s) %*% tmp_inv_sim %*% r1s)%*%t(r1s) %*% tmp_inv_sim %*% x1s 
  
  w_sim = r1%*%alpha_j + cov_xiall_sim[idx, -idx]%*%tmp_inv_sim%*%(x1s - r1s%*%alpha_j)
  W_sim = cbind(1, w_sim)
  
  A1_sim = solve(t(W_sim)%*%CovUBE_GLS_inv%*%W_sim)
  EVar_sim = cbind(EVar_sim, diag(A1_sim))
  
  fols_j = lm(y~W_sim-1)
  beta_ols_j = fols_j$coef
  
  pc_j = r1 - cov_xiall_sim[idx, -idx]%*%solve(cov_xiall_sim[-idx,-idx]) %*% r1s
  CovU_j = cov_xiall_sim[idx, idx] - cov_xiall_sim[idx, -idx]%*%solve(cov_xiall_sim[-idx,-idx]) %*% t(cov_xiall_sim[idx, -idx]) + (pc_j)%*%solve(t(r1s)%*%solve(cov_xiall_sim[-idx,-idx])%*%r1s)%*%t(pc_j)
  CovUB_j = CovU_j*beta_ols_j[2]^2
  
  out_gls = MLE.fit_GLS(y, W_sim, Dall[idx,idx], CovUB_j, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CovE_GLS_j = out_gls$theta[2]*cor.mat(Dall[idx,idx], out_gls$theta[1], cov_ep, nug = 0)
  
  CovUBE_GLS_j = CovUB_j + CovE_GLS_j 
  
  B2_e = t(W_sim)%*%CovUBE_GLS_inv%*%CovUBE_GLS_j%*%CovUBE_GLS_inv%*%W_sim
  varbeta_sim_gls = diag(A1_sim%*%B2_e%*%A1_sim)               
  EVarGLS_sim = cbind(EVarGLS_sim, varbeta_sim_gls)
  
  A2_sim = solve(t(W_sim)%*%solve(CovUBE_GLS_j)%*%W_sim)
  EVar_sim2 = cbind(EVar_sim2, diag(A2_sim))
  
}
varbeta_A1 = rowMeans(EVar_sim)
varbeta_A2 = rowMeans(EVar_sim2)

varbeta_gls_median = apply(EVarGLS_sim, 1, FUN = median)
varbeta_gls_mean = rowMeans(EVarGLS_sim)
varbeta_gls_max = apply(EVarGLS_sim, 1, FUN = max)
varbeta_gls_min = apply(EVarGLS_sim, 1, FUN = min)
varbeta_gls_ave_max_min = (varbeta_gls_max+varbeta_gls_min)/2

q1 = quantile(EVarGLS_sim[1,], probs=c(.025,.975))
q2 = quantile(EVarGLS_sim[2,], probs=c(.025,.975))
index1 = (EVarGLS_sim[1,]<q1[1]) | (EVarGLS_sim[1,]>q1[2])
index2 = (EVarGLS_sim[2,]<q2[1]) | (EVarGLS_sim[2,]>q2[2])
m1 = mean(EVarGLS_sim[1,][index1])
m2 = mean(EVarGLS_sim[2,][index2])
varbeta_gls_95 = as.numeric(cbind(m1,m2))

M1 = mean(sort(EVarGLS_sim[1,])[2:(dim(EVarGLS_sim)[2]-1)])
M2 = mean(sort(EVarGLS_sim[2,])[2:(dim(EVarGLS_sim)[2]-1)])
varbeta_gls_r2m = as.numeric(cbind(M1,M2))

#weighted.mean(0)
#t(as.matrix(EVarGLS_sim))
#sample_lib_out
ss1 = as.numeric(sort(EVarGLS_sim[1,]))
ss2 = as.numeric(sort(EVarGLS_sim[2,]))



weight_sample = seq(1:dim(EVarGLS_sim)[2])/sum(seq(1:dim(EVarGLS_sim)[2]))
wmean1 = weighted.mean(ss1, weight_sample)
wmean2 = weighted.mean(ss2, weight_sample)
varbeta_wmean = as.numeric(cbind(wmean1, wmean2))


exp_weight = exp(seq(1:dim(EVarGLS_sim)[2]))/sum(exp(seq(1:dim(EVarGLS_sim)[2])))
ewmean1 = weighted.mean(ss1, exp_weight)
ewmean2 = weighted.mean(ss2, exp_weight)
varbeta_ewmean = as.numeric(cbind(ewmean1, ewmean2))
#--------------------------------------------------------------------------------------------------------------
list(alpha_est = out1$beta,
     xi_est = as.matrix(out1$theta),
     
     beta_ols = as.matrix(as.numeric(beta_ols)), 
     beta_RBEGLS = beta_RBEGLS,  
     beta_kr = beta_kr,
     
     theta_est = as.matrix(theta_RBEGLS),
     theta_kr = as.matrix(theta_kr),
     
     samples_outliers_mean = colMeans(samples[outliers,]),
     num_outliers = dim(samples[outliers,])[1],
     theta_gls = as.matrix(out_gls$theta),
     
     varbeta_ols = as.matrix(as.numeric(varbeta_ols)), 
     varbeta_kr = as.matrix(varbeta_kr),
     varbeta_rb = as.matrix(varbeta_rb),
     varbeta_mc = as.matrix(varbeta_mc),
     
     varbeta_wmean = as.matrix(varbeta_wmean),
     varbeta_ewmean = as.matrix(varbeta_ewmean),
     varbeta_gls_mean = as.matrix(varbeta_gls_mean),
     varbeta_gls_r2m = as.matrix(varbeta_gls_r2m),
     varbeta_A1 = as.matrix(varbeta_A1),
     varbeta_A2 = as.matrix(varbeta_A2),
     
     varbeta_gls_median = as.matrix(varbeta_gls_median),
     varbeta_gls_max = as.matrix(varbeta_gls_max),
     varbeta_gls_min = as.matrix(varbeta_gls_min),
     varbeta_gls_ave_max_min = as.matrix(varbeta_gls_ave_max_min),
     varbeta_gls_95q = as.matrix(varbeta_gls_95)
)


#}

# seq(1:100)/sum(seq(1:100))
# 
# 
plot(seq(1:5)/sum(seq(1:5)), xlab = "Order index", ylab = "Weight")
plot(seq(1:10)/sum(seq(1:10)), xlab = "Order index", ylab = "Weight")
plot(seq(1:50)/sum(seq(1:50)), xlab = "Order index", ylab = "Weight")
plot(seq(1:25)/sum(seq(1:25)), xlab = "Order index", ylab = "Weight")
