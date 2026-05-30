# # # # # test use:
rm(list = ls(all = TRUE)) 
nall = 200
idx = 1:100
beta1 = 4
# # # # # 
MC_iter = 100
# # # # # 
R_address = "C:/Users/qihan/Desktop/QH_sim2_ab6"


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
W = cbind(1, w)
fols = lm(y~W-1)
beta_ols = fols$coef
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

sdout = SDCal(Dall[-idx,-idx], out1$theta, cov_xi, nug = 0)
num_sim_samples = 10000
samples <- rmvnorm(num_sim_samples, mean = out1$theta, sigma = sdout$covmat)  
samples_in = samples
for (q in 1:dim(samples)[1]){
  while ((samples[q,1] < 0)|(samples[q,2] < 0)) {
    samples[q,] = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  }
}
n_samples <- nrow(samples)
p_samples <- ncol(samples)
Beta <- n_samples/(n_samples-1)^2*mahalanobis(samples, center = colMeans(samples), cov = cov(samples))
F_stat <- ((n_samples-p_samples-1)/p_samples)*(Beta/(1-Beta))
outliers <- which(F_stat > qf(0.05, df1 = p_samples, df2 = n_samples-p_samples-1, lower.tail=FALSE))
n_out = length(outliers)

n_samples1 <- nrow(samples_in)
p_samples1 <- ncol(samples_in)
Beta1 <- n_samples1/(n_samples1-1)^2*mahalanobis(samples_in, center = colMeans(samples_in), cov = cov(samples_in))
F_stat1 <- ((n_samples1-p_samples1-1)/p_samples1)*(Beta1/(1-Beta1))
outliers1 <- which(F_stat1 > qf(0.05, df1 = p_samples1, df2 = n_samples1-p_samples1-1, lower.tail=FALSE))
n_out1 = length(outliers1)

# Visualization
plot(samples_in[,1], samples_in[,2])
points(x = samples_in[,1], y = samples_in[,2],col = "grey")
points(x = samples_in[outliers1,1], y = samples_in[outliers1,2],col = "orange")
points(x = samples[outliers,1], y = samples[outliers,2],col = "red")
points(x = out1$theta[1], y = out1$theta[2], col = "blue", pch = 19)
points(x = xi0[1], y = xi0[2], col = "green", pch = 19)
points(x = colMeans(samples[outliers,])[1], 
       y = colMeans(samples[outliers,])[2], col = "red", pch = 19)
points(x = colMeans(samples_in[outliers1,])[1], 
       y = colMeans(samples_in[outliers1,])[2], col = "yellow", pch = 19)

xi_out = colMeans(samples[outliers,])
xi_sim = colMeans(samples)

xi_out1 = colMeans(samples_in[outliers1,])
xi_sim1 = colMeans(samples_in)

varxi = diag(sdout$covmat)
varxi_out = c(var(samples[outliers,1]), var(samples[outliers,2]))
varxi_sim = c(var(samples[,1]), var(samples[,2]))
#--------------------------------------------------------------------------------------------------------------
list(beta_RBEGLS = beta_RBEGLS,  
     varbeta_rb = as.matrix(varbeta_rb),
     xi_est = as.matrix(out1$theta),
     xi_out = as.matrix(xi_out),
     xi_sim = as.matrix(xi_sim),
     varxi = as.matrix(varxi),
     varxi_out = as.matrix(varxi_out), # should not use this one since the nature of outliers are very far away from each other, therefore the variance will be very large.
     varxi_sim = as.matrix(varxi_sim),
     n_out = n_out
)


#}
hist(samples[outliers,1],breaks=20, freq = FALSE)
density_data <- density(samples[outliers,1])
plot(density_data, main = "Density Plot of Data", xlab = "Value")

hist(samples[outliers,2],breaks=20, freq = FALSE)
density_data <- density(samples[outliers,2])
plot(density_data, main = "Density Plot of Data", xlab = "Value")

hist(samples[,1],breaks=20, freq = FALSE)
hist(samples[,2],breaks=20, freq = FALSE)



