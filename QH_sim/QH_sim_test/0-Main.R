rm(list = ls(all = TRUE))
setwd("C:/Users/qihan/Desktop/QH_sim_test")

nug = 0

beta0 = 2
beta1 = 4
alpha0 = c(1,0.5)

MC_iter = 100
num_sim_samples = 200


xi0 = c(0.01, 0.01)  # initial values of (phi, sigma^2)
cov_xi = "Mat32"
xi.ini = c(0.008)
up.bound = c(1)
lo.bound = c(0.005)

cov_ep = "Exp"
ep_theta0 = c(2,1) 
theta_ep.ini = c(3.5,1.2)
up.bounde = c(20, 20)
lo.bounde = c(0.01, 0.01)



NALL = c(200) #, 400, 800, 1600)
for (ni in 1:1){
  nall = NALL[ni]
  idx = 1:(nall/2) 
  source("0-main_sim.R")
}







