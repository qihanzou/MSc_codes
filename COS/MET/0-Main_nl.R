rm(list = ls(all = TRUE))
setwd("C:/Users/qihan/Desktop/QH_Nonlinear")


beta1 = 4

MC_iter = 100


NALL = c(200) #, 400, 800, 1600)

for (ni in 1:1){
  nall = NALL[ni]
  idx = 1:(nall/2) 
  source("0-main_sim_nl.R")
}
