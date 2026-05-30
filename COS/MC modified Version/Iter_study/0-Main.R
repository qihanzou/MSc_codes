rm(list = ls(all = TRUE))
setwd("C:/Users/qihan/Desktop/QH_sim2_ab4_iter_study")


beta1 = 4


#for (MC_I in 1:1){
  #MC_iter = MC_I
  MC_iter = 20
  
  NALL = c(200)#, 400, 800, 1600)
  for (ni in 1:1){
    nall = NALL[ni]
    idx = 1:(nall/2) 
    source("0-main_sim.R")
  }
#}
  
  
  