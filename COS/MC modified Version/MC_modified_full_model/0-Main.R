rm(list = ls(all = TRUE))
setwd("C:/Users/qihan/Desktop/QH_new")


Nver = 1
beta1 = 4


NALL = c(200, 400, 800, 1600)
for (ni in 1:4){
  nall = NALL[ni]
  idx = 1:(nall/2) 
  source("0-main_TDS.R")
}
