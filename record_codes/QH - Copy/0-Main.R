rm(list = ls(all = TRUE))
setwd("C:/Users/qihan/Desktop/QH - Copy")

#nall = 400 # 1000, 1500, 2000
#idx = sample(1:nall, 100)
Nver = 1
beta1 = 4


#M = (1:6)*100
#N = 200
NALL = c(200, 400, 800, 1600)
for (ni in 1:4){
  nall = NALL[ni]
  idx = 1:(nall/2) #1:45
  source("0-main_TDS.R")
}




