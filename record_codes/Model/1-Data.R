rm(list = ls(all = TRUE)) # The line removes all variables from the current environment
setwd("C:/Users/qihan/Desktop/Model")
library("reticulate")
require("reticulate")
library(reticulate)
library(dplyr)
import("scipy")
pd <- import("pandas")
np <- import("numpy")
source("R_gas_functions.R")
library(pracma)
library(mvtnorm) 
library(geoR)

# ---------------------------------------------------------------------------- #
# Loading the data & Applying diagnostics
# ---------------------------------------------------------------------------- #
data_gal <- pd$read_pickle("MUSE_N1300_ver1_copt.pkl")
data_mgs    <- pd$read_pickle("N1300_mgsd.pkl")


#data_mgs$CO21_mgd = log(data_mgs$CO21_mgd)
data_mgs$CO21_mgsd = log(data_mgs$CO21_mgsd)
# ---------------------------------------------------------------------------- #

idx = is.finite((data_gal$Z_N2S2Ha))
data_gal = data_gal[(idx == 1), ]

# 3 diagnostics for DIG/HII regions

#data_gal = subset(data_gal, data_gal$N2_BPT < 1) 
#data_gal = subset(data_gal, data_gal$S2_BPT < 1) 
data_gal = subset(data_gal, data_gal$S2_DIG > 0.9) 




dim(data_gal)

# ---------------------------------------------------------------------------- #
# Create coordinates:
# ---------------------------------------------------------------------------- #
gal_coor = To_xy_coor(data_gal$RA, data_gal$DEC, RA_gal = 49.921167, 
                        DEC_gal = -19.411361, PA_gal = 278, i_gal = 61.1842,
                        D_gal = 18.99)
data_gal$X <- gal_coor[,1]
data_gal$Y <- gal_coor[,2]

# ---------------------------------------------------------------------------- #
# add an index
# ---------------------------------------------------------------------------- #
data_gal <- cbind(data_gal, index = 1:dim(data_gal)[1])
data_mgs <- cbind(data_mgs, index = 1:dim(data_mgs)[1]) 

# ---------------------------------------------------------------------------- #
# Visualization:
# ---------------------------------------------------------------------------- #
set.seed(1174116)
sample1 = sample.int(n = nrow(data_mgs), size = floor(1000), replace = F)
data_mgs = data_mgs[sample1, ]

set.seed(1174116)
sample2 = sample.int(n = nrow(data_gal), size = floor(1000), replace = F)
data_gal = data_gal[sample2, ]

#plot(data_gal$X, data_gal$Y, cex = 0.3, col = "red")
plot(data_mgs$X, data_mgs$Y, cex = 0.3)
points(data_gal$X, data_gal$Y, cex = 0.3, col = "red")


# a = 8
# b = 6
# c = -6
# d = -7.5
# Localdata = c()
# for (i in 1:dim(data_gal)[1]){
#   if (data_gal[i,]$X < a){
#     if (data_gal[i,]$X > b){
#       if (data_gal[i,]$Y < c){
#         if (data_gal[i,]$Y > d){
#           Localdata = rbind(Localdata, data_gal[i,])
#         }
#       }
#     }
#   }
# }
# 
# Localgas = c()
# for (i in 1:dim(data_mgs)[1]){
#   if (data_mgs[i,]$X < a){
#     if (data_mgs[i,]$X > b){
#       if (data_mgs[i,]$Y < c){
#         if (data_mgs[i,]$Y > d){
#           Localgas = rbind(Localgas, data_mgs[i,])
#         }
#       }
#     }
#   }
# }
# # ---------------------------------------------------------------------------- #
# # Visualization
# # ---------------------------------------------------------------------------- #
# plot(Localgas$X, Localgas$Y, cex = 0.3)
# points(Localdata$X, Localdata$Y, cex = 0.3, col = "red", pch = 16)
# 
# set.seed(2024)
# sample1 = sample.int(n = nrow(Localgas), size = floor(0.3*nrow(Localgas)), replace = F)
# Localgas2 = Localgas[sample1, ]
# 
# set.seed(2024)
# sample2 = sample.int(n = nrow(Localdata), size = floor(0.085*nrow(Localdata)), replace = F)
# Localdata2 = Localdata[sample2, ]
# 
# plot(Localgas2$X, Localgas2$Y, cex = 0.3)
# points(Localdata2$X, Localdata2$Y, cex = 0.3, col = "red", pch = 16)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
Local_gas  = data_mgs
Local_data = data_gal
dim(Local_gas)
dim(Local_data)

save(Local_gas, Local_data, file = "1-Data_N1300_N2S2Ha_S2DIG_rs1000.Rdata")





















