rm(list = ls(all = TRUE)) # The line removes all variables from the current environment
#setwd("C:/Users/qihan/Desktop/Model")
library("reticulate")
require("reticulate")
library(reticulate)
library(dplyr)
import("scipy")
pd <- import("pandas")
np <- import("numpy")
#source("R_gas_functions.R")
library(pracma)
library(mvtnorm) 
library(geoR)

# ---------------------------------------------------------------------------- #
# Loading the data & Applying diagnostics
# ---------------------------------------------------------------------------- #

data_mgs_2as    <- pd$read_pickle("C:/Users/qihan/Desktop/N1300_mgsd_2as.pkl")
data_mgs_7p5as    <- pd$read_pickle("C:/Users/qihan/Desktop/N1300_mgsd_7p5as.pkl")
data_mgs_11as    <- pd$read_pickle("C:/Users/qihan/Desktop/N1300_mgsd_11as.pkl")
data_mgs_15as    <- pd$read_pickle("C:/Users/qihan/Desktop/N1300_mgsd_15as.pkl")




plot(data_mgs_2as$X, data_mgs_2as$Y, cex = 0.01)
plot(data_mgs_7p5as$X, data_mgs_7p5as$Y, cex = 0.01)
plot(data_mgs_11as$X, data_mgs_11as$Y, cex = 0.01)
plot(data_mgs_15as$X, data_mgs_15as$Y, cex = 0.01)



m1 = mean(data_mgs_2as$CO21_mgsd)
m2 = mean(data_mgs_7p5as$CO21_mgsd)
m3 = mean(data_mgs_11as$CO21_mgsd)
m4 = mean(data_mgs_15as$CO21_mgsd)

plot(c(15,11,7.5,2), c(m4,m3,m2,m1))




hist(data_mgs_2as$CO21_mgsd, breaks = 100, probability = T)
hist(data_mgs_7p5as$CO21_mgsd, breaks = 100, probability = T)
hist(data_mgs_11as$CO21_mgsd, breaks = 100, probability = T)
hist(data_mgs_15as$CO21_mgsd, breaks = 100, probability = T)



data1 = cbind(data_mgs_2as["X"],data_mgs_2as["Y"],data_mgs_2as["CO21_mgsd"])
library(ggplot2)
ggplot(data1, aes(x = X, y = Y, color = CO21_mgsd)) +
  geom_point(size = 1) +
  scale_color_viridis_c() +
  theme_bw()

data2 = cbind(data_mgs_7p5as["X"],data_mgs_7p5as["Y"],data_mgs_7p5as["CO21_mgsd"])
library(ggplot2)
ggplot(data2, aes(x = X, y = Y, color = CO21_mgsd)) +
  geom_point(size = 1) +
  scale_color_viridis_c() +
  theme_bw()

data3 = cbind(data_mgs_11as["X"],data_mgs_11as["Y"],data_mgs_11as["CO21_mgsd"])
library(ggplot2)
ggplot(data3, aes(x = X, y = Y, color = CO21_mgsd)) +
  geom_point(size = 1) +
  scale_color_viridis_c() +
  theme_bw()


