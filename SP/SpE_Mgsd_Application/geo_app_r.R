
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

# ---------------------------------------------------------------------------- #
# Loading the data & Applying diagnostics
# ---------------------------------------------------------------------------- #
data_N5236 <- pd$read_pickle("N5236_25_1_2024.pkl")
data_mgs    <- pd$read_pickle("N5236_mgs_data.pkl")
data_mgs$CO21_mgd = log(data_mgs$CO21_mgd)

idx = is.finite(data_N5236$Z_N2S2Ha)
data_N5236 = data_N5236[(idx == 1), ]
data_N5236 = subset(data_N5236, data_N5236$N2_BPT < 1) 

# ---------------------------------------------------------------------------- #
# Data reduction: Choose 2% of whole data
# ---------------------------------------------------------------------------- #
set.seed(2024)
sample = sample.int(n = nrow(data_mgs), size = floor(0.02*nrow(data_mgs)), replace = F)
d_test = data_mgs[sample, ]

# ---------------------------------------------------------------------------- #
# Create coordinates:
# ---------------------------------------------------------------------------- #
N5236_coor = To_xy_coor(data_N5236$RA, data_N5236$DEC, RA_gal = 204.2516509, 
                        DEC_gal = -29.86547, PA_gal = 54, i_gal = 12.50496329,
                        D_gal = 4.89778819)
data_N5236$X <- N5236_coor[,1]
data_N5236$Y <- N5236_coor[,2]

# ---------------------------------------------------------------------------- #
# add an index
# ---------------------------------------------------------------------------- #
data_N5236 <- cbind(data_N5236, index = 1:dim(data_N5236)[1])
data_mgs <- cbind(data_mgs, index = 1:dim(data_mgs)[1]) 

# ---------------------------------------------------------------------------- #
# Visualization:
# ---------------------------------------------------------------------------- #
plot(data_N5236$X, data_N5236$Y, cex = 0.3, col = "red")
points(d_test$X, d_test$Y, cex = 0.3)


# ---------------------------------------------------------------------------- #
# Applying restriction
# ---------------------------------------------------------------------------- #
a = 5.5
b = -5.5
c = 6
d = -6
Localdata = c()
for (i in 1:dim(data_N5236)[1]){
  if (data_N5236[i,]$X < a){
    if (data_N5236[i,]$X > b){
      if (data_N5236[i,]$Y < c){
        if (data_N5236[i,]$Y > d){
          Localdata = rbind(Localdata, data_N5236[i,])
        }
      }
    }
  }
}

Localgas = c()
for (i in 1:dim(d_test)[1]){
  if (d_test[i,]$X < a){
    if (d_test[i,]$X > b){
      if (d_test[i,]$Y < c){
        if (d_test[i,]$Y > d){
          Localgas = rbind(Localgas, d_test[i,])
        }
      }
    }
  }
}


load("1-Data.Rdata")

# ---------------------------------------------------------------------------- #
# Visualization
# ---------------------------------------------------------------------------- #
plot(Localgas$X, Localgas$Y, cex = 0.3)
points(Localdata$X, Localdata$Y, cex = 0.3, col = "red", pch = 16)

plot(Local_gas$X, Local_gas$Y, cex = 0.3, pch = 16)
points(Local_data$X, Local_data$Y, cex = 0.3, col = "red", pch = 16)


# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
Local_gas  = Localgas
Local_data = Localdata
dim(Local_gas)
dim(Local_data)

D_gas = deprojected_distances(Local_gas$RA, Local_gas$DEC, Local_gas$RA, Local_gas$DEC)
D_data = deprojected_distances(Local_data$RA, Local_data$DEC, Local_data$RA, Local_data$DEC)
D_data_to_gas = deprojected_distances(Local_data$RA, Local_data$DEC, Local_gas$RA, Local_gas$DEC)
y_gas = as.matrix(select(Local_gas, c("CO21_mgd")))
X_gas = as.matrix(cbind(1, select(Local_gas, c("proj_dist"))))

gas_coor = cbind(Local_gas$X, Local_gas$Y)
v_gas = cbind(gas_coor, y_gas, X_gas)
V_gas = as.geodata(v_gas, coords.col = 1:2, data.col = 3, covar.col = 4:5, covar.names = c("int", "dist")) 
vario = variog(V_gas, max.dist=2, trend = ~ 1 + dist)
plot(vario,xlab = "h", ylab = "gamma(h)")


# ---------------------------------------------------------------------------- #
# Applying Universal Kriging
# ---------------------------------------------------------------------------- #
source_python("BFuns.py")

res1 = MLE_fit_gas(y = y_gas, X = X_gas, D = D_gas, cov_model = "Exp")
res1

thetahat1 = res1[[1]] 
betahat1  = res1[[4]] 
Sigma1 = thetahat1[2]*exp(-D_gas/thetahat1[1]) 
diag(Sigma1) = thetahat1[2] + thetahat1[3]
d0 = cbind(1, Local_data$proj_dist)
cmat1 = thetahat1[2]*exp(-D_data_to_gas/thetahat1[1]) 
pl1 = as.matrix(d0)%*%betahat1 
pe1 = cmat1%*% solve(Sigma1)%*%(y_gas - X_gas%*%betahat1)
w = pl1 + pe1
# ---------------------------------------------------------------------------- #
# Create data with the predicted w
# ---------------------------------------------------------------------------- #
new_data = Local_data
new_data$w = w




# ---------------------------------------------------------------------------- #
# Training set and testing set
# ---------------------------------------------------------------------------- #
set.seed(2024)
idx = sample.int(n = nrow(new_data), size = floor(0.2*nrow(new_data)), replace = F)
test= new_data[idx, ]
train = new_data[-idx, ] 
Wtrain = as.matrix(cbind(1, select(train, c("w", "proj_dist"))))
Ytrain = as.matrix(select(train, c("Z_N2S2Ha")))


# ---------------------------------------------------------------------------- #
# Beta OLS:
# ---------------------------------------------------------------------------- #
beta_ols = solve(t(Wtrain)%*%Wtrain)%*%t(Wtrain)%*%Ytrain



# ---------------------------------------------------------------------------- #
# KR method
# ---------------------------------------------------------------------------- #
res2 = MLE_fit(y = Ytrain, X = Wtrain, D = D_data[-idx,-idx], cov_model = "Exp")
res2

# ---------------------------------------------------------------------------- #
# KR RMSE and MAD
# ---------------------------------------------------------------------------- #
thetahat2 = res2[[1]] # theta_est
betahat2  = res2[[4]] # beta_est
Sigma2 = thetahat2[2]*exp(-D_data[-idx,-idx]/thetahat2[1]) 
d1 = cbind(1, test$w, test$proj_dist)
cmat2 = thetahat2[2]*exp(-D_data[idx,-idx]/thetahat2[1]) 
pl2 = as.matrix(d1)%*%betahat2       
pe2 = cmat2%*% solve(Sigma2)%*%(Ytrain - Wtrain%*%betahat2) 
pred_Y = pl2 + pe2
(RMSE_geo_dist_kr <- sqrt(mean((pred_Y - test$Z_N2S2Ha)^2)))
(MAD_geo_dist_kr <- mean(abs(pred_Y - test$Z_N2S2Ha)))

# ---------------------------------------------------------------------------- #
# OLS RMSE and MAD
# ---------------------------------------------------------------------------- #
d1 = cbind(1, test$w, test$proj_dist)
pols = as.matrix(d1)%*%beta_ols
(RMSE_geo_dist_ols <- sqrt(mean((pols - test$Z_N2S2Ha)^2)))
(MAD_geo_dist_ols <- mean(abs(pols - test$Z_N2S2Ha)))


# ---------------------------------------------------------------------------- #
# Variance for thetas
# ---------------------------------------------------------------------------- #
source("SDCal_modified.R")
sdout1 = SDCal(D_gas, thetahat1, 'Exp', nug = 1)
sdout2 = SDCal(D_data[-idx,-idx], thetahat2, 'Exp', nug = 1)
# ---------------------------------------------------------------------------- #



# ---------------------------------------------------------------------------- #
# Compute CovU
# ---------------------------------------------------------------------------- #
inv_Sigma = solve(Sigma1)
Cz0 = (thetahat1[2] + thetahat1[3])*exp(-D_data/thetahat1[1]) 
pc = as.matrix(d0) - cmat1%*%inv_Sigma%*%X_gas
CovU = Cz0 - cmat1%*%inv_Sigma%*%t(cmat1) + (pc)%*%solve(t(X_gas)%*%inv_Sigma%*%X_gas)%*%t(pc)

source_python("BFuns.py")
source_python("Internal_Fun.py")
# ---------------------------------------------------------------------------- #
# RBEGLS method and iterations
# ---------------------------------------------------------------------------- #
res3 = MLE_fit_GLS(y = Ytrain, X = Wtrain, D = D_data[-idx,-idx], kSigma = CovU[-idx,-idx]*beta_ols[2]^2, cov_model = "Exp")
res3
thetahat3 = res3[[1]] # theta_est
CovE_GLS = thetahat3[2]*exp(-D_data[-idx,-idx]/thetahat3[1]) 
CovUBE_GLS = CovU[-idx,-idx]*beta_ols[2]^2 + CovE_GLS
beta_RBEGLS_initial = solve(t(Wtrain)%*%solve(CovUBE_GLS)%*%Wtrain)%*%t(Wtrain)%*%solve(CovUBE_GLS)%*%Ytrain


# ---------------------------------------------------------------------------- #
# 1st iteration
# ---------------------------------------------------------------------------- #
res4 = MLE_fit_GLS(y = Ytrain, X = Wtrain, D = D_data[-idx,-idx], kSigma = CovU[-idx,-idx]*beta_RBEGLS_initial[2]^2, cov_model = "Exp")
res4

thetahat4 = res4[[1]] # theta_est
CovE_GLS4 = thetahat4[2]*exp(-D_data[-idx,-idx]/thetahat4[1]) 
CovUBE_GLS4 = CovU[-idx,-idx]*beta_RBEGLS_initial[2]^2 + CovE_GLS4
beta_RBEGLS_4 = solve(t(Wtrain)%*%solve(CovUBE_GLS4)%*%Wtrain)%*%t(Wtrain)%*%solve(CovUBE_GLS4)%*%Ytrain

# ---------------------------------------------------------------------------- #
# 2nd iteration
# ---------------------------------------------------------------------------- #
res5 = MLE_fit_GLS(y = Ytrain, X = Wtrain, D = D_data[-idx,-idx], kSigma = CovU[-idx,-idx]*beta_RBEGLS_4[2]^2, cov_model = "Exp")
res5

thetahat5 = res5[[1]] # theta_est
CovE_GLS5 = thetahat5[2]*exp(-D_data[-idx,-idx]/thetahat5[1]) 
CovUBE_GLS5 = CovU[-idx,-idx]*beta_RBEGLS_4[2]^2 + CovE_GLS5
beta_RBEGLS_5 = solve(t(Wtrain)%*%solve(CovUBE_GLS5)%*%Wtrain)%*%t(Wtrain)%*%solve(CovUBE_GLS5)%*%Ytrain

# ---------------------------------------------------------------------------- #
# 3rd iteration
# ---------------------------------------------------------------------------- #
res6 = MLE_fit_GLS(y = Ytrain, X = Wtrain, D = D_data[-idx,-idx], kSigma = CovU[-idx,-idx]*beta_RBEGLS_5[2]^2, cov_model = "Exp")
res6

thetahat6 = res6[[1]] # theta_est
CovE_GLS6 = thetahat6[2]*exp(-D_data[-idx,-idx]/thetahat6[1]) 
CovUBE_GLS6 = CovU[-idx,-idx]*beta_RBEGLS_5[2]^2 + CovE_GLS6
beta_RBEGLS_6 = solve(t(Wtrain)%*%solve(CovUBE_GLS6)%*%Wtrain)%*%t(Wtrain)%*%solve(CovUBE_GLS6)%*%Ytrain

# ---------------------------------------------------------------------------- #
# RBEGLS beta variance
# ---------------------------------------------------------------------------- #
CovUB6 = CovU[-idx,-idx]*beta_RBEGLS_6[2]^2 
CovUBE_GLS6 = CovUB6 + CovE_GLS6
A1 = solve(t(Wtrain)%*%solve(CovUBE_GLS6)%*%Wtrain)
varbeta_cl0 = diag(A1)

# ---------------------------------------------------------------------------- #
# RBEGLS RMSE and MAD
# ---------------------------------------------------------------------------- #
d1 = cbind(1, test$w, test$proj_dist)
cmat6 = thetahat6[2]*exp(-D_data[idx,-idx]/thetahat6[1]) 
pl6 = as.matrix(d1)%*%beta_RBEGLS_6
pe6 = cmat6%*% solve(CovUBE_GLS6)%*%(Ytrain - Wtrain%*%beta_RBEGLS_6) 
pred6 = pl6 + pe6
(RMSE_geo_dist_6 <- sqrt(mean((pred6 - test$Z_N2S2Ha)^2)))
(MAD_geo_dist_6 <- mean(abs(pred6 - test$Z_N2S2Ha)))

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #
inv_CovUBE_GLS6 = solve(CovUBE_GLS6)
Cz6 =  thetahat6[2]*exp(-D_data[idx,idx]/thetahat6[1]) 
pc6 = as.matrix(d1) - cmat6%*%inv_CovUBE_GLS6%*%Wtrain
Cov6 = Cz6 - cmat6%*%inv_CovUBE_GLS6%*%t(cmat6) + (pc6)%*%solve(t(Wtrain)%*%inv_CovUBE_GLS6%*%Wtrain)%*%t(pc6)
CI6 = cbind(pred6-2.576*sqrt(diag(Cov6)), pred6 + 2.576*sqrt(diag(Cov6)))
v = 0
for (s in 1:length(test$Z_N2S2Ha)){
  if (test$Z_N2S2Ha[s] > CI6[s,][1] ){
    if (test$Z_N2S2Ha[s] < CI6[s,][2]){
      v = v + 1
    }
  }else{
    v = v + 0
  }
}
print(length(test$Z_N2S2Ha))
print(v)
print(v/length(test$Z_N2S2Ha))








# ---------------------------------------------------------------------------- #
# Proposed method 1:
# ---------------------------------------------------------------------------- #

z0 = beta_RBEGLS_6[1]  + beta_RBEGLS_6[2]*betahat1[1]
z1 = beta_RBEGLS_6[2]*betahat1[2]
zeta = t(cbind(z0,z1))

Sigm = Cz0[-idx, -idx]*beta_RBEGLS_6[2]^2 + CovE_GLS6



source_python("BFuns.py")
source_python("Internal_Fun.py")

out_tau = MLE_fit_CS_gau(Sigm, D_data[-idx,-idx], zeta, Ytrain, d0[-idx,])

out_tau2 = MLE_fit_CS(Sigm, D_data[-idx,-idx], zeta, Ytrain, d0[-idx,])


CovT  = out_tau[[1]][2]*exp(-D_data[-idx,-idx]/out_tau[[1]][1])

CovT = out_tau2[[1]][1]*matrix(1, dim(D_data[-idx,-idx])[1], dim(D_data[-idx,-idx])[2])

CovUBE_GLS6T = CovUBE_GLS6 + CovT
A1_pro = solve(t(Wtrain)%*%solve(CovUBE_GLS6T)%*%Wtrain)
varbeta_pro = diag(A1_pro)



d1 = cbind(1, test$w, test$proj_dist)
cmat6T = thetahat6[2]*exp(-D_data[idx,-idx]/thetahat6[1])  + out_tau2[[1]][1]*matrix(1, dim(D_data[idx,-idx])[1], dim(D_data[idx,-idx])[2])
pl6T = as.matrix(d1)%*%beta_RBEGLS_6


pe6T = cmat6T%*% solve(CovUBE_GLS6T)%*%(Ytrain - Wtrain%*%beta_RBEGLS_6) 
pred6T = pl6T + pe6T
(RMSE_geo_dist_6T <- sqrt(mean((pred6T - test$Z_N2S2Ha)^2)))
(MAD_geo_dist_6T <- mean(abs(pred6T - test$Z_N2S2Ha)))

inv_CovUBE_GLS6T = solve(CovUBE_GLS6T)
Cz6T =  thetahat6[2]*exp(-D_data[idx,idx]/thetahat6[1]) + out_tau2[[1]][1]*matrix(1, dim(D_data[idx,idx])[1], dim(D_data[idx,idx])[2])

pc6T = as.matrix(d1) - cmat6T%*%inv_CovUBE_GLS6T%*%Wtrain
Cov6T = Cz6T - cmat6T%*%inv_CovUBE_GLS6T%*%t(cmat6T) + (pc6T)%*%solve(t(Wtrain)%*%inv_CovUBE_GLS6T%*%Wtrain)%*%t(pc6T)


CI6T = cbind(pred6T-2.576*sqrt(diag(Cov6T)), pred6T + 2.576*sqrt(diag(Cov6T)))
# 2.576, 1.96, 1.645
v = 0
for (s in 1:length(test$Z_N2S2Ha)){
  if (test$Z_N2S2Ha[s] > CI6T[s,][1] ){
    if (test$Z_N2S2Ha[s] < CI6T[s,][2]){
      v = v + 1
    }
  }else{
    v = v + 0
  }
}
print(length(test$Z_N2S2Ha))
print(v)
print(v/length(test$Z_N2S2Ha))

# ---------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------- #






























