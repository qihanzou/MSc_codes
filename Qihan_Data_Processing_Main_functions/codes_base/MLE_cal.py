# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 21:34:29 2024

@author: qihan
"""
from TYPHOON_wrangling import *
import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import pickle
from AST2 import *
from Internal_Fun import *
from BFuns import *
from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
from sklearn import metrics
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
from dec_ra_to_xy import *
import scipy.spatial as spatial
import itertools
import random
from sklearn.metrics import mean_squared_error 
import statistics
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from numba import cuda
import numba
from numba import jit, cuda, njit
from scipy import linalg
import matplotlib.pyplot as plt

meta = meta_getter('N5236')
data_galaxy = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_25_1_2024.pkl')
#gas_density_N2S2Ha = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_inference_gas_density_3664.pkl')
#gas_density_N2S2Ha_w2 = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_inference_gas_density_3664_w2.pkl')
##### BPT cut:
#index_below_1_BPT = data_galaxy["N2_BPT"] < 1 # BPT cut for Hii regions
#data_galaxy = data_galaxy[index_below_1_BPT]

idx_S2_BPT = data_galaxy["S2_BPT"] < 1
data_galaxy = data_galaxy[idx_S2_BPT]

##### CHii cut:
#index_above_CHii = data_galaxy["S2_DIG"] > 0.9 # CHii cut for Hii regions
#data_galaxy = data_galaxy[index_above_CHii]

index_above_0_Z = data_galaxy["Z_N2O2"] > 0  ### need to change here!!! Z_N2S2Ha
data_galaxy = data_galaxy[index_above_0_Z]

index_above_0_DIG = data_galaxy["S2_DIG"] > 0
data_galaxy = data_galaxy[index_above_0_DIG]

#data_galaxy.insert(0, "gas_density_N2S2Ha", gas_density_N2S2Ha.values, True)
#data_galaxy.insert(0, "gas_density_N2S2Ha_w2", gas_density_N2S2Ha.values, True)

#index_above_0_gas = data_galaxy["gas_density_N2S2Ha"] > 0 
#index_above_0_gas = data_galaxy["gas_density_N2S2Ha_w2"] > 0 
#data_galaxy = data_galaxy[index_above_0_gas] #3512



'''
index_below_1_BPT = data_galaxy["N2_BPT"] < 1 # BPT cut for Hii regions
data_galaxy = data_galaxy[index_below_1_BPT]
'''

'''
RA_galaxy_list = data_galaxy['RA']
DEC_galaxy_list = data_galaxy['DEC']
proj_dist_galaxy_list = data_galaxy['proj_dist']
coor_full_galaxy = RA_DEC_to_xy(RA_galaxy_list, DEC_galaxy_list, meta)
coor_full_galaxy = np.transpose(coor_full_galaxy)
X_galaxy = coor_full_galaxy[:,0]
Y_galaxy = coor_full_galaxy[:,1]

plt.scatter(X_galaxy, Y_galaxy, c = 'black', marker='.', s=0.1)
plt.grid()
plt.show()
'''

testset = list()
trainset = list()
random.seed(2026) # seed here!!!
idx = random.sample(range(1, len(data_galaxy)), math.ceil(0.2*len(data_galaxy)))
for kk in range(len(data_galaxy)):
    if kk in idx:
        testset.append(data_galaxy.iloc[[kk]])
    else:
        trainset.append(data_galaxy.iloc[[kk]])
testset = pd.concat(testset)
trainset = pd.concat(trainset)


### LOG SFR here!!!!
trainset[["SFR_density"]] = np.log(trainset[["SFR_density"]])
testset[["SFR_density"]] = np.log(testset[["SFR_density"]])


#X_train = trainset[["proj_dist", "SFR_density", "S2_DIG", "gas_density_N2S2Ha_w2"]]
#X_train = trainset[["proj_dist", "SFR_density", "S2_DIG", "gas_density_N2S2Ha"]]
#X_train = trainset[["proj_dist", "SFR_density", "S2_DIG", "stellar_mass_t16form"]]
X_train = trainset[["proj_dist", "SFR_density", "S2_DIG"]]
#X_train = trainset[["proj_dist", "SFR_density"]]
#X_train = trainset[["proj_dist"]]

#X_pred = testset[["proj_dist", "SFR_density", "S2_DIG", "gas_density_N2S2Ha_w2"]]
#X_pred = testset[["proj_dist", "SFR_density", "S2_DIG", "gas_density_N2S2Ha"]]
#X_pred = testset[["proj_dist", "SFR_density", "S2_DIG", "stellar_mass_t16form"]]
X_pred = testset[["proj_dist", "SFR_density", "S2_DIG"]]
#X_pred = testset[["proj_dist", "SFR_density"]]
#X_pred = testset[["proj_dist"]]




#Y_train = trainset["Z_N2S2Ha"] 
#Y_train = trainset["Z_O3N2"] 
#Y_train = trainset["Z_RS32"] 
Y_train = trainset["Z_N2O2"] 

#Y_true = testset['Z_N2S2Ha']
#Y_true = testset['Z_O3N2']
#Y_true = testset['Z_RS32']
Y_true = testset['Z_N2O2']







X_train.insert(0, "intersect", np.ones(X_train.shape[0]), True)
eta_ini_value = np.array([0.1, 0.0012])
lower_bound = np.array([1e-5, 1e-10])
upper_bound = np.array([1,1])
RAtrain = trainset['RA']
DECtrain = trainset['DEC']
RAtest = testset['RA']
DECtest = testset['DEC']
Dtrain = deprojected_distances(RAtrain, DECtrain)
Dtest = deprojected_distances(RAtest, DECtest, RA2 = RAtrain, DEC2 = DECtrain)


MLE_result = MLE_fit(y = Y_train, X = X_train, D = Dtrain, cov_model = "Exp", eta_ini = eta_ini_value, nug = True, opt = "LB", lo_bound = lower_bound, up_bound = upper_bound)
thetahat = MLE_result[0] # theta_est
betahat = MLE_result[3]  # beta_est
Sigma = thetahat[1]*np.exp(-Dtrain/thetahat[0])
Sigma[np.diag_indices_from(Sigma)] = thetahat[1] + thetahat[2]
X_pred.insert(0, "intersect", np.ones(X_pred.shape[0]), True)
cmat = thetahat[1]*np.exp(-Dtest/thetahat[0])
pl = X_pred @ betahat
pe = cmat @ linalg.inv(Sigma) @ (Y_train - X_train @ betahat)
Y_pred = pl + pe
MSE = mean_squared_error(Y_true.values,Y_pred.values) 
RMSE = np.sqrt(MSE)
MAD = statistics.mean(abs(Y_pred.values - Y_true.values))