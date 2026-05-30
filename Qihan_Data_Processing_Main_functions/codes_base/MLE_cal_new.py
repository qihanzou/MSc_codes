# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 17:54:54 2024

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
from dec_ra_to_xy import *
import scipy.spatial as spatial
import itertools
import random
from sklearn.metrics import mean_squared_error 
import statistics
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.stats import linregress
import time
start_time = time.time()


meta = meta_getter('N5236') #N5236_26_04
trainset = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_25_1_2024.pkl') #orginal!!!


##### BPT cut:
#index_below_1_BPT = trainset["N2_BPT"] < 1 # BPT cut for Hii regions
#trainset = trainset[index_below_1_BPT]




##### CHii cut:
index_above_CHii = trainset["S2_DIG"] > 0.9 # CHii cut for Hii regions
trainset = trainset[index_above_CHii]

index_above_0_Z = trainset["Z_N2S2Ha"] > 0  ### need to change here!!!
#index_above_0_Z = trainset["Z_O3N2"] > 0  ### need to change here!!!
#index_above_0_Z = trainset["Z_RS32"] > 0  ### need to change here!!!
#index_above_0_Z = trainset["Z_N2O2"] > 0  ### need to change here!!!
trainset = trainset[index_above_0_Z]


#index_above_0_DIG = trainset["S2_DIG"] > 0
#trainset = trainset[index_above_0_DIG]



testset = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_25_1_2024.pkl') # orginal !!!

#index_below_1_BPT = testset["N2_BPT"] > 0 # BPT cut for DIG regions
#testset = testset[index_below_1_BPT]



##### CHii cut:
index_above_CHii = testset["S2_DIG"] <= 0.9 # CHii cut for DIG regions
testset = testset[index_above_CHii]


#index_above_0_Z = testset["Z_O3S2_kumari_N2"] > 0 
#index_above_0_Z = testset["Z_O3N2_kumari_N2"] > 0 
index_above_0_Z = testset["Z_N2S2Ha"] > 0  ### need to change here!!!
#index_above_0_Z = testset["Z_O3N2"] > 0  ### need to change here!!!
#index_above_0_Z = testset["Z_RS32"] > 0  ### need to change here!!!
#index_above_0_Z = testset["Z_N2O2"] > 0  ### need to change here!!!
testset = testset[index_above_0_Z]

#index_above_0_DIG = testset["S2_DIG"] > 0
#testset = testset[index_above_0_DIG]






# for RS32 only.
#index_not_inf = trainset["Z_RS32"] < 15
#trainset = trainset[index_not_inf]
#index_not_inf = testset["Z_RS32"] < 15
#testset = testset[index_not_inf]
#index_not_inf = testset["Z_O3S2_kumari_N2"] < 15
#testset = testset[index_not_inf]




#X_train= trainset[["proj_dist", "SFR_density", "S2_DIG"]]
#X_train= trainset[["proj_dist", "SFR_density"]]
X_train = trainset[["proj_dist"]]



#X_pred = testset[["proj_dist", "SFR_density", "S2_DIG"]]
#X_pred = testset[["proj_dist", "SFR_density"]]
X_pred = testset[["proj_dist"]]


Y_train = trainset["Z_N2S2Ha"] 
#Y_train = trainset["Z_O3N2"] 
#Y_train = trainset["Z_RS32"] 
#Y_train = trainset["Z_N2O2"] 

Y_true = testset['Z_N2S2Ha']
#Y_true = testset['Z_O3N2']
#Y_true = testset['Z_RS32']
#Y_true = testset['Z_N2O2']
#Y_true = testset['Z_O3N2_kumari_N2']
#Y_true = testset['Z_O3S2_kumari_N2']



### For general:
RAtrain = trainset['RA']
DECtrain = trainset['DEC']
RAtest = testset['RA']
DECtest = testset['DEC']
Dtest = deprojected_distances(RAtest, DECtest, RA2 = RAtrain, DEC2 = DECtrain)
X_train.insert(0, "intersect", np.ones(X_train.shape[0]), True)
eta_ini_value = np.array([0.1, 0.0012])
lower_bound = np.array([1e-5, 1e-10])
upper_bound = np.array([1,1])
Dtrain = deprojected_distances(RAtrain, DECtrain)
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




print("--- %s seconds ---" % (time.time() - start_time))

'''
plt.axline((0, 0), slope=1, linestyle="--", color="grey")
plt.hist2d(Y_true.values, Y_pred.values, 
           bins = 100,  
           norm = colors.LogNorm(),  
           #cmap =plt.cm.jet, range = [[8.2, 9.6], [8.2, 9.6]]) # plt.cm.BuPu, plt.cm.jet
           cmap =plt.cm.jet, range = [[8.5, 8.9], [8.5, 8.9]])
plt.title('D+SFR Geo-model (DIG) + K19 ') # change name here. + Kumari 19 correction
plt.xlabel('log([O/H])+12 measured')
plt.ylabel('log([O/H])+12 predicted')
#plt.savefig('s2_k19_20_RS32_DIG_dist_sfr.eps', dpi=300)
'''
