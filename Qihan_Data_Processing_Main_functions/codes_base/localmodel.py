# -*- coding: utf-8 -*-
"""
Created on Wed Apr 24 22:39:12 2024

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


meta = meta_getter('N5236')
trainset = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_25_1_2024.pkl')

index_below_1_BPT = trainset["N2_BPT"] < 1 # BPT cut for Hii regions
trainset = trainset[index_below_1_BPT]
index_above_0_Z = trainset["Z_O3N2"] > 0  ### need to change here!!!
#index_above_0_Z = trainset["Z_RS32"] > 0  ### need to change here!!!
trainset = trainset[index_above_0_Z]
index_above_0_DIG = trainset["S2_DIG"] > 0
trainset = trainset[index_above_0_DIG]


testset = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_25_1_2024.pkl')
index_below_1_BPT = testset["N2_BPT"] > 0 # BPT cut for DIG regions
testset = testset[index_below_1_BPT]
index_above_0_Z = testset["Z_O3N2_kumari_N2"] > 0 
#index_above_0_Z = testset["Z_O3S2_kumari_N2"] > 0 
testset = testset[index_above_0_Z]
index_above_0_DIG = testset["S2_DIG"] > 0
testset = testset[index_above_0_DIG]

# for RS32 only.
#index_not_inf = trainset["Z_RS32"] < 15
#trainset = trainset[index_not_inf]
#index_not_inf = testset["Z_RS32"] < 15
#testset = testset[index_not_inf]
#index_not_inf = testset["Z_O3S2_kumari_N2"] < 15
#testset = testset[index_not_inf]

#plt.scatter(trainset['RA'], trainset['DEC'], s=1, c='red')
#plt.scatter(testset['RA'], testset['DEC'], s=1, c='blue')
#trainset[["SFR_density"]] = np.log(trainset[["SFR_density"]])
#testset[["SFR_density"]] = np.log(testset[["SFR_density"]])


RAtrain = trainset['RA']
DECtrain = trainset['DEC']
RAtest = testset['RA']
DECtest = testset['DEC']
Y_true = testset['Z_O3N2_kumari_N2']
#Y_true = testset['Z_O3S2_kumari_N2']




pred_test_set = list()
k = 50
for ii in range(testset.shape[0]):
    print(ii)
    testpoint = pd.DataFrame(testset).iloc[ii]
    RApoint = pd.DataFrame(RAtest).iloc[ii][0]
    DECpoint = pd.DataFrame(DECtest).iloc[ii][0]
    D_testpoint_to_train = deprojected_distances(RApoint, DECpoint, RA2 = RAtrain, DEC2 = DECtrain)    
    idx_first_k_smallest = sort(np.argpartition(D_testpoint_to_train,k)[0][0:k])
    sur_train = pd.DataFrame(trainset).iloc[idx_first_k_smallest] 
    RAsur_train = sur_train['RA']
    DECsur_train = sur_train['DEC']   
    Dtrain = deprojected_distances(RAsur_train, DECsur_train)
    Dtest = deprojected_distances(RApoint, DECpoint, RA2 = RAsur_train, DEC2 = DECsur_train)   
    X_train = sur_train[["proj_dist", "SFR_density", "S2_DIG"]]
    X_pred = testpoint[["proj_dist", "SFR_density", "S2_DIG"]]   
    Y_train = sur_train["Z_O3N2"] 
    #Y_train = sur_train["Z_RS32"]     
    X_train.insert(0, "intersect", np.ones(X_train.shape[0]), True)
    eta_ini_value = np.array([0.1, 0.0012])
    lower_bound = np.array([1e-5, 1e-10])
    upper_bound = np.array([1,1])  
    MLE_result = MLE_fit(y = Y_train, X = X_train, D = Dtrain, cov_model = "Exp", eta_ini = eta_ini_value, nug = True, opt = "LB", lo_bound = lower_bound, up_bound = upper_bound)  
    thetahat = MLE_result[0] 
    betahat = MLE_result[3]  
    Sigma = thetahat[1]*np.exp(-Dtrain/thetahat[0])    
    Sigma[np.diag_indices_from(Sigma)] = thetahat[1] + thetahat[2]  
    X_pred = insert(matrix(X_pred), 0,1, True)   
    cmat = thetahat[1]*np.exp(-Dtest/thetahat[0])
    pl = X_pred @ betahat
    pe = cmat @ linalg.inv(Sigma) @ (Y_train - X_train @ betahat)
    Y_pred = pl + pe
    pred_test_set.append(Y_pred)


new_set = list()
for j in range(testset.shape[0]):
    new_set.append(float(pd.DataFrame(pred_test_set[j]).values))
    

MSE = mean_squared_error(new_set,Y_true.values) 
RMSE = np.sqrt(MSE)
MAD = statistics.mean(abs(new_set- Y_true.values))

'''
plt.axline((0, 0), slope=1, linestyle="--", color="grey")
plt.hist2d(Y_true.values, new_set, 
           bins = 100,  
           norm = colors.LogNorm(),  
           #cmap =plt.cm.jet, range = [[8.2, 9.6], [8.2, 9.6]]) # plt.cm.BuPu, plt.cm.jet
           cmap =plt.cm.jet, range = [[8.5, 8.9], [8.5, 8.9]])
plt.title('RS32 local Geostatistical model (DIG)') # change name here. + Kumari 19 correction
plt.xlabel('log([O/H])+12 measured')
plt.ylabel('log([O/H])+12 predicted')
#plt.savefig('RS32_local_500.eps', dpi=300)

print(RMSE)
print(MAD)
'''
print("--- %s seconds ---" % (time.time() - start_time))
