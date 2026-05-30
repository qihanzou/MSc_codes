# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 19:51:56 2024

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
import pandas as pd
import numpy as np
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_score
from sklearn import linear_model, tree, ensemble
start_time = time.time()


meta = meta_getter('N5236') #N5236_26_04
trainset = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_25_1_2024.pkl') #orginal!!!


##### BPT cut:
#index_below_1_BPT = trainset["N2_BPT"] < 1 # BPT cut for Hii regions
#trainset = trainset[index_below_1_BPT]

##### CHii cut:
index_above_CHii = trainset["S2_DIG"] > 0.9 # CHii cut for Hii regions
trainset = trainset[index_above_CHii]

#index_above_0_Z = trainset["Z_N2S2Ha"] > 0  ### need to change here!!!
#trainset = trainset[index_above_0_Z]

#index_above_0_Z = trainset["Z_O3N2"] > 0  ### need to change here!!!
#trainset = trainset[index_above_0_Z]

#index_above_0_Z = trainset["Z_RS32"] > 0  ### need to change here!!!
#trainset = trainset[index_above_0_Z]

#index_above_0_DIG = trainset["S2_DIG"] > 0
#trainset = trainset[index_above_0_DIG]

#X_train= trainset[["proj_dist", "SFR_density", "S2_DIG"]]
X1 = trainset[["proj_dist"]]
Y1 = trainset["Z_N2S2Ha"] 



kf =KFold(n_splits=10, shuffle=True, random_state=2024)

cnt = 1
# split()  method generate indices to split data into training and test set.
for train_index, test_index in kf.split(X1, Y1):
    print(f'Fold:{cnt}, Train set: {len(train_index)}, Test set:{len(test_index)}')
    cnt += 1


trainset = trainset.reset_index(drop=True)
Yo_list = list()
Yf_list = list()
for train_index, test_index in kf.split(X1, Y1):
    X_train, Y_train = trainset["proj_dist"][train_index], trainset["Z_N2S2Ha"][train_index]
    X_test, Y_test = trainset["proj_dist"][test_index], trainset["Z_N2S2Ha"][test_index]
    
    RAtrain = trainset['RA'][train_index]
    DECtrain = trainset['DEC'][train_index]
    RAtest = trainset['RA'][test_index]
    DECtest = trainset['DEC'][test_index]
    Dtest = deprojected_distances(RAtest, DECtest, RA2 = RAtrain, DEC2 = DECtrain)
    
    X_train = pd.DataFrame(X_train)
    X_train.insert(0, "intersect", np.ones(X_train.shape[0]), True)
    
    X_test = pd.DataFrame(X_test)
    X_test.insert(0, "intersect", np.ones(X_test.shape[0]), True)
    
    eta_ini_value = np.array([0.1, 0.0012])
    lower_bound = np.array([1e-5, 1e-10])
    upper_bound = np.array([1,1])
    Dtrain = deprojected_distances(RAtrain, DECtrain)
    MLE_result = MLE_fit(y = Y_train, X = X_train, D = Dtrain, cov_model = "Exp", eta_ini = eta_ini_value, nug = True, opt = "LB", lo_bound = lower_bound, up_bound = upper_bound)
    thetahat = MLE_result[0] # theta_est
    betahat = MLE_result[3]  # beta_est
    Sigma = thetahat[1]*np.exp(-Dtrain/thetahat[0])
    Sigma[np.diag_indices_from(Sigma)] = thetahat[1] + thetahat[2]

    cmat = thetahat[1]*np.exp(-Dtest/thetahat[0])
    pl = X_test @ betahat
    pe = cmat @ linalg.inv(Sigma) @ (Y_train - X_train @ betahat)
    Y_pred = pl + pe
    Yf_list.extend(Y_pred)
    Yo_list.extend(Y_test)
    










plt.axline((0, 0), slope=1, linestyle="--", color="grey")
plt.hist2d(Yo_list, Yf_list, 
           bins = 100,  
           norm = colors.LogNorm(),  
           cmap = 'plasma', range = [[8.2, 9.7], [8.2, 9.7]])
plt.title('O3N2 Geo-model with Distance (HII)') # change name here. + Kumari 19 correction
plt.xlabel('log([O/H])+12 measured')
plt.ylabel('log([O/H])+12 predicted')

plt.savefig('pic_for_Ben_N2S2Ha_prediction_for_HII_S2DIG_O3N2.png', dpi=300)



#plt.scatter(Y_pred.values, Y_train.values - Y_pred.values)





