# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 13:17:02 2024

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



def calculate_gas_density_CO21(ICO21, meta):
    '''
    parameters:
        1. R21 = 0.65, 1/R21*alpha_Co(1-0) gives constant around 6.7. 
        adopted CO(2-1) to CO(1-0) line ratio (den Brok et al. 2021; Leroy et al. 2022)
        2. alpha_Co(1-0)=4.35, unit: M_sun pc-2 (K km s-1)^-1
        above are coefficents. 
        3. ICO21, Integrated CO(2-1) line intensity (moment-0). unit: K km s−1
        4. cosi, inclination correction.
    return:
        1. molecular gas surface density (unit: M_sun pc-2)
    '''
    i  =  np.radians(meta['i']) 
    gas_density = 6.7*ICO21.data*np.cos(i) 
    return gas_density # unit: M_sun pc-2


def create_regions_from_whole(num_parts_axis, min_x, min_y, max_x, max_y, total_length, data):
    '''
    *parameters:
    1. num_parts_axis: number of parts for both axis: eg: 1, 2, 4, 6, 12, 24
    2. min_x: minimum value for x axis
    3. min_y: minimum value for y axis
    4. total_length: total length of axis, need to be same for x y for now.
    
    *results:
        1. local_df: regions for the whole data set, if we choose num_parts_axis = 12, then we
        will get 144 regions. This df contains data for different variables.
        2. xy_boundary: the min x, max x, min y, max y values for every regions.
    
    *use the following codes:
    local_df, xy_boundary = create_regions_from_whole(num_parts_axis = 12, min_x = -6, min_y = -6, total_length = 12, data = data)
    '''
    M = num_parts_axis
    local_df = list()
    xy_boundary = list()
    xy_center_coor = list()
    add = total_length/M # add segments
    for jj in range(M):
        for ii in range(M):
           region_test11 = data[(data["X"] > min_x + ii*add)]
           region_test12 = region_test11[(region_test11["X"] < min_x + (ii + 1)*add)]
           region_test13 = region_test12[(region_test12["Y"] > min_y + jj*add)]
           region_test14 = region_test13[(region_test13["Y"] < min_y + (jj + 1)*add)]
           max_min_vec = np.array([min_x + ii*add, min_x + (ii + 1)*add, min_y + jj*add, min_y + (jj + 1)*add])  
           x_y_center_coor_vec = np.array([((min_x + ii*add) + (min_x + (ii + 1)*add))/2, ((min_y + jj*add) + (min_y + (jj + 1)*add))/2])
           local_df.append(region_test14)  # store "regions" data
           xy_boundary.append(max_min_vec) # store corresponding x and y boundaries.
           xy_center_coor.append(x_y_center_coor_vec) 
    print(len(local_df))    # length check
    print(len(xy_boundary)) # length check
    print(len(xy_center_coor))
    return local_df, xy_boundary, xy_center_coor



def test_performance_of_models_based_on_regions(local_df, precentage_for_test_set):
    '''
    parameters:
        1. local_df: data frames for regions
        2. precentage_for_test_set: eg: 0.2, cut 0.2 points as test set, others be tranning set.
        
    returns:
        1. MSE_list: a list for MSEs (for each regions, combine together)
        2. MAD_list: a list for MADs (for each regions, combine together)
        3. RMSE_list: a list for RMSEs (for each regions, combine together)
        
    test_performance_of_models_based_on_regions(local_df, 0.2)
    '''
    MSE_list = list()
    MAD_list = list()
    RMSE_list = list()
    for qq in range(len(local_df)):
        print(qq)
        local_now = local_df[qq]
        check = local_now.empty
        if check == False:
        # checking
           testset = list()
           trainset = list()
           idx = random.sample(range(1, len(local_now)), math.ceil(precentage_for_test_set*len(local_now)))
           for kk in range(len(local_now)):
               if kk in idx:
                   testset.append(local_now.iloc[[kk]])
               else:
                   trainset.append(local_now.iloc[[kk]])
           testset = pd.concat(testset)
           trainset = pd.concat(trainset)
           #X_train = trainset[["X", "Y", "proj_dist"]]
           X_train = trainset[["proj_dist"]]
           #X_train = trainset[["X", "Y"]]
           X_train.insert(0, "intersect", np.ones(X_train.shape[0]), True)
           Y_train = trainset["CO21_gas_density"] 
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
           #X_pred = testset[["X", "Y", "proj_dist"]]
           X_pred = testset[["proj_dist"]]
           #X_pred = testset[["X", "Y"]]
           X_pred.insert(0, "intersect", np.ones(X_pred.shape[0]), True)
           cmat = thetahat[1]*np.exp(-Dtest/thetahat[0])
           pl = X_pred @ betahat
           pe = cmat @ np.linalg.inv(Sigma) @ (Y_train - X_train @ betahat)
           Y_pred = pl + pe
           Y_true = testset['CO21_gas_density']
           MSE = mean_squared_error(Y_true.values,Y_pred.values) 
           RMSE = np.sqrt(MSE)
           MAD = statistics.mean(abs(Y_pred.values - Y_true.values))
           MSE_list.append(MSE)
           MAD_list.append(MAD)
           RMSE_list.append(RMSE)
        else:
           print('DataFrame is empty!')
    print(statistics.mean(np.sqrt(MSE_list)))
    print(statistics.mean(MSE_list))
    print(statistics.mean(MAD_list))
    return MSE_list, MAD_list, RMSE_list


def construct_geo_dist_models_for_regions(local_df):
    '''
    parameters:
        1. local_df
        
    returns:
        1. thetahat_list
        2. betahat_list
        3. Sigma_list
        4. error_list
        5. inv_Sigma_list
        * for faster.
        
    '''
    thetahat_list = list()      
    betahat_list = list()
    Sigma_list = list()  
    error_list = list() 
    inv_Sigma_list = list()
    for zz in range(len(local_df)):
        print(zz)
        local_now = local_df[zz]
        check = local_now.empty
        if check == False:
           #X_train = local_now[["X", "Y", "proj_dist"]]
           X_local_now = local_now[["proj_dist"]]
           #X_train = local_now[["X", "Y"]]
           X_local_now.insert(0, "intersect", np.ones(X_local_now.shape[0]), True)
           Y_local_now = local_now["CO21_gas_density"] 
           Y_local_now = np.log(Y_local_now) ############### take log here, may work? 
           # take log to solve the problem: some prediction values are nagative if we 
           # do not take log. after take log, do exp again, results should always positive.
           eta_ini_value = np.array([0.1, 0.0012])
           lower_bound = np.array([1e-5, 1e-10])
           upper_bound = np.array([1,1])
           RAlocal_now = local_now['RA']
           DEClocal_now = local_now['DEC']
           Dlocal_now = deprojected_distances(RAlocal_now, DEClocal_now)
           MLE_result = MLE_fit(y = Y_local_now, X = X_local_now, D = Dlocal_now, cov_model = "Exp", eta_ini = eta_ini_value, nug = True, opt = "LB", lo_bound = lower_bound, up_bound = upper_bound)
           thetahat = MLE_result[0] # theta_est
           betahat = MLE_result[3]  # beta_est
           Sigma = thetahat[1]*np.exp(-Dlocal_now/thetahat[0])
           Sigma[np.diag_indices_from(Sigma)] = thetahat[1] + thetahat[2]
           inv_Sigma = linalg.inv(Sigma)
           # store thetahat, betahat, Sigma, error
           error = Y_local_now - X_local_now @ betahat
           thetahat_list.append(thetahat)
           betahat_list.append(betahat)
           Sigma_list.append(Sigma)
           error_list.append(error)
           inv_Sigma_list.append(inv_Sigma)
        else:
           thetahat_list.append([])
           betahat_list.append([])
           Sigma_list.append([])
           error_list.append([])
           inv_Sigma_list.append([])
           print('DataFrame is empty!')
    return thetahat_list, betahat_list, Sigma_list, error_list, inv_Sigma_list



def pred_gas_density(error_list, Sigma_list, betahat_list, thetahat_list, xy_boundary, local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x, max_y, min_x, min_y):
    inference_gas_density = list()
    for gg in range(data_galaxy.shape[0]):
        print(gg)
        x_point = X_galaxy[gg]
        y_point = Y_galaxy[gg]
        RA_point = RA_galaxy_list.values[gg]
        DEC_point = DEC_galaxy_list.values[gg]
        X_infer = np.array(proj_dist_galaxy_list.values[gg])
        X_infer = insert(X_infer,0,1)
        if (x_point > max_x) or (x_point < min_x) or (y_point > max_y) or (y_point < min_y):
           gas_density_pred = np.nan
           inference_gas_density.append(gas_density_pred)
        else:
           for t in range(len(xy_boundary)):
               xy_region_now = xy_boundary[t]
               x_min_now = xy_region_now[0]
               x_max_now = xy_region_now[1]
               y_min_now = xy_region_now[2]
               y_max_now = xy_region_now[3]
               if (x_min_now <= x_point <= x_max_now) and (y_min_now <= y_point <= y_max_now):
                  index_xy_region = t
           local_now = local_df[index_xy_region]
           thetahat = thetahat_list[index_xy_region]
           betahat = betahat_list[index_xy_region]
           Sigma = Sigma_list[index_xy_region]
           error_values = np.matrix(error_list[index_xy_region]) 
           inv_Sigma = inv_Sigma_list[index_xy_region]
           RAlocal_now = local_now['RA']
           DEClocal_now = local_now['DEC']
           Dpoint = deprojected_distances(RA_point, DEC_point, RA2 = RAlocal_now, DEC2 = DEClocal_now)      
           cmat = thetahat[1]*np.exp(-Dpoint/thetahat[0])
           pl = X_infer @ betahat
           pe = cmat @ inv_Sigma @ error_values.T
           gas_density_pred = pl + pe.item()
           inference_gas_density.append(gas_density_pred)         
    print(len(inference_gas_density))
    return inference_gas_density








# 1 arcsecond = 4.84813681 × 10-6 radians
def arcsecond_to_rad(arcsecond):
    rad = arcsecond*(4.84813681*10**(-6))
    return rad

def arcseond_to_pc_resolusion(arcsecond, DMpc):
    # DMpc: Distance from Earth in Mpc， 4.89778819
    # arcsecond: the angular resolusion in as, we need to change to radian
    # TYPHOON 1.65 as = 39.12 pc
    # ALMA 2 as       = 47.41 pc
    # MUSE 0.2 as     = 4.749029 pc
    rad = arcsecond_to_rad(arcsecond)
    pc = DMpc*rad * 1000000
    return pc


def direct_convert_factor(as1, as2, DMpc):
    re_pc1 = arcseond_to_pc_resolusion(as1, DMpc)
    re_pc2 = arcseond_to_pc_resolusion(as2, DMpc)
    factor1over2 = (re_pc1**2)/(re_pc2**2)
    factor2over1 = 1/factor1over2
    return factor1over2, factor2over1



def Co21_data_processing(RA_grid, DEC_grid, gas_density):
    data_dict1 = {
    	'RA':                   RA_grid.flatten(),
    	'DEC':                  DEC_grid.flatten(),
        'CO21_gas_density':     gas_density.flatten()
    	}
    df1 = pd.DataFrame(data_dict1)
    index_above_0 = df1["CO21_gas_density"] > 0
    above_0 = df1[index_above_0]
    RA_full = above_0['RA']
    DEC_full = above_0['DEC']
    CO21_gas_density = above_0['CO21_gas_density']
    proj_dist = RA_DEC_to_radius(RA_full, DEC_full)
    X_full = np.transpose(RA_DEC_to_xy(RA_full, DEC_full, meta))
    data_dict_new = {
    	'RA':                   RA_full.values.flatten(),
    	'DEC':                  DEC_full.values.flatten(),
        'X':                    X_full[:,0].flatten(),
        'Y':                    X_full[:,1].flatten(),
        'CO21_gas_density':     CO21_gas_density.values.flatten(),
        'proj_dist':            proj_dist.flatten()
    	}
    data = pd.DataFrame(data_dict_new)
    return data



def pred_gas_density_weight(xy_center_coor, error_list, Sigma_list, betahat_list, thetahat_list, xy_boundary, p, local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x, max_y, min_x, min_y):
    inference_gas_density = list()
    for gg in range(data_galaxy.shape[0]):
        print(gg)
        x_point = X_galaxy[gg]
        y_point = Y_galaxy[gg]
        coor = [x_point,y_point]
        RA_point = RA_galaxy_list.values[gg]
        DEC_point = DEC_galaxy_list.values[gg]
        X_infer = np.array(proj_dist_galaxy_list.values[gg])
        X_infer = insert(X_infer,0,1)
        if (x_point > max_x) or (x_point < min_x) or (y_point > max_y) or (y_point < min_y):
           gas_density_pred = np.nan
           inference_gas_density.append(gas_density_pred)
        else:
           for t in range(len(xy_boundary)):
               xy_region_now = xy_boundary[t]
               x_min_now = xy_region_now[0]
               x_max_now = xy_region_now[1]
               y_min_now = xy_region_now[2]
               y_max_now = xy_region_now[3]
               if (x_min_now <= x_point <= x_max_now) and (y_min_now <= y_point <= y_max_now):
                  index_xy_region = t
                  if (y_max_now <= 5) and (y_min_now >= -5) and (x_max_now <= 5) and (x_min_now >= -5):
                       index_up = t+12
                       index_down = t-12
                       index_left = t-1
                       index_right = t+1
                       index_right_up = index_up +1
                       index_right_down = index_down +1
                       index_left_up = index_up-1
                       index_left_down = index_down-1                           
                       pred_self = MLE_local(index_xy_region, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_self = call_weight(xy_center_coor[index_xy_region], coor, p)               
                       pred_up = MLE_local(index_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_up = call_weight(xy_center_coor[index_up], coor, p)                
                       pred_down = MLE_local(index_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_down = call_weight(xy_center_coor[index_down], coor, p)               
                       pred_left = MLE_local(index_left, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left = call_weight(xy_center_coor[index_left], coor, p)               
                       pred_right = MLE_local(index_right, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right = call_weight(xy_center_coor[index_right], coor, p)                
                       pred_right_up = MLE_local(index_right_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right_up = call_weight(xy_center_coor[index_right_up], coor, p)               
                       pred_right_down = MLE_local(index_right_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right_down = call_weight(xy_center_coor[index_right_down], coor, p)               
                       pred_left_up = MLE_local(index_left_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left_up = call_weight(xy_center_coor[index_left_up], coor, p)               
                       pred_left_down = MLE_local(index_left_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left_down = call_weight(xy_center_coor[index_left_down], coor, p)               
                       gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_down*pred_down + w_left*pred_left + w_right*pred_right + w_right_up*pred_right_up + w_right_down*pred_right_down + w_left_up*pred_left_up + w_left_down*pred_left_down
                       gas_density_pred2 = w_self + w_up + w_down + w_left + w_right + w_right_up + w_right_down + w_left_up + w_left_down
                       gas_density_pred = gas_density_pred1/gas_density_pred2
                       inference_gas_density.append(gas_density_pred)                       
                  elif (y_min_now == 5) and (y_max_now == 6) and (x_max_now <= 5) and (x_min_now >= -5): # up line
                       index_down = t-12
                       index_left = t-1
                       index_right = t+1
                       index_right_down = index_down +1
                       index_left_down = index_down-1                             
                       pred_self = MLE_local(index_xy_region, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_self = call_weight(xy_center_coor[index_xy_region], coor, p)                                
                       pred_down = MLE_local(index_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_down = call_weight(xy_center_coor[index_down], coor, p)               
                       pred_left = MLE_local(index_left, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left = call_weight(xy_center_coor[index_left], coor, p)                
                       pred_right = MLE_local(index_right, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right = call_weight(xy_center_coor[index_right], coor, p)               
                       pred_right_down = MLE_local(index_right_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right_down = call_weight(xy_center_coor[index_right_down], coor, p)              
                       pred_left_down = MLE_local(index_left_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left_down = call_weight(xy_center_coor[index_left_down], coor, p)              
                       gas_density_pred1 = w_self*pred_self + w_down*pred_down + w_left*pred_left + w_right*pred_right + w_right_down*pred_right_down + w_left_down*pred_left_down
                       gas_density_pred2 = w_self + w_down + w_left + w_right + w_right_down + w_left_down
                       gas_density_pred = gas_density_pred1/gas_density_pred2
                       inference_gas_density.append(gas_density_pred)                      
                  elif (y_min_now == -6) and (y_max_now == -5) and (x_max_now <= 5) and (x_min_now >= -5): # down line
                       index_up = t+12
                       index_left = t-1
                       index_right = t+1
                       index_right_up = index_up +1
                       index_left_up = index_up-1                              
                       pred_self = MLE_local(index_xy_region, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_self = call_weight(xy_center_coor[index_xy_region], coor, p)                
                       pred_up = MLE_local(index_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_up = call_weight(xy_center_coor[index_up], coor, p)               
                       pred_left = MLE_local(index_left, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left = call_weight(xy_center_coor[index_left], coor, p)               
                       pred_right = MLE_local(index_right, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right = call_weight(xy_center_coor[index_right], coor, p)                
                       pred_right_up = MLE_local(index_right_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right_up = call_weight(xy_center_coor[index_right_up], coor, p)
                       pred_left_up = MLE_local(index_left_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left_up = call_weight(xy_center_coor[index_left_up], coor, p)                
                       gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_left*pred_left + w_right*pred_right + w_right_up*pred_right_up + w_left_up*pred_left_up 
                       gas_density_pred2 = w_self + w_up + w_left + w_right + w_right_up  + w_left_up 
                       gas_density_pred = gas_density_pred1/gas_density_pred2
                       inference_gas_density.append(gas_density_pred)                
                  elif (x_min_now ==-6) and (x_max_now == -5) and (y_max_now <= 5) and (y_min_now >= -5): # left line
                       index_up = t+12
                       index_down = t-12
                       index_right = t+1
                       index_right_up = index_up +1
                       index_right_down = index_down +1                             
                       pred_self = MLE_local(index_xy_region, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_self = call_weight(xy_center_coor[index_xy_region], coor, p)                
                       pred_up = MLE_local(index_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_up = call_weight(xy_center_coor[index_up], coor, p)                
                       pred_down = MLE_local(index_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_down = call_weight(xy_center_coor[index_down], coor, p)
                       pred_right = MLE_local(index_right, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right = call_weight(xy_center_coor[index_right], coor, p)                
                       pred_right_up = MLE_local(index_right_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right_up = call_weight(xy_center_coor[index_right_up], coor, p)                
                       pred_right_down = MLE_local(index_right_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right_down = call_weight(xy_center_coor[index_right_down], coor, p)
                       gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_down*pred_down  + w_right*pred_right + w_right_up*pred_right_up + w_right_down*pred_right_down 
                       gas_density_pred2 = w_self + w_up + w_down + w_right + w_right_up + w_right_down 
                       gas_density_pred = gas_density_pred1/gas_density_pred2
                       inference_gas_density.append(gas_density_pred)                      
                  elif (x_min_now == 5) and (x_max_now == 6) and (y_max_now <= 5) and (y_min_now >= -5): # right line
                       index_up = t+12
                       index_down = t-12
                       index_left = t-1
                       index_left_up = index_up-1
                       index_left_down = index_down-1                                
                       pred_self = MLE_local(index_xy_region, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_self = call_weight(xy_center_coor[index_xy_region], coor, p)              
                       pred_up = MLE_local(index_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_up = call_weight(xy_center_coor[index_up], coor, p)              
                       pred_down = MLE_local(index_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_down = call_weight(xy_center_coor[index_down], coor, p)               
                       pred_left = MLE_local(index_left, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left = call_weight(xy_center_coor[index_left], coor, p)
                       pred_left_up = MLE_local(index_left_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left_up = call_weight(xy_center_coor[index_left_up], coor, p)            
                       pred_left_down = MLE_local(index_left_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left_down = call_weight(xy_center_coor[index_left_down], coor, p)           
                       gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_down*pred_down + w_left*pred_left + w_left_up*pred_left_up + w_left_down*pred_left_down
                       gas_density_pred2 = w_self + w_up + w_down + w_left + w_left_up + w_left_down
                       gas_density_pred = gas_density_pred1/gas_density_pred2
                       inference_gas_density.append(gas_density_pred)                      
                  elif (y_max_now == 6) and (y_min_now == 5) and (x_max_now ==6) and (x_min_now==5): #up right
                       index_down = t-12
                       index_left = t-1
                       index_left_down = index_down-1                             
                       pred_self = MLE_local(index_xy_region, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_self = call_weight(xy_center_coor[index_xy_region], coor, p) 
                       pred_down = MLE_local(index_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_down = call_weight(xy_center_coor[index_down], coor, p)                
                       pred_left = MLE_local(index_left, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left = call_weight(xy_center_coor[index_left], coor, p)
                       pred_left_down = MLE_local(index_left_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left_down = call_weight(xy_center_coor[index_left_down], coor, p)      
                       gas_density_pred1 = w_self*pred_self + w_down*pred_down + w_left*pred_left + w_left_down*pred_left_down
                       gas_density_pred2 = w_self + w_down + w_left + w_left_down
                       gas_density_pred = gas_density_pred1/gas_density_pred2
                       inference_gas_density.append(gas_density_pred)                      
                  elif (y_max_now == 6) and (y_min_now == 5) and (x_min_now ==-6) and (x_max_now == -5): # up left
                       index_down = t-12
                       index_right = t+1
                       index_right_down = index_down +1                            
                       pred_self = MLE_local(index_xy_region, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_self = call_weight(xy_center_coor[index_xy_region], coor, p) 
                       pred_down = MLE_local(index_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_down = call_weight(xy_center_coor[index_down], coor, p) 
                       pred_right = MLE_local(index_right, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right = call_weight(xy_center_coor[index_right], coor, p)
                       pred_right_down = MLE_local(index_right_down, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right_down = call_weight(xy_center_coor[index_right_down], coor, p)
                       gas_density_pred1 = w_self*pred_self + w_down*pred_down + w_right*pred_right + w_right_down*pred_right_down
                       gas_density_pred2 = w_self + w_down + w_right + w_right_down 
                       gas_density_pred = gas_density_pred1/gas_density_pred2
                       inference_gas_density.append(gas_density_pred)                      
                  elif (y_min_now == -6) and (y_max_now == -5) and (x_min_now ==-6) and (x_max_now == -5): # down left
                       index_up = t+12
                       index_right = t+1
                       index_right_up = index_up +1                             
                       pred_self = MLE_local(index_xy_region, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_self = call_weight(xy_center_coor[index_xy_region], coor, p)                
                       pred_up = MLE_local(index_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_up = call_weight(xy_center_coor[index_up], coor, p)
                       pred_right = MLE_local(index_right, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right = call_weight(xy_center_coor[index_right], coor, p)               
                       pred_right_up = MLE_local(index_right_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_right_up = call_weight(xy_center_coor[index_right_up], coor, p)     
                       gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_right*pred_right + w_right_up*pred_right_up 
                       gas_density_pred2 = w_self + w_up + w_right + w_right_up 
                       gas_density_pred = gas_density_pred1/gas_density_pred2
                       inference_gas_density.append(gas_density_pred)  
                  elif (y_min_now == -6) and (y_max_now == -5) and (x_max_now ==6) and (x_min_now ==5): # conditions not right
                       index_up = t + 12
                       index_left = t-1
                       index_left_up = index_up-1                              
                       pred_self = MLE_local(index_xy_region, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_self = call_weight(xy_center_coor[index_xy_region], coor, p)               
                       pred_up = MLE_local(index_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_up = call_weight(xy_center_coor[index_up], coor, p)
                       pred_left = MLE_local(index_left, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left = call_weight(xy_center_coor[index_left], coor, p)
                       pred_left_up = MLE_local(index_left_up, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                       w_left_up = call_weight(xy_center_coor[index_left_up], coor, p)                    
                       gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_left*pred_left + w_left_up*pred_left_up 
                       gas_density_pred2 = w_self + w_up + w_left +  w_left_up
                       gas_density_pred = gas_density_pred1/gas_density_pred2
                       inference_gas_density.append(gas_density_pred)          
    print(len(inference_gas_density))
    return inference_gas_density


def MLE_local(index, X_infer, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list):
    local_now = local_df[index]
    thetahat = thetahat_list[index]
    betahat = betahat_list[index]
    error_values = np.matrix(error_list[index]) 
    inv_Sigma = inv_Sigma_list[index]
    RAlocal_now = local_now['RA']
    DEClocal_now = local_now['DEC']
    Dpoint = deprojected_distances(RA_point, DEC_point, RA2 = RAlocal_now, DEC2 = DEClocal_now)      
    cmat = thetahat[1]*np.exp(-Dpoint/thetahat[0])
    pl = X_infer @ betahat
    pe = cmat @ inv_Sigma @ error_values.T
    gas_density_pred = pl + pe.item()
    return gas_density_pred

def call_weight(coor1, coor2, p):
    weight = 1/((math.dist(coor1, coor2))**p)
    return weight

def plot_variograms_for_each_regions(local_df, model = 'exponential'):
    '''
    model = 'exponential'
    'spherical', 'exponential', 'gaussian', 'matern', 'stable', 'cubic'
    parameters:
        1. local_df
        
    returns:

    '''
    for zz in range(len(local_df)):
        print(zz)
        count = 0
        local_now = local_df[zz]
        check = local_now.empty
        if check == False:
           V1 = skg.Variogram(local_now[['X', 'Y']].values, local_now.CO21_gas_density.values, model = model)
           V1.plot()
           count = count + 1
        else:
           print('DataFrame is empty!')
    return count

