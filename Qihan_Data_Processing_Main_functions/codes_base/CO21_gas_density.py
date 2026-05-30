# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 21:30:01 2024

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

#cuda.detect()
'''
MUSE_N5068 = fits.open("C:/Users/qihan/Desktop/q/MUSE/NGC5068_maps_copt.fits")
MUSE_N7496 = fits.open("C:/Users/qihan/Desktop/q/MUSE/NGC7496_maps_copt.fits")
MUSE_N4535 = fits.open("C:/Users/qihan/Desktop/q/MUSE/NGC4535_maps_copt.fits")
MUSE_N4321 = fits.open("C:/Users/qihan/Desktop/q/MUSE/NGC4321_maps_copt.fits")
MUSE_N4303 = fits.open("C:/Users/qihan/Desktop/q/MUSE/NGC4303_maps_copt.fits")
'''
'''
N5236_co21_strict_mom0
N5236_co21_15as_strict_mom0
N5236_co21_11as_strict_mom0
N5236_co21_2as_strict_mom0
N5236_co21_7p5as_strict_mom0
'''
# In my opinion, we should use 2as (arcseconds) fixed angular resolution. 
data = fits.open(r"C:\Users\qihan\Desktop\q\ALMA_CO21\N5236_co21\N5236_co21_2as_strict_mom0.fits") # check path before use
ICO21 = data[0] # unit: K km s-1 from PHANGS-ALMA readme.
RA_grid, DEC_grid = make_RA_DEC_grid(data[0].header) # create RA and DEC grid.
meta = meta_getter('N5236') # get basic property form meta data.

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
    i  =  np.radians(meta['i']) # i should be in radian, right?
    gas_density = 6.7*ICO21.data*np.cos(i) 
    return gas_density # unit: M_sun pc-2

gas_density = calculate_gas_density_CO21(ICO21, meta)


with open('N5236_co21dist.pkl', 'rb') as f:
    data = pickle.load(f)
    print(data.columns.values)
    
df = pd.DataFrame(data)
RA = df['RA']
DEC = df['DEC']
proj_dist = RA_DEC_to_radius(RA, DEC)

out_path = 'C:/Users/qihan/Desktop/q/'
data_dict = {
	'RA':                   RA_grid.flatten(),
	'DEC':                  DEC_grid.flatten(),
    'CO21_gas_density':     gas_density.flatten(),
    'ICO21':                ICO21.data.flatten(),
    'proj_dist':            proj_dist.flatten()
	}
result_df = pd.DataFrame(data_dict)
result_df.to_pickle(out_path+'N5236_co21_gas_density'+'.pkl')

# plot it: 
df1 = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_co21_gas_density.pkl')
positions_are_NaNs = isnan(df1)
df1[positions_are_NaNs] = 0
index_above_0 = df1["CO21_gas_density"] > 0
above_0 = df1[index_above_0]
RA_full = above_0['RA']
DEC_full = above_0['DEC']
CO21_gas_density = above_0['CO21_gas_density']
proj_dist = above_0['proj_dist']
X_full = RA_DEC_to_xy(RA_full, DEC_full, meta)
X_full = np.transpose(X_full)
plt.scatter(X_full[:, 0], X_full[:, 1], c = 'black', marker='.', s=0.1)
plt.grid()
plt.show()



plt.axhline(y = -6, color = 'r', linestyle = '-')
plt.axhline(y = -5, color = 'r', linestyle = '-')
plt.axhline(y = -4, color = 'r', linestyle = '-')
plt.axhline(y = -3, color = 'r', linestyle = '-')
plt.axhline(y = -2, color = 'r', linestyle = '-')
plt.axhline(y = -1, color = 'r', linestyle = '-')
plt.axhline(y = 0, color = 'r', linestyle = '-')
plt.axhline(y = 1, color = 'r', linestyle = '-')
plt.axhline(y = 2, color = 'r', linestyle = '-')
plt.axhline(y = 3, color = 'r', linestyle = '-')
plt.axhline(y = 4, color = 'r', linestyle = '-')
plt.axhline(y = 5, color = 'r', linestyle = '-')
plt.axhline(y = 6, color = 'r', linestyle = '-')
plt.axvline(x = -6, color = 'r', linestyle = '-')
plt.axvline(x = -5, color = 'r', linestyle = '-')
plt.axvline(x = -4, color = 'r', linestyle = '-')
plt.axvline(x = -3, color = 'r', linestyle = '-')
plt.axvline(x = -2, color = 'r', linestyle = '-')
plt.axvline(x = -1, color = 'r', linestyle = '-')
plt.axvline(x = 0, color = 'r', linestyle = '-')
plt.axvline(x = 1, color = 'r', linestyle = '-')
plt.axvline(x = 2, color = 'r', linestyle = '-')
plt.axvline(x = 3, color = 'r', linestyle = '-')
plt.axvline(x = 4, color = 'r', linestyle = '-')
plt.axvline(x = 5, color = 'r', linestyle = '-')
plt.axvline(x = 6, color = 'r', linestyle = '-')
plt.scatter(X_full[:, 0], X_full[:, 1], marker='.', s=0.005, c = 'Black')
plt.xlim(-6.02,6)
plt.ylim(-6,6.01)
plt.grid()
plt.show()




data_dict_new = {
	'RA':                   RA_full.values.flatten(),
	'DEC':                  DEC_full.values.flatten(),
    'X':                    X_full[:,0].flatten(),
    'Y':                    X_full[:,1].flatten(),
    'CO21_gas_density':     CO21_gas_density.values.flatten(),
    'proj_dist':            proj_dist.values.flatten()
	}
	# save it
data = pd.DataFrame(data_dict_new)



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

local_df, xy_boundary, xy_center_coor = create_regions_from_whole(num_parts_axis = 24, min_x = -6, min_y = -6, max_x = 6, max_y = 6, total_length = 12, data = data)

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
    


def plot_variograms_for_each_regions(local_df, model):
    '''
    model = 'exponential'
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
           V1.plot(show = True)
           count = count + 1
        else:
           print('DataFrame is empty!')
    return count

def plot_df_gasdensity_vs_projdist(local_df):
    '''
    parameters:
        1. local_df
        
    returns:
        1. many plots
        * my computer can plots all points at once, I do not need to sepearate them.
    '''
    count = 0 # just for checking
    for zz in range(len(local_df)):
        print(zz)
        local_now = local_df[zz]
        check = local_now.empty
        if check == False:
           y = local_now["CO21_gas_density"]
           #y = np.log(y)
           x = local_now[["proj_dist"]]
           #x = np.log(x)
           plt.scatter(x, y)
           plt.ylabel('gas density')
           plt.xlabel('dist')
           plt.show()
           count = count + 1
        else:
           print('DataFrame is empty!')
    return count



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

thetahat_list, betahat_list, Sigma_list, error_list, inv_Sigma_list = construct_geo_dist_models_for_regions(local_df)

       
def open_pkl_and_data_processing(pkl_data):
    # data for galaxy, eg. NGC5236:
    data_galaxy = pd.read_pickle(pkl_data)
    index_above_0_Z = data_galaxy["Z_N2S2Ha"] > 0 # ignore NA
    data_galaxy = data_galaxy[index_above_0_Z]
    index_below_1_BPT = data_galaxy["N2_BPT"] < 1 # BPT cut for Hii regions
    data_galaxy = data_galaxy[index_below_1_BPT] # result: 7193 rows of data
    RA_galaxy_list = data_galaxy['RA']
    DEC_galaxy_list = data_galaxy['DEC']
    proj_dist_galaxy_list = data_galaxy['proj_dist']
    coor_full_galaxy = RA_DEC_to_xy(RA_galaxy_list, DEC_galaxy_list, meta)
    coor_full_galaxy = np.transpose(coor_full_galaxy)
    X_galaxy = coor_full_galaxy[:,0]
    Y_galaxy = coor_full_galaxy[:,1]
    return data_galaxy, RA_galaxy_list, DEC_galaxy_list, proj_dist_galaxy_list, coor_full_galaxy, X_galaxy, Y_galaxy
    
data_galaxy, RA_galaxy_list, DEC_galaxy_list, proj_dist_galaxy_list, coor_full_galaxy, X_galaxy, Y_galaxy = open_pkl_and_data_processing(pkl_data = 'C:/Users/qihan/Desktop/q/N5236_25_1_2024.pkl')



def pred_gas_density(local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x, max_y, min_x, min_y):
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
        
inference_gas_density = pred_gas_density(local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x = 6, max_y = 6, min_x = -6, min_y = -6)



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


def Inv_dist_weight_averaging(X,Y,local_df, xy_boundary, xy_center_coor, RA_point, DEC_point):
    x_point = X # x coor of pred point
    y_point = Y # y coor of pred point
    coor = [x_point, y_point]
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
                local_now = local_df[index_xy_region]
                local_up = local_df[index_up]
                local_down = local_df[index_down]
                local_left = local_df[index_left]
                local_right = local_df[index_right]
                local_right_up = local_df[index_right_up]
                local_right_down = local_df[index_right_down]
                local_left_up = local_df[index_left_up]
                local_left_down = local_df[index_left_down]               
                pred_self = MLE_local(index_xy_region, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_self = call_weight(xy_center_coor[index_xy_region], coor, p)               
                pred_up = MLE_local(index_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_up = call_weight(xy_center_coor[index_up], coor, p)                
                pred_down = MLE_local(index_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_down = call_weight(xy_center_coor[index_down], coor, p)               
                pred_left = MLE_local(index_left, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left = call_weight(xy_center_coor[index_left], coor, p)               
                pred_right = MLE_local(index_right, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right = call_weight(xy_center_coor[index_right], coor, p)                
                pred_right_up = MLE_local(index_right_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right_up = call_weight(xy_center_coor[index_right_up], coor, p)               
                pred_right_down = MLE_local(index_right_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right_down = call_weight(xy_center_coor[index_right_down], coor, p)               
                pred_left_up = MLE_local(index_left_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left_up = call_weight(xy_center_coor[index_left_up], coor, p)               
                pred_left_down = MLE_local(index_left_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left_down = call_weight(xy_center_coor[index_left_down], coor, p)               
                gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_down*pred_down + w_left*pred_left + w_right*pred_right + w_right_up*pred_right_up + w_right_down*pred_right_down + w_left_up*pred_left_up + w_left_down*pred_left_down
                gas_density_pred2 = w_self + w_up + w_down + w_left + w_right + w_right_up + w_right_down + w_left_up + w_left_down
                gas_density_pred = gas_density_pred1/gas_density_pred2
                inference_gas_density.append(gas_density_pred)
                
           elif (y_max_now > 5) and (x_max_now <= 5) and (x_min_now >= -5): # up line
                index_down = t-12
                index_left = t-1
                index_right = t+1
                index_right_down = index_down +1
                index_left_down = index_down-1               
                local_now = local_df[index_xy_region]
                local_down = local_df[index_down]
                local_left = local_df[index_left]
                local_right = local_df[index_right]                
                local_right_down = local_df[index_right_down]               
                local_left_down = local_df[index_left_down]               
                pred_self = MLE_local(index_xy_region, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_self = call_weight(xy_center_coor[index_xy_region], coor, 2)                                
                pred_down = MLE_local(index_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_down = call_weight(xy_center_coor[index_down], coor, 2)               
                pred_left = MLE_local(index_left, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left = call_weight(xy_center_coor[index_left], coor, 2)                
                pred_right = MLE_local(index_right, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right = call_weight(xy_center_coor[index_right], coor, 2)               
                pred_right_down = MLE_local(index_right_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right_down = call_weight(xy_center_coor[index_right_down], coor, 2)              
                pred_left_down = MLE_local(index_left_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left_down = call_weight(xy_center_coor[index_left_down], coor, 2)              
                gas_density_pred1 = w_self*pred_self + w_down*pred_down + w_left*pred_left + w_right*pred_right + w_right_down*pred_right_down + w_left_down*pred_left_down
                gas_density_pred2 = w_self + w_down + w_left + w_right + w_right_down + w_left_down
                gas_density_pred = gas_density_pred1/gas_density_pred2
                inference_gas_density.append(gas_density_pred)
                
           elif (y_min_now < -5) and (x_max_now <= 5) and (x_min_now >= -5): # down line
                index_up = t+12
                index_left = t-1
                index_right = t+1
                index_right_up = index_up +1
                index_left_up = index_up-1                
                local_now = local_df[index_xy_region]
                local_up = local_df[index_up]    
                local_left = local_df[index_left]
                local_right = local_df[index_right]
                local_right_up = local_df[index_right_up]
                local_left_up = local_df[index_left_up]                
                pred_self = MLE_local(index_xy_region, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_self = call_weight(xy_center_coor[index_xy_region], coor, 2)                
                pred_up = MLE_local(index_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_up = call_weight(xy_center_coor[index_up], coor, 2)               
                pred_left = MLE_local(index_left, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left = call_weight(xy_center_coor[index_left], coor, 2)               
                pred_right = MLE_local(index_right, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right = call_weight(xy_center_coor[index_right], coor, 2)                
                pred_right_up = MLE_local(index_right_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right_up = call_weight(xy_center_coor[index_right_up], coor, 2)
                pred_left_up = MLE_local(index_left_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left_up = call_weight(xy_center_coor[index_left_up], coor, 2)                
                gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_left*pred_left + w_right*pred_right + w_right_up*pred_right_up + w_left_up*pred_left_up 
                gas_density_pred2 = w_self + w_up + w_left + w_right + w_right_up  + w_left_up 
                gas_density_pred = gas_density_pred1/gas_density_pred2
                inference_gas_density.append(gas_density_pred)

                
           elif (x_min_now <-5) and (y_max_now <= 5) and (y_min_now >= -5): # left line
                index_up = t+12
                index_down = t-12
                index_right = t+1
                index_right_up = index_up +1
                index_right_down = index_down +1               
                local_now = local_df[index_xy_region]
                local_up = local_df[index_up]
                local_down = local_df[index_down]
                local_right = local_df[index_right]
                local_right_up = local_df[index_right_up]
                local_right_down = local_df[index_right_down]               
                pred_self = MLE_local(index_xy_region, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_self = call_weight(xy_center_coor[index_xy_region], coor, 2)                
                pred_up = MLE_local(index_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_up = call_weight(xy_center_coor[index_up], coor, 2)                
                pred_down = MLE_local(index_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_down = call_weight(xy_center_coor[index_down], coor, 2)
                pred_right = MLE_local(index_right, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right = call_weight(xy_center_coor[index_right], coor, 2)                
                pred_right_up = MLE_local(index_right_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right_up = call_weight(xy_center_coor[index_right_up], coor, 2)                
                pred_right_down = MLE_local(index_right_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right_down = call_weight(xy_center_coor[index_right_down], coor, 2)
                gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_down*pred_down  + w_right*pred_right + w_right_up*pred_right_up + w_right_down*pred_right_down 
                gas_density_pred2 = w_self + w_up + w_down + w_right + w_right_up + w_right_down 
                gas_density_pred = gas_density_pred1/gas_density_pred2
                inference_gas_density.append(gas_density_pred)
                
           elif (x_max_now > 5) and (y_max_now <= 5) and (y_min_now >= -5): # right line
                index_up = t+12
                index_down = t-12
                index_left = t-1
                index_left_up = index_up-1
                index_left_down = index_down-1                
                local_now = local_df[index_xy_region]
                local_up = local_df[index_up]
                local_down = local_df[index_down]
                local_left = local_df[index_left]
                local_left_up = local_df[index_left_up]
                local_left_down = local_df[index_left_down]                
                pred_self = MLE_local(index_xy_region, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_self = call_weight(xy_center_coor[index_xy_region], coor, 2)              
                pred_up = MLE_local(index_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_up = call_weight(xy_center_coor[index_up], coor, 2)              
                pred_down = MLE_local(index_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_down = call_weight(xy_center_coor[index_down], coor, 2)               
                pred_left = MLE_local(index_left, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left = call_weight(xy_center_coor[index_left], coor, 2)
                pred_left_up = MLE_local(index_left_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left_up = call_weight(xy_center_coor[index_left_up], coor, 2)            
                pred_left_down = MLE_local(index_left_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left_down = call_weight(xy_center_coor[index_left_down], coor, 2)           
                gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_down*pred_down + w_left*pred_left + w_left_up*pred_left_up + w_left_down*pred_left_down
                gas_density_pred2 = w_self + w_up + w_down + w_left + w_left_up + w_left_down
                gas_density_pred = gas_density_pred1/gas_density_pred2
                inference_gas_density.append(gas_density_pred)
                
           elif (y_max_now == 6) and (x_max_now ==6): #up right
                index_down = t-12
                index_left = t-1
                index_left_down = index_down-1               
                local_now = local_df[index_xy_region]
                local_down = local_df[index_down]
                local_left = local_df[index_left]
                local_left_down = local_df[index_left_down]               
                pred_self = MLE_local(index_xy_region, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_self = call_weight(xy_center_coor[index_xy_region], coor, 2) 
                pred_down = MLE_local(index_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_down = call_weight(xy_center_coor[index_down], coor, 2)                
                pred_left = MLE_local(index_left, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left = call_weight(xy_center_coor[index_left], coor, 2)
                pred_left_down = MLE_local(index_left_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left_down = call_weight(xy_center_coor[index_left_down], coor, 2)      
                gas_density_pred1 = w_self*pred_self + w_down*pred_down + w_left*pred_left + w_left_down*pred_left_down
                gas_density_pred2 = w_self + w_down + w_left + w_left_down
                gas_density_pred = gas_density_pred1/gas_density_pred2
                inference_gas_density.append(gas_density_pred)
                
           elif (y_max_now == 6) and (x_min_now ==-6): # up left
                index_down = t-12
                index_right = t+1
                index_right_down = index_down +1               
                local_now = local_df[index_xy_region]
                local_down = local_df[index_down]
                local_right = local_df[index_right]
                local_right_down = local_df[index_right_down]               
                pred_self = MLE_local(index_xy_region, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_self = call_weight(xy_center_coor[index_xy_region], coor, 2) 
                pred_down = MLE_local(index_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_down = call_weight(xy_center_coor[index_down], coor, 2) 
                pred_right = MLE_local(index_right, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right = call_weight(xy_center_coor[index_right], coor, 2)
                pred_right_down = MLE_local(index_right_down, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right_down = call_weight(xy_center_coor[index_right_down], coor, 2)
                gas_density_pred1 = w_self*pred_self + w_down*pred_down + w_right*pred_right + w_right_down*pred_right_down
                gas_density_pred2 = w_self + w_down + w_right + w_right_down 
                gas_density_pred = gas_density_pred1/gas_density_pred2
                inference_gas_density.append(gas_density_pred)
                
           elif (y_min_now == -6) and (x_min_now ==-6): # down left
                index_up = t+12
                index_right = t+1
                index_right_up = index_up +1              
                local_now = local_df[index_xy_region]
                local_up = local_df[index_up]
                local_right = local_df[index_right]
                local_right_up = local_df[index_right_up]                
                pred_self = MLE_local(index_xy_region, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_self = call_weight(xy_center_coor[index_xy_region], coor, 2)                
                pred_up = MLE_local(index_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_up = call_weight(xy_center_coor[index_up], coor, 2)
                pred_right = MLE_local(index_right, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right = call_weight(xy_center_coor[index_right], coor, 2)               
                pred_right_up = MLE_local(index_right_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_right_up = call_weight(xy_center_coor[index_right_up], coor, 2)     
                gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_right*pred_right + w_right_up*pred_right_up 
                gas_density_pred2 = w_self + w_up + w_right + w_right_up 
                gas_density_pred = gas_density_pred1/gas_density_pred2
                inference_gas_density.append(gas_density_pred)                
           elif (y_min_now == -6) and (x_max_now ==6):
                index_up = t + 12
                index_left = t-1
                index_left_up = index_up-1               
                local_now = local_df[index_xy_region]
                local_up = local_df[index_up]
                local_left = local_df[index_left]
                local_left_up = local_df[index_left_up]                
                pred_self = MLE_local(index_xy_region, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_self = call_weight(xy_center_coor[index_xy_region], coor, 2)               
                pred_up = MLE_local(index_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_up = call_weight(xy_center_coor[index_up], coor, 2)
                pred_left = MLE_local(index_left, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left = call_weight(xy_center_coor[index_left], coor, 2)
                pred_left_up = MLE_local(index_left_up, RA_point, DEC_point, local_df, thetahat_list, betahat_list, error_list, inv_Sigma_list)
                w_left_up = call_weight(xy_center_coor[index_left_up], coor, 2)                    
                gas_density_pred1 = w_self*pred_self + w_up*pred_up + w_left*pred_left + w_left_up*pred_left_up 
                gas_density_pred2 = w_self + w_up + w_left +  w_left_up
                gas_density_pred = gas_density_pred1/gas_density_pred2
                inference_gas_density.append(gas_density_pred)

             
p=2
             
def pred_gas_density_weight(p, local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x, max_y, min_x, min_y):
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
        
#inference_gas_density1 = pred_gas_density_weight(1, local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x = 6, max_y = 6, min_x = -6, min_y = -6)

#inference_gas_density2 = pred_gas_density_weight(2, local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x = 6, max_y = 6, min_x = -6, min_y = -6)

inference_gas_density4 = pred_gas_density(local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x = 6, max_y = 6, min_x = -6, min_y = -6)

           

'''
data_dict = {
	'gas_density_w_p1':  pd.DataFrame(inference_gas_density1).values.flatten(),
    'gas_density_w_p2':  pd.DataFrame(inference_gas_density2).values.flatten(),
    'infer_gas_density':  pd.DataFrame(inference_gas_density3).values.flatten() 
	}
result_df = pd.DataFrame(data_dict)
result_df.to_pickle(out_path+'N5236_three_gas_density'+'.pkl')


data_dict = {
	'gas_density_2424':  pd.DataFrame(inference_gas_density4).values.flatten()
	}
result_df = pd.DataFrame(data_dict)
result_df.to_pickle(out_path+'N5236_2424_gas_density'+'.pkl')


'''

'''
gas1 = pd.DataFrame(inference_gas_density1)
gas2 = pd.DataFrame(inference_gas_density2)
gas3 = pd.DataFrame(inference_gas_density3)
gas4 = pd.DataFrame(inference_gas_density4)
data_galaxy.insert(0, "gas_density_w_p1", gas1.values, True)
data_galaxy.insert(0, "gas_density_w_p2", gas2.values, True)
data_galaxy.insert(0, "infer_gas_density", gas3.values, True)
data_galaxy.insert(0, "infer_gas_density2424", gas4.values, True)

data_galaxy = data_galaxy[["Z_N2S2Ha", "RA", "DEC", "proj_dist", "SFR_density", "gas_density_w_p1", "gas_density_w_p2", "infer_gas_density", "infer_gas_density2424", "surface_gas_density"]]
data_galaxy = pd.DataFrame(data_galaxy).dropna()




####################################################################################################
####################################################################################################

avg_gas_co21 = statistics.mean(np.exp(gas_density1)) # 106.80725064428266
med_gas_co21 = statistics.median(np.exp(gas_density1)) # 38.523297210405964
avg_gas_sfr = statistics.mean(gas_density2) # 49.52473453900626
med_gas_sfr = statistics.median(gas_density2) # 36.443010330200195

min_gas_co21 = min(np.exp(gas_density1))
min_gas_sfr = min(gas_density2)
max_gas_co21 = max(np.exp(gas_density1))
max_gas_sfr = max(gas_density2)
print(avg_gas_co21)
print(max_gas_co21) # 6345.653259006231
print(min_gas_co21) # 0.7991455764765604
print(avg_gas_sfr)
print(max_gas_sfr) # 1241.4769287109375
print(min_gas_sfr) # 9.047096252441406
MSE_gas = mean_squared_error(gas_density1,gas_density2) 
RMSE_gas = np.sqrt(MSE_gas)
MAD_gas = statistics.mean(abs(gas_density1 - gas_density2))
print(MSE_gas)
print(RMSE_gas)
print(MAD_gas)

print(statistics.mean(data_galaxy['SFR_density'])) # 0.07082442741188968
print(statistics.mean(np.log(data_galaxy['SFR_density']))) # -3.159027255460935
print(statistics.mean(gas_density1)) # 3.488981483052133
print(statistics.mean(data_galaxy['proj_dist'])) # 3.428938716897325


## Do the geo for metallicity, parameters: dist, sfr, molecular gas density
testset = list()
trainset = list()
random.seed(2024)
idx = random.sample(range(1, len(data_galaxy)), math.ceil(0.2*len(data_galaxy)))
for kk in range(len(data_galaxy)):
    if kk in idx:
        testset.append(data_galaxy.iloc[[kk]])
    else:
        trainset.append(data_galaxy.iloc[[kk]])
testset = pd.concat(testset)
trainset = pd.concat(trainset)
###
trainset[["SFR_density"]] = np.log(trainset[["SFR_density"]])
trainset[["av_gas_density"]] = np.log(trainset[["av_gas_density"]])
###
X_train = trainset[["proj_dist", "SFR_density", "av_gas_density"]]
X_train.insert(0, "intersect", np.ones(X_train.shape[0]), True)
Y_train = trainset["Z_N2S2Ha"] 
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
###
testset[["SFR_density"]] = np.log(testset[["SFR_density"]])
testset[["av_gas_density"]] = np.log(testset[["av_gas_density"]])
###
X_pred = testset[["proj_dist", "SFR_density", "av_gas_density"]]
X_pred.insert(0, "intersect", np.ones(X_pred.shape[0]), True)
cmat = thetahat[1]*np.exp(-Dtest/thetahat[0])
pl = X_pred @ betahat
pe = cmat @ linalg.inv(Sigma) @ (Y_train - X_train @ betahat)
Y_pred = pl + pe
Y_true = testset['Z_N2S2Ha']
MSE = mean_squared_error(Y_true.values,Y_pred.values) 
RMSE = np.sqrt(MSE)
MAD = statistics.mean(abs(Y_pred.values - Y_true.values))
'''

'''
residuals_test = Y_true.values - Y_pred.values
plt.scatter(residuals_test,Y_pred.values)
plt.scatter(Y_true.values, Y_pred.values)
plt.scatter(testset[["SFR_density"]].values, Y_pred.values)
plt.scatter(testset[["inference_gas_density"]].values, Y_pred.values)
plt.scatter(testset[["proj_dist"]].values, Y_pred.values)
plt.show()
'''
'''
MSE
d:                   0.004107056841850173
d+sfr:               0.0038789263003132708
d + logsfr:          0.0028408531736398525
d+sfr+log(gas):      0.0037574180030229003
d + logsfr + loggas: 12x12: 0.002829007590246874
d + logsfr + loggas: 12x12: weight: 0.0028321661266113318
av: 0.0028292178995594237
d + logsfr + loggas: 24x24: 0.002829424096627136
factor               0.002829008061518731

RMSE
d:                   0.0640863233603721
d+sfr:               0.0612977813874442
d+logsfr:            0.05329965453583965
d+sfr+log(gas):      0.0612977813874442
d + logsfr + loggas: 12x12: 0.05318841594038005
d + logsfr + loggas: 12x12: weight: 0.053218099614805225
av: 0.053190392925409224
d + logsfr + loggas: 24x24: 0.05319233118248472
factor               0.053188420370591294

MAD
d:                   0.04641932159815537
d+sfr:               0.04502872682470499
d+logsfr             0.03958424815524499
d+sfr+log(gas):      0.04477661186817579
d+logsfr+loggas:     12x12: 0.039432098833867
d + logsfr + loggas: 12x12: weight: 0.03947540756202831
av: 0.03941834865304587
d + logsfr + loggas: 24x24: 0.03942651956680603
factor               0.03943210146796548
'''

'''
d+sfr+log(gas):
betahat:
array([8.6760006 , 0.00961281, 0.1763438 , 0.02395221])
thetahat:
array([1.17893185e-01, 9.40986450e-03, 9.40986450e-13])
thetahat[0] = phi (kpc)
thetahat[1] = sigma^2
thetahat[2] = tau^2 = sigma^2*theta3

Zc = 8.6760006
beta2 = 0.00961281
SFR = 0.1763438
gas_density = 0.02395221
'''


'''



######################################################################################
y1 = data_galaxy["gas_density_w_p1"]
y2 = data_galaxy["gas_density_w_p2"]
y3 = data_galaxy['infer_gas_density']
y4 = data_galaxy['surface_gas_density']
y5 = data_galaxy['infer_gas_density2424']
x = data_galaxy[["proj_dist"]]


mean_squared_error(gas_density1,gas_density2)
statistics.mean(abs(gas_density1 - gas_density2))

statistics.mean()
statistics.median()


plt.scatter(x, np.exp(y5), s=0.1)
plt.xlabel("Proj_dist (kpc)")
plt.ylabel("CO(2-1) molecular gas surface density (M_sun pc-2)")
plt.title("Without weight, 24x24 grid")

plt.scatter(x, y5, s=0.1)
plt.xlabel("Proj_dist (kpc)")
plt.ylabel("CO(2-1) molecular gas surface density (M_sun pc-2)")
plt.title("log scale, Without weight, 24x24 grid")



plt.scatter(x, y4, s=0.1)
plt.xlabel("Proj_dist (kpc)")
plt.ylabel("SFR molecular gas surface density (M_sun pc-2)")
plt.title("SFR based")

plt.scatter(x, np.log(y4), s=0.1)
plt.xlabel("Proj_dist (kpc)")
plt.ylabel("SFR molecular gas surface density (M_sun pc-2)")
plt.title("log scale, SFR based")

#####################################################################################
#####################################################################################
    '''



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
    







