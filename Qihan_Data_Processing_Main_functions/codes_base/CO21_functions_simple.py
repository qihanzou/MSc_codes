# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 16:41:15 2024

Modified on Fri Jul 12 14:08:12 2024

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
from dec_ra_to_xy import *
from sklearn.metrics import mean_squared_error 
import statistics
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform
import skgstat as skg

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

def gal_data_processing(data_galaxy, Z_train = "Z_N2S2Ha", DIG_CUT = "N2_BPT"):
    if DIG_CUT == "N2_BPT":
       index_below_1_BPT = data_galaxy["N2_BPT"] < 1 
       data_galaxy = data_galaxy[index_below_1_BPT]   
    if DIG_CUT == "S2_BPT":
       idx_S2_BPT1 = data_galaxy["S2_BPT"] < 1
       data_galaxy = data_galaxy[idx_S2_BPT1]
    if DIG_CUT == "S2_DIG":
       index_above_CHii = data_galaxy["S2_DIG"] > 0.9 
       data_galaxy = data_galaxy[index_above_CHii]
    index_above_0_Z1 = data_galaxy[Z_train] > 0  
    data_galaxy = data_galaxy[index_above_0_Z1]
    return data_galaxy


def obtain_useful_vars(data_galaxy, meta):
    RA_galaxy_list = data_galaxy['RA']
    DEC_galaxy_list = data_galaxy['DEC']
    proj_dist_galaxy_list = data_galaxy['proj_dist']
    coor_full_galaxy = np.transpose(RA_DEC_to_xy(RA_galaxy_list, DEC_galaxy_list, meta))
    X_galaxy = coor_full_galaxy[:,0]
    Y_galaxy = coor_full_galaxy[:,1]
    return RA_galaxy_list, DEC_galaxy_list, proj_dist_galaxy_list, coor_full_galaxy, X_galaxy, Y_galaxy




def create_regions_from_whole(num_parts_axis, min_x, min_y, max_x, max_y, total_length, data):
    '''
    *parameters:
    1. num_parts_axis: number of parts for both axis: eg: 1, 2, 3, 4, 5, 6, ..., 12, 24
    2. min_x: minimum value for x axis, check pkl to set this value
    3. min_y: minimum value for y axis, check pkl to set this value
    4. total_length: total length of axis, need to be same for x y for now. should be abs(min_x) + max_x
    
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
    plot_regions_from_whole(num_parts_axis = num_parts_axis, min_x = min_x, min_y = min_y, max_x = max_x, max_y = max_y, total_length = total_length)
    plt.scatter(data["X"], data["Y"], marker='.', s=0.005, c = 'Black')
    plt.show()
    return local_df, xy_boundary, xy_center_coor

def plot_regions_from_whole(num_parts_axis, min_x, min_y, max_x, max_y, total_length):
    add = total_length/num_parts_axis
    for vy in range(num_parts_axis+1):
        plt.axhline(y = min_y + vy*add, color = 'r', linestyle = '-')
    for hx in range(num_parts_axis+1):
        plt.axvline(x = min_x + hx*add, color = 'r', linestyle = '-')
        
        

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
           V1 = skg.Variogram(local_now[['X', 'Y']].values, local_now.CO21_gas_density.values, model = 'exponential')
           var_range = V1.parameters[0]
           var_sill = V1.parameters[1]
           var_nug = V1.parameters[2]
           eta_ini_value = np.array([var_sill-var_nug, var_range])
           lower_bound = np.array([1e-5, 1e-10])
           upper_bound = np.array([eta_ini_value[0]*5,eta_ini_value[1]*5])
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










def construct_geo_dist_models_for_regions(local_df):
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
           X_local_now = local_now[["proj_dist"]]
           X_local_now.insert(0, "intersect", np.ones(X_local_now.shape[0]), True)
           Y_local_now = local_now["CO21_gas_density"] 
           Y_local_now = np.log(Y_local_now)
           eta_ini_value = np.array([0.01, 0.1])
           lower_bound = np.array([1e-5, 1e-10])
           upper_bound = np.array(0,1)
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




