# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 23:33:56 2024

@author: Qihan Zou

Last updated: 20/11/2024
"""


import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import pickle
import matplotlib.pyplot as plt
from MUSE_data_processing import *


def compute_molecular_gas_surface_density_from_CO21(ICO21, inc_angle_deg):
    '''
    This function can be used to calculate molecular gas surface density from CO(2-1)
    emission line. 
    
    inc_angle is the inclination angle in galaxy info excel file. It changes based
    on different assumptions (from different resources).
    
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
    i  =  np.radians(inc_angle_deg) 
    gas_density = 6.7*ICO21.data*np.cos(i) 
    return gas_density # unit: M_sun pc-2


def gal_data_processing(data_galaxy, Z_train = "Z_N2S2Ha", DIG_CUT = "N2_BPT"):
    '''
    Parameters
    ----------
    data_galaxy : the data of gaalxy before applying any cut.
    Z_train : the metallicity diagnostic you want to use.
        DESCRIPTION. The default is "Z_N2S2Ha".
    DIG_CUT : the DIG cut you want to use.
        DESCRIPTION. The default is "N2_BPT".

    Returns
    -------
    data_galaxy : the data of galaxy after applying metalliicty diagnostic cut 
                  and DIG diagnostic cut.
    '''
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


def create_regions_from_whole(num_parts_axis, min_x, min_y, max_x, max_y, total_length, data):
    '''
    This function is created for cutting the whole spatial data set to local retangle regions.
    Please look the resulting regions carefully to make sure the function provides the outcomes
    that you want. you can plot them to see what happened.
    
    *parameters:
    1. num_parts_axis: number of parts for both axis: eg: 1, 2, 3, 4, 5, 6, ..., 12, 24;
       based on trials of same structure of code in R, it may suitable for non-interger
       values. However, I did not test it in PYTHON. 
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
    print(len(xy_center_coor)) # length check
    plot_regions_from_whole(num_parts_axis = num_parts_axis, min_x = min_x, min_y = min_y, max_x = max_x, max_y = max_y, total_length = total_length)
    plt.scatter(data["X"], data["Y"], marker='.', s=0.005, c = 'Black')
    plt.show()
    return local_df, xy_boundary, xy_center_coor


def plot_regions_from_whole(num_parts_axis, min_x, min_y, max_x, max_y, total_length):
    '''
    This function can plot the local regions, I add it inside the above function.
    Then, we can see what going on.
    '''
    add = total_length/num_parts_axis
    for vy in range(num_parts_axis+1):
        plt.axhline(y = min_y + vy*add, color = 'r', linestyle = '-')
    for hx in range(num_parts_axis+1):
        plt.axvline(x = min_x + hx*add, color = 'r', linestyle = '-')





def CO21_processing(RA_grid, DEC_grid, gas_density, Dist, PA_gal, inc_gal, RA_galaxy, DEC_galaxy):
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
    #proj_dist = RA_DEC_to_radius(RA_full, DEC_full)
    proj_dist = RA_DEC_to_radius_MUSE(Dist, PA_gal, inc_gal, RA_full, DEC_full, RA_galaxy, DEC_galaxy)
    X_full = np.transpose(RA_DEC_to_xy_MUSE(RA_full, DEC_full, RA_galaxy, DEC_galaxy, PA_gal, inc_gal, Dist))
    
    data_dict_new = {
    	'RA':                   RA_full.values.flatten(),
    	'DEC':                  DEC_full.values.flatten(),
        'X':                    X_full[:,0].flatten(),
        'Y':                    X_full[:,1].flatten(),
        'proj_dist':            proj_dist.flatten(),
        'CO21_mgsd':             CO21_gas_density.values.flatten()
    	}
    data = pd.DataFrame(data_dict_new)
    return data



