# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 16:40:05 2024

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
from CO21_functions import *
import skgstat as skg


data = fits.open(r"C:\Users\qihan\Desktop\q\ALMA_CO21\N5236_co21\N5236_co21_2as_strict_mom0.fits") 
# check path before use
ICO21 = data[0] # unit: K km s-1 from PHANGS-ALMA readme.
RA_grid, DEC_grid = make_RA_DEC_grid(data[0].header) 
meta = meta_getter('N5236') 
gas_density = calculate_gas_density_CO21(ICO21, meta)
# data processing:
data = Co21_data_processing(RA_grid, DEC_grid, gas_density)
local_df, xy_boundary, xy_center_coor = create_regions_from_whole(num_parts_axis = 12, min_x = -6, min_y = -6, max_x = 6, max_y = 6, total_length = 12, data = data)

# Local models
thetahat_list, betahat_list, Sigma_list, error_list, inv_Sigma_list = construct_geo_dist_models_for_regions(local_df)


data_galaxy = pd.read_pickle('C:/Users/qihan/Desktop/q/N5236_25_1_2024.pkl')
index_above_0_Z = data_galaxy["Z_O3S2_kumari_N2"] > 0 # ignore NA
data_galaxy = data_galaxy[index_above_0_Z]
index_below_1_BPT = data_galaxy["N2_BPT"] >0  # BPT cut for Hii regions
data_galaxy = data_galaxy[index_below_1_BPT] # result: 7193 rows of data
#index_above_09_CHII = data_galaxy["S2_DIG"] > 0.9
#data_galaxy = data_galaxy[index_above_09_CHII]
RA_galaxy_list = data_galaxy['RA']
DEC_galaxy_list = data_galaxy['DEC']
proj_dist_galaxy_list = data_galaxy['proj_dist']
coor_full_galaxy = RA_DEC_to_xy(RA_galaxy_list, DEC_galaxy_list, meta)
coor_full_galaxy = np.transpose(coor_full_galaxy)
X_galaxy = coor_full_galaxy[:,0]
Y_galaxy = coor_full_galaxy[:,1]


        
inference_gas_density = pred_gas_density(error_list, Sigma_list, betahat_list, thetahat_list, xy_boundary, local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x = 6, max_y = 6, min_x = -6, min_y = -6)
inference_gas_weight2 = pred_gas_density_weight(xy_center_coor, error_list, Sigma_list, betahat_list, thetahat_list, xy_boundary, 2, local_df, X_galaxy, Y_galaxy, coor_full_galaxy, proj_dist_galaxy_list,DEC_galaxy_list, RA_galaxy_list, inv_Sigma_list, data_galaxy, max_x = 6, max_y = 6, min_x = -6, min_y = -6)

'''
y = inference_gas_density
x = data_galaxy[["proj_dist"]]
plt.scatter(x, y, s=0.1)

y = exp(inference_gas_density)
x = data_galaxy[["proj_dist"]]
plt.scatter(x, y, s=0.1)  
'''  

#sum(math.isnan(x) for x in inference_gas_density)

'''
out_path = 'C:/Users/qihan/Desktop/'
data_dict = {
    'infer_gas_density':  pd.DataFrame(inference_gas_density).values.flatten(),
    'inference_gas_weight2':  pd.DataFrame(inference_gas_weight2).values.flatten() 
	}
result_df = pd.DataFrame(data_dict)
result_df.to_pickle(out_path+'N5236_CO21_RS32_test'+'.pkl')
'''



data_dict_CO21 = {
	'X':                 data['X'].values.flatten(),
	'Y':                 data['Y'].values.flatten(),
	'CO21_gas_density':  data['CO21_gas_density'].values.flatten(), 
	'RA':                data['RA'].values.flatten(),
	'DEC':               data['DEC'].values.flatten(),
	'proj_dist':         data['proj_dist'].values.flatten()
	}

result_df = pd.DataFrame(data_dict_CO21)
#result_df.to_pickle('C:/Users/qihan/Desktop/'+'N5236_CO21_orginal'+'.pkl')


