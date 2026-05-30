# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 10:21:25 2024

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
from dec_ra_to_xy import *

meta = meta_getter('N5236') 

def calculate_gas_density_CO21(ICO21, i_deg):
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
    i  =  np.radians(i_deg) 
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
        'proj_dist':            proj_dist.flatten(),
        'CO21_mgd':             CO21_gas_density.values.flatten()
    	}
    data = pd.DataFrame(data_dict_new)
    return data

data = fits.open(r"C:\Users\qihan\Desktop\q\ALMA_CO21\N5236_co21\N5236_co21_2as_strict_mom0.fits") 
# check path before use
ICO21 = data[0] # unit: K km s-1 from PHANGS-ALMA readme.
RA_grid, DEC_grid = make_RA_DEC_grid(data[0].header) 
i_deg = 12.50496329
gas_density = calculate_gas_density_CO21(ICO21, i_deg)
data = Co21_data_processing(RA_grid, DEC_grid, gas_density)

#data.to_pickle('C:/Users/qihan/Desktop/'+'N5236_CO21_molecular_gas_surface_density'+'.pkl')


import numpy as np
from scipy.spatial import KDTree


points = pd.concat([data['X'],data['Y']], axis=1)
tree = KDTree(points)
rows_to_fuse = tree.query_pairs(r=0.05)    


points[list(rows_to_fuse)]







