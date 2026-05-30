# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:02:32 2024

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
from ALMA_data_processing import *


gal_name = 'N1300' # need to change it! name: N + number, not NGC + number!
data = fits.open(r"C:\Users\qihan\Desktop\ngc1300_12m+7m+tp_co21_15as_strict_mom0.fits") 

# you need to change the path:
galaxy_info = pd.read_excel("C:/Users/qihan/Desktop/Data processing/main_files/Gal_info_data/galaxydata.xlsx")
idx = galaxy_info["Gal_ID"] == gal_name
galaxy_info_now = galaxy_info[idx]
inc_angle = float(galaxy_info_now['i'])
Distance = float(galaxy_info_now['D'])
RA_galaxy = float(galaxy_info_now['RA'])
DEC_galaxy = float(galaxy_info_now['DEC'])
PA_galaxy = float(galaxy_info_now['PA_MUSE'])


# check path before use
ICO21 = data[0] # unit: K km s-1 from PHANGS-ALMA readme.
x_dim = data[0].data.shape[1] # dim of x axis
y_dim = data[0].data.shape[0] # dim of y axis
RA_grid, DEC_grid = make_RA_DEC_grid_MUSE(data[0].header, x_dim, y_dim) 


gas_density = compute_molecular_gas_surface_density_from_CO21(ICO21, inc_angle)
data = CO21_processing(RA_grid, DEC_grid, gas_density, Distance, PA_galaxy, inc_angle, RA_galaxy, DEC_galaxy)


# use to save the file as pkl file:
data.to_pickle('C:/Users/qihan/Desktop/'+'N1300_mgsd_15as'+'.pkl')






