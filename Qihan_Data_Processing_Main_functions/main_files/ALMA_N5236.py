# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:06:10 2024

@author: qihan
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
import matplotlib.pyplot as plt

musk = fits.open(r"C:\Users\qihan\Desktop\Data processing\main_files\ALMA_fits_data\N5236_co21\ngc5236_12m+7m+tp_co21_2as_strictmask.fits")
plt.imshow(sum(musk[0].data))



musk2 = fits.open(r"C:\Users\qihan\Desktop\Data processing\main_files\ALMA_fits_data\N5236_co21\ngc5236_12m+7m+tp_co21_broadmask.fits")
plt.imshow(sum(musk2[0].data))


data_2as = fits.open(r"C:\Users\qihan\Desktop\Data processing\main_files\ALMA_fits_data\N5236_co21\N5236_co21_2as_strict_mom0.fits")
data_2as.info()
data_2as[0].header()

data_gal = data_2as[0].data
data_gal[np.isnan(data_gal)] = 0
plt.imshow(data_gal)
plt.colorbar()


i  =  np.radians(12.5) 
gas_density = 6.7*data_gal*np.cos(i) 
plt.imshow(np.log(gas_density), cmap ='hot_r')
plt.colorbar()
np.count_nonzero(gas_density)


gas_density2 = gas_density
gas_density2[gas_density2!=0] = 1
plt.imshow(gas_density2)
plt.colorbar()




data_15as = fits.open(r"C:\Users\qihan\Desktop\Data processing\main_files\ALMA_fits_data\N5236_co21\N5236_co21_15as_strict_mom0.fits")

data_15 = data_15as[0].data
data_15[np.isnan(data_15)] = 0
plt.imshow(data_gal)
plt.colorbar()

i  =  np.radians(12.5) 
gas_density_15 = 6.7*data_15*np.cos(i) 
plt.imshow(np.log(gas_density_15), cmap ='hot_r')
plt.colorbar()

gas_density2_15 = gas_density_15
gas_density2_15[gas_density2_15!=0] = 1
plt.imshow(gas_density2_15)
plt.colorbar()

np.count_nonzero(gas_density_15)





