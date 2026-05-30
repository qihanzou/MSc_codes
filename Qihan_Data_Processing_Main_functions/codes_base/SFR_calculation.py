# -*- coding: utf-8 -*-
"""
Created on Fri Oct  6 19:28:48 2023

@author: qihan
"""
import numpy as np
from numpy import *
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from TYPHOON_wrangling import *

import pandas as pd 
from astropy.wcs import WCS
import astropy.units as u
# from sklearn.metrics.pairwise import euclidean_distances 
# from statsmodels.regression.linear_model import GLS
# from scipy.linalg import cho_factor, cho_solve
from extinction import ccm89, apply
from Z_diags import *


def calculate_SFR_density(Ha_map, meta):
	'''
    fits_file = fits.open('C:/Users/qihan/Desktop/q/N5236_lowres_cal_1_comp_WCS.fits')
    # meta=meta_getter('N5236')
    # Ha_map = fits_file[8]
    # Ha_map = fits_file[18] check it in the new version
    # Let sfrd = calculate_SFR_density(Ha_map, meta)
    # plt.figure()
    # plt.imshow(sfrd,origin = 'lower', norm = LogNorm(), cmap= 'Greys') # camp='Greys' means Grey scale, if not use it, will be color.
    
	Determine the star formation rate using Kennicutt+Evans98:
	https://ui.adsabs.harvard.edu/abs/2012ARA%26A..50..531K/abstract
	
	Formula: log(SFR) = log(LHα) − 41.27
    LHa is in unit erg s−1
    SFR is in M_sun yr−1
	
	where SB_Ha is in units of erg s−1 kpc−2
	
	Parameters
	----------
	
	Ha_map: hdu
		the intensity of the Ha line data for this galaxy
	
	meta: dict-like object including:
	 D - distance to an object in Mpc.
	 i - inclination in degrees
	
	Returns
	-------
	
	SFR_density_map: np array
		SFR density of each spaxel, in units of solar masses/year/kpc^2
	'''
	# Finagle out the deprojected area (units: kpc)
	world = WCS(Ha_map.header)
	pix_solid_angle = world.proj_plane_pixel_area().to(u.steradian).value
	plane_area	= pix_solid_angle * (meta['D']*1000)**2 # in kpc^2
	i  = np.radians(meta['i'])
	deproj_area = plane_area / np.cos(i)
	# convert flux (units:1e-17 erg/s/cm2) to Luminosity (units: erg/s)
	log_Ha_luminosity_map = np.log10(Ha_map.data) + np.log10(4*np.pi*meta['D']**2 * 9.5234) + 31
	# convert to surface brightnesss (units: erg/s/kpc^2)
	log_Ha_SB_map = log_Ha_luminosity_map - np.log10(deproj_area)
	# Constant from Kennicutt+Evans12
	log_SFR_density_map   = log_Ha_SB_map - 41.27
	return 10**(log_SFR_density_map)


def calculate_deproj_area(Ha_map, meta):
	'''
	Parameters
	----------
	Ha_map: hdu
		the intensity of the Ha line data for this galaxy
	
	meta: dict-like object including:
	 D - distance to an object in Mpc.
	 i - inclination in degrees
	
	Returns
	-------
	SFR_density_map: np array
		SFR density of each spaxel, in units of solar masses/year/kpc^2
	'''
	# Finagle out the deprojected area (units: kpc)
	world = WCS(Ha_map.header)
	pix_solid_angle = world.proj_plane_pixel_area().to(u.steradian).value
	plane_area	= pix_solid_angle * (meta['D']*1000)**2 # in kpc^2
	i  = np.radians(meta['i'])
	deproj_area = plane_area / np.cos(i)
	return deproj_area


def calculate_total_SFR_of_galaxy(SFR_density,deproj_area):
    '''
    ----------
    Parameters:
    ----------
    SFR_density : result from the function "calculate_SFR_density(Ha_map, meta)"
    deproj_area : result from the function "calculate_deproj_area(Ha_map, meta)"
    -------
    Returns:
    -------
    Total_SFR_Of_Galaxy : A number, total SFR of target galaxy, 
                          in units of solar masses/year, we can 
                          use this number to compare the total SFR
                          from other resourses to check whether the 
                          local SFR we got is reasonable.
    '''
    # Find positions of NaNs and remove all of them to be 0.
    positions_are_NaNs_in_SFR_density = isnan(SFR_density)
    SFR_density[positions_are_NaNs_in_SFR_density] = 0
    # sum over the local SFR to get the Total SFR
    Total_SFR_Of_Galaxy = np.sum(SFR_density * deproj_area)
    return Total_SFR_Of_Galaxy

def calculate_SFR_density_from_SFR(SFR,):
    
    '''
    SFR_density.shape
    Out[162]: (666, 246)

    4.1/(deproj_area*666*246)
    Out[163]: 0.01591585852392809

    4.8/(deproj_area*666*246)
    Out[164]: 0.018633200223135325

    np.mean(SFR_density)
    Out[145]: 0.018758316
    '''

'''
line_df_1 = fits.open('C:/Users/qihan/Desktop/q/N5236_lowres_cal_1_comp_WCS.fits')
wavelengths = np.array([3726.0, 3729.0, 4861.3, 5007.0, 6562.8, 6583.0, 6716.0, 6731.0])
#Ha, 6562
line_df_2 = SN_cut(line_df_1, threshold=3)
line_df_3 = SN_cut_O2(line_df_2, threshold=1)
line_df_4= extinction_correction(line_df_3, wavelengths, R_V=3.1)

Ha_map = line_df_4[8]
meta = meta_getter('N5236')

SFR_density=calculate_SFR_density(Ha_map, meta)

#positions_are_NaNs_in_SFR_density = isnan(SFR_density)
#SFR_density[positions_are_NaNs_in_SFR_density] = 0
#SFR_density = SFR_density[~np.isnan(SFR_density)]
#np.mean(SFR_density)
#np.sum(SFR_density)

deproj_area=calculate_deproj_area(Ha_map, meta)


calculate_total_SFR_of_galaxy(SFR_density,deproj_area)
# the results for N5236 is 4.8322296 (0.6841475613010073) by our code. compare to 10**0.62, 4.168693834703354
'''
	

def calculate_SFR_density_from_fits_file(fits_file):
    # only work for N5236 now.
    wavelengths = np.array([3726.0, 3729.0, 4861.3, 5007.0, 6562.8, 6583.0, 6716.0, 6731.0])
    line_df_2 = SN_cut(fits_file, threshold=3)
    line_df_3 = SN_cut_O2(line_df_2, threshold=1)
    line_df_4= extinction_correction(line_df_3, wavelengths, R_V=3.1)
    Ha_map = line_df_4[8]
    meta = meta_getter('N5236')
    SFR_density=calculate_SFR_density(Ha_map, meta)
    positions_are_NaNs_in_SFR_density = isnan(SFR_density)
    SFR_density[positions_are_NaNs_in_SFR_density] = 0
    return SFR_density


def calculate_total_SFR_from_fits_file(fits_file):
    # only work for N5236 now.
    wavelengths = np.array([3726.0, 3729.0, 4861.3, 5007.0, 6562.8, 6583.0, 6716.0, 6731.0])
    line_df_2 = SN_cut(fits_file, threshold=3)
    line_df_3 = SN_cut_O2(line_df_2, threshold=1)
    line_df_4= extinction_correction(line_df_3, wavelengths, R_V=3.1)
    Ha_map = line_df_4[8]
    meta = meta_getter('N5236')
    SFR_density=calculate_SFR_density(Ha_map, meta)
    deproj_area=calculate_deproj_area(Ha_map, meta)
    Total_SFR=calculate_total_SFR_of_galaxy(SFR_density,deproj_area)
    return Total_SFR