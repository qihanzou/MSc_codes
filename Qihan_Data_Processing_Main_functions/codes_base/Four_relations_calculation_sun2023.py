# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 15:53:02 2023

@author: qihan
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 18:40:25 2023

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
    meta=meta_getter('N5236')
    Ha_map = fits_file[8]

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
		SFR density of each spaxel, in units of M_sun/year/kpc^2 !!!
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

def calculate_mKS_sun_2023(SFR_density, a_mks, b_mks):
    '''
    from sun2023:
    SFR density in unit M_sun yr^-1 kpc^-2
    molecular gas density in unit M_sun pc^-2
    
                   a_mks   b_mks    sigma
Fiducial:          -2.40    1.00     0.36
FUV+W4 SFR         -2.34    0.88     0.29
Av-corr Ha SFR     -2.23    0.93     0.29
MW a_co            -2.43    0.92     0.37
B13 a_co           -2.36    1.21     0.38
G20 a_co           -2.22    1.18     0.35
average            -2.33    1.0199999999999998
    '''
    log_molecular_gas_density = np.log10(10) + (np.log10(SFR_density)-a_mks)/b_mks
    molecular_gas_density = 10**(log_molecular_gas_density)
    return molecular_gas_density


def calculate_mES_sun_2023(SFR_density, a_mES, b_mES):
    '''
    from sun2023:
    SFR density in unit M_sun yr^-1 kpc^-2
    molecular gas density/t_orbit in unit M_sun yr^-1 kpc^-2
    
                   a_mES   b_mES    sigma
Fiducial:          -2.23    0.77     0.31
FUV+W4 SFR         -2.19    0.67     0.26
Av-corr Ha SFR     -2.06    0.78     0.28
MW a_co            -2.26    0.69     0.33
B13 a_co           -2.17    0.90     0.32
G20 a_co           -2.11    0.77     0.33
    '''
    log_molecular_gas_density_timescale_orbit = np.log10(0.1) + (np.log10(SFR_density)-a_mES)/b_mES
    molecular_gas_density_timescale_orbit = 10**(log_molecular_gas_density_timescale_orbit)
    return molecular_gas_density_timescale_orbit


def calculate_FFTR_sun_2023(SFR_density, a_FFTR, b_FFTR):
    '''
    from sun2023:
    SFR density in unit M_sun yr^-1 kpc^-2
    molecular gas density/t_ff in unit M_sun yr^-1 kpc^-2
    
                  a_FFTR   b_FFTR    sigma
Fiducial:          -2.32    0.65     0.34
FUV+W4 SFR         -2.28    0.57     0.28
Av-corr Ha SFR     -2.16    0.62     0.28
MW a_co            -2.34    0.62     0.34
B13 a_co           -2.29    0.75     0.36
G20 a_co           -2.20    0.76     0.33
    '''
    log_molecular_gas_density_timescale_ff = (np.log10(SFR_density)-a_FFTR)/b_FFTR
    molecular_gas_density_timescale_ff = 10**(log_molecular_gas_density_timescale_ff)
    return molecular_gas_density_timescale_ff


def calculate_PR_sun_2023(SFR_density, a_PR, b_PR):
    '''
    from sun2023:
    SFR density in unit M_sun yr^-1 kpc^-2
    dynamical equilibrium pressure P_DE in unit k_B K cm^-3
    
                   a_PR   b_PR    sigma
Fiducial:          -2.95    0.93     0.33
FUV+W4 SFR         -2.84    0.84     0.24
Av-corr Ha SFR     -2.72    0.85     0.25
MW a_co            -2.94    0.86     0.33
B13 a_co           -2.95    1.08     0.32
G20 a_co           -2.87    1.05     0.31
    '''
    log_pressure_de = np.log10(10^4) + (np.log10(SFR_density)-a_PR)/b_PR
    pressure_de = 10**(log_pressure_de)
    return pressure_de





'''
    fits_file = fits.open('C:/Users/qihan/Desktop/q/N5236_lowres_cal_1_comp_WCS.fits')
    meta=meta_getter('N5236')
    Ha_map = fits_file[8]
    SFR_density = calculate_SFR_density(Ha_map, meta)
    
    # Fiducial:
    mKS_Fiducial = calculate_mKS_sun_2023(SFR_density, -2.40, 1.00)
    mES_Fiducial = calculate_mES_sun_2023(SFR_density, -2.23, 0.77)
    FFTR_Fiducial = calculate_FFTR_sun_2023(SFR_density, -2.32, 0.65)
    PR_Fiducial = calculate_PR_sun_2023(SFR_density, -2.95, 0.93)
    
    # FUV+W4 SFR: FUV_W4_SFR
    mKS_FUV_W4_SFR = calculate_mKS_sun_2023(SFR_density, -2.34, 0.88)
    mES_FUV_W4_SFR = calculate_mES_sun_2023(SFR_density, -2.19, 0.67)
    FFTR_FUV_W4_SFR = calculate_FFTR_sun_2023(SFR_density, -2.28, 0.57)
    PR_FUV_W4_SFR = calculate_PR_sun_2023(SFR_density, -2.84, 0.84)
    
    # Av-corr Ha SFR: Av_corr_Ha_SFR
    mKS_Av_corr_Ha_SFR = calculate_mKS_sun_2023(SFR_density, -2.23, 0.93)
    mES_Av_corr_Ha_SFR = calculate_mES_sun_2023(SFR_density, -2.06, 0.78)
    FFTR_Av_corr_Ha_SFR = calculate_FFTR_sun_2023(SFR_density, -2.16, 0.62)
    PR_Av_corr_Ha_SFR = calculate_PR_sun_2023(SFR_density, -2.72, 0.85)
    
    # MW a_co: MW_a_co
    mKS_MW_a_co = calculate_mKS_sun_2023(SFR_density, -2.43, 0.92)
    mES_MW_a_co = calculate_mES_sun_2023(SFR_density, -2.26, 0.69)
    FFTR_MW_a_co = calculate_FFTR_sun_2023(SFR_density, -2.34, 0.62)
    PR_MW_a_co = calculate_PR_sun_2023(SFR_density, -2.94, 0.86)
    
    # B13 a_co: B13_a_co
    mKS_B13_a_co = calculate_mKS_sun_2023(SFR_density, -2.36, 1.21)
    mES_B13_a_co = calculate_mES_sun_2023(SFR_density, -2.17, 0.90)
    FFTR_B13_a_co = calculate_FFTR_sun_2023(SFR_density, -2.29, 0.75)
    PR_B13_a_co = calculate_PR_sun_2023(SFR_density, -2.95, 1.08)
    
    # G20 a_co: G20_a_co
    mKS_G20_a_co = calculate_mKS_sun_2023(SFR_density, -2.22, 1.18)
    mES_G20_a_co = calculate_mES_sun_2023(SFR_density, -2.11, 0.77)
    FFTR_G20_a_co = calculate_FFTR_sun_2023(SFR_density, -2.20, 0.76)
    PR_G20_a_co = calculate_PR_sun_2023(SFR_density, -2.87, 1.05)


'''






