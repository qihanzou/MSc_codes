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
		SFR density of each spaxel, in units of mass_sun/year/kpc^2
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


def calculate_gas_density(SFR_density):
    '''
    By the Kennicutt_Schmidt law, a empirical relation between 
    the surface gas density and the SFR:
    # SFR_denisty must in M_sun yr−1 kpc-2, the result of above part is in our target unit.
    # result gas density in M_sun pc-2
    '''
    gas_density = (SFR_density / (2.5*10**(-4)))**(1/(1.4))
    return gas_density

def molecular_gas_surface_density(SFR_density):
    '''
    # SFR_denisty must in M_sun yr−1 kpc-2, the result of above part is in our target unit.
    # result gas density in M_sun pc-2
    '''
    gas_density = (SFR_density / (10**(-2.1)))*10
    return gas_density

# gas_density.shape
# gas_density.size

def calculate_sum_of_gas_surface_density(gas_density):
    total_gas_surface_density = np.sum(gas_density)
    return total_gas_surface_density 