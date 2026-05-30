# -*- coding: utf-8 -*-
"""
Created on Tue Nov 28 17:44:41 2023

@author: qihan
"""

'''
Galaxies show a strong correlation between their SFR and stellar mass M.
The SFR-M relation.
'The star formation mass sequence out to z=2.5.'
'''

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

import scipy.io


def calculate_SFR_density(Ha_map, meta):
	'''
    fits_file = fits.open('C:/Users/qihan/Desktop/q/N5236_lowres_cal_1_comp_WCS.fits')
    # meta=meta_getter('N5236')
    # Ha_map = fits_file[8]

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


def calculate_SFR_of_galaxy(SFR_density,deproj_area):
    '''
    ----------
    Parameters:
    ----------
    SFR_density : result from the function "calculate_SFR_density(Ha_map, meta)"
    deproj_area : result from the function "calculate_deproj_area(Ha_map, meta)"
    -------
    Returns:
    -------
    SFR_Of_Galaxy : A number, SFR of target positions, in units of solar masses/year.
    '''
    # from the local SFR to get the SFR
    SFR_Of_Galaxy = SFR_density * deproj_area
    return SFR_Of_Galaxy

def calculate_stellar_mass_whitaker2012star(SFR_Of_Galaxy,z):
    '''
    By the SFR-M relation, a fitting relation between 
    the stellar mass and the SFR:
    # SFR must in M_sun yr−1.
    return: stellar_mass in unit M_sun
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    alpha_z = 0.70 + 0.13*z
    beta_z = 0.38 + 1.14*z - 0.19*(z**2)
    log_stellar_mass = 10.5 + ((log_SFR-beta_z)/alpha_z)
    stellar_mass = 10**(log_stellar_mass)
    return stellar_mass

def calculate_total_mass_of_galaxy(stellar_mass):
    '''
    Returns:
    -------
    Total_mass_Of_Galaxy : A number, total mass of target galaxy, 
                          in units of solar masses.
    '''
    # Find positions of NaNs and remove all of them to be 0.
    positions_are_NaNs_in_stellar_mass = isnan(stellar_mass)
    stellar_mass[positions_are_NaNs_in_stellar_mass] = 0
    Total_mass_Of_Galaxy = np.sum(stellar_mass)
    return Total_mass_Of_Galaxy




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
deproj_area=calculate_deproj_area(Ha_map, meta)
SFR_Of_Galaxy=calculate_SFR_of_galaxy(SFR_density,deproj_area)
# redshift z = 0.001721+-0.000013
stellar_mass=calculate_stellar_mass(SFR_Of_Galaxy,0.001721)
stellar_mass.shape
stellar_mass.size
Total_mass_Of_Galaxy=calculate_total_mass_of_galaxy(stellar_mass)

scipy.io.savemat('stellar_massdata.mat',{'name':stellar_mass})


import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
plt.figure()
plt.imshow(stellar_mass,origin = 'lower', norm = LogNorm(), cmap= 'Greys') # camp='Greys' means Grey scale, if not use it, will be color.
plt.colorbar()
plt.show()


log_SFR = np.log10(4.8322296)
alpha_z = 0.70 + 0.13*0.001721
beta_z = 0.38 + 1.14*0.001721 - 0.19*(0.001721**2)
log_stellar_mass = 10.5 + ((log_SFR-beta_z)/alpha_z)
stellar_mass = 10**(log_stellar_mass)
'''
'''
https://iopscience.iop.org/article/10.3847/0004-637X/817/2/118
'''
def calculate_stellar_mass_tomczak2016sfr_all_galaxies(SFR_Of_Galaxy,z):
    '''
    By the SFR-M relation, a fitting relation between 
    the stellar mass and the SFR:
    # SFR must in M_sun yr−1.
    return: stellar_mass in unit M_sun
    
    12405054000
    10.093598
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    s0 = 0.195 +1.157*z-0.143*(z**2)
    log_M0= 9.244+0.753*z-0.09*(z**2)
    para_gamma = 1.118
    stellar_mass = (10**log_M0)*(((10**(s0-log_SFR))-1)**(-1/para_gamma))
    return stellar_mass

def calculate_stellar_mass_tomczak2016sfr_star_forming_galaxies(SFR_Of_Galaxy,z):
    '''
    By the SFR-M relation, a fitting relation between 
    the stellar mass and the SFR:
    # SFR must in M_sun yr−1.
    return: stellar_mass in unit M_sun
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    s0 = 0.448 +1.220*z-0.174*(z**2)
    log_M0= 9.458+0.865*z-0.132*(z**2)
    para_gamma = 1.091
    stellar_mass = (10**log_M0)*(((10**(s0-log_SFR))-1)**(-1/para_gamma))
    return stellar_mass

def calculate_stellar_mass_2(SFR_Of_Galaxy,z):
    '''
    By the SFR-M relation, a fitting relation between 
    the stellar mass and the SFR:
    # SFR must in M_sun yr−1.
    return: stellar_mass in unit M_sun
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    s0 = 0.8
    log_M0= 10.03
    para_gamma = 0.92
    stellar_mass = (10**log_M0)*(((10**(s0-log_SFR))-1)**(-1/para_gamma))
    return stellar_mass
'''
log_SFR = 0.62
s0 = 0.195 +1.157*0.001721-0.143*(0.001721**2)
log_M0= 9.244+0.753*0.001721-0.09*(0.001721**2)
para_gamma = 1.118
stellar_mass = (10**log_M0)*(((10**(s0-log_SFR))-1)**(-1/para_gamma))
'''


'''
A HIGHLY CONSISTENT FRAMEWORK FOR THE EVOLUTION OF THE
STAR-FORMING “MAIN SEQUENCE” FROM z ∼ 0–6

speagle2014highly

https://ui.adsabs.harvard.edu/abs/2014ApJS..214...15S/abstract
'''

def calculate_stellar_mass_E11(SFR_Of_Galaxy,z):
    ''' 
    19237482000.0
    10.284148
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    log_stellar_mass=(log_SFR-(-9.6))/1
    stellar_mass = 10**(log_stellar_mass)
    return stellar_mass

def calculate_stellar_mass_C14_1(SFR_Of_Galaxy,z):
    '''
    56906440.0
    7.7551613
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    log_stellar_mass=(log_SFR-(-4.57))/0.477
    stellar_mass = 10**(log_stellar_mass)
    return stellar_mass

def calculate_stellar_mass_O10(SFR_Of_Galaxy,z):
    '''
    11371965000.0
    10.055835514373113
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    log_stellar_mass=(log_SFR-(-7.88))/0.77
    stellar_mass = 10**(log_stellar_mass)
    return stellar_mass

def calculate_stellar_mass_C09_1(SFR_Of_Galaxy,z):
    '''
    15655616.0
    7.194670160431863
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    log_stellar_mass=(log_SFR-(-3.56))/0.35
    stellar_mass = 10**(log_stellar_mass)
    return stellar_mass

def calculate_stellar_mass_S07(SFR_Of_Galaxy,z):
    '''
    952629400.0
    8.978923980556571
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    log_stellar_mass=(log_SFR-(-6.33))/0.65
    stellar_mass = 10**(log_stellar_mass)
    return stellar_mass

def calculate_stellar_mass_W12(SFR_Of_Galaxy,z):
    '''
    687408600.0
    8.837214961186689
    '''
    log_SFR = np.log10(SFR_Of_Galaxy)
    log_stellar_mass=(log_SFR-(-6.36))/0.67
    stellar_mass = 10**(log_stellar_mass)
    return stellar_mass



