# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 21:39:44 2024

@author: Qihan Zou

These codes are mainly written for PHANGES MUSE dataset.

Last updated: 25/06/2024 
"""
import numpy as np
from numpy import *
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pandas as pd 
from astropy.wcs import WCS
import astropy.units as u
from extinction import ccm89, apply
from Z_diags import *
from TYPHOON_wrangling import *
from AST2 import *
from Z_diags import *
import numpy as np 
from scipy import linalg
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
from sklearn.metrics.pairwise import euclidean_distances 
from scipy.special import kv
from scipy.special import gamma
from Internal_Fun import *
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from scipy import linalg
from matplotlib import colors
from scipy.stats import linregress
from sklearn.model_selection import ShuffleSplit
import skgstat as skg
from skgstat import models
from collections import defaultdict


STARBURST = 0
SEYFERT   = 1
LINER     = 2

##############################################################################
def make_RA_DEC_grid_MUSE(header,x_dim,y_dim):
    '''
    Given a header file, create a grid of RA//DEC for each pixel in that file.
    '''
    world = WCS(header)
    x = np.arange(x_dim)
    y = np.arange(y_dim)
    X, Y = np.meshgrid(x, y)
    RA_grid, DEC_grid = world.wcs_pix2world(X, Y, 0)
    return RA_grid, DEC_grid  
###############################################################################
#### This one should be OK now.
def SN_cut_MUSE(line_df, threshold=3):
    '''
    Replace all spaxels with SN<3 in a certain line with NANs.
    lines_df: hdu list, lines where S/N < threshold replaced with np.nan
    '''
    n_lines = int(len(line_df)/2)
    x_max, y_max = line_df[1].data.shape
    for l in range(8): # range(8): 0 1 2 3 4 5 6 7
        signal = line_df[6*(l+1)].data
        noise  = line_df[6*(l+1)+1].data
        too_low = signal <= threshold*noise
        for ii in range(x_max):
            for jj in range(y_max):
                # replace low signals/no signals with NANs.
                if too_low[ii,jj]:
                    signal[ii,jj] = np.nan
                    noise[ii,jj]  = np.nan
    return line_df

##############################################################################
def determine_DIG_Ha_Zhang17_MUSE(Ha_map, inc_angle, Distance):
    '''
    Determine whether a spaxel is Hii/DIG dominated, by applying a surface
    brightness cut in Ha.
    
    Formula: dig if log_10( SB_Ha ) < 39
    
    SB_Ha is in units of erg s−1 kpc−2  
    Ha_map: the intensity of the Ha line data for this galaxy
        
    DIG_map: np array
        0 if a spaxel is DIG dominated
        1 if it's a Hii region
        nan if SN too low to decide.
    '''
    # Finagle out the deprojected area (units: kpc)
    world = WCS(Ha_map.header)
    pix_solid_angle = world.proj_plane_pixel_area().to(u.steradian).value
    plane_area  = pix_solid_angle * (Distance*1000)**2 # in kpc^2
    i  = np.radians(inc_angle)
    deproj_area = plane_area / np.cos(i)
    # convert flux (units:1e-20 erg/s/cm2) to Luminosity (units: erg/s)
    log_Ha_luminosity_map = np.log10(Ha_map.data) + np.log10(4*np.pi*(Distance*3.086*10**24)**2) - 20 # Here been changed.
    # convert to surface brightnesss (units: erg/s/kpc^2)
    log_Ha_SB_map = log_Ha_luminosity_map - np.log10(deproj_area)
    is_Hii = log_Ha_SB_map > 39
    is_nan = np.isnan(log_Ha_SB_map)
    DIG_map = np.zeros(log_Ha_SB_map.shape)+ 1.0*is_Hii
    X_MAX, Y_MAX = log_Ha_SB_map.shape
    for ii in range(X_MAX):
        for jj in range(Y_MAX):
            if is_nan[ii,jj]:
                DIG_map[ii,jj] = np.nan
    return DIG_map
##############################################################################
def extinction_correction_MUSE(line_df, wavelengths, R_V=3.1):
	'''
	Parameters
	----------		
	wavelengths: np.array, Wavelength of each of the 8 lines in this data cube, in Angstroms.
    For MUSE data wavelength:
    wavelengths = np.array([4861.35, 4958.91, 5006.84, 6548.05, 6562.79, 6583.45, 6716.44, 6730.82])
	R_V: float, The free parameter in ccm89 extinction law. Set (kept) at 3.1.
	
	Returns
	-------
	corrected_lines_df: hdu list, Corrections for all lines using the calibration of ccm89.
	'''
	Ha_map = line_df[30].data
	Hb_map = line_df[6].data
	# To convert balmer decrement to extinction, need these...
	HA_EXT =  ccm89(np.array([6562.8]), 1.0, R_V)[0]
	HB_EXT =  ccm89(np.array([4861.3]), 1.0, R_V)[0]
	Ha_Hb_ratio	 = Ha_map/Hb_map
	balmer_decrement = 2.5*np.log10(Ha_Hb_ratio / 2.86)
	A_V = balmer_decrement/(HB_EXT - HA_EXT) 
	A_V_positive = A_V * (A_V > 0) # sets negatives to zero
	
	# Use this to correct obs and error for each wavelength
	for l in range(8):
		extinction_at_wav = ccm89(wavelengths[l:l+1], 1, R_V)[0]
		extinction_map = extinction_at_wav*A_V_positive
		# correct signal and noise
		line_df[6*(l+1)].data	 = line_df[6*(l+1)].data * 10**(0.4 * extinction_map)
		line_df[6*(l+1)+1].data	 = line_df[6*(l+1)+1].data * 10**(0.4 * extinction_map)
	return line_df
##############################################################################
def classify_S2_BPT_MUSE(line_df):
    ''' 
    For each spaxel, specify whether it is SEYFERT, LINER, or SF
    using the diagnostics of Kewley+01 and Kewley+06 and the S2-BPT diagram.

    S2_BPT_classification: np array
        NAN if line data is missing/has too low S/N
        0 if starburst
        1 if Seyfert
        2 if LINER
    
    1s/2s will be treated as DIG
    '''
    #line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    #O3Hb = np.log10( line_df[line_IDs.index('OIII5006_FLUX')].data /  line_df[line_IDs.index('HB4861_FLUX')].data )
    #S2Ha = np.log10( (line_df[line_IDs.index('SII6716_FLUX')].data+line_df[line_IDs.index('SII6730_FLUX')].data)/line_df[line_IDs.index('HA6562_FLUX')].data    )   
    O3Hb = np.log10( line_df[18].data /  line_df[6].data )
    S2Ha = np.log10( (line_df[42].data+line_df[48].data)/line_df[30].data)    
    is_starburst = O3Hb < ( 0.72/(S2Ha-0.32) + 1.3 )
    is_seyfert   = (1.89*S2Ha +0.76) < O3Hb
    is_liner     = (1.89*S2Ha +0.76) >= O3Hb # Neither implies it's a nan.
    S2_BPT_classification = np.ones(O3Hb.shape) * np.nan
    X_MAX, Y_MAX = O3Hb.shape
    for ii in range(X_MAX):
        for jj in range(Y_MAX):
            if is_starburst[ii,jj] and S2Ha[ii,jj] < 0.32:
                S2_BPT_classification[ii,jj] = STARBURST
            elif is_seyfert[ii,jj]:
                S2_BPT_classification[ii,jj] = SEYFERT
            elif is_liner[ii,jj]:
                S2_BPT_classification[ii,jj] = LINER
    return S2_BPT_classification

def classify_S2_BPT_MUSE_ver2(line_df):
    S2Ha = np.log10((line_df[42].data+line_df[48].data)/line_df[30].data) 
    O3Hb = np.log10( line_df[18].data /  line_df[6].data ) 
    is_starburst = O3Hb < ( 0.72/(S2Ha-0.32) + 1.3 ) 
    return is_starburst & (S2Ha < 0.32)
##############################################################################
def classify_N2_BPT_MUSE(line_df, rule="Kauffmann03"):
    '''
    For each spaxel, specify whether it is LINER or SF, using the diagnostic of Kewley+01, and the N2-BPT diagram.
        
    N2_BPT_classification: np array
        NAN if line data is missing/has too low S/N
        0 if starburst
        2 if LINER
    
    2s will be treated as DIG
    '''
    #line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    #O3Hb = np.log10( line_df[line_IDs.index('OIII5006_FLUX')].data/line_df[line_IDs.index('HB4861_FLUX')].data )
    #N2Ha = np.log10( line_df[line_IDs.index('NII6583_FLUX')].data/line_df[line_IDs.index('HA6562_FLUX')].data    )   
    O3Hb = np.log10( line_df[18].data/line_df[6].data )
    N2Ha = np.log10( line_df[36].data/line_df[30].data )   
    if rule=='Kewley01':
        is_starburst = O3Hb < 0.61/(N2Ha-0.47) + 1.19
        is_LINER     = O3Hb >= 0.61/(N2Ha-0.47) + 1.19 # otherwise it's a NAN
    elif rule=='Kauffmann03':
        is_starburst = (O3Hb < 0.61/(N2Ha-0.05) + 1.3) 
        is_LINER     = (O3Hb > 0.61/(N2Ha-0.05) + 1.3)
    else:
        print("Error: classsify_N2_BPT only works when 'rule' is either 'Kewley01' or 'Kauffmann03'.")
        exit(1)
        return None
    N2_BPT_classification = np.ones(O3Hb.shape) * np.nan
    X_MAX, Y_MAX = O3Hb.shape
    for ii in range(X_MAX):
        for jj in range(Y_MAX):
            if is_starburst[ii,jj] and N2Ha[ii,jj] < 0.05:
                N2_BPT_classification[ii,jj] = STARBURST
            elif is_LINER[ii,jj]:
                N2_BPT_classification[ii,jj] = LINER
    return N2_BPT_classification

def classify_N2_BPT_MUSE_ver2(line_df, rule="Kauffmann03"):
    O3Hb = np.log10( line_df[18].data/line_df[6].data )
    N2Ha = np.log10( line_df[36].data/line_df[30].data ) 
    if rule=='Kewley01': 
        is_starburst = O3Hb < 0.61/(N2Ha-0.47) + 1.19
        is_LINER	 = O3Hb >= 0.61/(N2Ha-0.47) + 1.19 # otherwise it's a NAN
    elif rule=='Kauffmann03':
        is_starburst = (O3Hb < 0.61/(N2Ha-0.05) + 1.3) 
        is_LINER	 = (O3Hb > 0.61/(N2Ha-0.05) + 1.3)
    else:
        print("Error: classsify_N2_BPT only works when 'rule' is either 'Kewley01' or 'Kauffmann03'.") 
        exit(1) 
        return None
    return is_starburst & (N2Ha < 0.05)
###############################################################################
def determine_DIG_S2_Kaplan16_MUSE(line_df, n_spaxels=100, max_prop=0.05):
    '''
    Assuming that:
    1. The Sii/Ha line ratio is significantly different for DIG/Hii regions;
    2. The intrinsic distributions of Sii/Ha are (infinitely) narrow for purely Hii/DIG regions
    
    Compute the fraction of Ha-light originating from Hii regions for each spaxel
    (C_Hii), using the formalism of Kaplan+2016:
    https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1642K
    
    The formula:
    
    [Sii/Ha] = C_Hii [Sii/Ha]_Hii + C_DIG [Sii/Ha]_DIG
    C_Hii + C_DIG = 1
    [Sii/Ha]_Hii = the median [Sii/Ha] of the brightest spaxels
    [Sii/Ha]_DIG = the median [Sii/Ha] of the faintest spaxels
    
    Parameters
    ----------
    
    lines_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
        
    n_spaxels: int
        A hyperparameter. Set to be 100 fiducually, it sets the number of spaxels
        used to compute intrinsic Hii or DIG Sii/Ha values.
        
    max_prop: float
        only use up to [max_prop] of spaxels to compute the intrinsic Hii or DIG 
        Sii/Ha values.
        
    Returns
    -------
    
    C_Hii: np. array
        Fraction of Ha light emanating from Hii regions, for each spaxel.
        
    Does not return errors associated with this calculation.
    '''
    # open wanted line data/get wanted line ratios
    #line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    #Ha = line_df[line_IDs.index('HA6562_FLUX')].data
    #S2Ha = np.log10( (line_df[line_IDs.index('SII6716_FLUX')].data+line_df[line_IDs.index('SII6730_FLUX')].data)/line_df[line_IDs.index('HA6562_FLUX')].data    )   
    Ha = line_df[30].data
    S2Ha = np.log10( (line_df[42].data+line_df[48].data)/line_df[30].data)   
    # only consider spaxels where both of these measures exceed the S/N theshold
    useful_spaxels = ~np.isnan(S2Ha)
    n_useful_spaxels = np.sum(useful_spaxels)
    useful_Ha = Ha[useful_spaxels]
    useful_S2Ha = S2Ha[useful_spaxels]
    n_per_group = np.min((n_spaxels, int(max_prop*n_useful_spaxels)))
    if n_per_group < n_spaxels:
        print("Computing median Sii/Ha for DIG/Hii regions using {0} of spaxels (5%) per group...".format(n_per_group))
    brightness_order = np.argsort(useful_Ha)
    # Assume the least bright spaxels are pure DIG
    # Compute their median S2Ha
    DIG_S2Ha = np.median(useful_S2Ha[brightness_order][:n_per_group])
    # Do the same for Hii
    Hii_S2Ha = np.median(useful_S2Ha[brightness_order][-1*n_per_group:])
    # Using these values, convert S2Ha to C_Hii
    C_Hii = (S2Ha - DIG_S2Ha)/(Hii_S2Ha - DIG_S2Ha)
    return C_Hii

##############################################################################
Z_sun = 8.69
O3N2_cal_data      = np.loadtxt('C:/Users/qihan/Desktop/q/Curti17_O3N2.txt')
O3S2_cal_data_2017 = np.loadtxt('C:/Users/qihan/Desktop/q/Curti17_O3S2.txt')
RS32_cal_data_2020 = np.loadtxt('C:/Users/qihan/Desktop/q/Curti20_RS32.txt')
##############################################################################
def compute_Z_N2S2Ha_Dop16_MUSE(line_df):
    '''
    Given a set of deredenned emission line maps+error, 
    compute metallicity maps+error, using the
    N2S2Ha diagnostic of Dopita+2016:
    
    
    Parameters
    ----------
    line_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
    https://ui.adsabs.harvard.edu/abs/2016Ap&SS.361...61D
    Returns
    -------
    Z: array
        Metallicity using this diagnostic
        
    e_Z: array
        Error in metallicity using this diagnostic, 
        computed via linear error propagation.
    
    # Unpack the wanted lines: f_NII, f_SII6717, f_SII6731, f_Ha,
    # and their errors.
    #line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_NII = line_df[line_IDs.index('NII6583_FLUX')].data
    e_NII = line_df[line_IDs.index('NII6583_FLUX_ERR')].data
    f_SII6716 = line_df[line_IDs.index('SII6716_FLUX')].data
    e_SII6716 = line_df[line_IDs.index('SII6716_FLUX_ERR')].data
    f_SII6731 = line_df[line_IDs.index('SII6730_FLUX')].data
    e_SII6731 = line_df[line_IDs.index('SII6730_FLUX_ERR')].data
    f_Ha = line_df[line_IDs.index('HA6562_FLUX')].data
    e_Ha = line_df[line_IDs.index('HA6562_FLUX_ERR')].data
    '''
    f_NII = line_df[36].data
    e_NII = line_df[37].data
    f_SII6716 = line_df[42].data
    e_SII6716 = line_df[43].data
    f_SII6731 = line_df[48].data
    e_SII6731 = line_df[49].data
    f_Ha = line_df[30].data
    e_Ha = line_df[31].data
    # compute relevant line ratios
    N2S2  = np.log10(f_NII/(f_SII6716+f_SII6731) )
    N2Ha  = np.log10(f_NII/f_Ha)
    N2S2Ha = N2S2 + 0.264*N2Ha
    Z_N2S2Ha_low = 8.77 + N2S2Ha
    Z_N2S2Ha_upper_correction = 0.45 * (N2S2Ha+0.3)**5
    Z_N2S2Ha = Z_N2S2Ha_low + Z_N2S2Ha_upper_correction*(Z_N2S2Ha_low > 9.05)
    # and errors
    dZ_dN2S2Ha = 1 + 2.25*(N2S2Ha+0.3)**4 *(Z_N2S2Ha_low > 9.05)
    dratio_dN2 = 1.264/(np.log(10)*f_NII)
    dratio_dS2 = 1.0/(np.log(10)*(f_SII6716+f_SII6731))
    dratio_dHa = 0.264/(np.log(10)*f_Ha)
    e_Z2 = dZ_dN2S2Ha**2 * ((dratio_dN2*e_NII)**2 + (dratio_dS2*e_SII6716)**2 \
                            + (dratio_dS2*e_SII6731)**2 + (dratio_dHa*e_Ha)**2)
    e_Z  = np.sqrt(e_Z2)
    return Z_N2S2Ha, e_Z
###############################################################################
def compute_Z_O3N2_Curti17_MUSE(line_df, kumari_correction = False):
    '''
    Given a set of deredenned emission line maps+error, 
    compute metallicity maps+error, using the
    N2S2Ha diagnostic of Curti+2017:
    https://ui.adsabs.harvard.edu/abs/2017MNRAS.465.1384C/abstract
    
    Optionally, add a DIG correction as devised by Kumari+19:
    https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..367K
    
    Parameters
    ----------
    line_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
        
    kumari_correction: bool or str
        if False; no correction
        if 'N2': add the published correction for DIG as selected using an
            N2-BPT diagram
        if 'S2': add the published correction for DIG as selected using an
            S2-BPT diagram
    
    Returns
    -------
    Z: array
        Metallicity using this diagnostic
        
    e_Z: array
        Error in metallicity using this diagnostic, 
        computed via linear error propagation.
    
    # Unpack the wanted lines: f_NII, f_OIII, f_Ha, f_Hb,
    # and their errors.
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_NII  = line_df[line_IDs.index('NII6583_FLUX')].data
    e_NII  = line_df[line_IDs.index('NII6583_FLUX_ERR')].data
    f_OIII = line_df[line_IDs.index('OIII5006_FLUX')].data
    e_OIII = line_df[line_IDs.index('OIII5006_FLUX_ERR')].data
    f_Ha   = line_df[line_IDs.index('HA6562_FLUX')].data
    e_Ha   = line_df[line_IDs.index('HA6562_FLUX_ERR')].data
    f_Hb   = line_df[line_IDs.index('HB4861_FLUX')].data
    e_Hb   = line_df[line_IDs.index('HB4861_FLUX_ERR')].data
    '''
    f_NII  = line_df[36].data
    e_NII  = line_df[37].data
    f_OIII = line_df[18].data
    e_OIII = line_df[19].data
    f_Ha   = line_df[30].data
    e_Ha   = line_df[31].data
    f_Hb   = line_df[6].data
    e_Hb   = line_df[7].data
    # define line ratios
    O3     = np.log10(f_OIII/f_Hb) 
    N2     = np.log10(f_NII/f_Ha)
    O3N2   = O3 - N2 
    # compute Z_O3N2 by interpolating the inverse function
    Z_o3n2 = np.interp(O3N2,O3N2_cal_data[:,0],O3N2_cal_data[:,1],left=np.nan, right=np.inf)
    x_o3n2 = Z_o3n2 - Z_sun
    # compute error function
    dZ_dratio= 1.0/(4.765+4.536*x_o3n2)
    if kumari_correction==False:
        e_Z2 = (dZ_dratio/np.log(10))**2 * ((e_NII/f_NII)**2 + (e_OIII/f_OIII)**2 + (e_Hb/f_Hb)**2+ (e_Ha/f_Ha)**2)
    if kumari_correction=='S2':
        # Then add a DIG correction based on the O3 line;
        # and our square error is a little different.
        Z_o3n2 = Z_o3n2 + 0.156*O3
        e_Z2 = (dZ_dratio/np.log(10))**2 * ((e_NII/f_NII)**2 + (1.156*e_OIII/f_OIII)**2 + (1.156*e_Hb/f_Hb)**2+ (e_Ha/f_Ha)**2)
    if kumari_correction=='N2':
        # same idea
        Z_o3n2 = Z_o3n2 + 0.033+0.127*O3
        e_Z2 = (dZ_dratio/np.log(10))**2 * ((e_NII/f_NII)**2 + (1.127*e_OIII/f_OIII)**2 + (1.127*e_Hb/f_Hb)**2+ (e_Ha/f_Ha)**2)
    e_Z = np.sqrt(e_Z2)
    return Z_o3n2, e_Z
################################################################################
def compute_Z_O3S2_Curti17_MUSE(line_df, kumari_correction=False):
    '''
    Given a set of deredenned emission line maps+error, 
    compute metallicity maps+error, using the
    N2S2Ha diagnostic of Curti+2017:
    https://ui.adsabs.harvard.edu/abs/2017MNRAS.465.1384C/abstract
    
    Optionally, add a DIG correction as devised by Kumari+19:
    https://ui.adsabs.harvard.edu/abs/2019MNRAS.485..367K
    
    Parameters
    ----------
    line_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
        
    kumari_correction: bool or str
        if False; no correction
        if 'N2': add the published correction for DIG as selected using an
            N2-BPT diagram
        if 'S2': add the published correction for DIG as selected using an
            S2-BPT diagram
    
    Returns
    -------
    Z: array
        Metallicity using this diagnostic
        
    e_Z: array
        Error in metallicity using this diagnostic, 
        computed via linear error propagation.
    '''
    # Unpack the wanted lines: f_SII6717, f_SII6731, f_OIII, f_Ha, f_Hb,
    # and their errors.
    f_SII6716 = line_df[42].data
    e_SII6716 = line_df[43].data
    f_SII6731 = line_df[48].data
    e_SII6731 = line_df[49].data   
    f_OIII = line_df[18].data
    e_OIII = line_df[19].data
    f_Ha   = line_df[30].data
    e_Ha   = line_df[31].data
    f_Hb   = line_df[6].data
    e_Hb   = line_df[7].data
    # define handy dandy ratios
    O3     = np.log10(f_OIII/f_Hb) 
    combo_flux = f_OIII/f_Hb + (f_SII6716+f_SII6731)/f_Ha
    combo_S2 = f_SII6716+f_SII6731
    O3S2   = np.log10(combo_flux)
    Z_o3s2 = np.interp(O3S2,O3S2_cal_data_2017[:,0],O3S2_cal_data_2017[:,1],left=np.nan, right=np.inf)
    x_o3s2 = Z_o3s2 - Z_sun
    # compute error function
    dZ_dratio= 1.0/(-2.223 - 2.146*x_o3s2 + 1.599* x_o3s2**2 )
    if kumari_correction==False:
        e_Z2 = (dZ_dratio/(np.log(10)*combo_flux))**2 * ( (e_OIII/f_Hb)**2 + (f_OIII*e_Hb/(f_Hb**2))**2+(combo_S2*e_Ha/(f_Ha**2))**2 +(e_SII6716/f_Hb)**2+ (e_SII6731/f_Hb)**2)
    if kumari_correction=='S2':
        # Then add a DIG correction based on the O3 line;
        # and the associated error.
        Z_o3s2 = Z_o3s2 + 0.075 + 0.309*O3 + 0.208*(O3**2)
        # first account for Ha,S2 contributions:
        e_Z2 = (dZ_dratio/(np.log(10)*combo_flux))**2 * ( (combo_S2*e_Ha/(f_Ha**2))**2 +(e_SII6716/f_Hb)**2+ (e_SII6731/f_Hb)**2)
        dZ_dO3 = 0.309+0.416*O3
        e_Z2_from_OIII = (dZ_dratio/(np.log(10)*combo_flux*f_Hb) +dZ_dO3/(np.log(10)*f_OIII)  )**2 * e_OIII**2
        e_Z2_from_Hb = (dZ_dratio*f_OIII/(np.log(10)*combo_flux*f_Hb*f_Hb) +dZ_dO3/(np.log(10)*f_Hb)  )**2 * e_Hb**2
        e_Z2 = e_Z2 + e_Z2_from_OIII + e_Z2_from_Hb
    if kumari_correction=='N2':
        # Then add a DIG correction based on the O3 line;
        # and the associated error
        Z_o3s2 = Z_o3s2 + 0.113 + 0.229*O3
        e_Z2 = (dZ_dratio/(np.log(10)*combo_flux))**2 * ( (combo_S2*e_Ha/(f_Ha**2))**2 +(e_SII6716/f_Hb)**2+ (e_SII6731/f_Hb)**2)
        dZ_dO3 = 0.229
        e_Z2_from_OIII = (dZ_dratio/(np.log(10)*combo_flux*f_Hb) +dZ_dO3/(np.log(10)*f_OIII)  )**2 * e_OIII**2
        e_Z2_from_Hb = (dZ_dratio*f_OIII/(np.log(10)*combo_flux*f_Hb*f_Hb) +dZ_dO3/(np.log(10)*f_Hb)  )**2 * e_Hb**2
        e_Z2 = e_Z2 + e_Z2_from_OIII + e_Z2_from_Hb
    e_Z = np.sqrt(e_Z2)
    return Z_o3s2, e_Z

################################################################################
def compute_Z_RS32_Curti20_MUSE(line_df):
    '''
    Given a set of deredenned emission line maps+error, 
    compute metallicity maps+error, using the
    N2S2Ha diagnostic of Curti+20:
    https://ui.adsabs.harvard.edu/abs/2020MNRAS.491..944C/abstract
    
    Parameters
    ----------
    line_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
    
    Returns
    -------
    Z: array
        Metallicity using this diagnostic
        
    e_Z: array
        Error in metallicity using this diagnostic, 
        computed via linear error propagation.
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    f_SII6716 = line_df[line_IDs.index('SII6716_FLUX')].data
    e_SII6716 = line_df[line_IDs.index('SII6716_FLUX_ERR')].data
    f_SII6731 = line_df[line_IDs.index('SII6730_FLUX')].data
    e_SII6731 = line_df[line_IDs.index('SII6730_FLUX_ERR')].data
    f_OIII = line_df[line_IDs.index('OIII5006_FLUX')].data
    e_OIII = line_df[line_IDs.index('OIII5006_FLUX_ERR')].data
    f_Ha   = line_df[line_IDs.index('HA6562_FLUX')].data
    e_Ha   = line_df[line_IDs.index('HA6562_FLUX_ERR')].data
    f_Hb   = line_df[line_IDs.index('HB4861_FLUX')].data
    e_Hb   = line_df[line_IDs.index('HB4861_FLUX_ERR')].data
    '''
    f_SII6716 = line_df[42].data
    e_SII6716 = line_df[43].data
    f_SII6731 = line_df[48].data
    e_SII6731 = line_df[49].data
    f_OIII = line_df[18].data
    e_OIII = line_df[19].data
    f_Ha   = line_df[30].data
    e_Ha   = line_df[31].data
    f_Hb   = line_df[6].data
    e_Hb   = line_df[7].data
    # define handy dandy ratios
    O3     = np.log10(f_OIII/f_Hb) 
    combo_flux = f_OIII/f_Hb + (f_SII6716+f_SII6731)/f_Ha
    combo_S2 = f_SII6716+f_SII6731
    O3S2   = np.log10(combo_flux)
    # Interpolate metallicities
    Z_o3s2 = np.interp(O3S2,RS32_cal_data_2020[:,0],RS32_cal_data_2020[:,1],left=np.nan, right=np.inf)
    x_o3s2 = Z_o3s2 - Z_sun
    # compute error function
    # 1/Derivative (d/dx) of:  -0.054 -2.546*x -1.970*(x**2) + 0.082* (x**3) + 0.222*(x**4)
    dZ_dratio= 1.0/(-2.546 - 3.940*x_o3s2 + 0.246*(x_o3s2**2) + 0.888*(x_o3s2**3) )
    # Compute error using linear error propagation
    e_Z2 = (dZ_dratio/(np.log(10)*combo_flux))**2 * ( (e_OIII/f_Hb)**2 + (f_OIII*e_Hb/(f_Hb**2))**2+(combo_S2*e_Ha/(f_Ha**2))**2 +(e_SII6716/f_Hb)**2+ (e_SII6731/f_Hb)**2)
    e_Z = np.sqrt(e_Z2)
    return Z_o3s2, e_Z
###############################################################################
def calculate_SFR_density_MUSE(Ha_map, D, inc):
    # Finagle out the deprojected area (units: kpc)
	world = WCS(Ha_map.header)
	pix_solid_angle = world.proj_plane_pixel_area().to(u.steradian).value
	plane_area	= pix_solid_angle * (D*1000)**2 # in kpc^2
	i  = np.radians(inc)
	deproj_area = plane_area / np.cos(i)
	# convert flux (units:1e-20 erg/s/cm2) to Luminosity (units: erg/s) L=4pi r^2 flux
	log_Ha_luminosity_map = np.log10(Ha_map.data) + np.log10(4*np.pi*(D*3.086*10**24)**2) - 20
	# convert to surface brightnesss (units: erg/s/kpc^2)
	log_Ha_SB_map = log_Ha_luminosity_map - np.log10(deproj_area)
	# Constant from Kennicutt+Evans12
	log_SFR_density_map   = log_Ha_SB_map - 41.27
	return 10**(log_SFR_density_map)
###############################################################################
def calculate_deproj_area_MUSE(Ha_map, inc_angle, Distance):
	# Finagle out the deprojected area (units: kpc)
	world = WCS(Ha_map.header)
	pix_solid_angle = world.proj_plane_pixel_area().to(u.steradian).value
	plane_area	= pix_solid_angle * (Distance*1000)**2 # in kpc^2
	i  = np.radians(inc_angle)
	deproj_area = plane_area / np.cos(i)
	return deproj_area
################################################################################
def calculate_total_SFR_of_galaxy_MUSE(SFR_density,deproj_area):
    # Find positions of NaNs and remove all of them to be 0.
    positions_are_NaNs_in_SFR_density = isnan(SFR_density)
    SFR_density[positions_are_NaNs_in_SFR_density] = 0
    # sum over the local SFR to get the Total SFR
    Total_SFR_Of_Galaxy = np.sum(SFR_density * deproj_area)
    return Total_SFR_Of_Galaxy
###############################################################################
def cal_Inclinations_from_axis_ratio(log_a_b_from_LEDA): # logr25 in LEDA website
    a_b_from_LEDA = 10**log_a_b_from_LEDA
    b_a = 1/a_b_from_LEDA
    Inclinations = np.arccos(b_a)*180/np.pi
    return Inclinations

def cal_Inclinations_from_axis_ratio_MUSE(log_a_b_from_LEDA, qz = 0.2): # logr25 in LEDA website, qz=0.2
    a_b_from_LEDA = 10**log_a_b_from_LEDA
    b_a = 1/a_b_from_LEDA
    Inclinations = np.arccos(np.sqrt(((b_a**2 - qz**2))/(1-qz**2)))*180/np.pi
    return Inclinations

def RA_DEC_to_radius_MUSE(Dist, PA_gal, inc_gal, RA, DEC, RA_galaxy, DEC_galaxy):
    return deprojected_distances_MUSE(Dist, PA_gal, inc_gal, RA, DEC, RA2 = RA_galaxy, DEC2 = DEC_galaxy).T[0]

def deprojected_distances_MUSE(Dist, PA_gal, inc_gal, RA1, DEC1, RA2 = None, DEC2 = None):
    '''
    Computes the deprojected distances between one set of RAs/DECs and
    another, for a known galaxy.
    
    Parameters
    ----------
    
    RA1: float, list, or np array-like
        List of (first) RA values. Must be in degrees.
        
    DEC1: float, list, or np array-like
        List of (first) DEC values. Must be in degrees.
        
    RA2: float, list, or np array-like
        (Optional) second list of RA values. Must be in degrees.
        If no argument is provided, then the first list will be used again.
        
    DEC2: float, list, or np array-like
        (Optional) second list of DEC values. Must be in degrees.
        If no argument is provided, then the first list will be used again.    
    
    meta: dict
        Metadata used to calculate the distances. Must contain:
        PA: float
            Principle Angle of the galaxy, degrees.
        i: float
            inclination of the galaxy along this principle axis, degrees.
        D: float
            Distance from this galaxy to Earth, Mpc.
        
    Returns
    -------
    dists: np array
        Array of distances between all RA, DEC pairs provided.
        Units: kpc.
    
    '''
    try:
        assert len(RA1) == len(DEC1), "Error: len of RA1 must match len of DEC1"
        RA1 = np.array(RA1)
        DEC1 = np.array(DEC1)
    except TypeError:
        assert type(RA1) == type(DEC1), "Error: type of RA1 must match type of DEC1"  
        # Then cast them to arrays
        RA1 = np.array([RA1])
        DEC1 = np.array([DEC1])
        
    if type(RA2) == type(None):
        RA2 = RA1
    if type(DEC2) == type(None):
        DEC2 = DEC1
    
    try:
        assert len(RA2) == len(DEC2), "Error: len of RA2 must match len of DEC2"
        RA2 = np.array(RA2)
        DEC2 = np.array(DEC2)
    except TypeError:
        assert type(RA2) == type(DEC2), "Error: type of RA2 must match type of DEC2" 
        RA2 = np.array([RA2])
        DEC2 = np.array([DEC2])
    
    PA = np.radians(PA_gal)
    i = np.radians(inc_gal)
    # 1: Rotate RA, DEC by PA to get y (major axis direction) and x (minor axis direction)
    x1 = RA1*np.cos(PA) - DEC1*np.sin(PA)
    y1 = DEC1*np.cos(PA) + RA1*np.sin(PA)
    x2 = RA2*np.cos(PA) - DEC2*np.sin(PA)
    y2 = DEC2*np.cos(PA) + RA2*np.sin(PA)
    # 2: Stretch x values to remove inclination effects
    long_x1 = x1 /np.cos(i)
    long_x2 = x2 /np.cos(i)
    # 3: Compute Euclidean Distances between x1,y1 and x2,y2 to get angular offsets (degrees).
    vec1 = np.stack((y1.flatten(), long_x1.flatten())).T
    vec2 = np.stack((y2, long_x2)).T
    deg_dists = euclidean_distances(vec1, np.array(np.matrix(vec2)))
    rad_dists = np.radians(deg_dists)
    # 4: Convert angular offsets to kpc distances using D, and the small-angle approximation.
    # Mpc_dists = rad_dists * meta['D']
    Mpc_dists = rad_dists * Dist
    kpc_dists = Mpc_dists * 1000
    
    return kpc_dists

def RA_DEC_to_xy_MUSE(RA, DEC, RA_gal, DEC_gal, PA, inc, D):
    '''
    Parameters
    ----------
    RA: np array-like
    List of RA values. Must be in degrees.
    DEC: np array-like
    List of DEC values. Must be in degrees.  
    meta: dict
    Metadata used to calculate the distances. Must contain:
    PA: float
    Principle Angle of the galaxy, degrees.
    i: float
    inclination of the galaxy along this principle axis, degrees.
    D: float
    Distance from this galaxy to Earth, Mpc.
    Returns
    -------
    x: np array
    Deprojected distances along the direction of the minor axis (kpc)
    y: np array
    Deprojected distances along the direction of the major axis (kpc)
    '''
    RA = np.array(RA) - RA_gal
    DEC = np.array(DEC) - DEC_gal
    # Now onto the maths
    PA = np.radians(PA)
    i  = np.radians(inc)
    # 1: Rotate RA, DEC by PA to get y (major axis direction) and x (minor axis direction)
    x = RA*np.cos(PA) - DEC*np.sin(PA)
    y = DEC*np.cos(PA) + RA*np.sin(PA)
    # 2: Stretch x values to remove inclination effects
    x = x /np.cos(i)
    # 3: Convert deg to kpc
    x = np.radians(x)*D*1000
    y = np.radians(y)*D*1000
    return x, y
##############################################################################################
def MUSE_plot_galaxy(gal_name, file_path, DIG_CUT="S2_DIG", Z_train = "Z_N2S2Ha", 
                     galaxy_info_path = "C:/Users/qihan/Desktop/q/galaxydata.xlsx"):
    '''
    DIG_CUT = "N2_BPT"
    DIG_CUT = "S2_BPT"
    DIG_CUT = "S2_DIG"
    #
    Z_train = "Z_N2S2Ha"
    Z_train = "Z_O3N2"
    Z_train = "Z_RS32"
    #
    eg. MUSE_plot_galaxy(gal_name = 'N4321', 
                         file_path = 'C:/Users/qihan/Desktop/q/MUSE_pkl/N4321_ver1_copt.pkl', 
                         DIG_CUT="S2_DIG", Z_train = "Z_N2S2Ha", galaxy_info_path = 
                         "C:/Users/qihan/Desktop/q/galaxydata.xlsx")
    '''
    trainset = pd.read_pickle(file_path) 
    ### Get information of galaxy. ###
    galaxy_info = pd.read_excel(galaxy_info_path)
    idx = galaxy_info["Gal_ID"] == gal_name
    galaxy_info_now = galaxy_info[idx]
    inc_angle = float(galaxy_info_now['i_leda'])
    Distance = float(galaxy_info_now['D'])
    RA_galaxy = float(galaxy_info_now['RA'])
    DEC_galaxy = float(galaxy_info_now['DEC'])
    PA_galaxy = float(galaxy_info_now['PA_MUSE'])
    ### 
    if DIG_CUT == "N2_BPT":
        index_below_1_BPT = trainset["N2_BPT"] < 1 
        trainset = trainset[index_below_1_BPT]   
    if DIG_CUT == "S2_BPT":
        idx_S2_BPT1 = trainset["S2_BPT"] < 1
        trainset = trainset[idx_S2_BPT1]
    if DIG_CUT == "S2_DIG":
        index_above_CHii = trainset["S2_DIG"] > 0.9 
        trainset = trainset[index_above_CHii]
    index_above_0_Z1 = trainset[Z_train] > 0  
    trainset = trainset[index_above_0_Z1]
    RA1 = trainset['RA']
    DEC1 = trainset['DEC']
    X1 = RA_DEC_to_xy_MUSE(RA1, DEC1, RA_galaxy, DEC_galaxy, PA_galaxy, inc_angle, Distance)
    X1 = np.transpose(X1)
    plt.scatter(X1[:, 0], X1[:, 1], c = trainset[Z_train], cmap='hot_r', marker='.', s=1)
    plt.title(f'{gal_name} {Z_train} {DIG_CUT} HII regions') 
    plt.xlabel('x (kpc)', color = "black")
    plt.ylabel('y (kpc)', color = "black")
    return None

##############################################################################################
def MLE_fit_MUSE(y, X, D, cov_model, eta_ini, nug, opt, lo_bound, up_bound):
    #bound = cov_bound(cov_model = cov_model, nug = nug)
    q = eta_ini.size + 1
    bound = []
    for i in range(q-1):
        bound.append((lo_bound[i], up_bound[i]))
    
    ## Optimization algorithms, L-BFGS-B
    if opt == "LB": 
        soln = minimize(profile_nll, eta_ini, bounds=bound, args=(y, X, D, cov_model, nug), method="L-BFGS-B")
    # Nelder-Mead
    elif opt == "NM":
        soln = minimize(profile_nll, eta_ini, bounds=bound, args=(y, X, D, cov_model, nug), method="Nelder-Mead")
    elif opt == "DE":
        soln = differential_evolution(profile_nll, bound, args = (y, X, D, cov_model, nug), seed = 2024, tol = 0.0001)
    else:
        soln = minimize(profile_nll, eta_ini, bounds=bound, args=(y, X, D, cov_model, nug), method = opt)
    eta = soln.x
    suc = soln.success
    nll = soln.fun
    
    n = y.size
    cormat = cor_mat(D = D, eta = eta, cov_model = cov_model, nug = nug)
    L = np.linalg.cholesky(cormat)
    white_X = np.linalg.solve(L, X)
    white_y = np.linalg.solve(L, y)
    
    beta_est = linalg.solve(white_X.T @ white_X, white_X.T @ white_y, assume_a='sym')
    white_resids = white_y - white_X @ beta_est
    sill = white_resids.T @ white_resids/n
    beta_var = linalg.inv(white_X.T @ white_X)*sill

    ## reparametrization
    theta_est = np.copy(eta)
    if nug == True: 
        #theta_est [ntheta-1] = theta[ntheta-1]*sill    # nugget effect
        nug_effect = sill*eta[q-2]
        psill = sill*(1-eta[q-2])
        theta_est = np.append(theta_est[0:(q-2)], [psill, nug_effect])
    elif nug == False:
        psill = sill
        theta_est = np.append(theta_est,  psill)
    return theta_est, suc, nll, beta_est, beta_var
#####################################################################################
def MUSE_compute_geo_local_models(file_path, gal_name, test_set_percentage = 0.2, 
                             n_split=10, seed = 2024, k=200, DIG_CUT = "S2_DIG", 
                             Z_train = "Z_N2S2Ha", X_var = ["proj_dist"], 
                             galaxy_info_path = "C:/Users/qihan/Desktop/q/galaxydata.xlsx"):
    '''
    MLE used EXP model and include the beta0 in all calculation.
    
    # Choose different DIG cut
    DIG_CUT = "N2_BPT"
    DIG_CUT = "S2_BPT"
    DIG_CUT = "S2_DIG"
    
    # Choose different metallicity diagnostic
    Z_train = "Z_N2S2Ha"
    Z_train = "Z_O3N2"
    Z_train = "Z_RS32"
    
    X_var: the variables be used in the MLE local models. please see the galaxy data file to
           see what variables avaible in the dataset. By default for MUSE, we can choose 
           X_var = ["proj_dist"] or X_var = ["proj_dist", "SFR_density", "S2_DIG"]. 
    
    k: the number of nearby points for local model, default = 200. k should keep 
       small for efficiency. k=100, 150, 200, 250, 300 should be good enough.
       
    seed: The seed number for n-fold CV. 
    
    n_split: The number of n in n-fold CV, noramlly, choose n = 10. However, it will make the 
             program run very slow, hence, if we use the code for test purpose, always choose n=1.
    
    test_set_percentage: the percentage of the test set from full dataset. normally, choose
                         20% of test data and 80% of train data. Or, choose 30% of test data.
                         By default, we use 0.2 (20%).
                    
    return:
        RMSE_list: the RMSE list of n-fold CV
        
        MAD_list: the MAD list of n-fold CV
        
        Pred_lib: the prediction values of local models
        
        True_lib: the true values of the test set
        
    example:
        RMSE_list, MAD_list, Pred_lib, True_lib = MUSE_compute_geo_local_models(
                                 file_path = 'C:/Users/qihan/Desktop/q/MUSE_pkl/N4321_ver1_copt.pkl', 
                                 gal_name = 'N4321', test_set_percentage = 0.2, n_split=1, seed = 2024, 
                                 k=50, DIG_CUT = "S2_DIG", Z_train = "Z_N2S2Ha", X_var = ["proj_dist"], 
                                 galaxy_info_path = "C:/Users/qihan/Desktop/q/galaxydata.xlsx")
    '''
    fullset = pd.read_pickle(file_path) 
    ### Get information of galaxy. ###
    galaxy_info = pd.read_excel(galaxy_info_path)
    idx = galaxy_info["Gal_ID"] == gal_name
    galaxy_info_now = galaxy_info[idx]
    inc_angle = float(galaxy_info_now['i_leda'])
    Distance = float(galaxy_info_now['D'])
    RA_galaxy = float(galaxy_info_now['RA'])
    DEC_galaxy = float(galaxy_info_now['DEC'])
    PA_galaxy = float(galaxy_info_now['PA_MUSE'])

    # remove DIG regions and NANs:
    if DIG_CUT == "N2_BPT":
        index_below_1_BPT = fullset["N2_BPT"] < 1 
        fullset = fullset[index_below_1_BPT]   
    if DIG_CUT == "S2_BPT":
        idx_S2_BPT1 = fullset["S2_BPT"] < 1
        fullset = fullset[idx_S2_BPT1]
    if DIG_CUT == "S2_DIG":
        index_above_CHii = fullset["S2_DIG"] > 0.9 
        fullset = fullset[index_above_CHii]
    index_above_0_Z1 = fullset[Z_train] > 0  
    fullset = fullset[index_above_0_Z1]

    # by default, we use logSFR instead of SFR:
    fullset[["SFR_density"]] = np.log(fullset[["SFR_density"]])
    # initializtion:
    RMSE_list = list()
    MAD_list = list()
    MAD_lib = list()
    Pred_lib = defaultdict(list)
    True_lib = defaultdict(list)

    rs = ShuffleSplit(n_splits=1, test_size=test_set_percentage, random_state=seed)
    for i, (train_index, test_index) in enumerate(rs.split(fullset)):
        print(f"Fold {i}:")
        trainset = fullset.iloc[train_index]
        testset = fullset.iloc[test_index]
        RAtrain = trainset['RA']
        DECtrain = trainset['DEC']
        RAtest = testset['RA']
        DECtest = testset['DEC']
        Y_true = testset[Z_train]

        pred_test_set = list()
        for ii in range(testset.shape[0]):
            print(ii)
            testpoint = pd.DataFrame(testset).iloc[ii]
            RApoint = pd.DataFrame(RAtest).iloc[ii][0]
            DECpoint = pd.DataFrame(DECtest).iloc[ii][0]
            D_testpoint_to_train = deprojected_distances_MUSE(Distance, PA_galaxy, inc_angle, RApoint, DECpoint, RA2 = RAtrain, DEC2 = DECtrain)        
            idx_first_k_smallest = sort(np.argpartition(D_testpoint_to_train,k)[0][0:k])
            sur_train = pd.DataFrame(trainset).iloc[idx_first_k_smallest] 
            RAsur_train = sur_train['RA']
            DECsur_train = sur_train['DEC']          
            Dtrain = deprojected_distances_MUSE(Distance, PA_galaxy, inc_angle, RAsur_train, DECsur_train)
            Dtest = deprojected_distances_MUSE(Distance, PA_galaxy, inc_angle, RApoint, DECpoint, RA2 = RAsur_train, DEC2 = DECsur_train)
            coords = RA_DEC_to_xy_MUSE(RAsur_train, DECsur_train, RA_galaxy, DEC_galaxy, PA_galaxy, inc_angle, Distance)
            # get the initial values, however, it should be not important for most case since the MLE will convergent.
            Vg = skg.Variogram(np.transpose(np.array(coords)), sur_train[Z_train], model = "exponential")
            coef = Vg.parameters
            sigma2_ini = coef[1]
            phi_ini = coef[0]  
            # 
            X_train = sur_train[X_var]
            X_pred = testpoint[X_var]
            Y_train = sur_train[Z_train]     
            X_train.insert(0, "intersect", np.ones(X_train.shape[0]), True)
            eta_ini_value = np.array([phi_ini/np.sqrt(3), sigma2_ini])
            lower_bound = np.array([1e-5, 1e-10])
            upper_bound = np.array([1,1])  
            MLE_result = MLE_fit_MUSE(y = Y_train, X = X_train, D = Dtrain, cov_model = "Exp", eta_ini = eta_ini_value, nug = True, opt = "LB", lo_bound = lower_bound, up_bound = upper_bound)  
            thetahat = MLE_result[0] 
            betahat = MLE_result[3]  
            Sigma = thetahat[1]*np.exp(-Dtrain/thetahat[0])    
            Sigma[np.diag_indices_from(Sigma)] = thetahat[1] + thetahat[2]  
            X_pred = insert(matrix(X_pred), 0,1, True)   
            cmat = thetahat[1]*np.exp(-Dtest/thetahat[0])
            pl = X_pred @ betahat
            pe = cmat @ linalg.inv(Sigma) @ (Y_train - X_train @ betahat)
            Y_pred = pl + pe
            pred_test_set.append(Y_pred)
            
        new_set = list()
        for j in range(testset.shape[0]):
            new_set.append(float(pd.DataFrame(pred_test_set[j]).values))
        MSE = mean_squared_error(new_set,Y_true.values) 
        RMSE = np.sqrt(MSE)
        MAD = statistics.mean(abs(new_set- Y_true.values))
        RMSE_list.append(RMSE)
        MAD_list.append(MAD)
        Pred_lib[i].append(new_set)
        True_lib[i].append(Y_true.values)
        return RMSE_list, MAD_list, Pred_lib, True_lib
    
##################################################################################
def plot_corr_matrix(maxdist, file_path, DIG_cut = "S2_DIG", 
                     var = ["proj_dist", "SFR_density", "S2_DIG"], save_eps = False):
    '''
    This code use to create correlation matrix with variables of the trainning set.
    Note that we should always draw corr matrix for train set, do not use test set
    since we should not "see" the test set before we compare our prediction results
    to the true values of the test set.
    '''
    trainset = pd.read_pickle(file_path) 
    index_above_0_DIG = trainset[DIG_cut] > 0
    trainset = trainset[index_above_0_DIG]
    if "proj_dist" in var:
       index_above_0_dist = trainset["proj_dist"] > 0
       trainset = trainset[index_above_0_dist]
    if "SFR_density" in var:
       index_above_0_SFR = trainset["SFR_density"] > 0
       trainset = trainset[index_above_0_SFR]
    df = trainset[var]
    f = plt.figure(figsize=(19, 15))
    plt.matshow(df.corr(), fignum=f.number, cmap='Blues')
    plt.xticks(range(df.select_dtypes(['number']).shape[1]), df.select_dtypes(['number']).columns, fontsize=20)
    plt.yticks(range(df.select_dtypes(['number']).shape[1]), df.select_dtypes(['number']).columns, fontsize=20)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=14)
    plt.title('Correlation Matrix for some of Variables', fontsize=16)
    if save_eps == True:
        plt.savefig('coor_matrix.eps', dpi=300)
    return None
    
    
def url_downloader(target_url, output_path, download = True):
    '''
    This function used to download url fits file from website. However, this method
    is not important for now since we can download the fits file from PHANGES website
    directly now. 
    eg: 
    target_url = 'https://ws-cadc.canfar.net/......'
    output_path = 'N5236_annulus_0p5kpc.fits'
    '''
    data_url = fits.open(target_url)
    data_url.writeto(output_path)  

# 1 arcsecond = 4.84813681 × 10-6 radians
def arcsecond_to_rad(arcsecond):
    rad = arcsecond*(4.84813681*10**(-6))
    return rad

def arcseond_to_pc_resolusion(arcsecond, DMpc):
    # DMpc: Distance from Earth in Mpc， 4.89778819
    # arcsecond: the angular resolusion in as, we need to change to radian
    # TYPHOON 1.65 as = 39.12 pc
    # ALMA 2 as       = 47.41 pc
    rad = arcsecond_to_rad(arcsecond)
    pc = DMpc*rad * 1000000
    return pc

def asec_fast(gal_name, galaxy_info_path = "C:/Users/qihan/Desktop/q/galaxydata.xlsx"):
    galaxy_info = pd.read_excel(galaxy_info_path)
    idx = galaxy_info["Gal_ID"] == gal_name
    galaxy_info_now = galaxy_info[idx]
    DMpc = float(galaxy_info_now['D'])
    arcsecond = float(galaxy_info_now['asec'])
    pc = arcseond_to_pc_resolusion(arcsecond, DMpc)
    return pc
    


    