'''
TYPHOON_wrangling.py

Set of helper functions designed to help process TYPHOON data. 

Created by: Benjamin Metha
Last updated: Dec 07, 2022 (trimmed to remove unneeded functionality)
'''

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

data_path = 'C:/Users/qihan/Desktop/q/'

# BPT classification codes
STARBURST = 0
SEYFERT   = 1
LINER     = 2

ASEC_PER_RAD = 206265.0

class InputError(Exception):
    pass

#############################################
#                                           #   
#             input/output                  #
#                                           #
#############################################

def check_diag(diagnostic):
    '''
    Tells you whether a diagnostic is valid; throws an error if not.
    '''
    diagnostic = diagnostic.lower()
    ok_diags=['n2o2', 'o3n2', 'n2s2ha', 'o3s2_old','rs32']
    if diagnostic not in ok_diags:
        raise ValueError("Error: "+ diagnostic +" is not a valid diagnostic.\n Choose from:"+str(ok_diags)[1:-1])
    return 0

def open_metallicity_map(diagnostic='o3n2', gal_name='N5236', error=False):
    '''
    Opens metallicity map for a given galaxy, using a given metallicity diagnostic
    
    Parameters
    ----------
    diagnostic: str
        Tells us which metallicity diagnostic to use.
        Options are ['n2ha', 'n2o2', 'n2s2', 'o3n2', 'r23'].
        Upper case values also accepted.
        Defaults to 'o3n2'
    
    gal_name: str
        Name of the galaxy. 
        Currently only N5236 is downloaded/has data, so it defaults to that.
        But it's a parameter you can change in the future.
        
    error: bool
        Do you want the error map for this diagnostic, too?

    Returns
    -------
    
    Z_map: np.array
       Map for the metallicity of the galaxy at each location
    
    e_Z_map: np.array
        if error==True, I'll also give you the error map for this diagnostic.
    '''
    diagnostic = diagnostic.lower()
    check_diag(diagnostic)
    Z_hdu = fits.open(data_path+'{0}/{0}_full_41pc_MetIon_{1}.fits'.format(gal_name, diagnostic))
    if error:
        e_Z_hdu = fits.open(data_path+'{0}/{0}_full_41pc_meterror_{1}.fits'.format(gal_name, diagnostic))
        return Z_hdu[0].data, e_Z_hdu[0].data
    # Otherwise just give me the metallicity map, no error
    return Z_hdu[0].data

def open_Hii_map(gal_name='N5236'):
    """
    Parameters
    ----------
    gal_name: str
        Name of the galaxy. 
        Currently only N5236 is downloaded/has data, so it defaults to that.
        But it's a parameter you can change in the future.
        
    Returns
    -------
    
    Hii_map: np.array
       Value of 1 if we're in a Hii region, 
       0 if we are DIG dominated,
       and NaN if there is no data.
    """
    Hii_hdu = fits.open(data_path+gal_name+'/Hii_region.fits')
    return Hii_hdu[0].data

def M83_metadata():
    '''Yes I should save this as a data object I can read. 
    No, I won't do that.
    '''
    M83_meta = {'D': 4.47,       #Mpc (Tully et al. 2008)
                'i': 32.5,       #Lauberts & Valentijn (1989)
               'PA': 44.9,       #Lauberts & Valentijn (1989)
               'RA':204.2539583, #Dıaz et al. (2006)
              'DEC':-29.8654167  #Dıaz et al. (2006)
                   }
    return M83_meta

def open_line_df(gal_ID):
    return fits.open(data_path + '{0}_lowres_cal_1_comp_WCS.fits'.format(gal_ID))

def open_Hii_df(gal_name='N5236', by='S2_Kaplan16'):
    return pd.read_pickle(data_path +'Handmade/{0}/Hii_by_{1}.pkl'.format(gal_name,by))

def open_DIG_df(gal_name='N5236', by='S2_Kaplan16'):
    return pd.read_pickle('../Data/Handmade/{0}/DIG_by_{1}.pkl'.format(gal_name, by))
    
def Re_kpc(meta):
    '''Takes in a galaxy's metadata. Returns a conversion factor that converts
    kpc to Re.'''
    Re_radians = meta['Re (arcmin)']*60/ASEC_PER_RAD
    Re_kpc = Re_radians*1000*meta['D']
    return Re_kpc
    
#############################################
#                                           #
#             data wrangling                #
#                                           #
#############################################

def make_RA_DEC_grid(header):
    '''
    Given a header file, create a grid of RA//DEC for each pixel in that file.
    '''
    world = WCS(header)
    x = np.arange(header['NAXIS1'])
    y = np.arange(header['NAXIS2'])
    X, Y = np.meshgrid(x, y)
    RA_grid, DEC_grid = world.wcs_pix2world(X, Y, 0)
    return RA_grid, DEC_grid        

def unpack_and_trim(Hii_df, diag, dtype='f4'):
    '''
    Trim nans. Return Z, e_Z, RA, and DEC associated with the non-nan values for
    a specified diagnostic.
    '''
    check_diag(diag)
    wanted_spaxels = ~np.isnan(Hii_df['Z_'+diag])
    Z   = np.array(Hii_df['Z_'+diag][wanted_spaxels],dtype=dtype)
    e_Z = np.array(Hii_df['e_Z_'+diag][wanted_spaxels],dtype=dtype)
    RA  = np.array(Hii_df['RA'][wanted_spaxels],dtype=dtype)
    DEC = np.array(Hii_df['DEC'][wanted_spaxels],dtype=dtype)
    return RA, DEC, Z, e_Z

def count_non_nans(A):
    '''Handy debugginator'''
    return np.sum(~np.isnan(A))

#############################################
#                                           #
#     processing emission line data into    #
#        metallicity and DIG/Hii maps       #
#                                           #
#############################################

def SN_cut(line_df, threshold=3):
    '''
    Replace all spaxels with SN<3 in a certain line with NANs.
    
    Parameters
    ----------
    
    lines_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
        
    threshold: float
        At what S/N do we cut a line? (Defaulted to 3)
        
    Returns
    -------
    lines_df: hdu list
        The same hdu list, but with lines where S/N < threshold
        replaced with np.nan
    '''
    n_lines = int(len(line_df)/2)
    x_max, y_max = line_df[0].data.shape
    for l in range(n_lines):
        signal = line_df[2*l].data
        noise  = line_df[2*l+1].data
        too_low = signal <= threshold*noise
        # Skip this step if the line is O2 - there's an error in the estimates
        # in this, so S/N cuts ought not be trusted. (Battisti, private 
        # communication, Oct 1 2021)
        # In this case, just exclude spaxels with S/N=0.
        if line_df[2*l].header['EXTNAME'] == 'OII3726' or line_df[2*l].header['EXTNAME'] == 'OII3729':
            too_low = (signal <= 0)
        for ii in range(x_max):
            for jj in range(y_max):
                # replace low signals/no signals with NANs.
                if too_low[ii,jj]:
                    signal[ii,jj] = np.nan
                    noise[ii,jj]  = np.nan
    
    return line_df


    
def SN_cut_O2(line_df, threshold=1):
    '''
    Separate program for the O2 line because it's special.
    '''
    n_lines = int(len(line_df)/2)
    x_max, y_max = line_df[0].data.shape
    for l in range(n_lines):
        signal = line_df[2*l].data
        noise  = line_df[2*l+1].data
        if line_df[2*l].header['EXTNAME'] == 'OII3726' or line_df[2*l].header['EXTNAME'] == 'OII3729':
            too_low = signal <= threshold*noise
        else:
            too_low = (signal <= 0)
        for ii in range(x_max):
            for jj in range(y_max):
                # replace low signals/no signals with NANs.
                if too_low[ii,jj]:
                    signal[ii,jj] = np.nan
                    noise[ii,jj]  = np.nan
    
    return line_df
    

def determine_DIG_S2_Kaplan16(line_df, n_spaxels=100, max_prop=0.05):
    '''
    Assuming that:
    1. The Sii/Ha line ratio is significantly different for DIG/Hii regions;
    2. The intrinsic distributions of Sii/Ha are (infinitely) narrow for purely Hii/DIG regions
    
    Compute the fraction of Ha-light originating from Hii regions for each spaxel
    (C_Hii), using the formalism of Kaplan+2016:
    https://ui.adsabs.harvard.edu/abs/2016MNRAS.462.1642K
    
    The formula:
    
    [Sii/Ha] = C_Hii [Sii/Ha]_Hii + C_DIG [Sii/Ha]_DIG
    
    Solvable when you use:
    
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
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    Ha = line_df[line_IDs.index('HALPHA')].data
    S2Ha = np.log10( (line_df[line_IDs.index('SII6716')].data+line_df[line_IDs.index('SII6731')].data)/line_df[line_IDs.index('HALPHA')].data    )
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

def determine_DIG_Ha_Zhang17(Ha_map, meta):
    '''
    Determine whether a spaxel is Hii/DIG dominated, by applying a surface
    brightness cut in Ha.
    
    
    Formula: dig if log_10( SB_Ha ) < 39
    
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
    
    DIG_map: np array
        0 if a spaxel is DIG dominated
        1 if it's a Hii region
        nan if SN too low to decide.
    '''
    # Finagle out the deprojected area (units: kpc)
    world = WCS(Ha_map.header)
    pix_solid_angle = world.proj_plane_pixel_area().to(u.steradian).value
    plane_area  = pix_solid_angle * (meta['D']*1000)**2 # in kpc^2
    i  = np.radians(meta['i'])
    deproj_area = plane_area / np.cos(i)
    # convert flux (units:1e-17 erg/s/cm2) to Luminosity (units: erg/s)
    log_Ha_luminosity_map = np.log10(Ha_map.data) + np.log10(4*np.pi*meta['D']**2 * 9.5234) + 31
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

def classify_S2_BPT(line_df):
    '''
    For each spaxel
    specify whether it is SEYFERT, LINER, or SF
    using the diagnostics of Kewley+01 and Kewley+06
    and the S2-BPT diagram.
    
    Parameters
    ----------
    
    lines_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
        
    Returns
    -------
    
    S2_BPT_classification: np array
        NAN if line data is missing/has too low S/N
        0 if starburst
        1 if Seyfert
        2 if LINER
    
    1s/2s will be treated as DIG
    '''
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    O3Hb = np.log10( line_df[line_IDs.index('OIII5007')].data /  line_df[line_IDs.index('HBETA')].data )
    S2Ha = np.log10( (line_df[line_IDs.index('SII6716')].data+line_df[line_IDs.index('SII6731')].data)/line_df[line_IDs.index('HALPHA')].data    )
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

def classify_N2_BPT(line_df, rule="Kauffmann03"):
    '''
    For each spaxel
    specify whether it is LINER or SF
    using the diagnostic of Kewley+01
    and the N2-BPT diagram.
    
    Parameters
    ----------
    
    lines_df: hdu list
        A big guy containing all the different emission line data reduced
        from TYPHOON data cubes
        
    Returns
    -------
    
    N2_BPT_classification: np array
        NAN if line data is missing/has too low S/N
        0 if starburst
        2 if LINER
    
    2s will be treated as DIG
    '''
    line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))]
    O3Hb = np.log10( line_df[line_IDs.index('OIII5007')].data/line_df[line_IDs.index('HBETA')].data )
    N2Ha = np.log10( line_df[line_IDs.index('NII6583')].data/line_df[line_IDs.index('HALPHA')].data    )
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



def calculate_SFR_density(Ha_map, meta):
	'''
    fits_file = fits.open('C:/Users/qihan/Desktop/q/N5236_lowres_cal_1_comp_WCS.fits')
    # meta=meta_getter('N5236')
    # Ha_map = fits_file[8]
    # Let sfrd = calculate_SFR_density(Ha_map, meta)
    # plt.figure()
    # plt.imshow(sfrd,origin = 'lower', norm = LogNorm(), cmap= 'Greys') # camp='Greys' means Grey scale, if not use it, will be color.
    
	Determine the star formation rate using Kennicutt+Evans98:
	https://ui.adsabs.harvard.edu/abs/2012ARA%26A..50..531K/abstract
	
	Formula: .... 
	
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
	# convert flux (units:1e-17 erg/s/cm2) to Luminosity (units: erg/s) L=4pi r^2 flux
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


'''
line_df_1 = fits.open('C:/Users/qihan/Desktop/q/N5236_lowres_cal_1_comp_WCS.fits')
wavelengths = np.array([3726.0, 3729.0, 4861.3, 5007.0, 6562.8, 6583.0, 6716.0, 6731.0])
line_df_2 = SN_cut(line_df_1, threshold=3)
line_df_3 = SN_cut_O2(line_df_2, threshold=1)
line_df_4= extinction_correction(line_df_3, wavelengths, R_V=3.1)

Ha_map = line_df_4[8]
meta = meta_getter('N5236')

SFR_density=calculate_SFR_density(Ha_map, meta)
SFR =calculate_SFR_MUSE(Ha_map, meta)
deproj_area=calculate_deproj_area(Ha_map, meta)

calculate_total_SFR_of_galaxy(SFR_density,deproj_area)
calculate_total_SFR_of_galaxy_MUSE(SFR)
# the results for N5236 is 4.8322296.
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


def extinction_correction(line_df, wavelengths, R_V=3.1):
	'''
	Parameters
	----------
	lines_df: hdu list, A big guy containing all the different emission line data reduced from TYPHOON data cubes
		
	wavelengths: np.array, Wavelength of each of the 8 lines in this data cube, in Angstroms.
        
	wavelengths = np.array([3726.0, 3729.0, 4861.3, 5007.0, 6562.8, 6583.0, 6716.0, 6731.0])
          
	R_V: float, The free parameter in ccm89 extinction law. Set (kept) at 3.1.
	
	Returns
	-------
	corrected_lines_df: hdu list, Corrections for all lines using the calibration of ccm89.
	'''
	line_IDs = [line_df[x].header['EXTNAME'] for x in range(len(line_df))] # the who's who of line data
	Ha_map = line_df[line_IDs.index('HALPHA')].data
	Hb_map = line_df[line_IDs.index('HBETA')].data
	# To convert balmer decrement to extinction, need these...
	HA_EXT =  ccm89(np.array([6562.8]), 1.0, R_V)[0]
	HB_EXT =  ccm89(np.array([4861.3]), 1.0, R_V)[0]
	Ha_Hb_ratio	 = Ha_map/Hb_map
	balmer_decrement = 2.5*np.log10(Ha_Hb_ratio / 2.86)
	A_V = balmer_decrement/(HB_EXT - HA_EXT) 
	A_V_positive = A_V * (A_V > 0) # sets negatives to zero
	
	# Use this to correct obs and error for each wavelength
	for l in range(len(wavelengths)):
		extinction_at_wav = ccm89(wavelengths[l:l+1], 1, R_V)[0]
		extinction_map = extinction_at_wav*A_V_positive
		# correct signal and noise
		line_df[2*l].data	 = line_df[2*l].data * 10**(0.4 * extinction_map)
		line_df[2*l+1].data	 = line_df[2*l+1].data * 10**(0.4 * extinction_map)
	return line_df



def open_metadata():
	'''
	Open the metadata file; saves from re-writing this code in every script,
	makes scripts more readable, and future-proofs code if I need to update
	or move the metadata file.
	'''
	meta_df = pd.read_csv(data_path+'metadata.csv')
	metadata = meta_df.to_dict(orient='records')
	return metadata
	
def meta_getter(gal_ID):
	metadata = open_metadata()
	meta   = [x for x in metadata if str(gal_ID) in x['Gal_ID']][0]
	return meta

    


