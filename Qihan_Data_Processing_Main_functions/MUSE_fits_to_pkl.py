# -*- coding: utf-8 -*-
"""
##############################################################################
Created on Fri Apr  5 21:23:43 2024

@author: Qihan Zou 

Last updated: 18/11/2024

Based on original codes from Benjamin Metha (Version: 26/10/2021) and his helps. 

##############################################################################
This file is mainly constructed for MUSE dataset.
1. Data cleaning and data processing
2. Create simple potential useful variables
3. Create pkl file from fits file

Main functions for:
1. Apply SN cut
2. apply an extinction correction
3. Create two BPT diagrams DIG cut variables S2BPT, N2BPT
4. Create two further DIG cut variables: Ha DIG, S2 DIG
4. Compute 3 different metallicity diagnostics: N2S2Ha, O3N2, RS32
5. Compute 4 different K19 corrected metallicity diagnostics 
##############################################################################
"""



import pandas as pd
import numpy as np
import logging
from astropy.io import fits
import pickle
from scipy import linalg
from sklearn.metrics.pairwise import euclidean_distances 
from MUSE_data_processing import *

'''
### ### ### ### ### ### ### ### ### 
###    Input fits file here.    ###
### ### ### ### ### ### ### ### ### 
'''
#plt.imshow(gal_df[30].data,origin='lower')
# change here for your target galaxy!
gal_df = fits.open(r"C:\Users\qihan\Desktop\Data processing\main_files\MUSE_fits_data\NGC7496_MAPS_copt_0.89asec.fits")
gal_name = 'N7496' # need to change it! name: N + number, not NGC + number!
output_version = '_ver1' # need to be changed, this is the version for the name of output file.
types = '_copt' # please do not change it, normally, we should use copt rather than native.
#types = '_native' # If you want to use native raw data, please use this option.

# you need to change the path:
galaxy_info = pd.read_excel("C:/Users/qihan/Desktop/Data processing/main_files/Gal_info_data/galaxydata.xlsx")



# :
out_name = 'MUSE_'
#####################################################################################################
#####################################################################################################
'''
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###                          Get information of galaxy.                     ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
'''
idx = galaxy_info["Gal_ID"] == gal_name
galaxy_info_now = galaxy_info[idx]
inc_angle = float(galaxy_info_now['i'])
Distance = float(galaxy_info_now['D'])
RA_galaxy = float(galaxy_info_now['RA'])
DEC_galaxy = float(galaxy_info_now['DEC'])
PA_galaxy = float(galaxy_info_now['PA_MUSE'])

'''
# Another way to compute the Inclinations for galaxy. qz = 0.2 based on assumptions.
# Please read the function directly if you want to use it. Otherwise, please use the
# inc_angle from above xlsx file. 
# inc_angle = cal_Inclinations_from_axis_ratio_MUSE(log_a_b_from_LEDA, qz = 0.2)
'''
'''
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###             Basic information from original raw data                    ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

There are some baisc information from the raw fits data. 
Please do not change this part unless you have to change them.
'''
x_dim = gal_df[30].data.shape[1] # dim of x axis
y_dim = gal_df[30].data.shape[0] # dim of y axis
# Wavelengths for MUSE dataset. 
wavelengths = np.array([4861.35, 4958.91, 5006.84, 6548.05, 6562.79, 6583.45, 6716.44, 6730.82]) 
# index of H alpha in emission lines fits.file for MUSE. index 30 is Ha emission lines.
Ha_map = gal_df[30] 

###############################################################################
###############################################################################
'''
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
###                      Calculation start here:                            ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
'''
'''
Step one:
    1. Create RA and DEC grid
    2. Apply SN cut
    3. Create Ha DIG cut
    4. Apply extinction correction
'''
RA_grid, DEC_grid = make_RA_DEC_grid_MUSE(gal_df[1].header, x_dim, y_dim) 
gal_df = SN_cut_MUSE(gal_df, threshold=3) # Already changed.
Ha_DIG = determine_DIG_Ha_Zhang17_MUSE(Ha_map, inc_angle, Distance) # Already changed 
gal_df = extinction_correction_MUSE(gal_df, wavelengths)

'''
Step two:
    1. Create S2 BPT cut
    2. Create N2 BPT cut
    3. Create S2 DIG cut
'''
S2_BPT_classification = classify_S2_BPT_MUSE(gal_df)
N2_BPT_classification = classify_N2_BPT_MUSE(gal_df) 
#S2_BPT_classification = classify_S2_BPT_MUSE_ver2(gal_df)
#N2_BPT_classification = classify_N2_BPT_MUSE_ver2(gal_df)
S2_DIG = determine_DIG_S2_Kaplan16_MUSE(gal_df)

'''
Step three:
    1. N2S2Ha Metallicity Diagnostic
    2. O3N2 Metallicity Diagnostic
    3. RS32 Metallicity Diagnostic
'''
Z_N2S2Ha    = compute_Z_N2S2Ha_Dop16_MUSE(gal_df)
Z_O3N2      = compute_Z_O3N2_Curti17_MUSE(gal_df)
Z_RS32      = compute_Z_RS32_Curti20_MUSE(gal_df)

'''
Step four:
    1. O3N2 Metallicity Diagnostic with K19 correction (N2)
    2. O3N2 Metallicity Diagnostic with K19 correction (S2)
    3. RS32 Metallicity Diagnostic with K19 correction (N2)
    4. RS32 Metallicity Diagnostic with K19 correction (S2)
'''
Z_O3N2_kumari_N2 = compute_Z_O3N2_Curti17_MUSE(gal_df, kumari_correction='N2')
Z_O3S2_kumari_N2 = compute_Z_O3S2_Curti17_MUSE(gal_df, kumari_correction='N2')
Z_O3N2_kumari_S2 = compute_Z_O3N2_Curti17_MUSE(gal_df, kumari_correction='S2')
Z_O3S2_kumari_S2 = compute_Z_O3S2_Curti17_MUSE(gal_df, kumari_correction='S2')


'''
Step five:
    1. Create Variable SFR density
    2. Create Variable Projection distance from galaxy center
    3. 4 relation from sun2023
    
If you want to create variable mass, you can add the functions from MUSE_data_processing.
I already created these functions. Make sure you use SFR rather than SFR density and you
need to include variable z which is the redshift. For example z = 0.001721 for NGC5236.
But for MUSE, every galaxies should have different z, please cheak them!
'''
Ha_map = gal_df[30]
SFR_density = calculate_SFR_density_MUSE(Ha_map, Distance, inc_angle)
proj_dist = RA_DEC_to_radius_MUSE(Distance, PA_galaxy, inc_angle, RA_grid, DEC_grid, RA_galaxy, DEC_galaxy)

surface_gas_density = calculate_gas_density_MUSE(SFR_density)


# 4 relation from sun2023:
# Fiducial:
mKS_Fiducial = calculate_mKS_sun_2023_MUSE(SFR_density, -2.40, 1.00)
mES_Fiducial = calculate_mES_sun_2023_MUSE(SFR_density, -2.23, 0.77)
FFTR_Fiducial = calculate_FFTR_sun_2023_MUSE(SFR_density, -2.32, 0.65)
PR_Fiducial = calculate_PR_sun_2023_MUSE(SFR_density, -2.95, 0.93)
# FUV+W4 SFR: FUV_W4_SFR
mKS_FUV_W4_SFR = calculate_mKS_sun_2023_MUSE(SFR_density, -2.34, 0.88)
mES_FUV_W4_SFR = calculate_mES_sun_2023_MUSE(SFR_density, -2.19, 0.67)
FFTR_FUV_W4_SFR = calculate_FFTR_sun_2023_MUSE(SFR_density, -2.28, 0.57)
PR_FUV_W4_SFR = calculate_PR_sun_2023_MUSE(SFR_density, -2.84, 0.84)
# Av-corr Ha SFR: Av_corr_Ha_SFR
mKS_Av_corr_Ha_SFR = calculate_mKS_sun_2023_MUSE(SFR_density, -2.23, 0.93)
mES_Av_corr_Ha_SFR = calculate_mES_sun_2023_MUSE(SFR_density, -2.06, 0.78)
FFTR_Av_corr_Ha_SFR = calculate_FFTR_sun_2023_MUSE(SFR_density, -2.16, 0.62)
PR_Av_corr_Ha_SFR = calculate_PR_sun_2023_MUSE(SFR_density, -2.72, 0.85)
# MW a_co: MW_a_co
mKS_MW_a_co = calculate_mKS_sun_2023_MUSE(SFR_density, -2.43, 0.92)
mES_MW_a_co = calculate_mES_sun_2023_MUSE(SFR_density, -2.26, 0.69)
FFTR_MW_a_co = calculate_FFTR_sun_2023_MUSE(SFR_density, -2.34, 0.62)
PR_MW_a_co = calculate_PR_sun_2023_MUSE(SFR_density, -2.94, 0.86)
# B13 a_co: B13_a_co
mKS_B13_a_co = calculate_mKS_sun_2023_MUSE(SFR_density, -2.36, 1.21)
mES_B13_a_co = calculate_mES_sun_2023_MUSE(SFR_density, -2.17, 0.90)
FFTR_B13_a_co = calculate_FFTR_sun_2023_MUSE(SFR_density, -2.29, 0.75)
PR_B13_a_co = calculate_PR_sun_2023_MUSE(SFR_density, -2.95, 1.08)
# G20 a_co: G20_a_co
mKS_G20_a_co = calculate_mKS_sun_2023_MUSE(SFR_density, -2.22, 1.18)
mES_G20_a_co = calculate_mES_sun_2023_MUSE(SFR_density, -2.11, 0.77)
FFTR_G20_a_co = calculate_FFTR_sun_2023_MUSE(SFR_density, -2.20, 0.76)
PR_G20_a_co = calculate_PR_sun_2023_MUSE(SFR_density, -2.87, 1.05)


'''
Step six:
    1. Create dictionary
    2. Output pkl file
'''
data_dict = {
	'RA':                   RA_grid.flatten(),
	'DEC':                  DEC_grid.flatten(),
    'proj_dist':            proj_dist.flatten(),
	'S2_BPT':               S2_BPT_classification.flatten(),
	'N2_BPT':               N2_BPT_classification.flatten(), 
	'S2_DIG':               S2_DIG.flatten(),
	'Ha_DIG':               Ha_DIG.flatten(),
	'Z_N2S2Ha':             Z_N2S2Ha.flatten(),
	'Z_O3N2':               Z_O3N2.flatten(),
	'Z_RS32':               Z_RS32.flatten(),
    'Z_O3N2_kumari_N2':     Z_O3N2_kumari_N2.flatten(),
    'Z_O3S2_kumari_N2':     Z_O3S2_kumari_N2.flatten(),  
    'Z_O3N2_kumari_S2':     Z_O3N2_kumari_S2.flatten(),
    'Z_O3S2_kumari_S2':     Z_O3S2_kumari_S2.flatten(),
    'SFR_density':          SFR_density.flatten(),
    'surface_gas_density':  surface_gas_density.flatten(),
    'mKS_Fiducial':         mKS_Fiducial.flatten(),
    'mES_Fiducial':         mES_Fiducial.flatten(),
    'FFTR_Fiducial':        FFTR_Fiducial.flatten(),
    'PR_Fiducial':          PR_Fiducial.flatten(),
    'mKS_FUV_W4_SFR':       mKS_FUV_W4_SFR.flatten(),
    'mES_FUV_W4_SFR':       mES_FUV_W4_SFR.flatten(),
    'FFTR_FUV_W4_SFR':      FFTR_FUV_W4_SFR.flatten(),
    'PR_FUV_W4_SFR':        PR_FUV_W4_SFR.flatten(),
    'mKS_Av_corr_Ha_SFR':   mKS_Av_corr_Ha_SFR.flatten(),
    'mES_Av_corr_Ha_SFR':   mES_Av_corr_Ha_SFR.flatten(),
    'FFTR_Av_corr_Ha_SFR':  FFTR_Av_corr_Ha_SFR.flatten(),
    'PR_Av_corr_Ha_SFR':    PR_Av_corr_Ha_SFR.flatten(),
    'mKS_MW_a_co':          mKS_MW_a_co.flatten(),
    'mES_MW_a_co':          mES_MW_a_co.flatten(),
    'FFTR_MW_a_co':         FFTR_MW_a_co.flatten(),
    'PR_MW_a_co':           PR_MW_a_co.flatten(),
    'mKS_B13_a_co':         mKS_B13_a_co.flatten(),
    'mES_B13_a_co':         mES_B13_a_co.flatten(),
    'FFTR_B13_a_co':        FFTR_B13_a_co.flatten(),
    'PR_B13_a_co':          PR_B13_a_co.flatten(),
    'mKS_G20_a_co':         mKS_G20_a_co.flatten(),
    'mES_G20_a_co':         mES_G20_a_co.flatten(),
    'FFTR_G20_a_co':        FFTR_G20_a_co.flatten(),
    'PR_G20_a_co':          PR_G20_a_co.flatten()
	}

result_df = pd.DataFrame(data_dict)
result_df.to_pickle(out_name + gal_name + output_version + types + '.pkl')
logging.info("OK.\n")

'''
This part used to check total SFR of the galaxy based on local SFR from above 
calculation. It is important to check the assumption that Global SFR fomula is
suitable to apply to local regions for a specific galaxy. 
'''
deproj_area = calculate_deproj_area_MUSE(Ha_map, inc_angle, Distance)
Total_SFR_Of_Galaxy = calculate_total_SFR_of_galaxy_MUSE(SFR_density,deproj_area)




