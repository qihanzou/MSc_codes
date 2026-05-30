# -*- coding: utf-8 -*-
"""
Created on Fri Dec  8 12:13:48 2023

@author: qihan
"""

'''
classify_DIG_and_apply_Z_diagnostics.py

Data flow to convert emission line data into metallicity data.
For each galaxy in our subsample:

* apply SN cut
* apply an extinction correction
* Using two BPT diagrams, classify spaxels into star forming, Seyfert, or AGN
* Compute 8 different metallicity diagnostics and their associated errors
* Store all of this data in a pickled pandas df

Created by: Benjamin Metha
Last updated: Oct 26, 2021
'''

from TYPHOON_wrangling import *
import pandas as pd
import numpy as np
import logging
from astropy.io import fits

from SFR_calculation import *
from surface_gas_density_calculation import *
from stellar_mass_calculation import *
from Four_relations_calculation_sun2023 import *

import pickle
from AST2 import *
from Internal_Fun import *
from BFuns import *

out_path = 'C:/Users/qihan/Desktop/q/'
wavelengths = np.array([3726.0, 3729.0, 4861.3, 5007.0, 6562.8, 6583.0,6716.0,6731.0])

HALPHA = 8 # index of H alpha in line dfs
gal_name = 'N5236'


if __name__=='__main__':
    meta = M83_metadata()
logging.info('------------------------------------------\n\n')
    
gal_df=open_line_df('N5236')
    
RA_grid, DEC_grid = make_RA_DEC_grid(gal_df[0].header)
	# preprocessing
	# replace lines other than O2 where S/N < 3 with nans
SN_cut(gal_df, 3) 
	# replace O2 where S/N < 1 with nans
SN_cut_O2(gal_df, 1)
Ha_DIG = determine_DIG_Ha_Zhang17(gal_df[HALPHA], meta) # do this before extinction correction for some reason?
extinction_correction(gal_df, wavelengths)
	# BPT diagnostics
logging.info("Running BPT diagnostics...")
S2_BPT_classification = classify_S2_BPT(gal_df)
N2_BPT_classification = classify_N2_BPT(gal_df)
logging.info("Isolating HII regions...")
S2_DIG = determine_DIG_S2_Kaplan16(gal_df)
	# Some standard metallicity diagnostics
logging.info("Crafting metallicity maps...")
Z_N2S2Ha, e_Z_N2S2Ha  = compute_Z_N2S2Ha_Dop16(gal_df)
Z_O3N2, e_Z_O3N2      = compute_Z_O3N2_Curti17(gal_df)
Z_RS32, e_Z_RS32      = compute_Z_RS32_Curti20(gal_df)
Z_N2O2, e_Z_N2O2      = compute_Z_N2O2_Dop13(gal_df)
	# and some weird ones
Z_O3S2_old, e_Z_O3S2_old = compute_Z_O3S2_Curti17(gal_df)
Z_O3N2_kumari_N2, e_Z_O3N2_kumari_N2 = compute_Z_O3N2_Curti17(gal_df, kumari_correction='N2')
Z_O3S2_kumari_N2, e_Z_O3S2_kumari_N2 = compute_Z_O3S2_Curti17(gal_df, kumari_correction='N2')
Z_O3N2_kumari_S2, e_Z_O3N2_kumari_S2 = compute_Z_O3N2_Curti17(gal_df, kumari_correction='S2')
Z_O3S2_kumari_S2, e_Z_O3S2_kumari_S2 = compute_Z_O3S2_Curti17(gal_df, kumari_correction='S2')
	# maybe we will want to make our own?
O3, e_O3 = compute_O3(gal_df)
	# package this data into a pandas df
logging.info("Packaging data products...")
f_SII6716 = gal_df[('SII6716')].data
e_SII6716 = gal_df[('SII6716_ERR')].data
f_SII6731 = gal_df[('SII6731')].data
e_SII6731 = gal_df[('SII6731_ERR')].data
f_NII6583 = gal_df[('NII6583')].data
e_NII6583 = gal_df[('NII6583_ERR')].data
f_OIII = gal_df[('OIII5007')].data
e_OIII = gal_df[('OIII5007_ERR')].data
f_Ha   = gal_df[('HALPHA')].data
e_Ha   = gal_df[('HALPHA_ERR')].data
f_Hb   = gal_df[('HBETA')].data
e_Hb   = gal_df[('HBETA_ERR')].data


z = 0.001721

fits_file = fits.open('C:/Users/qihan/Desktop/q/N5236_lowres_cal_1_comp_WCS.fits')
meta=meta_getter('N5236')
Ha_map = fits_file[8]
# SFR density:
SFR_density = calculate_SFR_density(Ha_map, meta)

# Surface gas density: original KS formula.
surface_gas_density = calculate_gas_density(SFR_density)
#positions_are_NaNs_in_surface_gas_density = isnan(surface_gas_density)
#surface_gas_density[positions_are_NaNs_in_surface_gas_density] = 0

# MASS from E11...
line_df_2 = SN_cut(fits_file, threshold=3)
line_df_3 = SN_cut_O2(line_df_2, threshold=1)
line_df_4= extinction_correction(line_df_3, wavelengths, R_V=3.1)
Ha_map = line_df_4[8]
deproj_area = calculate_deproj_area(Ha_map, meta)
SFR_Of_Galaxy = calculate_SFR_of_galaxy(SFR_density,deproj_area)
# mass:
stellar_mass_E11 = calculate_stellar_mass_E11(SFR_Of_Galaxy,z)
stellar_mass_C14_1 = calculate_stellar_mass_C14_1(SFR_Of_Galaxy,z)
stellar_mass_O10 = calculate_stellar_mass_O10(SFR_Of_Galaxy,z)
stellar_mass_C09_1 = calculate_stellar_mass_C09_1(SFR_Of_Galaxy,z)
stellar_mass_S07 = calculate_stellar_mass_S07(SFR_Of_Galaxy,z)
stellar_mass_W12 = calculate_stellar_mass_W12(SFR_Of_Galaxy,z)
stellar_mass_whitaker2012star = calculate_stellar_mass_whitaker2012star(SFR_Of_Galaxy,z)
stellar_mass_tomczak2016sfr_all_galaxies = calculate_stellar_mass_tomczak2016sfr_all_galaxies(SFR_Of_Galaxy,z)
stellar_mass_tomczak2016sfr_star_forming_galaxies = calculate_stellar_mass_tomczak2016sfr_star_forming_galaxies(SFR_Of_Galaxy,z)
#positions_are_NaNs_in_stellar_mass_E11 = isnan(stellar_mass_E11)
#stellar_mass_E11[positions_are_NaNs_in_stellar_mass_E11] = 0
deproj_area = calculate_deproj_area(Ha_map, meta)
SFR1 = SFR_density * deproj_area
MASS_ALMA = 10**((log10(SFR1)+10.17)/(-0.32) +10)


# 4 relation from sun2023:
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


with open('N5236_projdist.pkl', 'rb') as f:
    data = pickle.load(f)
    print(data.columns.values)
    
df = pd.DataFrame(data)

RA = df['RA']
DEC = df['DEC']
proj_dist = RA_DEC_to_radius(RA, DEC)

    
data_dict = {
	'RA':                   RA_grid.flatten(),
	'DEC':                  DEC_grid.flatten(),
    'Ha_map':               Ha_map.data.flatten(),
    'proj_dist':            proj_dist.flatten(),
	'S2_BPT':               S2_BPT_classification.flatten(),
	'N2_BPT':               N2_BPT_classification.flatten(),
	'S2_DIG':               S2_DIG.flatten(),
	'Ha_DIG':               Ha_DIG.flatten(),
    'S16':                  f_SII6716.flatten(),
    'S16e':                 e_SII6716.flatten(),
    'S31':                  f_SII6731.flatten(),
    'S31e':                 e_SII6731.flatten(),
    'NII':                  f_NII6583.flatten(),
    'NIIe':                 e_NII6583.flatten(),
	'Z_N2S2Ha':             Z_N2S2Ha.flatten(),
	'e_Z_N2S2Ha':           e_Z_N2S2Ha.flatten(),
	'Z_O3N2':               Z_O3N2.flatten(),
	'e_Z_O3N2':             e_Z_O3N2.flatten(),
	'Z_RS32':               Z_RS32.flatten(),
	'e_Z_RS32':             e_Z_RS32.flatten(),
	'Z_N2O2':               Z_N2O2.flatten(),
	'e_Z_N2O2':             e_Z_N2O2.flatten(),
	'Z_O3S2_old':           Z_O3S2_old.flatten(),
	'e_Z_O3S2_old':         e_Z_O3S2_old.flatten(),
	'Z_O3N2_kumari_N2':     Z_O3N2_kumari_N2.flatten(),
	'e_Z_O3N2_kumari_N2':   e_Z_O3N2_kumari_N2.flatten(),
	'Z_O3N2_kumari_S2':     Z_O3N2_kumari_S2.flatten(),
	'e_Z_O3N2_kumari_S2':   e_Z_O3N2_kumari_S2.flatten(),
	'Z_O3S2_kumari_N2':     Z_O3S2_kumari_N2.flatten(),
	'e_Z_O3S2_kumari_N2':   e_Z_O3S2_kumari_N2.flatten(),
	'Z_O3S2_kumari_S2':     Z_O3S2_kumari_S2.flatten(),
	'e_Z_O3S2_kumari_S2':   e_Z_O3S2_kumari_S2.flatten(),
	'O3':                   O3.flatten(), 
	'e_O3':                 e_O3.flatten(),
    'SFR_density':          SFR_density.flatten(),
    'surface_gas_density':  surface_gas_density.flatten(),
    'stellar_mass_E11':     stellar_mass_E11.flatten(),
    'stellar_mass_C14_1':   stellar_mass_C14_1.flatten(),
    'stellar_mass_O10':     stellar_mass_O10.flatten(),
    'stellar_mass_C09_1':   stellar_mass_C09_1.flatten(),
    'stellar_mass_S07':     stellar_mass_S07.flatten(),
    'stellar_mass_W12':     stellar_mass_W12.flatten(),
    'stellar_mass_w12s':    stellar_mass_whitaker2012star.flatten(),
    'stellar_mass_t16all':  stellar_mass_tomczak2016sfr_all_galaxies.flatten(),
    'stellar_mass_t16form': stellar_mass_tomczak2016sfr_star_forming_galaxies.flatten(),
    'MASS_ALMA':            MASS_ALMA.flatten(),
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
result_df.to_pickle(out_path+'N5236_26_04'+'.pkl')
logging.info("OK.\n")


