# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 21:23:43 2024

@author: qihan
"""

'''
classify_DIG_and_apply_Z_diagnostics.py
* apply SN cut
* apply an extinction correction
* Using two BPT diagrams, classify spaxels into star forming, Seyfert, or AGN
* Compute 8 different metallicity diagnostics and their associated errors
'''
from TYPHOON_wrangling import *
import pandas as pd
import numpy as np
import logging
from astropy.io import fits
from SFR_calculation import *
import pickle
from AST2 import *
from Z_diags import *
import numpy as np 
from scipy import linalg
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
from sklearn.metrics.pairwise import euclidean_distances 
from scipy.optimize import differential_evolution
from scipy.special import kv
from scipy.special import gamma
from Internal_Fun import *
from MUSE_data_processing import *


#plt.imshow(gal_df[30].data,origin='lower')
gal_df = fits.open(r"C:\Users\qihan\Desktop\q\MUSE\NGC4321_MAPS_copt_1.16asec.fits")
gal_name = 'N4321'
version = '_ver2' # need to be changed
types = '_copt'
#types = '_native'



### Do not need to changse this part ###
### Get information of galaxy. ###
galaxy_info = pd.read_excel("C:/Users/qihan/Desktop/q/galaxydata.xlsx")
idx = galaxy_info["Gal_ID"] == gal_name
galaxy_info_now = galaxy_info[idx]
inc_angle = float(galaxy_info_now['i'])
Distance = float(galaxy_info_now['D'])
RA_galaxy = float(galaxy_info_now['RA'])
DEC_galaxy = float(galaxy_info_now['DEC'])
PA_galaxy = float(galaxy_info_now['PA_MUSE'])




x_dim = gal_df[30].data.shape[1] # dim of x axis
y_dim = gal_df[30].data.shape[0] # dim of y axis
out_path = 'C:/Users/qihan/Desktop/q/MUSE_pkl/'
#inc_angle = cal_Inclinations_from_axis_ratio_MUSE(log_a_b_from_LEDA, qz = 0.2)
wavelengths = np.array([4861.35, 4958.91, 5006.84, 6548.05, 6562.79, 6583.45, 6716.44, 6730.82]) # for MUSE
Ha_map = gal_df[30] # index of H alpha in emission lines fits.file for MUSE

### Calculation start here:
RA_grid, DEC_grid = make_RA_DEC_grid_MUSE(gal_df[1].header, x_dim, y_dim) 
gal_df = SN_cut_MUSE(gal_df, threshold=3) # Already changed.
Ha_DIG = determine_DIG_Ha_Zhang17_MUSE(Ha_map, inc_angle, Distance) # Already changed # do this before extinction correction for some reason?
gal_df = extinction_correction_MUSE(gal_df, wavelengths)

### BPT:
S2_BPT_classification = classify_S2_BPT_MUSE(gal_df)
N2_BPT_classification = classify_N2_BPT_MUSE(gal_df) 
#S2_BPT_classification = classify_S2_BPT_MUSE_ver2(gal_df)
#N2_BPT_classification = classify_N2_BPT_MUSE_ver2(gal_df)
S2_DIG = determine_DIG_S2_Kaplan16_MUSE(gal_df)

### Metallicity Diagnostic:
Z_N2S2Ha, e_Z_N2S2Ha  = compute_Z_N2S2Ha_Dop16_MUSE(gal_df)
Z_O3N2, e_Z_O3N2      = compute_Z_O3N2_Curti17_MUSE(gal_df)
Z_RS32, e_Z_RS32      = compute_Z_RS32_Curti20_MUSE(gal_df)


Z_O3N2_kumari_N2, e_Z_O3N2_kumari_N2 = compute_Z_O3N2_Curti17_MUSE(gal_df, kumari_correction='N2')
Z_O3S2_kumari_N2, e_Z_O3S2_kumari_N2 = compute_Z_O3S2_Curti17_MUSE(gal_df, kumari_correction='N2')
Z_O3N2_kumari_S2, e_Z_O3N2_kumari_S2 = compute_Z_O3N2_Curti17_MUSE(gal_df, kumari_correction='S2')
Z_O3S2_kumari_S2, e_Z_O3S2_kumari_S2 = compute_Z_O3S2_Curti17_MUSE(gal_df, kumari_correction='S2')

### SFR and proj_dist:
Ha_map = gal_df[30]
SFR_density = calculate_SFR_density_MUSE(Ha_map, Distance, inc_angle)
proj_dist = RA_DEC_to_radius_MUSE(Distance, PA_galaxy, inc_angle, RA_grid, DEC_grid, RA_galaxy, DEC_galaxy)
    
data_dict = {
	'RA':                   RA_grid.flatten(),
	'DEC':                  DEC_grid.flatten(),
    'proj_dist':            proj_dist.flatten(),
	'S2_BPT':               S2_BPT_classification.flatten(),
	'N2_BPT':               N2_BPT_classification.flatten(), # We only nedd this cut right?
	'S2_DIG':               S2_DIG.flatten(),
	'Ha_DIG':               Ha_DIG.flatten(),
	'Z_N2S2Ha':             Z_N2S2Ha.flatten(),
	'e_Z_N2S2Ha':           e_Z_N2S2Ha.flatten(),
	'Z_O3N2':               Z_O3N2.flatten(),
	'e_Z_O3N2':             e_Z_O3N2.flatten(),
	'Z_RS32':               Z_RS32.flatten(),
	'e_Z_RS32':             e_Z_RS32.flatten(),
    'Z_O3N2_kumari_N2':     Z_O3N2_kumari_N2.flatten(),
    'e_Z_O3N2_kumari_N2':   e_Z_O3N2_kumari_N2.flatten(),
    'Z_O3S2_kumari_N2':     Z_O3S2_kumari_N2.flatten(),
    'e_Z_O3S2_kumari_N2':   e_Z_O3S2_kumari_N2.flatten(),   
    'Z_O3N2_kumari_S2':     Z_O3N2_kumari_S2.flatten(),
    'e_Z_O3N2_kumari_S2':   e_Z_O3N2_kumari_S2.flatten(),
    'Z_O3S2_kumari_S2':     Z_O3S2_kumari_S2.flatten(),
    'e_Z_O3S2_kumari_S2':   e_Z_O3S2_kumari_S2.flatten(),
    'SFR_density':          SFR_density.flatten()
	}

result_df = pd.DataFrame(data_dict)
result_df.to_pickle(out_path + gal_name + version + types + '.pkl')
logging.info("OK.\n")


deproj_area = calculate_deproj_area_MUSE(Ha_map, inc_angle, Distance)
Total_SFR_Of_Galaxy = calculate_total_SFR_of_galaxy_MUSE(SFR_density,deproj_area)