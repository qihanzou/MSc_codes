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
import logging

out_path = '../Data/Handmade/Full_galaxy_Z_dfs/'
wavelengths = np.array([3726.0, 3729.0, 4861.3, 5007.0, 6562.8, 6583.0,6716.0,6731.0])

HALPHA = 8 # index of H alpha in line dfs

if __name__=='__main__':
    meta = M83_metadata()
	logging.info("              "+meta['Gal_ID']+"              ")
	logging.info('------------------------------------------\n\n')
	gal_df = open_line_df(meta['Gal_ID'])
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
	data_dict = {
	'RA':         RA_grid.flatten(),
	'DEC':        DEC_grid.flatten(),
	'S2_BPT':     S2_BPT_classification.flatten(),
	'N2_BPT':     N2_BPT_classification.flatten(),
	'S2_DIG':     S2_DIG.flatten(),
	'Ha_DIG':     Ha_DIG.flatten(),
	'Z_N2S2Ha':   Z_N2S2Ha.flatten(),
	'e_Z_N2S2Ha': e_Z_N2S2Ha.flatten(),
	'Z_O3N2':     Z_O3N2.flatten(),
	'e_Z_O3N2':   e_Z_O3N2.flatten(),
	'Z_RS32':     Z_RS32.flatten(),
	'e_Z_RS32':   e_Z_RS32.flatten(),
	'Z_N2O2':     Z_N2O2.flatten(),
	'e_Z_N2O2':   e_Z_N2O2.flatten(),
	'Z_O3S2_old':         Z_O3S2_old.flatten(),
	'e_Z_O3S2_old':       e_Z_O3S2_old.flatten(),
	'Z_O3N2_kumari_N2':   Z_O3N2_kumari_N2.flatten(),
	'e_Z_O3N2_kumari_N2': e_Z_O3N2_kumari_N2.flatten(),
	'Z_O3N2_kumari_S2':   Z_O3N2_kumari_S2.flatten(),
	'e_Z_O3N2_kumari_S2': e_Z_O3N2_kumari_S2.flatten(),
	'Z_O3S2_kumari_N2':   Z_O3S2_kumari_N2.flatten(),
	'e_Z_O3S2_kumari_N2': e_Z_O3S2_kumari_N2.flatten(),
	'Z_O3S2_kumari_S2':   Z_O3S2_kumari_S2.flatten(),
	'e_Z_O3S2_kumari_S2': e_Z_O3S2_kumari_S2.flatten(),
	'O3':   O3.flatten(), 
	'e_O3': e_O3.flatten()
	}
	# save it
	result_df = pd.DataFrame(data_dict)
	result_df.to_pickle(out_path+meta['Gal_ID']+'.pkl')
	logging.info("OK.\n")
