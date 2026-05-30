# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:18:47 2024

@author: Qihan Zou

This code design to create pkl file for spiral arms musk from fits file. 

Intend to examine relationship with arms or 
test clustering methods/classification methods and so on.

Data based on Querejeta et al. from PHANGES website.
"""
import numpy as np
import logging
from astropy.io import fits
import pickle
from MUSE_data_processing import *


#types = '_board_Arm_Musk'

arm_df = fits.open(r"C:\Users\qihan\Desktop\Data processing\main_files\Querejeta_spiral_arm_musk\NGC4321_spiral_mask_narrow.fits")

out_path = 'MUSE_'
gal_name = 'N4321'
types = '_narrow_Arm_Musk'
#types = '_board_Arm_Musk'


x_dim = arm_df[0].data.shape[1] # dim of x axis
y_dim = arm_df[0].data.shape[0] # dim of y axis
RA_arm, DEC_arm = make_RA_DEC_grid_MUSE(arm_df[0].header, x_dim, y_dim) 


data_dict = {
    'RA_arm':               RA_arm.flatten(),
    'DEC_arm':              DEC_arm.flatten(),
    'Arm_Musk':             arm_df[0].data.flatten()
	}

result_df = pd.DataFrame(data_dict)
result_df.to_pickle(out_path + gal_name + types + '.pkl')
logging.info("OK.\n")


