# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:18:47 2024

@author: qihan
"""
import numpy as np
import logging
from astropy.io import fits
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


out_path = 'C:/Users/qihan/Desktop/q/MUSE_pkl/'
gal_name = 'N1385'
version = '_ver1'
types = '_copt_Arm_Musk'

arm_df = fits.open(r"C:\Users\qihan\Desktop\q\Querejeta_spiral_arm_musk\NGC1385_spiral_mask_narrow.fits")
x_dim = arm_df[0].data.shape[1] # dim of x axis
y_dim = arm_df[0].data.shape[0] # dim of y axis
RA_arm, DEC_arm = make_RA_DEC_grid_MUSE(arm_df[0].header, x_dim, y_dim) 


    
data_dict = {
    'RA_arm':               RA_arm.flatten(),
    'DEC_arm':              DEC_arm.flatten(),
    'Arm_Musk':             arm_df[0].data.flatten()
	}

result_df = pd.DataFrame(data_dict)
result_df.to_pickle(out_path + gal_name + version + types + '.pkl')
logging.info("OK.\n")


