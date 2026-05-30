# -*- coding: utf-8 -*-
"""
Created on Wed May 22 19:28:14 2024

@author: qihan
"""


import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import pickle

import matplotlib.pyplot as plt
import scipy.spatial as spatial
import itertools
import random
from sklearn.metrics import mean_squared_error 
import statistics
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt
from scipy import linalg
import matplotlib.pyplot as plt
from matplotlib import colors
from scipy.stats import linregress
from MUSE_data_processing import *


trainset = pd.read_pickle(r"C:\Users\qihan\Desktop\Data processing\main_files\MUSE_pkl_data\MUSE_N4321_ver1_copt.pkl") 
arm_musk = pd.read_pickle(r"C:\Users\qihan\Desktop\Data processing\main_files\MUSE_pkl_data\MUSE_N4321_narrow_Arm_Musk.pkl")
gal_name = 'N4321'


musk1 = arm_musk

RA_arm = list()
DEC_arm = list()
Arm_musk = list()
for i in range(len(musk1)):
    print(i)
    if musk1['Arm_Musk'][i] > 0:
        RA_arm.append(musk1['RA_arm'][i])
        DEC_arm.append(musk1['DEC_arm'][i])
        Arm_musk.append(musk1['Arm_Musk'][i])

arm_data = {
  "RA_arm": RA_arm,
  "DEC_arm": DEC_arm,
  "Arm_Musk": Arm_musk
}

#load data into a DataFrame object:
arm_data = pd.DataFrame(arm_data)










### Do not need to changse this part ###
### Get information of galaxy. ###
galaxy_info = pd.read_excel("C:/Users/qihan/Desktop/Data processing/main_files/Gal_info_data/galaxydata.xlsx")
idx = galaxy_info["Gal_ID"] == gal_name
galaxy_info_now = galaxy_info[idx]
inc_angle = float(galaxy_info_now['i'])
Distance = float(galaxy_info_now['D'])
RA_galaxy = float(galaxy_info_now['RA'])
DEC_galaxy = float(galaxy_info_now['DEC'])
PA_galaxy = float(galaxy_info_now['PA_MUSE'])



### choose DIG diagnostic ###
#DIG_CUT = "N2_BPT"
#DIG_CUT = "S2_BPT"
DIG_CUT = "S2_DIG"
#############################
### choose HII Train Z diagnostic ###
#Z_train = "Z_N2S2Ha"
Z_train = "Z_O3N2"
#Z_train = "Z_RS32"





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
RA2 = arm_data['RA_arm']
DEC2 = arm_data['DEC_arm']
X2 = RA_DEC_to_xy_MUSE(RA2, DEC2, RA_galaxy, DEC_galaxy, PA_galaxy, inc_angle, Distance)
X2 = np.transpose(X2)
plt.scatter(X2[:, 0], X2[:, 1], c = arm_data['Arm_Musk'], marker='.', s=0.5)
plt.colorbar()
plt.title(f'{gal_name} arm regions') 
plt.xlabel('x (kpc)', color = "black")
plt.ylabel('y (kpc)', color = "black")
#plt.savefig(f'{gal_name}_{Z_train}_{DIG_CUT}_MUSK.eps', dpi=300)
#plt.savefig(f'{gal_name}_{Z_train}_{DIG_CUT}_MUSK.png', dpi=300)
plt.close()



around_arm_or_not = list(itertools.repeat(0, trainset.shape[0]))
RA_arm = arm_data['RA_arm']
DEC_arm = arm_data['DEC_arm']
min_dist = 0.01
for j in range(len(trainset)):
    print(j)
    check_point = pd.DataFrame(trainset).iloc[j]
    RApoint = check_point['RA']
    DECpoint = check_point['DEC']
    #D_point_to_arm = deprojected_distances(RApoint, DECpoint, RA2 = RA_arm, DEC2 = DEC_arm) 
    
    D_point_to_arm = deprojected_distances_MUSE(Distance, PA_galaxy, inc_angle, RApoint, DECpoint, RA2 = RA_arm, DEC2 = DEC_arm)
    for k in range(D_point_to_arm.shape[1]):
        if D_point_to_arm[0][k] < min_dist:
            around_arm_or_not[j] = arm_data['Arm_Musk'][k]
            #print(k)
        else:
            continue 
            
            
        
    


              
I_arm = around_arm_or_not

data1 = trainset
data1['I_arm'] = I_arm
arm1 = pd.DataFrame({"I_arm": I_arm})
idx = arm1['I_arm'] > 0 
data1 = trainset
data1 = data1.reset_index(drop=True)
data1 = data1[idx]


RA1 = trainset['RA']
DEC1 = trainset['DEC']
X1 = RA_DEC_to_xy_MUSE(RA1, DEC1, RA_galaxy, DEC_galaxy, PA_galaxy, inc_angle, Distance)
X1 = np.transpose(X1)
plt.scatter(X1[:, 0], X1[:, 1], c = trainset[Z_train], cmap='hot_r', marker='.', s=1)
#plt.colorbar()
#plt.title(f'{gal_name} {Z_train} {DIG_CUT} HII regions') 
#plt.xlabel('x (kpc)', color = "black")
#plt.ylabel('y (kpc)', color = "black")
#plt.savefig(f'{gal_name}_{Z_train}_{DIG_CUT}.eps', dpi=300)
#plt.savefig(f'{gal_name}_{Z_train}_{DIG_CUT}.png', dpi=300)
#plt.close()

#RA2 = arm_data['RA_arm']
#DEC2 = arm_data['DEC_arm']
#X2 = RA_DEC_to_xy_MUSE(RA2, DEC2, RA_galaxy, DEC_galaxy, PA_galaxy, inc_angle, Distance)
#X2 = np.transpose(X2)
#plt.scatter(X2[:, 0], X2[:, 1], c = arm_data['Arm_Musk'], marker='.', s=0.5)
#plt.title(f'{gal_name} arm regions') 
#plt.xlabel('x (kpc)', color = "black")
#plt.ylabel('y (kpc)', color = "black")

RA4 = data1['RA']
DEC4 = data1['DEC']
X4 = RA_DEC_to_xy_MUSE(RA4, DEC4, RA_galaxy, DEC_galaxy, PA_galaxy, inc_angle, Distance)
X4 = np.transpose(X4)
plt.scatter(X4[:, 0], X4[:, 1], c = data1['I_arm'], marker='.', s=1)
plt.colorbar()
plt.title(f'{gal_name} {Z_train} {DIG_CUT} arm regions') 
plt.xlabel('x (kpc)', color = "black")
plt.ylabel('y (kpc)', color = "black")
#plt.savefig(f'{gal_name}_{Z_train}_{DIG_CUT}_arms.eps', dpi=300)
#plt.savefig(f'{gal_name}_{Z_train}_{DIG_CUT}_arms.png', dpi=300)
#plt.close()


#names = ['RA', 'DEC', 'proj_dist', 'S2_BPT', 'N2_BPT', 'S2_DIG', 'Ha_DIG', 'Z_N2S2Ha', 'e_Z_N2S2Ha', 'Z_O3N2', 'e_Z_O3N2', 'Z_RS32', 'e_Z_RS32', 'Z_O3N2_kumari_N2', 'e_Z_O3N2_kumari_N2', 'Z_O3S2_kumari_N2', 'e_Z_O3S2_kumari_N2', 'Z_O3N2_kumari_S2', 'e_Z_O3N2_kumari_S2', 'Z_O3S2_kumari_S2', 'e_Z_O3S2_kumari_S2', 'SFR_density', 'I_arm']

#np.savetxt(f"{gal_name}_{Z_train}_{DIG_CUT}_arm.csv", data1, header=','.join(names) , delimiter=",", comments='')


