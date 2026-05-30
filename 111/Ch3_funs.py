# -*- coding: utf-8 -*-
"""
Created on Thu May 22 22:54:57 2025

@author: qihan
"""
import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

def make_RA_DEC_grid_general(header,x_dim,y_dim):
    # Given a header file, create a grid of RA//DEC for each pixel in that file.
    world = WCS(header)
    x = np.arange(x_dim)
    y = np.arange(y_dim)
    X, Y = np.meshgrid(x, y)
    RA_grid, DEC_grid = world.wcs_pix2world(X, Y, 0)
    return RA_grid, DEC_grid 

def compute_molecular_gas_surface_density_from_CO21(ICO21, inc_angle_deg):
#    This function can be used to calculate molecular gas surface density from CO(2-1) emission line. 
#    inc_angle is the inclination angle in the galaxy info Excel file. It changes based on different assumptions (from different resources).
#    parameters:
#        1. R21 = 0.65, 1/R21*alpha_Co(1-0) gives the constant around 6.7. adopted CO(2-1) to CO(1-0) line ratio (den Brok et al. 2021; Leroy et al. 2022)
#        2. alpha_Co(1-0)=4.35, unit: M_sun pc-2 (K km s-1)^-1 above are coefficients. 
#        3. ICO21, Integrated CO(2-1) line intensity (moment-0). unit: K km s−1
#        4. cosi, inclination correction.
#    return:
#        1. molecular gas surface density (unit: M_sun pc-2)
    i  =  np.radians(inc_angle_deg) 
    gas_density = 6.7*ICO21.data*np.cos(i)  
    return gas_density # unit: M_sun pc-2

def RA_DEC_to_xy_general(RA, DEC, RA_gal, DEC_gal, PA, inc, D):
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

def euclidean_distances_general(X, Y=None):
    X = np.atleast_2d(X)
    if Y is None:
        Y = X
    else:
        Y = np.atleast_2d(Y)

    dists_squared = (np.sum(X ** 2, axis=1)[:, np.newaxis] +
                     np.sum(Y ** 2, axis=1)[np.newaxis, :] -
                     2 * np.dot(X, Y.T))
    
    dists_squared = np.maximum(dists_squared, 0.0)
    return np.sqrt(dists_squared)

def deprojected_distances_general(Dist, PA_gal, inc_gal, RA1, DEC1, RA2 = None, DEC2 = None):
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
    #deg_dists = euclidean_distances(vec1, np.array(np.matrix(vec2))) # original
    deg_dists = euclidean_distances_general(vec1, np.array(np.matrix(vec2)))
    rad_dists = np.radians(deg_dists)
    # 4: Convert angular offsets to kpc distances using D, and the small-angle approximation.
    # Mpc_dists = rad_dists * meta['D']
    Mpc_dists = rad_dists * Dist
    kpc_dists = Mpc_dists * 1000
    return kpc_dists

def RA_DEC_to_radius_general(Dist, PA_gal, inc_gal, RA, DEC, RA_galaxy, DEC_galaxy):
    return deprojected_distances_general(Dist, PA_gal, inc_gal, RA, DEC, RA2 = RA_galaxy, DEC2 = DEC_galaxy).T[0]


