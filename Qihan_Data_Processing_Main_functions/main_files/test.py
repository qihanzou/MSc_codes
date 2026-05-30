# -*- coding: utf-8 -*-
"""
Created on Thu May 22 20:33:57 2025

@author: qihan
"""


import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


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

def compute_molecular_gas_surface_density_from_CO21(ICO21, inc_angle_deg):
    '''
    This function can be used to calculate molecular gas surface density from CO(2-1)
    emission line. 
    
    inc_angle is the inclination angle in galaxy info excel file. It changes based
    on different assumptions (from different resources).
    
    parameters:
        1. R21 = 0.65, 1/R21*alpha_Co(1-0) gives constant around 6.7. 
        adopted CO(2-1) to CO(1-0) line ratio (den Brok et al. 2021; Leroy et al. 2022)
        2. alpha_Co(1-0)=4.35, unit: M_sun pc-2 (K km s-1)^-1
        above are coefficents. 
        3. ICO21, Integrated CO(2-1) line intensity (moment-0). unit: K km s−1
        4. cosi, inclination correction.
    return:
        1. molecular gas surface density (unit: M_sun pc-2)
    '''
    i  =  np.radians(inc_angle_deg) 
    gas_density = 6.7*ICO21.data*np.cos(i) 
    return gas_density # unit: M_sun pc-2

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

def euclidean_distances_MUSE(X, Y=None):
    """
    Compute the Euclidean distance matrix between each row of X and each row of Y.
    
    Parameters:
    -----------
    X : array-like of shape (n_samples_X, n_features)
    Y : array-like of shape (n_samples_Y, n_features), optional (default=None)
        If None, uses X as Y.

    Returns:
    --------
    distances : ndarray of shape (n_samples_X, n_samples_Y)
        The distance matrix.
    """
    X = np.atleast_2d(X)
    if Y is None:
        Y = X
    else:
        Y = np.atleast_2d(Y)

    dists_squared = (
        np.sum(X ** 2, axis=1)[:, np.newaxis] +
        np.sum(Y ** 2, axis=1)[np.newaxis, :] -
        2 * np.dot(X, Y.T)
    )
    
    dists_squared = np.maximum(dists_squared, 0.0)
    return np.sqrt(dists_squared)

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
    #deg_dists = euclidean_distances(vec1, np.array(np.matrix(vec2))) # original
    deg_dists = euclidean_distances_MUSE(vec1, np.array(np.matrix(vec2)))
    rad_dists = np.radians(deg_dists)
    # 4: Convert angular offsets to kpc distances using D, and the small-angle approximation.
    # Mpc_dists = rad_dists * meta['D']
    Mpc_dists = rad_dists * Dist
    kpc_dists = Mpc_dists * 1000
    
    return kpc_dists

def RA_DEC_to_radius_MUSE(Dist, PA_gal, inc_gal, RA, DEC, RA_galaxy, DEC_galaxy):
    return deprojected_distances_MUSE(Dist, PA_gal, inc_gal, RA, DEC, RA2 = RA_galaxy, DEC2 = DEC_galaxy).T[0]


def CO21_processing(RA_grid, DEC_grid, gas_density, Dist, PA_gal, inc_gal, RA_galaxy, DEC_galaxy):
    data_dict1 = {
    	'RA':                   RA_grid.flatten(),
    	'DEC':                  DEC_grid.flatten(),
        'CO21_gas_density':     gas_density.flatten()
    	}
    df1 = pd.DataFrame(data_dict1)
    index_above_0 = df1["CO21_gas_density"] > 0
    above_0 = df1[index_above_0]
    RA_full = above_0['RA']
    DEC_full = above_0['DEC']
    CO21_gas_density = above_0['CO21_gas_density']
    #proj_dist = RA_DEC_to_radius(RA_full, DEC_full)
    proj_dist = RA_DEC_to_radius_MUSE(Dist, PA_gal, inc_gal, RA_full, DEC_full, RA_galaxy, DEC_galaxy)
    X_full = np.transpose(RA_DEC_to_xy_MUSE(RA_full, DEC_full, RA_galaxy, DEC_galaxy, PA_gal, inc_gal, Dist))
    
    data_dict_new = {
    	'RA':                   RA_full.values.flatten(),
    	'DEC':                  DEC_full.values.flatten(),
        'X':                    X_full[:,0].flatten(),
        'Y':                    X_full[:,1].flatten(),
        'proj_dist':            proj_dist.flatten(),
        'CO21_mgsd':             CO21_gas_density.values.flatten()
    	}
    data = pd.DataFrame(data_dict_new)
    return data


data_mgs = fits.open(r"C:\Users\qihan\Desktop\ngc1300_12m+7m+tp_co21_15as_strict_mom0.fits") 
galaxy_info = pd.read_excel("C:/Users/qihan/Desktop/Data processing/main_files/Gal_info_data/galaxydata.xlsx")
idx = galaxy_info["Gal_ID"] == 'N1300'
galaxy_info_now = galaxy_info[idx]
inc_angle = float(galaxy_info_now['i'])
Distance = float(galaxy_info_now['D'])
RA_galaxy = float(galaxy_info_now['RA'])
DEC_galaxy = float(galaxy_info_now['DEC'])
PA_galaxy = float(galaxy_info_now['PA_MUSE'])


# check path before use
ICO21 = data_mgs[0] # unit: K km s-1 from PHANGS-ALMA readme.
x_dim = data_mgs[0].data.shape[1] # dim of x axis
y_dim = data_mgs[0].data.shape[0] # dim of y axis
RA_grid, DEC_grid = make_RA_DEC_grid_MUSE(data_mgs[0].header, x_dim, y_dim) 


gas_density = compute_molecular_gas_surface_density_from_CO21(ICO21, inc_angle)
data = CO21_processing(RA_grid, DEC_grid, gas_density, Distance, PA_galaxy, inc_angle, RA_galaxy, DEC_galaxy)










