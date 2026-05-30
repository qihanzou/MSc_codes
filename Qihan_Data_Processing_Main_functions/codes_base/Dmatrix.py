# -*- coding: utf-8 -*-
"""
Created on Mon Sep 16 15:44:43 2024

@author: qihan
"""
import numpy as np 
from scipy import linalg
from sklearn.metrics.pairwise import euclidean_distances 

def deprojected_distances(RA1, DEC1, RA2 = None, DEC2 = None):
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
    # Check parameters
    #try:
    #    meta['PA'] 
    #except KeyError:
    #    assert False, "Error: PA not defined for metadata"
    #try:
    #    meta['i'] 
    #except KeyError:
    #    assert False, "Error: i not defined for metadata"
    #try:
    #    meta['D'] 
    #except KeyError:
    #    assert False, "Error: D not defined for metadata"
    
    # If RA1 and DEC1 are arrays, they must have the same length.
    # If one of them is a float, they must both be floats.
    # You can't supply only one of RA2 and DEC2
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
    
    # Now onto the maths
    PA = np.radians(54)
    i = np.radians(15.3)
    # 1: Rotate RA, DEC by PA to get y (major axis direction) and x (minor axis direction)
    x1 = RA1*np.cos(PA) - DEC1*np.sin(PA)
    y1 = DEC1*np.cos(PA) + RA1*np.sin(PA)
    x2 = RA2*np.cos(PA) - DEC2*np.sin(PA)
    y2 = DEC2*np.cos(PA) + RA2*np.sin(PA)
    # 2: Stretch x values to remove inclination effects
    long_x1 = x1 /np.cos(i)
    long_x2 = x2 /np.cos(i)
    # 3: Compute Euclidean Distances between x1,y1 and x2,y2 to get angular offsets (degrees).
    vec1 = np.stack((y1, long_x1)).T
    vec2 = np.stack((y2, long_x2)).T
    deg_dists = euclidean_distances(vec1, vec2)
    rad_dists = np.radians(deg_dists)
    # 4: Convert angular offsets to kpc distances using D, and the small-angle approximation.
    # Mpc_dists = rad_dists * meta['D']
    Mpc_dists = rad_dists * 4.66
    kpc_dists = Mpc_dists * 1000
    
    return kpc_dists