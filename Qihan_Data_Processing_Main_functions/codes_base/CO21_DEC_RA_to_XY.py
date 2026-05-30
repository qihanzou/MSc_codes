# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 10:32:20 2024

@author: qihan
"""

import numpy as np
from numpy import inf
from numpy import nan
from numpy import *
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from TYPHOON_wrangling import *
import pandas as pd 
from astropy.wcs import WCS
import astropy.units as u
from extinction import ccm89, apply
from Z_diags import *

def RA_DEC_to_xy(RA, DEC, meta):
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
    
    # Check parameters
    try:
        meta['PA']
    except KeyError:
        assert False, "Error: PA not defined for metadata"
    try:
        meta['i']
    except KeyError:
        assert False, "Error: i not defined for metadata"
    try:
        meta['D']
    except KeyError:
        assert False, "Error: D not defined for metadata"
        assert len(RA) == len(DEC), "Error: len of RA1 must match len of DEC1"
    RA = np.array(RA) - meta['RA']
    DEC = np.array(DEC) - meta['DEC']
    # Now onto the maths
    PA = np.radians(meta['PA'])
    i  = np.radians(meta['i'])
    # 1: Rotate RA, DEC by PA to get y (major axis direction) and x (minor axis direction)
    x = RA*np.cos(PA) - DEC*np.sin(PA)
    y = DEC*np.cos(PA) + RA*np.sin(PA)
    # 2: Stretch x values to remove inclination effects
    x = x /np.cos(i)
    # 3: Convert deg to kpc
    x = np.radians(x)*meta['D']*1000
    y = np.radians(y)*meta['D']*1000
    return x, y