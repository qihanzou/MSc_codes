# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 13:45:31 2022

@author: tingjinc
"""

import numpy as np 
from scipy import linalg
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
from sklearn.metrics.pairwise import euclidean_distances 

from scipy.optimize import differential_evolution
from scipy.special import kv
from scipy.special import gamma
from Internal_Fun import *




def RA_DEC_to_radius(RA, DEC):
    return deprojected_distances(RA, DEC, 204.253792, -29.86575).T[0]



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
    #PA = np.radians(meta['PA'])
    PA = np.radians(54)
    #i  = np.radians(meta['i'])
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




"""
MLE estimation: with covariates X. If not, use the function MLE_fitE
y: the response
X: covariates
D: distance matrix
cov_model: types of covariance functions. See cor_mat function for more details. 
eta_ini: initial value for L-BFGS-B 
nug: nugget effect.  if == False: nugget free model; if == True, a model with nugget
opt: optimization algorithm: LB(L_BFGS-B), DE (differential evolution)
"""
def MLE_fit_1(y, X, coords, cov_model, eta_ini, nug, opt, lo_bound, up_bound):
    
    D = deprojected_distances(coords["RA"], coords["DEC"], RA2 = None, DEC2 = None)
    #bound = cov_bound(cov_model = cov_model, nug = nug)
    q = eta_ini.size + 1
    bound = []
    for i in range(q-1):
        bound.append((lo_bound[i], up_bound[i]))
    
    ## Optimization algorithms
    if opt == "LB": 
        soln = minimize(profile_nll, eta_ini, bounds=bound, args=(y, X, D, cov_model, nug), method="L-BFGS-B")
    # Nelder-Mead, L-BFGS-B
    if opt == "DE":
        soln = differential_evolution(profile_nll, bound, args = (y, X, D, cov_model, nug), seed = 2019, tol = 0.0001)
    eta = soln.x
    suc = soln.success
    nll = soln.fun
    
    n = y.size
    cormat = cor_mat(D = D, eta = eta, cov_model = cov_model, nug = nug)
    L = np.linalg.cholesky(cormat)
    white_X = np.linalg.solve(L, X)
    white_y = np.linalg.solve(L, y)
    
    beta_est = linalg.solve(white_X.T @ white_X, white_X.T @ white_y, assume_a='sym')
    white_resids = white_y - white_X @ beta_est
    sill = white_resids.T @ white_resids/n
    beta_var = linalg.inv(white_X.T @ white_X)*sill

    ## reparametrization
    theta_est = np.copy(eta)
    if nug == True: 
        #theta_est [ntheta-1] = theta[ntheta-1]*sill    # nugget effect
        nug_effect = sill*eta[q-2]
        psill = sill*(1-eta[q-2])
        theta_est = np.append(theta_est[0:(q-2)], [psill, nug_effect])
    elif nug == False:
        psill = sill
        theta_est = np.append(theta_est,  psill)
    return beta_est, theta_est, suc, nll, beta_var



