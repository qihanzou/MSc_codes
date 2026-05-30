# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 10:44:54 2022

@author: tingjinc
"""

import numpy as np
from numpy.linalg import inv
import numpy as np 
from scipy import linalg
from scipy.optimize import minimize
from scipy.optimize import differential_evolution
from scipy.special import kv
from scipy.special import gamma


"""
Negative log-likelihood function
D: distance matrix
"""

def profile_nll(eta, y, X, D, cov_model, nug):
    n = y.size
    cormat = cor_mat(D = D, eta = eta, cov_model = cov_model, nug = nug)   
    L = np.linalg.cholesky(cormat)
    log_det_L = np.sum(np.log(np.diag(L)))
    log_det_cormat = 2*log_det_L
    white_X = np.linalg.solve(L, X)
    white_y = np.linalg.solve(L, y)

    beta = linalg.solve(white_X.T @ white_X, white_X.T @ white_y, assume_a='sym') # singular happen here
    #resids =  y - X @ best
    #white_resids = np.linalg.solve(L, resids)
    white_resids = white_y - white_X @ beta
    nll2 = n*np.log(2*np.pi) + n + log_det_cormat + n*np.log(white_resids.T @ white_resids/n)
    return nll2/n




def profile_nllE(eta, y, D, cov_model, nug):
    n = y.size
    cormat = cor_mat(D = D, eta = eta, cov_model = cov_model, nug = nug)   
    L = np.linalg.cholesky(cormat)
    log_det_L = np.sum(np.log(np.diag(L)))
    log_det_cormat = 2*log_det_L
    white_resids = np.linalg.solve(L, y)
    nll2 = n*np.log(2*np.pi) + n + log_det_cormat + n*np.log(white_resids.T @ white_resids/n)
    return nll2/n


"""
    return L[1:3, 1:5]

Types of covariance functions.
Mainly following Chapter 2.7 of Handbook of Spatial Statistics or geoR package "cov.spatial"
Inside the code, c = nugget/sill
"""
def cor_mat(D, eta, cov_model, nug):
    dphi = D/eta[0]
    if cov_model == "Cau":       # Cauchy covariance function, (2.15) with alpha = 2
        spatial_cor = (1+dphi**2)**(-eta[1])
    elif cov_model == "CauF":      # Cauchy family (2.15). eta[1] in (0,2], eta[2]>0 
        spatial_cor = (1+dphi**(eta[1]))**(-eta[2]/eta[1])   
    elif cov_model == "Cub":       # Cubic 
        spatial_cor = (1-7*dphi**2 + 8.75*dphi**3 - 3.5*dphi**5 + 0.75* dphi**7)*(dphi<1)
    elif cov_model == "Edc":       # Exponentially damped cosine (2.16), eta[1] >= 1/tan(pi/2d)
        spatial_cor = np.exp(-eta[1]*dphi)*np.cos(dphi)
    elif cov_model == "Exp":       # Exponential
        spatial_cor = np.exp(-dphi)
    elif cov_model == "Gau":       # Gaussian
        spatial_cor = np.exp(-dphi**2)
    elif cov_model == "MatOld":       # Matern class (2.13)
        np.fill_diagonal(dphi, 1)
        spatial_cor = (dphi/2)**eta[1]*2*kv(eta[1], dphi)/gamma(eta[1])
    elif cov_model == "Mat":
        dtmp = (2*eta[1])**0.5*dphi
        spatial_cor = 2**(1-eta[1])/gamma(eta[1])*(dtmp)**eta[1]*kv(eta[1],dtmp)
    elif cov_model == "Mat32":
        spatial_cor = (1 + 3**0.5*dphi)*np.exp(-3**0.5*dphi)
    elif cov_model == "Mat52":
        spatial_cor = (1 + 5**0.5*dphi + 5/3*dphi**2)*np.exp(-5**0.5*dphi)
    elif cov_model == "Sph":       # Spherical (2.17)
        spatial_cor = (1-1.5*dphi+0.5*dphi**3)*(dphi<1)
    elif cov_model == "PExp":      # Power Exponential (2.14), eta[1] in (0,2]
        spatial_cor = np.exp(-dphi**eta[1])
    elif cov_model == "Wave":      # Wave: can be negative
        np.fill_diagonal(dphi, 1)
        spatial_cor = np.sin(dphi)/dphi 
    elif cov_model == "GW0":       # GW with kappa = 0, mu \geq 1 + d/2 + kappa
        spatial_cor = (1 - dphi)**eta[1]*(dphi<1)
    elif cov_model == "GW02":       # GW ith kappa = 0 and mu = 2
        spatial_cor = (1 - dphi)**2*(dphi<1)
    elif cov_model == "GW1":       # GW with kappa = 1, mu \geq 1 + d/2 + kappa
        spatial_cor = (1 - dphi)**(eta[1]+1)*(dphi<1)*(1+dphi*(eta[1]+1))
    else: 
        print("Error: Covariance function not suppported!")
        return None
    
    if nug == False:
        cormat = spatial_cor
    elif nug == True:
        q = eta.size + 1
        c = eta[q-2]
        cormat = (1-c)*spatial_cor
    np.fill_diagonal(cormat, 1+1e-8) 
    return cormat


#######################################################################

"""
User specify the bound
"""

def cov_bound(cov_model, nug):
    if cov_model in ("Cub", "Exp", "Gau", "Sph", "Wave"): 
        bound = [(0.1, 1000)]
    elif cov_model in ("Cau", "Mat", ):
        bound = [(0.001, 15), (0.1, 5)]
    elif cov_model == "CauF":
        bound = [(0.001, 100), (0.001, 2), (0.001, 1000)]
    elif cov_model == "Edc":
        bound = [(0.001, 100), (1, 10)]
    elif cov_model == "PExp":
        bound = [(0.001, 100), (0.001, 2)]
    else: 
        print("Error: Covariance function not suppported!")
        return None
    if nug == True:
        bound.append((0, 0.9999))
    return bound