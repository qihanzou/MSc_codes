#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 09:47:18 2020

@author: tingjin.chu@unimelb.edu.au
"""

import numpy as np
from numpy.linalg import inv
import numpy as np 
from scipy import linalg
from scipy.optimize import minimize
from scipy.optimize import differential_evolution, dual_annealing
from scipy.special import kv
from scipy.special import gamma
from Internal_Fun import *

""" Parameters and Reparameters
theta = (range, partial sill, other parameters, nugget effect), size q
eta = (range, other parameters, nugget ratio), size q-1
n: sample size
p: dimension of X (including intercept)
q: dimension of theta
"""


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


def MLE_fit(y, X, D, cov_model, eta_ini, nug, opt, lo_bound, up_bound):
    #bound = cov_bound(cov_model = cov_model, nug = nug)
    q = eta_ini.size + 1
    bound = []
    for i in range(q-1):
        bound.append((lo_bound[i], up_bound[i]))
    
    ## Optimization algorithms, L-BFGS-B
    if opt == "LB": 
        soln = minimize(profile_nll, eta_ini, bounds=bound, args=(y, X, D, cov_model, nug), method="L-BFGS-B")
    # Nelder-Mead
    elif opt == "NM":
        soln = minimize(profile_nll, eta_ini, bounds=bound, args=(y, X, D, cov_model, nug), method="Nelder-Mead")
    elif opt == "DE":
        soln = differential_evolution(profile_nll, bound, args = (y, X, D, cov_model, nug), seed = 2019, tol = 0.0001)
    else:
        soln = minimize(profile_nll, eta_ini, bounds=bound, args=(y, X, D, cov_model, nug), method = opt)
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
    return theta_est, suc, nll, beta_est, beta_var




##### MLE fitting without covariates X
def MLE_fitE(y, D, cov_model, eta_ini, nug, opt, lo_bound, up_bound):
    q = eta_ini.size + 1
    bound = []
    for i in range(q-1):
        bound.append((lo_bound[i], up_bound[i]))
    
    ## Optimization algorithms
    if opt == "LB": 
        soln = minimize(profile_nllE, eta_ini, bounds = bound, args=(y, D, cov_model, nug), method="L-BFGS-B")
    # Nelder-Mead
    elif opt == "NM":
        soln = minimize(profile_nllE, eta_ini, bounds = bound, args=(y, D, cov_model, nug), method="Nelder-Mead")        
    # Nelder-Mead, L-BFGS-B
    elif opt == "DE":
        soln = differential_evolution(profile_nllE, bound, args = (y, D, cov_model, nug), seed = 2019, tol = 0.0001)
    elif opt == "DA":
        soln = dual_annealing(profile_nllE, bound, args = (y, D, cov_model, nug))        
    else:
       soln = minimize(profile_nllE, eta_ini, bounds = bound, args=(y, D, cov_model, nug), method = opt)
    eta = soln.x
    suc = soln.success
    nll = soln.fun
    
    n = y.size
    cormat = cor_mat(D = D, eta = eta, cov_model = cov_model, nug = nug)
    L = np.linalg.cholesky(cormat)
    white_resids = np.linalg.solve(L, y)
    sill = white_resids.T @ white_resids/n

    ## reparametrization
    q = eta.size + 1
    theta_est = np.copy(eta)
    if nug == True: 
        nug_effect = sill*eta[q-2]
        psill = sill*(1-eta[q-2])
        theta_est = np.append(theta_est[0:(q-2)], [psill, nug_effect])
    elif nug == False:
        psill = sill
        theta_est = np.append(theta_est,  psill)
    return theta_est, suc, nll, eta, sill



#### GLS estimator of beta given eta 
def beta_gls(eta, y, X, D, cov_model, nug, sill):
    cormat = cor_mat(D = D, eta = eta, cov_model = cov_model, nug = nug) 
    L = np.linalg.cholesky(cormat)
    white_X = np.linalg.solve(L, X)
    white_y = np.linalg.solve(L, y)
    beta_est = linalg.solve(white_X.T @ white_X, white_X.T @ white_y, assume_a='sym')
    beta_var = linalg.inv(white_X.T @ white_X)*sill          # not real variance, need sill part
    return beta_est, beta_var


