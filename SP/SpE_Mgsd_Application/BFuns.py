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



def MLE_fit(y, X, D, cov_model):
    opt = "LB"
    nug = False
    eta_ini = np.array([0.1])
    lo_bound = np.array([1e-10])
    up_bound = np.array([2,2])
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
    ##
    
    n = y.size # the number of data in total
    cormat = cor_mat(D = D, eta = eta, cov_model = cov_model, nug = nug) # the spatial cor, Cz(sigma^2, phi)
    L = np.linalg.cholesky(cormat) # return a matrix. A =LL* (^T for real)
    # Cholesky decomposition solve Ax=b, 
    # A is both Hermitian/symmetric and positive-definite.
    # First, we solve for y in Ly=b.
    # then solve for x in L*x=y for x.
    white_X = np.linalg.solve(L, X) # Coefficient matrix L. solve(a,b) == solve ax = b return "x"
    # Lx = X
    white_y = np.linalg.solve(L, y)
    
    beta_est = linalg.solve(white_X.T @ white_X, white_X.T @ white_y, assume_a='sym') # equ(7) in geoII
    # scipy.linalg.solve
    # scipy.linalg.solve(a, b, lower=False, overwrite_a=False, overwrite_b=False, check_finite=True, assume_a='gen', transposed=False)
    # Solves the linear equation set a @ x == b for the unknown x for square a matrix.
    # assume_a: generic matrix, gen, default. symmetric, sym. hermitian, her. positive definite, pos.
    white_resids = white_y - white_X @ beta_est # calculate the residual, Res = y - XB = y - mu
    
    # Why we can calculate sill as below???
    sill = white_resids.T @ white_resids/n 
    beta_var = linalg.inv(white_X.T @ white_X)*sill  # Why need to times sill here??? equ(8)

    ## reparametrization
    theta_est = np.copy(eta) # copy
    if nug == True: 
        #theta_est [ntheta-1] = theta[ntheta-1]*sill    # nugget effect
        nug_effect = sill*eta[q-2]
        psill = sill*(1-eta[q-2])
        theta_est = np.append(theta_est[0:(q-2)], [psill, nug_effect])
    elif nug == False:
        psill = sill
        theta_est = np.append(theta_est,  psill) #np.append, make result together
    return theta_est, suc, nll, beta_est, beta_var


def MLE_fit_gas(y, X, D, cov_model):
    opt = "LB"
    nug = True
    eta_ini = np.array([0.5, 0.001])
    lo_bound = np.array([1e-5, 1e-10])
    up_bound = np.array([2,2])
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
    ##
    
    n = y.size # the number of data in total
    cormat = cor_mat(D = D, eta = eta, cov_model = cov_model, nug = nug) # the spatial cor, Cz(sigma^2, phi)
    L = np.linalg.cholesky(cormat) # return a matrix. A =LL* (^T for real)
    # Cholesky decomposition solve Ax=b, 
    # A is both Hermitian/symmetric and positive-definite.
    # First, we solve for y in Ly=b.
    # then solve for x in L*x=y for x.
    white_X = np.linalg.solve(L, X) # Coefficient matrix L. solve(a,b) == solve ax = b return "x"
    # Lx = X
    white_y = np.linalg.solve(L, y)
    
    beta_est = linalg.solve(white_X.T @ white_X, white_X.T @ white_y, assume_a='sym') # equ(7) in geoII
    # scipy.linalg.solve
    # scipy.linalg.solve(a, b, lower=False, overwrite_a=False, overwrite_b=False, check_finite=True, assume_a='gen', transposed=False)
    # Solves the linear equation set a @ x == b for the unknown x for square a matrix.
    # assume_a: generic matrix, gen, default. symmetric, sym. hermitian, her. positive definite, pos.
    white_resids = white_y - white_X @ beta_est # calculate the residual, Res = y - XB = y - mu
    
    # Why we can calculate sill as below???
    sill = white_resids.T @ white_resids/n 
    beta_var = linalg.inv(white_X.T @ white_X)*sill  

    ## reparametrization
    theta_est = np.copy(eta) # copy
    if nug == True: 
        #theta_est [ntheta-1] = theta[ntheta-1]*sill    # nugget effect
        nug_effect = sill*eta[q-2]
        psill = sill*(1-eta[q-2])
        theta_est = np.append(theta_est[0:(q-2)], [psill, nug_effect])
    elif nug == False:
        psill = sill
        theta_est = np.append(theta_est,  psill) 
    return theta_est, suc, nll, beta_est, beta_var



def MLE_fit_GLS(y, X, D, kSigma, cov_model):
    opt = "LB"
    nug = False
    theta_ini = np.array([0.5, 0.001])
    lo_bound = np.array([1e-5, 1e-10])
    up_bound = np.array([2,2])
    q = theta_ini.size
    bound = []
    for i in range(q):
        bound.append((lo_bound[i], up_bound[i]))
        
    if opt == "LB": 
        soln = minimize(profile_nll_GLS, theta_ini, bounds=bound, args=(y, X, D, kSigma, cov_model, nug), method="L-BFGS-B")
    # Nelder-Mead
    elif opt == "NM":
        soln = minimize(profile_nll_GLS, theta_ini, bounds=bound, args=(y, X, D, kSigma, cov_model, nug), method="Nelder-Mead")
    elif opt == "DE":
        soln = differential_evolution(profile_nll_GLS, bound, args = (y, X, D, kSigma, cov_model, nug), seed = 2019, tol = 0.0001)
    else:
        soln = minimize(profile_nll_GLS, theta_ini, bounds=bound, args=(y, X, D, kSigma, cov_model, nug), method = opt)
    theta = soln.x
    suc = soln.success
    nll = soln.fun
    
    return theta, suc, nll


def MLE_fit_CS(covUBE, D, zeta, y, S):
    opt = "LB"
    tau_ini = np.array([0.0015])
    lo_bound = np.array([0.00000001])
    up_bound = np.array([5])
    q = tau_ini.size
    bound = []
    for i in range(q):
        bound.append((lo_bound[i], up_bound[i]))
        
    if opt == "LB": 
        soln = minimize(profile_nll_CS, tau_ini, bounds=bound, args=(covUBE, D, zeta, y, S), method="L-BFGS-B")
   
    tau = soln.x
    suc = soln.success
    nll = soln.fun
    
    return tau, suc, nll


def MLE_fit_CS_gau(covUBE, D, zeta, y, S):
    opt = "LB"
    tau_ini = np.array([100, 0.0015])
    lo_bound = np.array([0.0001, 0.00001])
    up_bound = np.array([2000, 100])
    q = tau_ini.size
    bound = []
    for i in range(q):
        bound.append((lo_bound[i], up_bound[i]))
        
    if opt == "LB": 
        #soln = minimize(profile_nll_c_gau, tau_ini, bounds=bound, args=(covUBE, D, zeta, y, S), method="L-BFGS-B")
        soln = differential_evolution(profile_nll_c_gau, bound, args = (covUBE, D, zeta, y, S), seed = 2024, tol = 0.0001)
    tau = soln.x
    suc = soln.success
    nll = soln.fun
    
    return tau, suc, nll






def MLE_fit_GLS_nug(y, X, D, kSigma, cov_model):
    opt = "LB"
    nug = True
    theta_ini = np.array([0.5, 0.001, 0.0001])
    lo_bound = np.array([1e-5, 1e-10, 1e-20])
    up_bound = np.array([2, 2, 1])
    q = theta_ini.size
    bound = []
    for i in range(q):
        bound.append((lo_bound[i], up_bound[i]))
        
    if opt == "LB": 
        soln = minimize(profile_nll_GLS_nug, theta_ini, bounds=bound, args=(y, X, D, kSigma, cov_model, nug), method="L-BFGS-B")
    # Nelder-Mead
    elif opt == "NM":
        soln = minimize(profile_nll_GLS_nug, theta_ini, bounds=bound, args=(y, X, D, kSigma, cov_model, nug), method="Nelder-Mead")
    elif opt == "DE":
        soln = differential_evolution(profile_nll_GLS_nug, bound, args = (y, X, D, kSigma, cov_model, nug), seed = 2019, tol = 0.0001)
    else:
        soln = minimize(profile_nll_GLS_nug, theta_ini, bounds=bound, args=(y, X, D, kSigma, cov_model, nug), method = opt)
    theta = soln.x
    suc = soln.success
    nll = soln.fun
    
    return theta, suc, nll




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


