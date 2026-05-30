# This function corresponds to MLE_fit and MLE_fitE in Python

MLE.fit_tr = function(eet, D, cov_model, theta.ini, nug, opt, lo_bound, up_bound){
  # Preliminaries:
  n = dim(eet)[1]
  q = length(theta.ini)
  # minimze negative concentrated log-likelihood function
  
  if (opt == "LB"){
    prof.min = optim(theta.ini, profile.nll_tr, method = "L-BFGS-B", lower = lo_bound, upper = up_bound, 
                     eet = eet, cov_model = cov_model, D = D, nug = nug)
  } else if (opt == "NM"){
    prof.min = optim(theta.ini, profile.nll_tr, method = "Nelder-Mead",
                     eet = eet, cov_model = cov_model, D = D, nug = nug)
  } else if (opt == "CG"){
    prof.min = optim(theta.ini, profile.nll_tr, method = "CG",
                     eet = eet, cov_model = cov_model, D = D, nug = nug)
  }
  
  # Return theta
  theta = prof.min$par
  suc = prof.min$convergence + 1
  nll = prof.min$value
  if (nug == TRUE){
    cratio = theta[q]/(theta[q-1]+theta[q])
    eta = c(theta[1:(q-2)], cratio)
    #nugget_effect = 
    #psill = sill - nugget_effect
    #theta_est = c(eta[1:(q-2)], psill, nugget_effect)
  } else if (nug == FALSE){
    eta = theta[1:(q-2)]
  }
  
  # suc = 1 means convergence of optimization
  rlist = list(theta = theta, suc, nll = nll, eta = eta)
  
  invisible(rlist)
}


# This function correpsonds to profile_nll and profile_nllE in Python. 
profile.nll_tr = function(theta, eet, D, cov_model, nug){
  n = dim(eet)[1]
  q = length(theta)
  
  if (nug == TRUE){
    sill = theta[q-1]+theta[q]
    cratio = theta[q]/sill
    eta = c(theta[1:(q-2)], cratio)
  } else if (nug == FALSE){
    eta = theta[1:(q-1)]
    sill = theta[q]
  }
  
  cormat = cor.mat(D, eta, cov_model, nug)
  covmat = cormat*sill
  
  L = t(chol(covmat))
  log_det_covmat = 2*sum(log(diag(L)))
  
  nll2 = log_det_covmat + sum(diag(solve(covmat)%*%eet)) + n*log(2*pi)
  
  #return(n*log(mean(white_resids^2)))
  return(nll2/n)
}


