MLE.fit_GLS = function(y, X, D, kSigma, cov_model, theta.ini, nug, opt, lo_bound, up_bound){
  n = length(y)
  q = length(theta.ini) 

  if (opt == "LB"){
    prof.min = optim(theta.ini, profile.nll_GLS, method = "L-BFGS-B", lower = lo_bound, upper = up_bound, 
                     y = y, X = X, kSigma = kSigma, cov_model = cov_model, D = D, nug = nug)
  } else if (opt == "NM"){
    prof.min = optim(theta.ini, profile.nll_GLS, method = "Nelder-Mead",
                     y = y, X = X, kSigma = kSigma, cov_model = cov_model, D = D, nug = nug)
  } else if (opt == "CG"){
    prof.min = optim(theta.ini, profile.nll_GLS, method = "CG",
                     y = y, X = X, kSigma = kSigma, cov_model = cov_model, D = D, nug = nug)
  }

  # Return theta
  theta = prof.min$par
  suc = prof.min$convergence + 1
  nll = prof.min$value
  
  if (nug == TRUE){
    cratio = theta[q]/(theta[q-1]+theta[q])
    eta = c(theta[1:(q-2)], cratio)
  } else if (nug == FALSE){
    eta = theta[1:(q-2)]
  }
  
  rlist = list(theta = theta, suc, nll = nll, eta = eta)
  
  invisible(rlist)
}


 
profile.nll_GLS = function(theta, y, X, D, kSigma, cov_model, nug){
  n = length(y)
  p = dim(X)[2]
  q = length(theta) 
  
  if (nug == TRUE){
    sill = theta[q-1]+theta[q]
    cratio = theta[q]/sill
    eta = c(theta[1:(q-2)], cratio)
  } else if (nug == FALSE){
    eta = theta[1:(q-1)]
    sill = theta[q]
  }
  
  
  # Full covariance matrix:
  cormat_e = cor.mat(D, eta, cov_model, nug) 
  covmat = sill*cormat_e + kSigma
  
  # term 1
  L = t(chol(covmat))
  log_det_covmat = 2*sum(log(diag(L))) 
  cov_inv = solve(covmat)
  
  # term 2
  L1 = t(chol(t(X)%*%cov_inv%*%X))
  log_det_X_covmat_inv_X = 2*sum(log(diag(L1)))
  
  # term 3
  r = y - X%*%solve(t(X)%*%cov_inv%*%X)%*%t(X)%*%cov_inv%*%y
  part3 = t(r)%*%cov_inv%*%r
  
  # term 4
  log_det_XX = 2*sum(log(diag(t(chol(t(X)%*%X))))) 

  
  nll2 = (log_det_covmat + log_det_X_covmat_inv_X + part3 - log_det_XX + (n-p)*log(2*pi))
  return(nll2/2)
}



