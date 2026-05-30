MLE.fit_tr = function(kSigma, eet, D, cov_model, theta.ini, nug, opt, lo_bound, up_bound){
  # Preliminaries:
  n = dim(eet)[1]
  q = length(theta.ini)
  
  if (opt == "LB"){
    prof.min = optim(theta.ini, profile.nll_tr_mod, method = "L-BFGS-B", lower = lo_bound, upper = up_bound, 
                     eet = eet, kSigma = kSigma, cov_model = cov_model, D = D, nug = nug)
  } else if (opt == "NM"){
    prof.min = optim(theta.ini, profile.nll_tr_mod, method = "Nelder-Mead",
                     eet = eet, kSigma = kSigma, cov_model = cov_model, D = D, nug = nug)
  } else if (opt == "CG"){
    prof.min = optim(theta.ini, profile.nll_tr_mod, method = "CG",
                     eet = eet, kSigma = kSigma, cov_model = cov_model, D = D, nug = nug)
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



profile.nll_tr_mod = function(theta, kSigma, eet, D, cov_model, nug){
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
  
  cormat_e = cor.mat(D, eta, cov_model, nug) 
  covmat = sill*cormat_e + kSigma
  
  L = t(chol(covmat))
  log_det_covmat = 2*sum(log(diag(L)))
  
  # The log-likelihood function:
  nll2 = log_det_covmat + sum(diag(solve(covmat)%*%eet)) + n*log(2*pi) 
  
  return(nll2/n)
}

