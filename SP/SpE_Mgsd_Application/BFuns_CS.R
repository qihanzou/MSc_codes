MLE.fit_CS = function(covUBE, D, zeta, y, S, tau.ini, opt, lo_bound, up_bound){

  if (opt == "LB"){
    prof.min = optim(tau.ini, profile.nll_c, method = "L-BFGS-B", lower = lo_bound, upper = up_bound, 
                     zeta = zeta, y = y, S = S, covUBE = covUBE, D=D)
  } else if (opt == "NM"){
    prof.min = optim(tau.ini, profile.nll_c, method = "Nelder-Mead",
                     zeta = zeta, y = y, S = S, covUBE = covUBE, D=D)
  } else if (opt == "CG"){
    prof.min = optim(tau.ini, profile.nll_c, method = "CG",
                     zeta = zeta, y = y, S = S, covUBE = covUBE, D=D)
  }
  
  # Return tau
  tau = prof.min$par
  suc = prof.min$convergence + 1
  nll = prof.min$value
  
  rlist = list(tau = tau, suc, nll = nll)
  
  invisible(rlist)
}



profile.nll_c = function(tau, zeta, y, S, covUBE, D){
  n = length(y)
  p = dim(S)[2]
  
  # Full covariance matrix: 
  covmat = covUBE + tau#*diag(dim(covUBE)[1])
  
  # term 1
  L = t(chol(covmat))
  log_det_covmat = 2*sum(log(diag(L))) 
  
  # term 2
  part2 = t(y-S%*%zeta)%*%solve(covmat)%*%(y-S%*%zeta)
 
  
  nll2 = (n*log(2*pi)/2 + log_det_covmat/2 + part2/2)
  return(nll2)
}



