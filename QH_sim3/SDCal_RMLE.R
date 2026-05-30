## Calculate the cov matrix (cov.gau) and standard error of each estiamtor (s.gau)
SDCal_RMLE <- function (D,theta, cov.model,nug, CovUBboot, W){
  
  if (cov.model != "PExp"){
    Jout = Jmat_RMLE(theta, cov.model=cov.model, D = D, nug=nug, CovUBboot=CovUBboot, W=W)
  } else if (cov.model == "PExp"){
    #Jout = JmatPExp(theta, cov.model=cov.model, D = D)
    print("PExp has not been updated")
  }
  
  cov.gau = solve(Jout$J.theta)
  s.gau = diag(cov.gau)^.5     ## version under Gaussian Assumption
  object <- list (se = s.gau, covmat = cov.gau)
  invisible(object)
}




## Derivatives of covariance matrices
Jmat_RMLE = function (theta, cov.model, D, nug, CovUBboot, W){
  # The covariance matrix
  n = dim(D)[1]
  if (cov.model == "Exp"){
    psi = theta[2] *exp( - D/theta[1])      
  } else if (cov.model=="Gau"){
    psi = theta[2] * exp(-D^2/theta[1]^2)
  } else if (cov.model=="Sph"){
    psi = theta[2] *(1-1.5*D/theta[1]+0.5*(D/theta[1])^3)*(D < theta[1])
  } else if (cov.model=="Cau"){
    psi = theta[2]*(1+(D/theta[1])^2)^{-pkap}  
  } else if (cov.model=="Cau2"){
    psi = theta[2]*(kcau+(D/theta[1])^2)^{-pkap}  
  } else if (cov.model == "Mat32"){
    psi = theta[2]*(1 + 3^.5*D/theta[1])*exp(-3^0.5*D/theta[1])
  } else {print("The input covariance function is not supported")}
  if (nug == 1){
    diag(psi) = theta[2]+theta[3]
  } else{
    diag(psi) = theta[2]
  }
  
  CovE = psi #???
  
  Sigma_boot = CovE + CovUBboot
  Sigma_boot_inv = solve(Sigma_boot)
    
  Pi = Sigma_boot_inv - Sigma_boot_inv%*%W%*%solve(t(W)%*%Sigma_boot_inv%*%W)%*%t(W)%*%Sigma_boot_inv
  psi.inv = Pi
  
  
  # Matrix Derivative with respect to theta
  if (cov.model == "Exp"){
    Gamma1 = theta[2]*exp(-D/theta[1])*D/theta[1]^2
    Gamma2 = exp(-D/theta[1])      
  } else if (cov.model=="Gau"){
    Gamma1 = 2*theta[2]*exp(-D^2/theta[1]^2)*D^2/theta[1]^3
    Gamma2 = exp(-D^2/theta[1]^2)    
  } else if (cov.model=="Sph"){
    Gamma1 = 1.5*theta[2]*(D/theta[1]^2-D^3/theta[1]^4)*(D < theta[1])
    Gamma2 = (1-1.5*D/theta[1]+0.5*(D/theta[1])^3)*(D < theta[1])
  } else if (cov.model=="Cau"){
    Gamma1 = 2*theta[2]*pkap*D^2*(1+(D/theta[1])^2)^{-pkap-1}/theta[1]^3
    Gamma2 = (1+(D/theta[1])^2)^{-pkap} 
  } else if (cov.model=="Cau2"){
    Gamma1 = 2*theta[2]*pkap*D^2*(kcau+(D/theta[1])^2)^{-pkap-1}/theta[1]^3
    Gamma2 = (kcau+(D/theta[1])^2)^{-pkap} 
  } else if (cov.model == "Mat32"){
    Gamma1 = ((3*(D^2)*theta[2])/(theta[1]^3))*exp(-sqrt(3)*D/theta[1])
    Gamma2 = exp(-sqrt(3)*D/theta[1]) + (sqrt(3)*D/theta[1])*exp(-sqrt(3)*D/theta[1])
  } else {print("The input covariance function is not supported")}
  Gamma3= diag(n)
  
  if (nug == 1){
    J.theta=matrix(0,3,3)
    J.theta[1,1]=sum(diag(psi.inv%*%Gamma1%*%psi.inv%*%Gamma1))/2
    J.theta[1,2]=J.theta[2,1]=sum(diag(psi.inv%*%Gamma1%*%psi.inv%*%Gamma2))/2
    J.theta[1,3]=J.theta[3,1]=sum(diag(psi.inv%*%Gamma1%*%psi.inv%*%Gamma3))/2
    J.theta[2,2]=sum(diag(psi.inv%*%Gamma2%*%psi.inv%*%Gamma2))/2
    J.theta[2,3]=J.theta[3,2]=sum(diag(psi.inv%*%Gamma2%*%psi.inv%*%Gamma3))/2
    J.theta[3,3]=sum(diag(psi.inv%*%Gamma3%*%psi.inv%*%Gamma3))/2
  } else{
    J.theta=matrix(0,2,2)
    J.theta[1,1]=sum(diag(psi.inv%*%Gamma1%*%psi.inv%*%Gamma1))/2
    J.theta[1,2]=J.theta[2,1]=sum(diag(psi.inv%*%Gamma1%*%psi.inv%*%Gamma2))/2
    J.theta[2,2]=sum(diag(psi.inv%*%Gamma2%*%psi.inv%*%Gamma2))/2
  }
  
  
  
  object = list(J.theta = J.theta, Pi = Pi, 
                psi.inv=psi.inv, Gam1=Gamma1, Gam2=Gamma2, Gam3=Gamma3)
  invisible(object)
}


