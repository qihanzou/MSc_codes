est_rho = function(N, beta_est, TP, Var, B){
  rho_set = seq(0, 1, length = N)
  CP90set = matrix(0, N)
  CP95set = matrix(0, N)
  D = matrix(0, N)
  for (i in 1:N){
    rho = rho_set[i]
    TSD = (rho*Var + (Var/B)*(1 - rho))^.5
    Betas = matrix(beta_est,  dim(TP)[1], 1, byrow = T)
    UP = TP + qnorm(0.95)*TSD
    LO = TP - qnorm(0.95)*TSD
    CPT =  (UP > Betas & LO < Betas)
    CP90 = mean(CPT)
    UP = TP + qnorm(0.975)*TSD
    LO = TP - qnorm(0.975)*TSD
    CPT =  (UP > Betas & LO < Betas)
    CP95 = mean(CPT)
    CP90set[i] = CP90
    CP95set[i] = CP95
    dist = sqrt((CP90 - 0.9)^2 + (CP95 - 0.95)^2)
    D[i] = dist
  }
  idx1 = which.min(D)
  
  # idx = which(D == min(D))
  # r_set = matrix(0, length(idx))
  # for (j in 1:length(idx)){
  #   C90n = CP90set[idx[j]]
  #   C95n = CP95set[idx[j]]
  #   if (C90n <= 90){
  #     if (C95n <= 95){
  #       r_set[j] = rho_set[idx[j]]
  #     }
  #   }
  # }
  # 
  # if (max(r_set) > 0){
  #   return(min(r_set))
  # }else{
  #   return(rho_set[idx1])
  # }
  return(rho_set[idx1])
}




