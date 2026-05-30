# test use:
nall = 200
idx = 1:100
beta1 = 4


R_address = "C:/Users/qihan/Desktop/QH_rw"

setwd(R_address)
library(mvtnorm)               
source("SimuData.R")
source("BFuns.R")
source("SDCal_modified.R")
source("SDCal_RMLE.R")
source("BFuns_RBEGLS.R")
source("RBEGLS_functions_rw.R")
source("loop_functions.R")

ci_level = 0.95
learning_rate = 0.1
max_epoches = 1000

remove(update_est)
# -------------------------------------------------------------------------- #
locRdata=strwrap(paste("locs/loc", nall, ".Rdata",sep=""))
load(locRdata)

Dall = matrix(0, nrow = nall, ncol = nall)            
for(i in 1:nall) {
  Dall[i, ] <- sqrt((locall$x[i] - locall$x)^2 + (locall$y[i] - locall$y)^2)
}
#------------------------------------------------------------------------------#

alpha0 = c(1,0.5)
xi0 = c(3, 4)  # initial values of (phi, sigma^2)
cov_xi = "Mat32"
xi.ini = c(4)
up.bound = c(20)
lo.bound = c(0.01)

cov_ep = "Exp"
ep_theta0 = c(2,1) 
theta_ep.ini = c(3.5,1.2)
up.bounde = c(20, 20)
lo.bounde = c(0.01, 0.01)


beta0 = 2

nug = 0


#run.sim <- function(){
  
  
  etaall = SimuData(xi0, rate, locall, cov_xi, nug, 1)
  epall = SimuData(ep_theta0, rate, locall, cov_ep, nug, 1)
  
  k=1
  
  distall = sqrt((locall$x - 0)^2 + (locall$y - 0)^2)
  
  # -------------------------------------------------------------------------- #
  r1all = as.matrix(cbind(1, distall)) 
  r1s = as.matrix(cbind(1, distall[-idx]))
  r1 = as.matrix(cbind(1,distall[idx]))
  
  x1all = r1all%*%alpha0 + etaall$Ymat[k,]  
  x1 = x1all[idx] 
  x1s = x1all[-idx] 
  N = length(x1) 
  M = length(x1s)
  
  
  #  simple model
  ep_error = epall$Ymat[k,idx]        # error of s1,...,sN
  y = beta0 + beta1*x1 + ep_error     # create y values for s1,...,sN
  
  

  # -------------------------------------------------------------------------- #
  out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
  alpha_est = out1$beta
  cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
  w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 
  
  W = cbind(1, w)
  fols = lm(y~W-1)
  beta_ols = fols$coef
  varbeta_ols = diag(vcov(fols))
  
  pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
  CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)
  CovUB = CovU*beta_ols[2]^2 
  
  out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
  CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
  CovE_GLS = out2$theta[2]*CorE_GLS
  CovUBE_GLS = CovUB + CovE_GLS
  beta_RBEGLS = solve(t(W)%*%solve(CovUBE_GLS)%*%W)%*%t(W)%*%solve(CovUBE_GLS)%*%y
  
  
  # loop start here:
  res1 = RBEGLS_loop_mod2(50, 0.001, 0.001, beta_RBEGLS, out2$theta, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde)
  
  beta_RBEGLS = res1$beta_est
  theta_RBEGLS  = res1$theta_est
  
  
  CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
  CovE_GLS = theta_RBEGLS[2]*CorE_GLS
  CovUB = CovU*beta_RBEGLS[2]^2
  CovUBE_GLS = CovUB + CovE_GLS
  CovUBE_GLS_inv = solve(CovUBE_GLS)
  A1 = solve(t(W)%*%CovUBE_GLS_inv%*%W)
  varbeta_rb = diag(A1)
  
  #--------------------------------------------------------------------------------------------------------------
  

  # initialization
  Var_aplha_est = out1$beta_var
  dalpha0 = qnorm(ci_level)*sqrt(diag(Var_aplha_est)[1])
  dalpha1 = qnorm(ci_level)*sqrt(diag(Var_aplha_est)[2])
  dalpha = cbind(dalpha0, dalpha1)
  
  sdout = SDCal(Dall[-idx,-idx], out1$theta, cov_xi, nug = 0)
  dtheta0 = qnorm(ci_level)*sdout$se[1]
  dtheta1 = qnorm(ci_level)*sdout$se[2]
  dtheta = cbind(dtheta0, dtheta1)
  
 
  nll1 = out2$nll
  nll2 = res1$nll
  
  
  alpha0_est_i = alpha_est[1]
  alpha1_est_i = alpha_est[2]
  theta0_est_i = out1$theta[1]
  theta1_est_i = out1$theta[2]
  
  alpha0_max = alpha0_est_i + 1*dalpha0
  alpha1_max = alpha1_est_i + 1*dalpha1
  theta0_max = theta0_est_i + 1*dtheta0
  theta1_max = theta1_est_i + 1*dtheta1
  alpha0_min = alpha0_est_i - 1*dalpha0
  alpha1_min = alpha1_est_i - 1*dalpha1
  theta0_min = theta0_est_i - 1*dtheta0
  theta1_min = theta1_est_i - 1*dtheta1
  
  num_epoches = 0
  
  for (i_epoches in 1:max_epoches){
  
    if (alpha0_est_i + learning_rate*dalpha0 <= alpha0_max){
      alpha0_up = alpha0_est_i + learning_rate*dalpha0
    }else{
      alpha0_up = alpha0_max
    }
    
    if (alpha1_est_i + learning_rate*dalpha1 <= alpha1_max){
      alpha1_up = alpha1_est_i + learning_rate*dalpha1
    }else{
      alpha1_up = alpha1_max
    }
    
    if (theta0_est_i + learning_rate*dtheta0 <= theta0_max){
      theta0_up = theta0_est_i + learning_rate*dtheta0
    }else{
      theta0_up = theta0_max
    }
    
    if (theta1_est_i + learning_rate*dtheta1 <= theta1_max){
      theta1_up = theta1_est_i + learning_rate*dtheta1
    }else{
      theta1_up = theta1_max
    }
    
    
    if (alpha0_est_i - learning_rate*dalpha0 >= alpha0_min){
      alpha0_down = alpha0_est_i - learning_rate*dalpha0
    }else{
      alpha0_down = alpha0_min
    }
    
    if (alpha1_est_i - learning_rate*dalpha1 >= alpha1_min){
      alpha1_down = alpha1_est_i - learning_rate*dalpha1
    }else{
      alpha1_down = alpha1_min
    }
    
    if (theta0_est_i - learning_rate*dtheta0 >= theta0_min){
      theta0_down = theta0_est_i - learning_rate*dtheta0
    }else{
      theta0_down = theta0_min
    }
    
    if (theta1_est_i - learning_rate*dtheta1 >= theta1_min){
      theta1_down = theta1_est_i - learning_rate*dtheta1
    }else{
      theta1_down = theta1_min
    }
    
    
    
    
    w0000 = Comp_W(cbind(alpha0_down, alpha1_down), cbind(theta0_down, theta1_down))
    w0001 = Comp_W(cbind(alpha0_down, alpha1_down), cbind(theta0_down, theta1_up))
    w0010 = Comp_W(cbind(alpha0_down, alpha1_down), cbind(theta0_up, theta1_down))
    w0011 = Comp_W(cbind(alpha0_down, alpha1_down), cbind(theta0_up, theta1_up))
    w0100 = Comp_W(cbind(alpha0_down, alpha1_up), cbind(theta0_down, theta1_down))
    w0101 = Comp_W(cbind(alpha0_down, alpha1_up), cbind(theta0_down, theta1_up))
    w0110 = Comp_W(cbind(alpha0_down, alpha1_up), cbind(theta0_up, theta1_down))
    w0111 = Comp_W(cbind(alpha0_down, alpha1_up), cbind(theta0_up, theta1_up))
    w1000 = Comp_W(cbind(alpha0_up, alpha1_down), cbind(theta0_down, theta1_down))
    w1001 = Comp_W(cbind(alpha0_up, alpha1_down), cbind(theta0_down, theta1_up))
    w1010 = Comp_W(cbind(alpha0_up, alpha1_down), cbind(theta0_up, theta1_down))
    w1011 = Comp_W(cbind(alpha0_up, alpha1_down), cbind(theta0_up, theta1_up))
    w1100 = Comp_W(cbind(alpha0_up, alpha1_up), cbind(theta0_down, theta1_down))
    w1101 = Comp_W(cbind(alpha0_up, alpha1_up), cbind(theta0_down, theta1_up))
    w1110 = Comp_W(cbind(alpha0_up, alpha1_up), cbind(theta0_up, theta1_down))
    w1111 = Comp_W(cbind(alpha0_up, alpha1_up), cbind(theta0_up, theta1_up))
    
    res0000 = loop_funs(w0000, cbind(alpha0_down, alpha1_down), cbind(theta0_down, theta1_down))
    res0001 = loop_funs(w0001, cbind(alpha0_down, alpha1_down), cbind(theta0_down, theta1_up))
    res0010 = loop_funs(w0010, cbind(alpha0_down, alpha1_down), cbind(theta0_up, theta1_down))
    res0011 = loop_funs(w0011, cbind(alpha0_down, alpha1_down), cbind(theta0_up, theta1_up))
    res0100 = loop_funs(w0100, cbind(alpha0_down, alpha1_up), cbind(theta0_down, theta1_down))
    res0101 = loop_funs(w0101, cbind(alpha0_down, alpha1_up), cbind(theta0_down, theta1_up))
    res0110 = loop_funs(w0110, cbind(alpha0_down, alpha1_up), cbind(theta0_up, theta1_down))
    res0111 = loop_funs(w0111, cbind(alpha0_down, alpha1_up), cbind(theta0_up, theta1_up))
    res1000 = loop_funs(w1000, cbind(alpha0_up, alpha1_down), cbind(theta0_down, theta1_down))
    res1001 = loop_funs(w1001, cbind(alpha0_up, alpha1_down), cbind(theta0_down, theta1_up))
    res1010 = loop_funs(w1010, cbind(alpha0_up, alpha1_down), cbind(theta0_up, theta1_down))
    res1011 = loop_funs(w1011, cbind(alpha0_up, alpha1_down), cbind(theta0_up, theta1_up))
    res1100 = loop_funs(w1100, cbind(alpha0_up, alpha1_up), cbind(theta0_down, theta1_down))
    res1101 = loop_funs(w1101, cbind(alpha0_up, alpha1_up), cbind(theta0_down, theta1_up))
    res1110 = loop_funs(w1110, cbind(alpha0_up, alpha1_up), cbind(theta0_up, theta1_down))
    res1111 = loop_funs(w1111, cbind(alpha0_up, alpha1_up), cbind(theta0_up, theta1_up))
    
    res_diff_w = list(res0000, res0001, res0010, res0011, res0100, res0101, res0110, res0111, res1000, res1001, res1010, res1011, res1100, res1101, res1110,res1111)
    res_nll_w = cbind(res0000$nll, res0001$nll, res0010$nll, res0011$nll, res0100$nll, res0101$nll, res0110$nll, res0111$nll, res1000$nll, res1001$nll, res1010$nll, res1011$nll, res1100$nll, res1101$nll, res1110$nll, res1111$nll)
    index_min = which.min(res_nll_w)
    
    if (res_nll_w[index_min] < nll2){
      nll2 = res_nll_w[index_min]
      update_est = res_diff_w[index_min]
      
      alpha0_est_i = update_est[[1]]$alpha[1]
      alpha1_est_i = update_est[[1]]$alpha[2]
      theta0_est_i = update_est[[1]]$xi[1]
      theta1_est_i = update_est[[1]]$xi[2]
      
      num_epoches = num_epoches + 1
    } else{
      break
      num_epoches = num_epoches
    }
  }
  
  
  if (as.numeric(exists("update_est")) == 1  ){
    out_rw = list((as.matrix(as.numeric(update_est[[1]]$alpha))), as.matrix(as.numeric(update_est[[1]]$xi)), update_est[[1]]$beta, as.matrix(update_est[[1]]$theta), nll2)
  
  } else {
    out_rw = list(as.matrix(alpha_est), (as.matrix(out1$theta)), as.matrix(beta_RBEGLS), as.matrix(res1$theta), res1$nll)
  }
  
  
  cov_xiall_rw = out_rw[[2]][2]*cor.mat(Dall, out_rw[[2]][1], cov_xi, nug = 0)
  w_rw = r1%*%out_rw[[1]] + cov_xiall_rw[idx, -idx]%*%solve(cov_xiall_rw[-idx,-idx])%*%(x1s - r1s%*%out_rw[[1]]) 
  W_rw = cbind(1, w_rw)
  
  
  pc_rw = r1 - cov_xiall_rw[idx, -idx]%*%solve(cov_xiall_rw[-idx,-idx]) %*% r1s
  CovU_rw = cov_xiall_rw[idx, idx] - cov_xiall_rw[idx, -idx]%*%solve(cov_xiall_rw[-idx,-idx]) %*% t(cov_xiall_rw[idx, -idx]) + (pc_rw)%*%solve(t(r1s)%*%solve(cov_xiall_rw[-idx,-idx])%*%r1s)%*%t(pc_rw)
  
  
  
  CorE_GLS_rw = cor.mat(Dall[idx,idx], out_rw[[4]][1], cov_ep, nug = 0)
  CovE_GLS_rw = out_rw[[4]][2]*CorE_GLS_rw
  
  CovUB_rw = CovU_rw*out_rw[[3]][2]^2
  CovUBE_GLS_rw = CovUB_rw + CovE_GLS_rw
  CovUBE_GLS_rw_inv = solve(CovUBE_GLS_rw)
  A1_rw = solve(t(W_rw)%*%CovUBE_GLS_rw_inv%*%W_rw)
  varbeta_rw = diag(A1_rw)
  
  
  


  list(alpha_est = out1$beta,
       xi_est = as.matrix(out1$theta),
       alpha_rw = out_rw[[1]],
       xi_rw = out_rw[[2]],
       
       beta_ols = as.matrix(as.numeric(beta_ols)), 
       beta_RBEGLS = beta_RBEGLS,  
       theta_est = as.matrix(theta_RBEGLS),
       beta_rw = out_rw[[3]],
       theta_rw = out_rw[[4]],
       
       varbeta_ols = as.matrix(as.numeric(varbeta_ols)), 
       varbeta_rb = as.matrix(varbeta_rb),
       varbeta_rw = as.matrix(varbeta_rw),
       
       num_epoches = num_epoches,
       i_epoches = i_epoches,
       nll_rb = out2$nll,
       nll_rb_full = res1$nll,
       nll_rw = out_rw[[5]],
       rw_update_exist = as.numeric(exists("update_est")) == 1  )
  
  
#}

  