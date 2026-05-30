
RBEGLS_loop = function(iter, btor1, btor2, ttor1, ttor2, betaest, theta_ep.ini, CovU, y, W, Dall, cov_ep, lo.bounde,up.bounde){
  
  for (i in 1:iter){
    CovUB = CovU*betaest[2]^2
    out2 = MLE.fit_GLS(y, W, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
    CorE = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
    CovE = out2$theta[2]*CorE
    CovUBE = CovUB + CovE
    beta_update = solve(t(W)%*%solve(CovUBE)%*%W)%*%t(W)%*%solve(CovUBE)%*%y
    
    b_diff1 = beta_update[1] - betaest[1]
    b_diff2 = beta_update[2] - betaest[2]
    t_diff1 = out2$theta[1] - theta_ep.ini[1]
    t_diff2 = out2$theta[2] - theta_ep.ini[2]
    
    if (b_diff1 < btor1){
      if (b_diff2 < btor2){
         if (t_diff1 < ttor1){
           if (b_diff2 < ttor2){
               return (c(beta_update, out2$theta, i, round(b_diff1,digits = 3), round(b_diff2,digits = 3), round(t_diff1,digits = 3), round(t_diff2,digits = 3)))
               break
           }
         }
      }
    } 
    betaest = beta_update
    theta_ep.ini = out2$theta
  }
  
  return (c(beta_update, out2$theta, i, b_diff1, b_diff2, t_diff1, t_diff2))
}