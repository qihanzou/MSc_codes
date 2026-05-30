# # # # test use:
rm(list = ls(all = TRUE)) 
nall = 200
idx = 1:100
beta1 = 4
MC_iter = 100
R_address = "C:/Users/qihan/Desktop/LT"


setwd(R_address)
library(mvtnorm)  
library(moments)
library(revss)
library(asbio)
library(MASS)
library(r2spss)
source("SimuData.R")
source("BFuns.R")
source("BFuns_RBEGLS.R")
source("SDCal_modified.R")
source("RBEGLS_functions_final.R")
source("Matrix_FUns.R")
library(numDeriv)


# -------------------------------------------------------------------------- #
locRdata=strwrap(paste("locs/loc", nall, ".Rdata",sep=""))
load(locRdata)

Dall = matrix(0, nrow = nall, ncol = nall)            
for(i in 1:nall) {
  Dall[i, ] <- sqrt((locall$x[i] - locall$x)^2 + (locall$y[i] - locall$y)^2)
}

plot(locall$x[idx], locall$y[idx], col = "red")
#points(locall$x[-idx], locall$y[-idx], col = "blue")



assign_and_plot_regions <- function(x, y, unit_size = 1, xlim = NULL, ylim = NULL) {
  # Assign region indices
  x_region <- floor(x / unit_size)
  y_region <- floor(y / unit_size)
  region_id <- paste0("(", x_region, ",", y_region, ")")
  region_center_x <- x_region * unit_size + unit_size / 2
  region_center_y <- y_region * unit_size + unit_size / 2
  region_info <- data.frame(x = x, y = y, x_region, y_region, region_id, region_center_x, region_center_y)
  region_ids <- unique(region_info$region_id)
  region_colors <- rainbow(length(region_ids))
  names(region_colors) <- region_ids
  plot(NA, NA, xlim = xlim, ylim = ylim, xlab = "x", ylab = "y", main = paste0("Data Points Assigned to ", unit_size, " x ", unit_size, " Regions"))
  abline(v = seq(floor(xlim[1] / unit_size) * unit_size, ceiling(xlim[2] / unit_size) * unit_size, by = unit_size), col = "grey90")
  abline(h = seq(floor(ylim[1] / unit_size) * unit_size, ceiling(ylim[2] / unit_size) * unit_size, by = unit_size), col = "grey90")
  points(region_info$x, region_info$y, col = region_colors[region_info$region_id], pch = 16)
  points(region_info$region_center_x, region_info$region_center_y, pch = 4, col = "black", lwd = 2, cex = 1.2)
  return(invisible(region_info))
}
region_info <- assign_and_plot_regions(locall$x[idx], locall$y[idx], unit_size = 3, xlim = c(-6, 6), ylim = c(-6, 6))
#------------------------------------------------------------------------------#

alpha0 = c(1,0.5)
xi0 = c(3, 4)  # initial values of (phi, sigma^2) (3, 4)
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
ep_error = epall$Ymat[k,idx]         
y = beta0 + beta1*x1 + ep_error 

region_info <- data.frame(region_info, y_pval = y)
region_mean <- aggregate(cbind(x, y) ~ region_id, data = region_info, FUN = mean)
region_count <- aggregate(x ~ region_id, data = region_info, FUN = length)
names(region_count)[2] <- "n_points"
region_y_pval <- aggregate(y_pval ~ region_id, data = region_info, FUN = mean)
names(region_y_pval)[2] <- "y_bval"
region_centers <- unique(region_info[, c("region_id", "region_center_x", "region_center_y")])
region_summary <- merge(region_mean, region_count, by = "region_id")
region_summary <- merge(region_summary, region_y_pval, by = "region_id")
region_summary <- merge(region_summary, region_centers, by = "region_id")
region_summary





# Assume region_summary contains region_center_x, region_center_y, unit_size known
unit_size <- 3  # adjust as needed
spacing <- unit_size / 2  # spacing between grid points

# Relative positions (3x3 grid centered at 0,0)
offsets <- expand.grid(dx = c(-spacing, 0, spacing), dy = c(-spacing, 0, spacing))

# Generate expanded points for each region
expanded_points <- do.call(rbind, lapply(1:nrow(region_summary), function(i) {
  center_x <- region_summary$region_center_x[i]
  center_y <- region_summary$region_center_y[i]
  region_id <- region_summary$region_id[i]
  
  data.frame(
    x = center_x + offsets$dx,
    y = center_y + offsets$dy,
    region_id = region_id
  )
}))
plot(expanded_points$x, expanded_points$y, pch = 16, col = "black",
     xlab = "x", ylab = "y", main = "Evenly Distributed 9 Points per Region")




generate_region_grid_points <- function(region_summary, unit_size = 1, n_per_side = 3) {
  if (n_per_side < 1 || n_per_side %% 1 != 0) {
    stop("n_per_side must be a positive integer >= 1")
  }
  
  # Generate evenly spaced offsets centered at 0
  offset_vals <- seq(-unit_size / 2, unit_size / 2, length.out = n_per_side)
  offset_grid <- expand.grid(dx = offset_vals, dy = offset_vals)
  
  # Generate grid points for each region
  grid_points <- do.call(rbind, lapply(1:nrow(region_summary), function(i) {
    cx <- region_summary$region_center_x[i]
    cy <- region_summary$region_center_y[i]
    region_id <- region_summary$region_id[i]
    
    data.frame(
      x = cx + offset_grid$dx,
      y = cy + offset_grid$dy,
      region_id = region_id
    )
  }))
  
  return(grid_points)
}

# Suppose region_summary already exists and includes region_id, region_center_x, region_center_y
grid_points <- generate_region_grid_points(region_summary, unit_size = 3, n_per_side = 3)

# Plot
plot(grid_points$x, grid_points$y, pch = 16, col = "black", xlim = c(-6,6), ylim=c(-6,6),
     main = "Even Grid Points per Region", xlab = "x", ylab = "y")



#--------------------------------------------------------------------------------------------------------------
# w
out1 = MLE.fit(x1s, r1s, Dall[-idx,-idx], cov_xi, xi.ini, nug, "LB", lo.bound,up.bound)
alpha_est = out1$beta
cov_xiall = out1$theta[2]*cor.mat(Dall, out1$eta, cov_xi, nug = 0)
w = r1%*%alpha_est + cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx])%*%(x1s - r1s%*%alpha_est) 
pc = r1 - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% r1s
CovU = cov_xiall[idx, idx] - cov_xiall[idx, -idx]%*%solve(cov_xiall[-idx,-idx]) %*% t(cov_xiall[idx, -idx]) + (pc)%*%solve(t(r1s)%*%solve(cov_xiall[-idx,-idx])%*%r1s)%*%t(pc)

gdist_mc = DT_MC(w, CovU, g, 1000) # if use MCT method
CovU_nl_mc = gdist_mc$cov_z
z_mc = gdist_mc$z

gdist = DT_Delta_1st_jacobianfun(w, CovU, g) # 1st order Delta method
CovU_nl = gdist$cov_z
z = gdist$z

# KR
Z = cbind(1, z)
kr = MLE.fit(y, Z, Dall[idx,idx], cov_ep, c(2.5), nug, "LB", lo.bound,up.bound)
theta_kr = kr[[1]] 
beta_kr  = kr[[4]] 
varbeta_kr = diag(kr[[5]])
#--------------------------------------------------------------------------------------------------------------
# OLS dm
fols = lm(y~Z-1)
beta_ols = fols$coef
varbeta_ols = diag(vcov(fols))

# RBEGLS dm
CovUB = CovU_nl*beta_ols[2]^2 
out2 = MLE.fit_GLS(y, Z, Dall[idx,idx], CovUB, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
CorE_GLS = cor.mat(Dall[idx,idx], out2$eta, cov_ep, nug = 0)
CovE_GLS = out2$theta[2]*CorE_GLS
CovUBE_GLS = CovUB + CovE_GLS
beta_RBEGLS = solve(t(Z)%*%solve(CovUBE_GLS)%*%Z)%*%t(Z)%*%solve(CovUBE_GLS)%*%y
# loop:
res1 = RBEGLS_loop_final(50, 0.001, 0.001, beta_RBEGLS, out2$theta, theta_ep.ini, CovU_nl, y, Z, Dall, cov_ep, lo.bounde,up.bounde)
beta_RBEGLS = res1$beta_est
theta_RBEGLS  =  res1$theta_est 
CorE_GLS = cor.mat(Dall[idx,idx], theta_RBEGLS[1], cov_ep, nug = 0)
CovE_GLS = theta_RBEGLS[2]*CorE_GLS
CovUB = CovU_nl*beta_RBEGLS[2]^2
CovUBE_GLS = CovUB + CovE_GLS
CovUBE_GLS_inv = solve(CovUBE_GLS)
A1 = solve(t(Z)%*%CovUBE_GLS_inv%*%Z)
varbeta_rb_dm = diag(A1)

#--------------------------------------------------------------------------------------------------------------
# OLS mc
Z_mc = cbind(1, z_mc)
fols_mc = lm(y~Z_mc-1)
beta_ols_mc = fols_mc$coef

# RBEGLS mc
CovUB_mc = CovU_nl_mc*beta_ols_mc[2]^2 
out2_mc = MLE.fit_GLS(y, Z_mc, Dall[idx,idx], CovUB_mc, cov_ep, theta_ep.ini, nug, "LB", lo.bounde,up.bounde)
CorE_GLS_mc = cor.mat(Dall[idx,idx], out2_mc$eta, cov_ep, nug = 0)
CovE_GLS_mc = out2_mc$theta[2]*CorE_GLS_mc
CovUBE_GLS_mc = CovUB_mc + CovE_GLS_mc
beta_RBEGLS_mc = solve(t(Z_mc)%*%solve(CovUBE_GLS_mc)%*%Z_mc)%*%t(Z_mc)%*%solve(CovUBE_GLS_mc)%*%y
# loop:
res1_mc = RBEGLS_loop_final(50, 0.001, 0.001, beta_RBEGLS_mc, out2_mc$theta, theta_ep.ini, CovU_nl_mc, y, Z_mc, Dall, cov_ep, lo.bounde,up.bounde)
beta_RBEGLS_mc = res1_mc$beta_est
theta_RBEGLS_mc  =  res1_mc$theta_est 
CorE_GLS_mc = cor.mat(Dall[idx,idx], theta_RBEGLS_mc[1], cov_ep, nug = 0)
CovE_GLS_mc = theta_RBEGLS_mc[2]*CorE_GLS_mc
CovUB_mc = CovU_nl_mc*beta_RBEGLS_mc[2]^2
CovUBE_GLS_mc = CovUB_mc + CovE_GLS_mc
CovUBE_GLS_inv_mc = solve(CovUBE_GLS_mc)
A1_mc = solve(t(Z_mc)%*%CovUBE_GLS_inv_mc%*%Z_mc)
varbeta_rb_mc = diag(A1_mc)

#--------------------------------------------------------------------------------------------------------------
# MC1
sdout = SDCal(Dall[-idx,-idx], out1$theta, cov_xi, nug = 0)
ETAallboot_all = SimuData(out1$theta, rate, locall, cov_xi, nug, MC_iter)
EVar_dm = EVar_mc = matrix(0, 2, 2)
EV_dm = EV_mc = NULL
seq_mean_dm = seq_mean_mc = seq_median_dm = seq_median_mc = NULL
for (v in 1:MC_iter){
  x1allboot = r1all%*%alpha_est + ETAallboot_all$Ymat[v,]
  x1sboot = x1allboot[-idx]
  xie = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  while ((xie[1] < 0)|(xie[2] < 0)) {
    xie = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  }
  cov_xialle = xie[2]*cor.mat(Dall, xie[1], cov_xi, nug = 0)
  tmp_inv = solve(cov_xialle[-idx,-idx])
  alphaboot = solve(t(r1s) %*% tmp_inv %*% r1s)%*%t(r1s) %*% tmp_inv %*% x1sboot 
  wboot_e = r1%*%alphaboot + cov_xialle[idx, -idx]%*%tmp_inv%*%(x1sboot - r1s%*%alphaboot)
  Zboot_e = cbind(1, g(wboot_e))
  A1_e = solve(t(Zboot_e)%*%CovUBE_GLS_inv%*%Zboot_e)
  B1_e = solve(t(Zboot_e)%*%CovUBE_GLS_inv_mc%*%Zboot_e)
  EVar_dm = EVar_dm + A1_e   
  EV_dm = cbind(EV_dm, diag(A1_e))
  EVar_mc = EVar_mc + B1_e   
  EV_mc = cbind(EV_mc, diag(B1_e))
  
  seq_mean_dm = cbind(seq_mean_dm, c(mean(EV_dm[1, ]), mean(EV_dm[2, ])))
  seq_mean_mc = cbind(seq_mean_mc, c(mean(EV_mc[1, ]), mean(EV_mc[2, ])))
  seq_median_dm = cbind(seq_median_dm, c(median(EV_dm[1, ]), median(EV_dm[2, ])))
  seq_median_mc = cbind(seq_median_mc, c(median(EV_mc[1, ]), median(EV_mc[2, ])))
}
varbeta_mc_dm = diag(EVar_dm/MC_iter)
median_pair_dm <- c(median(EV_dm[1, ]), median(EV_dm[2, ]))
varbeta_mc_mc = diag(EVar_mc/MC_iter)
median_pair_mc <- c(median(EV_mc[1, ]), median(EV_mc[2, ]))

nout_dm_b1_mc1 = length(boxplot(EV_dm[1, ], plot = F)$out)
nout_dm_b2_mc1 = length(boxplot(EV_dm[2, ], plot = F)$out)
nout_mc_b1_mc1 = length(boxplot(EV_mc[1, ], plot = F)$out)
nout_mc_b2_mc1 = length(boxplot(EV_mc[2, ], plot = F)$out)

# ---------------------------------------------------------------------------- #
# MC2 modified
EVar_dm2 = EVar_mc2 = matrix(0, 2, 2)
EV_dm2 = EV_mc2 = NULL
seq_mean_dm2 = seq_mean_mc2 = seq_median_dm2 = seq_median_mc2 = NULL
for (v2 in 1:MC_iter){
  xie2 = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  while ((xie2[1] < 0)|(xie2[2] < 0)) {
    xie2 = out1$theta + rmvnorm(1, sigma=sdout$covmat)
  }
  ETAallboot_all2 = SimuData(xie2, rate, locall, cov_xi, nug, 1)
  x1allboot2 = r1all%*%alpha_est + ETAallboot_all2$Ymat[1,]
  x1sboot2 = x1allboot2[-idx]
  cov_xialle2 = xie2[2]*cor.mat(Dall, xie2[1], cov_xi, nug = 0)
  tmp_inv2 = solve(cov_xialle2[-idx,-idx])
  alphaboot2 = solve(t(r1s) %*% tmp_inv2 %*% r1s)%*%t(r1s) %*% tmp_inv2 %*% x1sboot2 
  wboot_e2 = r1%*%alphaboot2 + cov_xialle2[idx, -idx]%*%tmp_inv2%*%(x1sboot2 - r1s%*%alphaboot2)
  Wboot_e2 = cbind(1, wboot_e2)
  
  Zboot_e2 = cbind(1, g(wboot_e2))
  A1_e2 = solve(t(Zboot_e2)%*%CovUBE_GLS_inv%*%Zboot_e2)
  B1_e2 = solve(t(Zboot_e2)%*%CovUBE_GLS_inv_mc%*%Zboot_e2)
  EVar_dm2 = EVar_dm2 + A1_e2   
  EV_dm2 = cbind(EV_dm2, diag(A1_e2))
  EVar_mc2 = EVar_mc2 + B1_e2 
  EV_mc2 = cbind(EV_mc2, diag(B1_e2))       
  
  seq_mean_dm2 = cbind(seq_mean_dm2, c(mean(EV_dm2[1, ]), mean(EV_dm2[2, ])))
  seq_mean_mc2 = cbind(seq_mean_mc2, c(mean(EV_mc2[1, ]), mean(EV_mc2[2, ])))
  seq_median_dm2 = cbind(seq_median_dm2, c(median(EV_dm2[1, ]), median(EV_dm2[2, ])))
  seq_median_mc2 = cbind(seq_median_mc2, c(median(EV_mc2[1, ]), median(EV_mc2[2, ])))
}
varbeta_mc_dm2 = diag(EVar_dm2/MC_iter)
median_pair_dm2 <- c(median(EV_dm2[1, ]), median(EV_dm2[2, ]))
varbeta_mc_mc2 = diag(EVar_mc2/MC_iter)
median_pair_mc2 <- c(median(EV_mc2[1, ]), median(EV_mc2[2, ]))

nout_dm_b1_mc2 = length(boxplot(EV_dm2[1, ], plot = F)$out)
nout_dm_b2_mc2 = length(boxplot(EV_dm2[2, ], plot = F)$out)
nout_mc_b1_mc2 = length(boxplot(EV_mc2[1, ], plot = F)$out)
nout_mc_b2_mc2 = length(boxplot(EV_mc2[2, ], plot = F)$out)

#--------------------------------------------------------------------------------------------------------------
list(alpha_est = as.matrix(out1$beta),
     xi_est = as.matrix(out1$theta),
     
     beta_ols = as.matrix(as.numeric(beta_ols)), 
     beta_RBEGLS = as.matrix(as.numeric(beta_RBEGLS)),  
     beta_RBEGLS_mc = as.matrix(as.numeric(beta_RBEGLS_mc)), 
     beta_kr = as.matrix(as.numeric(beta_kr)),
     
     theta_est = as.matrix(as.numeric(theta_RBEGLS)),
     theta_est_mc = as.matrix(as.numeric(theta_RBEGLS_mc)),
     theta_kr = as.matrix(as.numeric(theta_kr)),
     
     varbeta_ols = as.matrix(as.numeric(varbeta_ols)), 
     varbeta_kr = as.matrix(as.numeric(varbeta_kr)),
     varbeta_rb_dm = as.matrix(as.numeric(varbeta_rb_dm)),
     varbeta_rb_mct = as.matrix(as.numeric(varbeta_rb_mc)),
     
     varbeta_mc_mean_dm = as.matrix(as.numeric(varbeta_mc_dm)),
     varbeta_mc_median_dm = as.matrix(as.numeric(median_pair_dm)),
     varbeta_mc_mean_mc = as.matrix(as.numeric(varbeta_mc_mc)),
     varbeta_mc_median_mc = as.matrix(as.numeric(median_pair_mc)),
     
     varbeta_mc_mean_dm2 = as.matrix(as.numeric(varbeta_mc_dm2)),
     varbeta_mc_median_dm2 = as.matrix(as.numeric(median_pair_dm2)),
     varbeta_mc_mean_mc2 = as.matrix(as.numeric(varbeta_mc_mc2)),
     varbeta_mc_median_mc2 = as.matrix(as.numeric(median_pair_mc2)),
     
     seq_mean_dm = seq_mean_dm,
     seq_mean_mc = seq_mean_mc,
     seq_median_dm = seq_median_dm,
     seq_median_mc = seq_median_mc,
     
     seq_mean_dm2 = seq_mean_dm2,
     seq_mean_mc2 = seq_mean_mc2,
     seq_median_dm2 = seq_median_dm2,
     seq_median_mc2 = seq_median_mc2,
     
     nout_dm_b1_mc1 = nout_dm_b1_mc2,
     nout_dm_b2_mc1 = nout_dm_b2_mc2,
     nout_mc_b1_mc1 = nout_mc_b1_mc2,
     nout_mc_b2_mc1 = nout_mc_b2_mc2,
     
     nout_dm_b1_mc2 = nout_dm_b1_mc2,
     nout_dm_b2_mc2 = nout_dm_b2_mc2,
     nout_mc_b1_mc2 = nout_mc_b1_mc2,
     nout_mc_b2_mc2 = nout_mc_b2_mc2
)


#}








