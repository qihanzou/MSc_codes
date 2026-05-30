#rm(list = ls(all = TRUE)) # The line removes all variables from the current environment
time1 = Sys.time() # record starting time
list.of.packages <- c("doSNOW","foreach","doRNG","iterators","data.table","parallel",
                      "proxy", "pls","mvtnorm")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(list.of.packages, require, character.only = TRUE))
# -------------------------------------------------- #
#             Set working Directory                  #
# -------------------------------------------------- #
R_address = "C:/Users/qihan/Desktop/Qihan"
w1_address= "C:/Users/qihan/Desktop/Qihan"
# -------------------------------------------------- #
# choose parameters
# -------------------------------------------------- #

#nall = 400 # 1000, 1500, 2000
#idx = 1:(nall/2) #1:450
#idx = sample(1:nall, 100)
#Nver = 1
#beta1 = 4


keepN = length(idx)

nreps = 1000   #  iteration use 400
iter0 = 1:nreps  # any subset of (1:400) 
boot_iter = 100


library(parallel)
worker.script = '1-Parameter_New2.R'
numCores <- detectCores()
cl <- makeCluster(min(numCores-1,nreps, 30))
clusterExport(cl, ls())
setwd(w1_address)
source(worker.script)
clusterEvalQ(cl, setwd(w1_address))
clusterEvalQ(cl, source(worker.script))
registerDoSNOW(cl)

res = foreach(iter = iter0, .combine=list,.maxcombine=max(nreps,2),.options.RNG=2024+beta1,
              .multicombine=TRUE, .errorhandling = 'pass') %dorng% {
                res =  run.sim()
                return(res)              
              }

# save results
rdata_name = strwrap(paste("X",nreps,"n_", nall,"_", keepN, ".Rdata",sep=""))

save.image(rdata_name)
stopCluster(cl)
time2 = Sys.time()

# return the estimated results and the running time
print(time2 - time1)

