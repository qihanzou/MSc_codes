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
R_address = "C:/Users/qihan/Desktop/rebound"
w1_address= "C:/Users/qihan/Desktop/rebound"
# -------------------------------------------------- #


keepN = length(idx)

nreps = 1000   
iter0 = 1:nreps   


library(parallel)
worker.script = '1-Parameter_rebound.R'
numCores <- detectCores()
cl <- makeCluster(min(numCores-1,nreps, 30))
clusterExport(cl, ls())
setwd(w1_address)
source(worker.script)
clusterEvalQ(cl, setwd(w1_address))
clusterEvalQ(cl, source(worker.script))
registerDoSNOW(cl)


seed2024 = as.matrix(read.csv("C:/Users/qihan/Desktop/rebound/seed2024.csv", header = TRUE))

res = foreach(iter = iter0, .combine=list,.maxcombine=max(nreps,2),.options.RNG = seed2024,
              .multicombine=TRUE, .errorhandling = 'pass') %dorng% {
                res =  run.sim()
                return(res)              
              }

rdata_name = strwrap(paste("rebound_",nreps,"n_", nall,"_", keepN, ".Rdata",sep=""))

save.image(rdata_name)
stopCluster(cl)
time2 = Sys.time()

print(time2 - time1)

