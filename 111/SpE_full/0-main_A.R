    # -------------------------------------------------- #
    #        Install and load needed R packages          #
    # -------------------------------------------------- #
    rm(list = ls(all = TRUE)) # The line removes all variables from the current environment
    time1 = Sys.time() # record starting time
    list.of.packages <- c("doSNOW","foreach","doRNG","iterators","data.table","parallel",
                          "proxy", "pls","mvtnorm")
    new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
    if(length(new.packages)) install.packages(new.packages)
    invisible(lapply(list.of.packages, require, character.only = TRUE))
    # -------------------------------------------------- #
    # choose parameters
    # -------------------------------------------------- #
    nall = 500 # 1000, 1500, 2000
    idx = 1:100
    keepN = length(idx)
    ver = 1
  
    # -------------------------------------------------- #
    #  nreps: number of iterations                       #
    #  1. To see result of $i$th iteration, set          #
    #     nreps = 1                                      #
    #     iter0 = i                                      #
    #  2. To see all results (time consuming) , set      #
    #     nreps = 400                                    #
    #     iter0 = 1:nreps                                #
    # -------------------------------------------------- #
    
    nreps = 400  #  iteration use 400
    iter0 = 1:nreps  # any subset of (1:400) 
    # -------------------------------------------------- #
    #             Set working Directory                  #
    # -------------------------------------------------- #
    R_address = "C:/Users/qihan/Desktop/SpE_full"
    w1_address = "C:/Users/qihan/Desktop/SpE_full"
    # -------------------------------------------------- #
    
    library(parallel)
    worker.script = '1-fittingNF_FULL.R'
    numCores <- detectCores()
    cl <- makeCluster(min(numCores-4,nreps, 30))
    clusterExport(cl, ls())
    setwd(w1_address)
    source(worker.script)
    clusterEvalQ(cl, setwd(w1_address))
    clusterEvalQ(cl, source(worker.script))
    registerDoSNOW(cl)
    
    seedM = rbind(rep(401,nreps),c(iter0),c(iter0))
    res = foreach(iter = iter0, .combine=list,.maxcombine=max(nreps,2),.options.RNG=seedM,
                  .multicombine=TRUE, .errorhandling = 'pass') %dorng% {
                    res =  run.sim()
                    return(res)              
                  }
  
    rdata_name = strwrap(paste("n_", nall,"_", keepN,"_ver",ver, "_FULL",".Rdata",sep=""))
    save.image(rdata_name)
    stopCluster(cl)
    time2 = Sys.time()
    
    # return the estimated results and the running time
    print(res[[1]])
    print(time2 - time1)
  
