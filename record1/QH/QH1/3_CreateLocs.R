
#setwd("C:/Users/qihan/Desktop/SpE_ncl/locs")

library(spatstat)

n = 1600  # 1000 1500 2000 2500 3000

rate  = 0.5                            # sample size
lbase = 5
l=lbase*(n/50)^rate                     # spatial domain of interest [-l/2,l/2]^2
lsize = (lbase/10)*(n/50)^(rate-1/2)    # 0.2*lsize: the minimum distance of two locations 

set.seed(2019+1000*rate)
locall = rSSI(0.2*lsize, n, win=as.owin(c(-l/2,l/2,-l/2,l/2)))

nameRdata=strwrap(paste("loc", n, ".Rdata", sep=""))
save(locall, file=nameRdata)


#x1all = rnorm(n)
#save(x1all, file = "x1.Rdata")

# 3_CreateLocs.R is for creating (x,y) locations