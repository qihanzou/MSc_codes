n <- 5000;
m <- array(rnorm(n^2),c(n,n));
m2 <- (m+t(m))/sqrt(2*n);# Make m symmetric
lambda <- eigen(m2, symmetric=T, only.values = T);
e <- lambda$values;
hist(e,breaks=seq(-2.01,2.01,.02),
     main=NA, xlab="Eigenvalues",freq=F)