# (1) Kriging prediction of sample means of statistics
# (2) Estimation of variance matrix by sample average approximation
library(qle)
data(normal)

# points to estimate statistics
p <- multiDimLHS(N=10,qsd$lower,qsd$upper,
		method="randomLHS",type="list")

# kriging mean
X <- as.matrix(qsd$qldata[,1:2])
# values of statistics
Tstat <- qsd$qldata[grep("mean.",names(qsd$qldata))]     
# Cholesky decompostions of variance matrices
Lstat <- qsd$qldata[grep("L+",names(qsd$qldata))]		   

# kriging prediction (low level functions)
(est.m1 <- estim(qsd$covT,p,X,Tstat,krig.type="var"))
(est.m2 <- estim(qsd$covT[[1]],p,X,Tstat[1],krig.type="var"))

stopifnot(
  all.equal(as.numeric(extract(est.m1,type="mean")[,1]),
			as.numeric(extract(est.m2,type="mean")))
)

# estimate derivatives
jacobian(qsd$covT,p,X,Tstat)

# average variance matrix interpolation
covarTx(qsd,theta=X,useVar=TRUE,doInvert=TRUE)

# predict and extract values 
predictKM(qsd$covT,p,X,Tstat,krig.type="var")

# calculate kriging prediction variances
varKM(qsd$covT,p,X,Tstat)
