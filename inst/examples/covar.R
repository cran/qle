# Compare prediction variances based on different covariance functions
# and reml parameter estimations and using cross-validation

library(qle)
data(normal)

## fit first covariance model
X <- as.matrix(qsd$qldata[,1:2])
T <- qsd$qldata["mean.T1"]
V <- qsd$qldata[["var.T1"]]

# get sample means of simulated statistics
Tstat <- qsd$qldata[grep("mean.",names(qsd$qldata))]

nsim <- attr(qsd$qldata,"nsim")

# Setup SIRF-2 (default) covariance model manually using simulation variances `V`
# 
# Covariance parameter 'alpha' is fixed and not estimated by REML (see below),
# however, any parameter of the covariance model can be estimated
# or excluded from the covariance fitting procedure
cvm <- setCovModel(model="sirfk",  param = c("scale"=0.001,"alpha"=2.0),
		fixed.param = c("alpha"), nugget=1e-4, npoints = nrow(X), trend = 2,
		  var.sim = V/nsim)

print(cvm)

# and fitt by REML	
fit <- fitCov(list(cvm),X,T,verbose=TRUE)[[1]]
# compare with original fit
(pnew <- fit[[1]]$param)
qsd$covT[[1]]$param

# test reml evaluation
val <- reml(list(cvm),pnew,T[1],X)[[1]]
stopifnot(attr(fit,"optres")$objective == val)

# Setup kriging (covariance) model for the second statistic
T2 <- qsd$qldata["mean.T2"]
V2 <- qsd$qldata[["var.T2"]]
# Use `matern` model and fit
cvm2 <- setCovModel(model="matern",
	     param = c("scale"=0.1,"nu"=2.5,"rho"=2),
		 nugget = 1e-4, npoints = nrow(X),  trend = 2,
		 var.sim = 1e-6)
# fitting by REML
(fit2 <- fitCov(list(cvm2),X,T2,verbose=TRUE)[[1]])

## all at once fit
# fitCov(
#    list(cvm,cvm2),         # both models
#	 X,
#    Tstat,					 # both statistics
#    verbose=TRUE)

# merge covartiance models 
qsd2 <- qsd
qsd2$covT <- structure(list(fit$model,fit2$model),class="krige")

## Grid for MSE estimation
x <- seq(qsd2$lower[1],qsd2$upper[1],by=0.05)
y <- seq(qsd2$lower[2],qsd2$upper[2],by=0.05)
p <- as.matrix(expand.grid(x,y))
 
## Kriging MSE
## old fit and new fit (second with `matern`)
kvar <- list(varKM(qsd$covT,p,X,Tstat),		
			 varKM(qsd2$covT,p,X,Tstat))

## Empirically integrated MSE (by estimated kriging variances):
colMeans(kvar[[1]]) # original
colMeans(kvar[[2]]) # 2nd is fit by `matern` covariance

## show prediction variances
## Prediction variances using SIRF-2 covariance

## first row: original fit, both statistics
## second row: new fit, both statistics
op <- par(mfrow=c(2,2))
for(j in 1:2){
 for(i in 1:2) {
	z1 <- matrix(kvar[[j]][,i],ncol=length(y))
	plot(x = 0, y = 0, type = "n", xlim=,range(x), ylim=range(y),xlab = "", ylab = "")
	contour(x, y, z1, col = "black", lty = "solid",
			nlevels = 50, add = TRUE,vfont = c("sans serif", "plain"))
	try(points(X,pch=23,cex=0.8,bg="black"),silent=TRUE) 
 }
}
par(op)


## Estimation of prediction errors by cross-validation
## using the original covariance fit with covariance function SIRF-k 
cv1 <- prefitCV(qsd)
cv2 <- prefitCV(qsd2)
cve <- list(crossValTx(qsd, cv1, p, type = "cve"),
			crossValTx(qsd2, cv2, p, type = "cve"))

# first row: original fit, both statistics
# second row: new fit, both statistics
dev.new()
op <- par(mfrow=c(2,2))
for(j in 1:2){	
	for(i in 1:2){	
	z3 <- matrix(cve[[j]][,i],ncol=length(y))
	plot(x = 0, y = 0, type = "n", xlim=,range(x), ylim=range(y),xlab = "", ylab = "")
	contour(x, y, z3, col = "black", lty = "solid",
			nlevels = 50, add = TRUE,vfont = c("sans serif", "plain"))
	try(points(X,pch=23,cex=0.8,bg="black"),silent=TRUE)
 }
}
par(op)
