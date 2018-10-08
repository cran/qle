## Quasi-likelihood simulation based estimation
##
## 1. use criterion `qle` for estimation
## 2. then `mahal` as (generalized) least squares
library(qle)
data(normal)

# setting number of local cores
options(mc.cores=8L)

## one step minimization
## no sammpling
# x0 <- c(2.5,1.5)
# qscoring(qsd, x0, pl=10,verbose=TRUE)

## alternatively use minimization by `nloptr`
# searchMinimizer(x0, qsd, method = c("bobyqa"), verbose=TRUE)

# main estimation with new evaluations
# (simulations of the statistical model)
OPT <- qle(qsd,qsd$simfn,nsim=20,
		global.opts=list("maxeval"=10),
		pl=10)

## intermediate results
# attr(OPT,"tracklist")

# check final result with quasi-deviance at solution
local <- OPT$final
info <- attr(OPT,"optInfo")
QD <- quasiDeviance(OPT$par,OPT$qsd,W=info$W,theta=info$theta,verbose=TRUE)[[1]]
QD$value
local$value

## here no quasi-scoring, therefore search on criterion function only
#OPT <- qle(qsd,qsd$simfn,nsim=10,		
#		global.opts=list("maxeval"=10),
#		local.opts=list("ftol_abs"=1e-10,"weights"=c(0.5),"useWeights"=FALSE,"test"=TRUE),
#		method="bobyqa",pl=10, plot=TRUE)

## only global search, no testing for roots
#OPT <- qle(qsd,qsd$simfn,nsim=10,
#		local.opts=list("ftol_abs"=0),
#		global.opts=list("maxeval"=10),
#		method="bobyqa",pl=10)

#OPT$final 
#attr(OPT,"tracklist")
#attr(OPT,"optInfo")

## restart estimation and do a pure global search (setting `ftol_abs=0`), 
## sample additional points for evaluation and select new candidates by
## criterion `var`
#GL <- qle(OPT$qsd, qsd$simfn, nsim=10,		
#		global.opts = list("maxiter"=10, "stopval"=0),
#		local.opts = list("nextSample"="var","ftol_abs"=0),
#		pl=10, iseed=1234)

## show final results
#library(rgl)
#cvm <- NULL
#qsd <- OPT$qsd
#
#x <- seq(qsd$lower[1],qsd$upper[1],by=0.05)
#y <- seq(qsd$lower[2],qsd$upper[2],by=0.05)
#points <- as.matrix(expand.grid(x,y))
#X <- as.matrix(qsd$qldata[,1:2])
#Xp <- quasiDeviance(X,qsd,cvm=cvm,value.only=TRUE)
#D <- quasiDeviance(points,qsd,cvm=cvm,value.only=TRUE)
#
#z <- matrix(D,ncol=length(y))
#persp3d(x,y,z,col="red", alpha=0.3, axes=TRUE)
#cnt <- contourLines(x,y,z,
#	      lev=seq(range(z)[1],range(z)[2],
#		  by=dist(range(z))/100))
#for (i in 1:length(cnt))
# with(cnt[[i]], lines3d(x, y, level, col="darkred"))
#points3d(X[,1],X[,2], Xp, size=3, col="red",add=TRUE)
#
## contour plot
#dev.new()
#z1 <- matrix(D,ncol=length(y))
#plot(x = 0, y = 0, type = "n", xlim=,range(x), ylim=range(y),xlab = "", ylab = "")
#contour(x, y, z1, col = "black", lty = "solid",
#		nlevels = 50, add = TRUE,vfont = c("sans serif", "plain"))
#try(points(X,pch=23,cex=0.8,bg="black"),silent=TRUE)

# sampling
#S <- attr(OPT,"optInfo")$W
#Y <- nextLOCsample(S, OPT$par,100,qsd$lower,qsd$upper, invert=TRUE)
#points(Y,cex=0.8,col="red")

