## Computing the quasi-deviance: 
## QD using a (weighted) average approximation of variance matrix
library(qle)
data(normal)

## at some parameter
# x <- c("mu"=2,"sigma"=1)		 # numeric vector 
# x <- cbind("mu"=2,"sigma"=1)   # as a matrix 
x <- list(c("mu"=2,"sigma"=1))   # list of numeric vectors

## compare both values:
## must be equal since number of statistics
## equals number of parameters to estimate
mahalDist(x,qsd)
quasiDeviance(x,qsd)

## alternatively fit cross-validation models
# cvm <- prefitCV(qsd, verbose=TRUE)
# use prediction errors based on these
cvm <- NULL
qsd$var.type <- "cholMean"
(QD1 <- quasiDeviance(x,qsd,cvm=cvm))

## (weighted) average variance matrix approximation
## by `wcholMean` option using a weighting matrix `W`
## and prediction variances calculated at `theta`
W <- QD1[[1]]$I
qsd$var.type <- "wcholMean"
(QD3 <- quasiDeviance(x,qsd,W=W,theta=c(2,1),cvm=cvm))

############################################################
## 3D plot of quasi-deviance, with kriged variance matrix
## and added CV prediction variances, requires package `rgl`
############################################################

#library(rgl)
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
