### R code from vignette source 'qle_with_R.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: qle_with_R.Rnw:1085-1087
###################################################
options(useFancyQuotes="UTF-8")
options(digits=4, prompt="R> ")


###################################################
### code chunk number 2: qle_with_R.Rnw:1154-1158
###################################################
library(qle)
library(graphics)
RNGkind("L'Ecuyer-CMRG")
set.seed(1356)


###################################################
### code chunk number 3: qle_with_R.Rnw:1164-1168
###################################################
cond <- list("n"=25)
simfn <- function(tet,cond){
	mean(rgeom(cond$n,prob=1-tet[1]))
}


###################################################
### code chunk number 4: qle_with_R.Rnw:1172-1174
###################################################
lb <- c("rho"=0.05)
ub <- c("rho"=0.95)


###################################################
### code chunk number 5: qle_with_R.Rnw:1180-1185
###################################################
nsim <- 10
X <- multiDimLHS(N=10,lb=lb,ub=ub,
      method="maximinLHS",type="matrix")

sim <- simQLdata(sim=simfn,cond=cond,nsim=nsim,X=X)


###################################################
### code chunk number 6: qle_with_R.Rnw:1191-1193
###################################################
qsd <- getQLmodel(sim, lb, ub, obs=c("N"=1),
		var.type="wlogMean",verbose=TRUE)


###################################################
### code chunk number 7: qle_with_R.Rnw:1197-1199
###################################################
S0 <- qscoring(qsd,x0=c("rho"=0.8),verbose=TRUE)
print(S0)


###################################################
### code chunk number 8: qle_with_R.Rnw:1210-1215
###################################################
OPT <- qle(qsd,simfn,cond=cond,	     	
  global.opts = list("maxeval"=5, "NmaxLam"=5),
  local.opts = list("nextSample"="score","weights"=0.5,"ftol_abs"=1e-4,
                    "lam_max"=1e-5,"useWeights"=TRUE),
  method = c("qscoring","bobyqa","direct"), iseed=1356) 


###################################################
### code chunk number 9: qle_with_R.Rnw:1220-1238
###################################################
## statistics
op <- par(xaxs='i', yaxs='i')
rho <- as.matrix(seq(0.1,0.9,by=0.001))
y <- as.numeric(unlist(simQLdata(sim=simfn,cond=cond,nsim=nsim,X=rho,mode="mean")))
T <- qsd$qldata[grep("mean.",names(qsd$qldata))]
Y <- predictKM(qsd$covT,rho,X,T,krig.type="var")
# steady state values
y0 <- rho/(1-rho)
plot(NULL, type="n", xlab=expression(rho),
		ylab="y",xlim=c(0,1), ylim=c(0,10))
lines(as.numeric(rho),y,col="black",lt=2,lwd=0.3)
lines(as.numeric(rho),Y,col="blue",lwd=0.3)
lines(as.numeric(rho),y0,col="red",lwd=0.3)
legend("topleft", c("Number of customers in the system",
	"Expected number at steady state","Kriging approximation"),
		lty=c(2,1,1),col=c("black","red","blue"),
		xpd=TRUE,pt.cex=1,cex=1)
par(op)


###################################################
### code chunk number 10: qle_with_R.Rnw:1240-1262
###################################################
op <- par(xaxs='i', yaxs='i')
p <- seq(lb,ub,by=0.0001)
QD <- quasiDeviance(X,qsd,value.only=TRUE)
qd <- quasiDeviance(as.matrix(p),qsd)
y <- sapply(qd,"[[","value")
score <- sapply(qd,"[[","score")
## plot quasi-deviance and quasi-score function
plot(NULL, type="n", xlab=expression(rho),
		ylab="",xlim=c(0,1), ylim=c(-10,50))
abline(h=0,col="gray")
points(X,QD,pch=3,cex=1)
lines(p,score, type='l',col="blue",lwd=1.5) 
lines(p,y,col="black",lwd=0.8)
legend("topleft", c("quasi-deviance","quasi-score","sample points", "approximate root","additional samples"),
		lty=c(1,1),lwd=c(1.5,1.5,NA,NA,NA),pch=c(NA,NA,3,5,8),
		col=c("black","blue","black","magenta","green"),pt.cex=1,cex=1)
points(S0$par,S0$val,col="magenta",pch=5,cex=1)
nmax <- OPT$ctls["maxeval","val"]
X <- as.matrix(qsd$qldata[,1])
Xnew <- OPT$qsd$qldata[(nrow(X)+1):(nrow(X)+nmax),1]
points(cbind(Xnew,0),pch=8,cex=2,col="green")
par(op)


###################################################
### code chunk number 11: qle_with_R.Rnw:1270-1271
###################################################
OPT


###################################################
### code chunk number 12: qle_with_R.Rnw:1276-1295
###################################################
op <-par(xaxs='i', yaxs='i')
qd <- quasiDeviance(as.matrix(p),OPT$qsd)
y <- sapply(qd,"[[","value")
score <- sapply(qd,"[[","score")
## plot quasi-deviance and quasi-score function
plot(NULL, type="n", xlab=expression(rho),
		ylab="",xlim=c(0,1), ylim=c(-10,50))
abline(h=0,col="gray")
lines(p,score, type='l',col="blue",lwd=1.5) 
lines(p,y,col="black",lwd=0.8)
legend("topleft", c("quasi-deviance","quasi-score","sample points", "QL estimate"),
		lty=c(1,1),lwd=c(1,1,NA,NA,NA),pch=c(NA,NA,3,5,8),
		col=c("black","blue","black","magenta","green"),pt.cex=1,cex=1)

X <- as.matrix(OPT$qsd$qldata[,1])
QD <- quasiDeviance(X,OPT$qsd,value.only=TRUE)
points(X,QD,pch=3,cex=1)
points(OPT$par,OPT$val,col="magenta",pch=5)
par(op)


###################################################
### code chunk number 13: qle_with_R.Rnw:1305-1306
###################################################
checkMultRoot(OPT,verbose = TRUE)


###################################################
### code chunk number 14: qle_with_R.Rnw:1316-1319
###################################################
X <- as.matrix(OPT$qsd$qldata[,1])
Tstat <- OPT$qsd$qldata[grep("mean.",names(qsd$qldata))]   
predictKM(OPT$qsd$covT,c("rho"=0.5),X,Tstat)


###################################################
### code chunk number 15: qle_with_R.Rnw:1354-1356
###################################################
tet0 <- c("rho"=0.5)
obs0 <- simQLdata(sim=simfn,cond=cond,nsim=100,X=tet0)


###################################################
### code chunk number 16: qle_with_R.Rnw:1361-1368
###################################################
mle <- do.call(rbind,
		lapply(obs0[[1]],function(y,n){
               tet <- 1-1/(1+y[[1]])
               c("mle.rho"=tet,"mle.var"=(tet*(1-tet)^2)/n)
            }, n=cond$n))
x <- mle[,1]-tet0
mle.var <- c(sum(x^2)/length(x),mean(mle[,2]))


###################################################
### code chunk number 17: qle_with_R.Rnw:1372-1390
###################################################
cl <- makeCluster(2)
clusterSetRNGStream(cl)
clusterExport(cl,list("qle"))
OPTS <- parLapplyLB(cl,obs0[[1]],
          function(obs,...) qle(...,obs=obs),
         qsd=qsd,
         sim=simfn, 
         cond=cond,
         global.opts = list("maxeval"=5,"NmaxLam"=5),
         local.opts = list("nextSample"="score","weights"=0.5,
				"ftol_abs"=1e-4,"lam_max"=1e-5,
				"useWeights"=TRUE),
         method = c("qscoring","bobyqa","direct"))
stopCluster(cl)
QLE <- do.call(rbind,lapply(OPTS,
 function(x) c("qle"=x$par,"qle.var"=1/as.numeric(x$final$I))))	
y <- QLE[,1]-tet0
qle.var <- c(sum(y^2)/length(y),mean(QLE[,2]))


###################################################
### code chunk number 18: qle_with_R.Rnw:1422-1424
###################################################
Stest <- qleTest(OPT,sim=simfn,cond=cond,obs=obs0)
Stest


###################################################
### code chunk number 19: qle_with_R.Rnw:1433-1435
###################################################
rbind("Var"=c("test"=attr(Stest,"aiqm"),"study"=qle.var[2]),
	  "MSE"=c(as.vector(Stest$param["RMSE"]^2),qle.var[1]))


###################################################
### code chunk number 20: qle_with_R.Rnw:1454-1462
###################################################
set.seed(123)
# setting the number of cores>1
# uses parallel simulations and computations
options(mc.cores=2)
simfunc <- function(pars) {	
	x <- rnorm(10,mean=pars["mu"],sd=pars["sigma"])    
	c("T1"=median(x),"T2"=mad(x))	
}


###################################################
### code chunk number 21: qle_with_R.Rnw:1467-1469
###################################################
lb <- c("mu"=0.5,"sigma"=0.1)
ub <- c("mu"=8.0,"sigma"=5.0)


###################################################
### code chunk number 22: qle_with_R.Rnw:1475-1478
###################################################
sim <- simQLdata(sim=simfunc,
           nsim=10,N=8,lb=lb,ub=ub,
              method="maximinLHS")


###################################################
### code chunk number 23: qle_with_R.Rnw:1484-1485
###################################################
obs <- structure(c("T1"=2,"T2"=1),class="simQL")


###################################################
### code chunk number 24: qle_with_R.Rnw:1488-1489
###################################################
qsd <- getQLmodel(sim,lb,ub,obs,var.type="wcholMean")


###################################################
### code chunk number 25: qle_with_R.Rnw:1493-1495
###################################################
QS <- qscoring(qsd, x0=c("mu"=5,"sigma"=3.0), verbose=TRUE)
print(QS)


###################################################
### code chunk number 26: qle_with_R.Rnw:1503-1510
###################################################
OPT <- qle(qsd,
        simfunc,		
        nsim=100,
        global.opts=list("maxeval"=50),
        local.opts=list("lam_max"=1e-4,"weights"=0.5),
		iseed=123)
print(OPT)


###################################################
### code chunk number 27: qle_with_R.Rnw:1514-1556
###################################################
op <- par(mfrow=c(1, 2), mar=c(5.1, 4.1, 1.1, 1.1),
		oma=c(5,4,1,1),xaxs='i', yaxs='i',
		cex=2.2, cex.axis=2.2, cex.lab=2.2)
# get points for plotting
theta0 <- c("T1"=2,"T2"=1)
x <- seq(qsd$lower[1],qsd$upper[1],by=0.05)
y <- seq(qsd$lower[2],qsd$upper[2],by=0.05)
p <- as.matrix(expand.grid(x,y))
X <- as.matrix(qsd$qldata[,1:2])
Tstat <- qsd$qldata[grep("mean.",names(qsd$qldata))]
Xp <- quasiDeviance(X,qsd,value.only=TRUE)
D <- quasiDeviance(p,qsd,value.only=TRUE)
z <- matrix(D,ncol=length(y))
Xnext <- as.matrix(OPT$qsd$qldata[,1:2])
Dnext <- quasiDeviance(p,OPT$qsd,value.only=TRUE)
znext <- matrix(Dnext,ncol=length(y))
nmax <- OPT$ctls["maxeval","val"]
Xnew <- OPT$qsd$qldata[(nrow(X)+1):(nrow(X)+nmax),c(1,2)]
# left
plot(x=0,y=0,type="n", xlim=range(x),ylim=range(y),
	 xlab=expression(mu),ylab=expression(sigma))
contour(x,y,z,col="black",lty="solid",nlevels=50,add=TRUE)
#
points(X,pch=23,cex=2,bg="black")
points(Xnew,pch=8,cex=2,col="green")
# right
plot(x=0,y=0,type="n", xlim=range(x),ylim=range(y),
	 xlab=expression(mu),ylab=expression(sigma))
contour(x,y,znext,col="black",lty="solid",nlevels=50,add=TRUE)
points(Xnext,pch=23,cex=2,bg="black")
points(rbind(OPT$par),pch=18,cex=2.5,col="magenta")
points(rbind(unlist(theta0)),pch=17,cex=2.5,col="red")
# legend 
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
cols <- c("black","green","magenta","red")
legend("bottomleft",text.width=c(0.45,0.45,0.45,0.45),
		legend=c("initial design points", "new sample points",
				"estimated parameter","true parameter"),
		pch=c(23,8,18,17),col=cols,pt.bg=cols,bty='n',
		horiz=TRUE,xpd=TRUE,pt.cex=2.5,cex=2.5)
par(op)


###################################################
### code chunk number 28: qle_with_R.Rnw:1565-1566
###################################################
checkMultRoot(OPT,verbose=TRUE)


###################################################
### code chunk number 29: qle_with_R.Rnw:1574-1577
###################################################
obs0 <- simQLdata(simfunc,X=OPT$par,nsim=1000,mode="matrix")[[1]]
var(obs0)
attr(OPT$final,"Sigma")


###################################################
### code chunk number 30: qle_with_R.Rnw:1580-1583
###################################################
Stest <- qleTest(OPT,method=c("qscoring","bobyqa"),
		   sim=simfunc,nsim=1000)
Stest


###################################################
### code chunk number 31: qle_with_R.Rnw:1607-1613
###################################################
data(matclust)
OPT <- matclust$OPT
qsd <- matclust$qsd
cvm <- matclust$cvm
Stest <- matclust$Stest
library(spatstat)


###################################################
### code chunk number 32: qle_with_R.Rnw:1617-1619
###################################################
data(redwood)
fitMat <- kppm(redwood, ~1, "MatClust")


###################################################
### code chunk number 33: qle_with_R.Rnw:1622-1623
###################################################
fitMat$modelpar


###################################################
### code chunk number 34: qle_with_R.Rnw:1625-1627 (eval = FALSE)
###################################################
## RNGkind("L'Ecuyer-CMRG")
## set.seed(297)


###################################################
### code chunk number 35: qle_with_R.Rnw:1631-1640
###################################################
simStat <- function(X,cond){
 x <- Kest(X,r=cond$rr,correction="best")
 x <- x[[attr(x,"valu")]]
 x <- x[x>0]
 if(anyNA(x) || any(!is.finite(x))) {
   warning(.makeMessage("`NA`, `NaN` or `Inf` detected.","\n"))
   x <- x[!is.nan(x) & is.finite(x)]}
 return(c(intensity(X),x))	
}


###################################################
### code chunk number 36: qle_with_R.Rnw:1644-1648
###################################################
simClust <- function(theta,cond){
 X <- rMatClust(theta["kappa"],theta["R"],theta["mu"],win=cond$win)	
 simStat(X,cond)
}


###################################################
### code chunk number 37: qle_with_R.Rnw:1653-1657
###################################################
nsim <- 50
Nsample <- 12
cond <- list(win=owin(c(0, 2),c(0, 2)),
             rr=seq(0,0.3,by=0.05))


###################################################
### code chunk number 38: qle_with_R.Rnw:1660-1662
###################################################
lb <- c("kappa"=20,"R"=0.01,"mu"=1)
ub <- c("kappa"=30,"R"=0.25,"mu"=5)	 


###################################################
### code chunk number 39: qle_with_R.Rnw:1666-1670 (eval = FALSE)
###################################################
## cl <- makeCluster(8)
## clusterSetRNGStream(cl)
## clusterCall(cl,fun=function(x) library("spatstat", character.only=TRUE))
## clusterExport(cl=cl,varlist=c("simStat"), envir=environment())


###################################################
### code chunk number 40: qle_with_R.Rnw:1674-1675 (eval = FALSE)
###################################################
## obs0 <- simStat(redwood,cond)


###################################################
### code chunk number 41: qle_with_R.Rnw:1678-1680 (eval = FALSE)
###################################################
## sim <- simQLdata(sim=simClust,cond=cond,nsim=nsim,
##           method="randomLHS",lb=lb,ub=ub,N=Nsample,cl=cl)


###################################################
### code chunk number 42: qle_with_R.Rnw:1683-1685 (eval = FALSE)
###################################################
## qsd <- getQLmodel(sim,lb,ub,obs0,criterion="qle",
##          var.type="kriging",verbose=TRUE)


###################################################
### code chunk number 43: qle_with_R.Rnw:1695-1696 (eval = FALSE)
###################################################
## cvm <- prefitCV(qsd, reduce=FALSE, verbose=TRUE)


###################################################
### code chunk number 44: qle_with_R.Rnw:1702-1703
###################################################
crossValTx(qsd, cvm, type = "acve")


###################################################
### code chunk number 45: qle_with_R.Rnw:1708-1709
###################################################
crossValTx(qsd, cvm, type = "mse")


###################################################
### code chunk number 46: qle_with_R.Rnw:1720-1721
###################################################
crossValTx(qsd, cvm, type = "ascve")


###################################################
### code chunk number 47: qle_with_R.Rnw:1731-1732
###################################################
attr(cvm,"type") <- "max"


###################################################
### code chunk number 48: qle_with_R.Rnw:1736-1739
###################################################
x0 <- c("kappa"=24,"R"=0.08,"mu"=2.5)
searchMinimizer(x0,qsd,info=TRUE,
		method="direct",cvm=cvm,verbose=TRUE)


###################################################
### code chunk number 49: qle_with_R.Rnw:1742-1745
###################################################
qscoring(qsd,x0,
  opts=list("ftol_rel"=1e-6,"slope_tol"=1e-4),
  cvm=cvm,verbose=TRUE)


###################################################
### code chunk number 50: qle_with_R.Rnw:1760-1772 (eval = FALSE)
###################################################
## OPT <- qle(qsd, simClust, cond=cond,  
## 		global.opts = list("maxiter"=10,
##                            "maxeval" = 20,
##                            "weights"=c(1,5,10)),
## 		local.opts = list("lam_max"=1e-2,
##                           "nextSample"="score",
##                           "ftol_abs"=1e-2,
##                           "weights"=c(0.55),
##                           "eta"=c(0.025,0.075),
##                           "test"=TRUE),
## 		method = c("qscoring","bobyqa","direct"),
## 		pl=2, errType="max", iseed=297, cl=cl)


###################################################
### code chunk number 51: qle_with_R.Rnw:1775-1776
###################################################
print(OPT)


###################################################
### code chunk number 52: qle_with_R.Rnw:1779-1780
###################################################
attr(OPT,"optInfo")


###################################################
### code chunk number 53: qle_with_R.Rnw:1785-1786
###################################################
OPT$final


###################################################
### code chunk number 54: qle_with_R.Rnw:1790-1792
###################################################
S0 <- searchMinimizer(OPT$par,OPT$qsd,
       method="bobyqa",cvm=OPT$cvm,verbose=TRUE)


###################################################
### code chunk number 55: qle_with_R.Rnw:1797-1800
###################################################
QS <- qscoring(OPT$qsd,OPT$par,
       opts=list("slope_tol"=1e-4,"score_tol"=1e-3),
       cvm=OPT$cvm)


###################################################
### code chunk number 56: qle_with_R.Rnw:1804-1806
###################################################
par <- rbind("QS"=QS$par,"S0"=S0$par)
checkMultRoot(OPT,par=par)


###################################################
### code chunk number 57: qle_with_R.Rnw:1809-1810
###################################################
OPT$par


###################################################
### code chunk number 58: qle_with_R.Rnw:1815-1818 (eval = FALSE)
###################################################
## Stest <- qleTest(OPT,sim=simClust,cond=cond, 
##   nsim=1000,method=c("qscoring","bobyqa","direct"),  
##   cl=cl,verbose=TRUE)


###################################################
### code chunk number 59: qle_with_R.Rnw:1820-1821
###################################################
print(Stest)


###################################################
### code chunk number 60: qle_with_R.Rnw:1824-1825 (eval = FALSE)
###################################################
## stopCluster(cl)


###################################################
### code chunk number 61: qle_with_R.Rnw:1834-1835
###################################################
diag(attr(Stest,"qi"))^0.5


###################################################
### code chunk number 62: qle_with_R.Rnw:1838-1839
###################################################
sqrt(diag(attr(Stest,"msem")))


###################################################
### code chunk number 63: qle_with_R.Rnw:1843-1844
###################################################
attr(Stest,"msem")


