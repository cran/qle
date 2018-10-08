# Example: apply `qle` to normal model with criterion `score`
# The following code is also part of the vignette
library(qle)

## a local cluster
cl <- makeCluster(8L)
clusterSetRNGStream(cl,1234)

## Multicore parallel processing:
# options(qle.multicore="mclapply")
# options(mc.cores=2L) 

# define a statistical model bysimulation function
simfunc <- function(pars) {	
	x <- rnorm(10,mean=pars["mu"],sd=pars["sigma"])    
	c("T1"=median(x),"T2"=mad(x))	
}

# box contraints defining the parameter space
lb <- c("mu"=0.5,"sigma"=0.1)
ub <- c("mu"=8.0,"sigma"=5.0)	   

## the (unknown) true parameter
theta0 <- c("mu"=2,"sigma"=1)

# simulate model at a minimum of required design points
sim <- simQLdata(sim=simfunc,nsim=10,N=8,
		method="maximinLHS",lb=lb,ub=ub)	 

# set number of simulations manually
# since otherwise only `nsim` would be used to 
# calculate sample average variance
attr(sim,"nsim") <- 100

# true and error-free observation
obs <- structure(c("T1"=2,"T2"=1), class="simQL")

# construct QL approximation model
qsd <- getQLmodel(sim,lb,ub,obs,var.type="wcholMean")

# quasi scoring first try
QS <- qscoring(qsd, x0=c("mu"=5,"sigma"=3.0))
print(QS)

# force only global searches and testing
options(mc.cores=8L)
options(qle.multicore="mclapply")

OPT <- qle(qsd,
		simfunc,		
		nsim=20,
		global.opts=list("maxeval"=50),
		local.opts=list("lam_max"=1e-3,"weights"=0.5,
				"useWeights"=FALSE,"test"=TRUE),cl=cl)

print(OPT)

OPT$final
OPT$why

## testing with criterion `mahal`
## here: no Iobs for best root selection
S0 <- multiSearch(theta0, qsd=OPT$qsd, method=c("bobyqa","cobyla","direct"),
		 nstart=25,	multi.start=TRUE, optInfo=TRUE, pl=10, verbose=TRUE)
 
## found roots are all the same up to 
## numerical precision
attr(S0,"roots")

# compare estimated variance matrix of statistics
obs0 <- simQLdata(simfunc,X=OPT$par,nsim=1000,mode="matrix")[[1]]
var(obs0)
attr(OPT$final,"Sigma")

stopCluster(cl)

#qsd$QS <- QS
#qsd$OPT <- OPT
#qsd$sim <- sim
#qsd$simfn <- simfunc

#save(qsd,file="normal.rda")
