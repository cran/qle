# One-Step Estimation Approach:
# find a suitable starting point
# and then search for a root of the
# quasi-score by fisher scoring.
library(qle)
data(normal)

# starting point
x0 <- c("mu"=3,"sigma"=2.5)

# box constraints for parameters
lower <- qsd$lower
upper <- qsd$upper

## use log average approximation 
## of variance matrix as an interpolation
qsd$var.type <- "logMean"

# direct minimization of Mahalanobis distance
ctls <- list("stopval"=1e-10,
			 "ftol_rel"=1e-6,
			 "maxeval"=1000)

## Using `nloptr`directly,
## though quite slow but possible	 
S0 <- nloptr::direct("mahalDist",
					 lower=lower,
					 upper=upper,
					 control=ctls,
			    qsd=qsd, value.only=TRUE)

print(S0) 

## A least squares approach to find suitable
## starting point, if possible, by a  global
# search strategy  called `direct` (see nloptr::direct). 
W <- diag(1,2)
qsd$criterion <- "mahal"

# one-step minimization
S1 <- searchMinimizer(x0, qsd, W=W, method="direct",
		control=list("stopval"=1e-3), verbose=TRUE)

# results
print(S1)

## now apply quasi-likelihood with quasi-scoring   
qsd$criterion <- "qle"
qscoring(qsd, S1$par)

