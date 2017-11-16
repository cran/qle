# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the L-GPL.
#
# File: 	test_scoring.R
# Date:  	12/04/2017
# Author: 	Markus Baaske
# 
# Testing the construction of the QL model for a normal distribution 

library(qle)

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
sim <- simQLdata(sim=simfunc,nsim=5,N=8,
    method="maximinLHS",lb=lb,ub=ub)	 

# true and error-free observation
obs <- structure(c("T1"=2,"T2"=1), class="simQL")

# construct QL approximation model
qsd <- getQLmodel(sim,lb,ub,obs,var.type="wcholMean")
