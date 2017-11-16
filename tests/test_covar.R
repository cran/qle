# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the L-GPL.
#
# File: 	test_covar.R
# Date:  	12/04/2017
# Author: 	Markus Baaske
# 
# Testing REML estimation 

library(qle)
data(normal)

# default reml optimization controls, see nloptr
attr(qsd,"opts")

# first covariance model
covT <- qsd$covT[1]

# sampled parameters
X <- as.matrix(qsd$qldata[1:2])

# 1st statistic
T <- qsd$qldata[c("mean.T1")]

# reml estimation
fit <- fitCov(covT,X,T,verbose=TRUE)[[1]]

# reml value at fitted parameters
p <- attr(fit,"optres")$solution
reml(covT,p,T,X)
attr(fit,"optres")$objective
