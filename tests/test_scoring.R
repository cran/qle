# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the L-GPL.
#
# File: 	test_scoring.R
# Date:  	12/04/2017
# Author: 	Markus Baaske
# 
# Testing quasi Fisher-scoring iteration

library(qle)
data(normal)

# Scoring with average variance approximation
qscoring(qsd,
	x0=c("mu"=3.5,"sigma"=1.5), 
	opts=list("pl"=10,
			  "ftol_stop"=1e-9,					# stopping value
			  "ftol_abs"=1e-3,
			  "score_tol"=1e-6,
			  "grad_tol"=1e-3),					# test for local minimum by < ftol_abs
	W=diag(2), theta=c(2,1), verbose=TRUE)

