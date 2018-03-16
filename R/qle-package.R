# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the GPL (>=3).
#
# File: 	qle-package.R
# Date:  	27/10/2017
# Author: 	Markus Baaske
# 
# General description of the package and data sets

#'  Simulation-Based Quasi-Likelihood Estimation
#' 
#'  We provide a method for parameter estimation of parametric statistical models which can be at least
#'  simulated and where standard methods, such as maximum likelihood, least squares or Bayesian
#'  algorithms (including MCMC) are not applicable. We follow the \emph{quasi-likelihood} theory [3]
#'  to estimate the unknown model parameter by finding a root of the so-called \dfn{quasi-score} estimating
#'  function. For an overview of our method and further in-depth examples please see the vignette.
#' 
#'  The basic idea is to transform the general parameter estimation problem into a global (black box) optimization problem
#'  (see [1]) with an expensive to evaluate objective function. This function can only be evaluated with substantial random
#'  errors due to the Monte Carlo simulation approach of the statistical model and the interpolation error of the involved
#'  approximating functions. The algorithm sequentially selects new evaluation points (which are the model parameters) for
#'  simulating the statistical model and aims on efficiently exploring the parameter space towards a root of the quasi-score
#'  vector as an estimate of the unknown model parameter by some weighted distance space-filling selection criteria of randomly
#'  generated candidate points.
#' 
#'  The main estimation process can be started by the function \code{\link{qle}} where other functions like, for example,
#'  \code{\link{qscoring}} or \code{\link{searchMinimizer}} search for a root or a local and global minimizer (without sampling new
#'  candidates) of some monitor function to control the estimation procedure.
#'  
#' @docType package
#' @name qle-package
#' 
#' @references  
#'  \enumerate{
#'   \item Baaske, M., Ballani, F., v.d. Boogaart,K.G. (2014). A quasi-likelihood
#'				 approach to parameter estimation for simulatable statistical models.
#'	 			 \emph{Image Analysis & Stereology}, 33(2):107-119.  
#' 	 \item Chiles, J. P., Delfiner, P. (1999). Geostatistics: modelling spatial uncertainty.
#'    		  \emph{J. Wiley & Sons}, New York.
#' 	 \item Heyde, C. C. (1997). Quasi-likelihood and its applications: a general approach
#' 		     to optimal parameter estimation. \emph{Springer}
#'   \item Kleijnen, J. P. C. & Beers, W. C. M. v. (2004). Application-driven sequential designs for simulation experiments:
#'        	Kriging metamodelling. \emph{Journal of the Operational Research Society}, 55(8), 876-883
#'   \item Mardia, K. V. (1996). Kriging and splines with derivative information. \emph{Biometrika}, 83, 207-221
#'   \item McFadden, D. (1989). A Method of Simulated Moments for Estimation of Discrete Response
#' 			 Models without Numerical Integration. \emph{Econometrica}, 57(5), 995-1026.
#'   \item Regis R. G., Shoemaker C. A. (2007). A stochastic radial basis function method for the global
#' 			optimization of expensive functions. \emph{INFORMS Journal on Computing}, 19(4), 497-509.  
#' 	 \item Wackernagel, H. (2003). Multivariate geostatistics. \emph{Springer}, Berlin.
#'   \item Zimmermann, D. L. (1989). Computationally efficient restricted maximum likelihood estimation
#' 			 of generalized covariance functions. \emph{Math. Geol.}. 21, 655-672
#' 	 \item Efron, B. and Tibshirani, R. J. (1993). An Introduction to the Bootstrap, Chapman & Hall, New York.
#'  }
#' 	
#' 
NULL

#' A normal model
#'
#' A statistical model of random numbers
#' 
#' This is a pedagogic example of a simulated data set for quasi-likelihood estimation using
#' normally distributed random numbers. The model outcome is a vector of summary statistics, that is,
#' simply the median and mean average deviation of \code{n=10} random numbers, which is evaluated at the
#' model parameter \eqn{\theta=(\mu,\sigma)} with mean \eqn{\mu} and standard deviation \eqn{\sigma} as
#' the parameters of the normal distribution. We estimate the model parameter given a specific
#' "observation" of those summary statistics. Clearly, maximum likelihood estimation would be the 
#' method of first choice if we had a real sample of observations. However, this example is used to demonstrate
#' the basic workflow of estimating the model parameter. We use this model as a standard example in the package
#' documentation.   
#' 
#' @docType data
#' @keywords datasets
#' @name qsd
#' @usage data(normal)
#' @format A list object named `\code{qsd}` of class \code{\link{QLmodel}} with additional elements
#'    \itemize{
#' 	   \item{simfn}{ simulation function }
#'     \item{sim}{ simulation results at design points, class `\code{simQL}`}
#' 	  }
#' @author M. Baaske
NULL

#' Matern cluster process data 
#' 
#' A data set of quasi-likelihood estimation results of estimating the parameters of a Matern cluster
#' point process model. In the vignette we apply our method to the `\code{redwood}` data set from the
#' package \code{spatstat}.
#' 
#' @docType data
#' @keywords datasets
#' @name matclust
#' @usage data(matclust)
#' @format A list object named `\code{matclust}` which consists of
#'   \itemize{
#' 	   \item{qsd}{ initial quasi-likelihood approximation model}
#'     \item{OPT}{ the results of estimation by \code{\link{qle}}}
#'     \item{Stest}{ score test results }
#'   }
#' 
#' @author M. Baaske
NULL
