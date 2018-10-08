# Copyright (C) 2018 Markus Baaske. All Rights Reserved.
# This code is published under the GPL (>=3).
#
# File: 	covariance.R
# Date:  	14/03/2018
# Author: 	Markus Baaske
#
# Define covariance structures: SIRF-k and Matern

#' @name
#'  setCovModel
#' 
#' @title 	
#' 	Set a covariance model
#' 
#' @description 
#' 	Set a covariance model for kriging the sample means of the involved statistics or for the variance matrix of the statistics.
#' 
#' @param model			name of covariance model: `\code{sirfk}` (default), `\code{matern}`, `\code{powexp}`, `\code{exp}`
#' @param param			numeric vector, \code{NULL} (default), starting values of covariance parameters for estimation 
#' @param npoints		number of sample points already evaluated for covariance parameter estimation
#' @param var.sim		numeric vector, \code{NULL} (default), local simulation variances (as local nugget variances)
#' @param nugget		starting value for (global nugget) variance estimation
#' @param trend			integer, \code{=2} (default) number of polynomial trend order: either set to linear (=1) or quadratic (=2) 					    
#' @param fixed.param	vector of names, corresponding to `\code{param}` of covariance parameters, which will be hold fixed for covariance parameter estimation by REML
#' @param lower			lower bounds of covariance parameters for REML estimation 
#' @param upper			upper bounds of covariance parameters for REML estimation
#' @param ...			additional arguments which can be stored
#' 
#' @return Object of class \code{covModel} as a list of the following objects
#'  \item{model}{ integer number of covariance function}
#'  \item{param}{ estimtated covariance parameter vector}
#' 	\item{start}{ start point for REML estimation of covariance parameters}
#'  \item{trend}{ trend order number}
#'  \item{fix.nugget}{ vector of (fixed) values used as local nugget variances}
#'  \item{free}{ index of free parameters for REML estimation}
#' 	\item{lower}{ lower bounds of covariance parameters for REML estimation}
#'  \item{upper}{ upper bounds of covariance parameters for REML estimation}
#'  \item{...}{ additional objects which can be stored}
#' 
#' @details The function defines a covariance model for the kriging approximation of the sample mean values of a summary statistic. The covariance model
#'  (which might include a polynomial trend) defines the spatial dependence between different locations (points) of the parameter space. Currently, the
#'  function provides the generalized covariance models (`\code{sirfk}`, see \code{\link{fitSIRFk}}) of order \eqn{k=1,2} 
#'  and the Mat\eqn{\textrm{\'{e}}}rn covariance model with scale (i.e. sill) parameter `\code{scale}`, smoothness parameter
#'  `\code{alpha}`, respectively, `\code{nu}`, and the range parameter `\code{rho}` defined only for the latter and the (power)
#'  exponential covariance model `\code{powexp}`.  
#'  
#'  \subsection{Use of simulation variance}{
#'  If a vector of simulation variances is statically set by `\code{var.sim}` for each location these are used as (local)
#'  nugget variance estimations which account for the sampling variability due to the repeated measurements of the statistics
#'  by simulations. The length should match the number of locations `\code{npoints}` otherwise the given vector components
#'  are recycled to the number of `\code{npoints}`. A global nugget value, which captures the variance of the underlying random function,
#'  could be set by `\code{nugget}` as a starting value for the REML estimation procedure for covariance estimation. Clearly, both types
#'  of nugget variances have a direct influence on the REML estimates in terms of smoothness and goodness-of-fit. 
#'  } 
#' 
#'  \subsection{Default parameters}{ 
#'  The default starting parameters are set to \code{("scale"=0.001,"alpha"=1.5)} for the `\code{sirfk}` model. The
#'  Mat\eqn{\textrm{\'{e}}}rn model uses the following parameters \code{("scale"=1.0,"nu"=2.5,"rho"=3.0)}. The default
#'  parameters for the power exponential covariance model are \code{"scale"=1.0}, (isotropic) range parameter
#'  \code{"phi"=1.0} and power \code{"kappa"=1.5} with \eqn{0<\kappa\leq 2}. 
#'  The corresponding lower and upper bounds are chosen such that the underlying random function remains
#'  twice continuously differentiable. Further, setting the names of the covariance parameters in `\code{fixed.param}`,
#'  excludes these parameters from subsequent REML estimations such that these are hold fixed and used as given in the
#'  starting parameter. 
#' 
#'  The above settings are applicable for a wide range of statistics but, however, generally depend on the kind of statistics to be interpolated
#'  and thus have to be chosen carefully. Note that a valid (generalized) covariance model for kriging requires at least \eqn{q+2} design points
#'  for the trend order \eqn{k=1} and \eqn{1+(q+1)(q+2)/2} for \eqn{k=2} where \eqn{q} is the dimension of the unknown model parameter `\code{param}`.
#'  }   
#'  
#' @examples 
#'  # set the standards sirf-2 covariance model
#'  setCovModel("sirfk",npoints=12)
#' 
#' @author M. Baaske
#' @rdname setCovModel
#' @export
setCovModel <- function(model = "sirfk", param = NULL, npoints = 0, var.sim = NULL,
						  nugget = 1.5e-4, trend = 2, fixed.param = NULL, 
						   lower = NULL, upper = NULL,...)
{	
	stopifnot(npoints>0)
	trend.nr <- pmatch(trend,c(1,2))
	if(anyNA(trend.nr)) {
	  msg <- paste0("Invalid trend order. Either choose `trend` as 1 (linear) or 2 (quadratic).")
	  message(msg)
	  stop(msg)
  	}
 	
    cov.list <- c("sirfk","matern","powexp","exp")
 	cov.nr <- pmatch(model, cov.list)
	if (is.na(cov.nr)) {
		msg <- paste("Unknown covariance model. Use one of ", paste(cov.list,collapse = ","))
		message(msg)
		stop(msg)
	}		
	
	fix.nugget <- NULL	
	if(!is.null(var.sim)) {		
	  if(!is.numeric(var.sim) || is.matrix(var.sim))
		 stop("`var.sim` has to be a numeric vector.")
	  fix.nugget <-
		# treated as local/global nugget
		if(length(var.sim) != npoints)
			as.numeric(unlist(lapply(list(var.sim), rep, length.out = npoints)))			
		else var.sim	 				
	}
	
	switch(model,
		"sirfk" = {			
			i <- min(which(trend < c(1,2,3), arr.ind=TRUE))
			param <-
			 if(is.null(param))
			  c("scale"=0.001,"alpha"=1.5,"nugget"=nugget)
			 else {			   
			   if(anyNA(match(names(param),c("scale","alpha"))))
				 stop("Invalid parameter vector for model `sirfk`.")
			   p <- c("scale"=0.001,"alpha"=1.5)  
			   p[names(param)] <- param 
			   c(p,"nugget"=nugget)
			 }
		 	
	        if(is.null(lower)) {
			  lower <- rep(-Inf,length(param))	
			}
			lower <-
			  if(anyNA(lower) | any(!is.finite(lower)))
				  c(1e-6,1.0,1e-4) else lower
			upper <-
			  if(is.null(upper)) {			  
				 c(5,i-1e-4,3)						
			  } else c(upper[1],min(upper[2],i-1e-4),upper[3])			 				
		},		
		"matern" = {
			param <-
			 if(is.null(param))
				c("scale"=1.0,"nu"=2.5,"rho"=3.0,"nugget"=nugget)
			 else {
				 if(anyNA(match(names(param),c("scale","nu","rho"))))
					 stop("Invalid parameter vector for model `matern`.")
				 p <- c("scale"=1.0,"nu"=2.5,"rho"=3.0)  
				 p[names(param)] <- param
			     c(p,"nugget"=nugget)
			 }
			if(is.null(lower) || is.null(upper)) {
				lower <- c(1e-6,2+1e-10,0.1,1e-4)
				upper <- c(10,3,10,10)
			}			
		},
		"powexp" = {
			param <-
			 if(is.null(param))
				param <- c("scale"=1.0,"phi"=1.0,"kappa"=1.5,"nugget"=nugget)
			 else {
				 if(anyNA(match(names(param),c("scale","phi","kappa"))))
					 stop("Invalid parameter vector for model `powexp`.")
				 p <- c("scale"=1.0,"phi"=1.0,"kappa"=1.5)  
				 p[names(param)] <- param				 
				 c(p,"nugget"=nugget)
			 }
	 
			if(is.null(lower) || is.null(upper)) {
				lower <- c(1e-6,1e-6,1e-4,1e-4)
				upper <- c(3,3,1.9999,10)
			}			
		},
		"exp" = {
			param <-
			if(is.null(param))
				param <- c("scale"=1.0,"phi"=1.0,"nugget"=nugget)
			else {
				if(anyNA(match(names(param),c("scale","phi"))))
					stop("Invalid parameter vector for model `exp`.")
				p <- c("scale"=1.0,"phi"=1.0)  
				p[names(param)] <- param				 
				c(p,"nugget"=nugget)
			}
			
			if(is.null(lower) || is.null(upper)) {
				lower <- c(1e-6,1e-6,1e-4)
				upper <- c(3,3,10)
			}			
		}
	)	
	if(length(param) != length(lower) ||
	   length(param) != length(upper)) {
	  stop("`lower` or `upper` bound lengths must be equal to the length of the covariance parameter vector.")
	}
	start <- param
	if(any(start<lower) || any(start>upper))
		stop("At least one parameter does not match constraints. Check `lower` and `upper`.")
	
			
	free <-  
	 if(is.null(fixed.param)){
		 seq(length(start))
	 } else {
		id <- which(is.na(match(names(start),fixed.param)))
		if(length(id) == 0) {
		  cat("All parameters are excluded estimation.")
		} else {
		  start <- start[id]	    # new start
		  lower <- lower[id]		# new bounds
		  upper <- upper[id]		  
		} 
		id
	 }  	
		
	structure(list("model"= cov.nr,
				   "param"= param,
				   "start"= start,				  				   				   
				   "trend"= trend,
				   "fix.nugget"= fix.nugget,	
				   "free"= free,
				   "lower"= lower,
				   "upper"= upper,...),
			class="covModel"
	)
	
}

#' @name 		reml
#' 
#' @title 		Restricted maximum likelihood (REML)
#' 
#' @description Compute the value of the REML function (without constant term) for a given covariance model and data 
#' 
#' @param models   	 object of class \code{krige} (list of covariance models) or class
#'					 	 \code{covModel} (a single covariance model), see \code{\link{setCovModel}}
#' @param pars 	     covariance parameter vector (including global scalar nugget value)
#' @param data	  	 data frame of simulated statistics, each column corresponds to a 
#' 					 single covariance model in `\code{models}`
#' @param Xs	  	 matrix of sample points
#' @param verbose	 logical, if \code{TRUE}, print intermediate output
#' 
#' @return List of REML function values.
#' 
#' @details Given a list of covariance models the function calculates the REML function values at
#'   the covariance parameter `\code{pars}`.
#' 
#' @examples
#' 	data(normal)
#' 
#'  # extract the sample points (the design)
#'  X <- as.matrix(qsd$qldata[1:2])
#' 
#'  # get the observed statistic
#'  T <- qsd$qldata[c("mean.T1")]
#'  reml(qsd$covT[1],pars=c(1e-4,1.5,0.1),T,X)
#' 
#' @author M. Baaske
#' @rdname reml
#' @export
reml <- function(models, pars, data, Xs, verbose = FALSE) {	
	stopifnot(is.data.frame(data))
	stopifnot(length(models) == ncol(data))
	
	lapply(1:length(models),
		function(i) {
			tryCatch({		
				 F <- .Call(C_Fmat,Xs,models[[i]]$trend)	
				 P <- .Call(C_Pmat,F)	
				 y <- crossprod(P,data[[i]]) # t(P)%*%z
				 fnREML(pars,y,Xs,P,models[[i]],verbose=verbose)
				}
			  ,error = function(e) {
				  msg <- .makeMessage("REML function evaluation failed: ",conditionMessage(e))
			  	  message(msg) 
				  structure( list( message=msg,  call = sys.call() ), class=c("error","condition"), error=e )
			 	}	  
			)			
		})
}

# intern, covModel has to be initialized! 
fnREML <- function(p, y, Xs, P, model, free = seq(length(p)), verbose = FALSE)
{	
	 tryCatch({			
		model$param[free] <- p		
		# covariance matrix
		Cmat <- .Call(C_covMatrix,Xs,model)		
		if (!is.matrix(Cmat) || anyNA(Cmat) || !is.numeric(Cmat))
		  stop("Error in covariance matrix: non-numeric values.")								
			
	 	msg <- NULL	
	 	W <- crossprod(P, Cmat %*% P)
		rc <- rcond(W)
				
	    if(rc > 1e-3){
			Wc <- try(chol(W),silent=TRUE)
			if(!inherits(Wc,"try-error") && all(diag(Wc)>0) ){
				w  <- backsolve(Wc,y,transpose=TRUE)	
				if(!inherits(w,"try-error")){
				  return (
				   structure(
					as.numeric(.5*(sum(w^2)) + sum(log(diag(Wc)))),  # beware of brackets: 0.5*2*sum(log(diag(Wc)))
				  	 "info"=list("rcond"=rc,"msg"=msg, "p"=p)) 
	              )
			  	} else {
					if(verbose) 
					 warning("Could not backsolve in `fnREML`.")
				}
			} 
		} else {		  
			msg <- paste0(c("reciprocal condition number (", rc,") of projected covariance matrix `W` is near zero at covariance parameter \n\t ",
					format(p, digits=6, justify="right"),"\n"," which might be unreliable."), collapse = " ")
			if(verbose)
				warning(msg)
		}				
		z <- try(gsiSolve(W,y,use.solve=FALSE),silent=TRUE)
		if(inherits(z,"try-error"))
		  stop("`gsiSolve` failed. Cannot continue REML estimation.")
	  	detW <- try(det(W),silent=TRUE)
		if(inherits(detW,"try-error") || detW < 0)
		 stop(.makeMessage("Could not compute determinant: ",detW," or negative value. \n"))
	 	if(verbose && detW < 1e-17)
		 warning("Determinant of projection matrix `W` (REML) is near zero.\n")
		
		return(	
		    structure(
				as.numeric( .5*(t(y) %*% z + log(detW)) ),
				 "info"=list("rcond"=rc,"msg"=msg, "p"=p) )	
		)		
	  
	} ,error = function(e) {
		 stop(.makeMessage("Error in function 'fnREML': ",	conditionMessage(e)))		  	
		}
	)		
	
}  

## TODO: Numerical derivatives of log likelihood
##  -> use gradient of loglik
## REML ->alpha
## alpha <- 1.5
## phi <- 1.0
## D <- as.matrix(dist(X))
## diag(D)<- (fix.nugget+nugget)
## H <- log(phi*D)
## exp(2*alpha*H)*2*H
#' @importFrom nloptr nl.grad
fnGradREML <- function(p, y, Xs, P, model, free = NULL, verbose = FALSE) {
	list("objective"=fnREML(p,y,Xs,P,model,free,verbose),
		 "gradient"=nloptr::nl.grad(p, fnREML, heps = .Machine$double.eps^(1/3),y,Xs,P,model,free)) 
}  

## TODO add data as parameter
# arguments '...' manipulate global options for nloptr 
#' @importFrom nloptr nl.opts
doREMLfit <- function(model, Xs, opts, verbose = FALSE )
{
	# return if all parameters are fixed
	if(!is.null(model$free) && length(model$free)==0L) {
	  return(
		structure(
		  list(model = model,convergence = 1L),
		optres=NULL, class = "reml") )
	}	
	err <- NULL	
	fn <- fnREML
	if(!is.null(opts$algorithm) &&
	   opts$algorithm == "NLOPT_LD_LBFGS") {
	   message("Caution: using `LBFGS` in REML function might not be stable.")
	   fn <- fnGradREML
 	}
	 
 	tryCatch({		
        nms <- names(model$param)				
		Fmat <- .Call(C_Fmat,Xs,model$trend)		
		P <- .Call(C_Pmat,Fmat)	
		
		# transform to mean zero data				
		y <- crossprod(P,as.numeric(model$dataT))
		
		model$dataT <- NULL		
		p0 <- .PROJMED(model$start,model$lower,model$upper)
		
		res <- nloptr::nloptr(p0, fn, lb = model$lower, ub = model$upper, opts = opts,
						y = y, Xs = Xs, P = P, model = model, free = model$free,
						 verbose = verbose)
		msg <- "Normal convergence."
		if(inherits(res,"error") || is.null(res) || anyNA(res$solution)){
			msg <- .makeMessage("Function call to 'nloptr' failed.")				
			message(msg)
			return(.qleError(message=msg,call=match.call(),error=res))
		}
		# do a final local search
		if(!is.null(opts$local_opts)){
			if(length(opts$local_opts) > 0L) {		
				locopts <- nloptr::nl.opts()
				locopts[names(opts$local_opts)] <- opts$local_opts 
			} else {
				locopts <- list("algorithm" = "NLOPT_LN_COBYLA","ftol_rel" = 1.0e-7,
								"xtol_rel" = 1.0e-6, "maxeval" = 100)
			}	
			if(verbose)
			  message("Do a final local search of covariance parameters.")
			res0 <- nloptr::nloptr(res$solution, fn, lb = model$lower, ub = model$upper, opts = locopts,
					 y = y, Xs = Xs, P = P, model = model, free = model$free,
					  verbose = verbose)
			if(inherits(res0,"error") || is.null(res0) || anyNA(res0$solution)){
			   warning(.makeMessage("Local function call to 'nloptr' failed after global optimization."))
			   res$final <- res0
		    } else res <- res0
		}
		
		converged <- FALSE
		sol <- res$solution			
		if(res$status >= 0L) {
			converged <- TRUE
			model$param[model$free] <- sol						   			   
	    } else {
		   verbose <- TRUE
		   msg <- .makeMessage("Estimation of covariance parameters did not converge.")		   
		   message(msg)
	    }		
	    structure(
		    list(model=model,
				 convergence=converged,
				 message=msg),
		  optres = if(verbose) res else NULL, class = "reml") 
				 
	 }, error = function(e){
			 msg <- .makeMessage("Nloptr error fitting covariance parameters: ",
					  conditionMessage(e))
			 message(msg)
			 structure(
				list(model=model,
					 convergence=FALSE,
					 message=msg,
					 call=sys.call()),
			  error=e)			  		
		}
	) 	 
}


#' @name fitCov
#' 
#' @title Covariance parameter estimation
#' 
#' @description The function estimates the (hyper)parameters of a list of covariance models `\code{models}` by
#' 	  the \emph{Restricted Maximum Likelihood} (REML) estimation method.
#' 
#' @param models  	 object either of class \code{krige}, a list of covariance models or an object of
#' 	 				 class \code{covModel} (a single covariance model)
#' @param Xs	 	 matrix of sample points (design points)
#' @param data		 data frame of simulated sample means of statistics
#' 					 first column corresponds to the first model in the list `\code{models}` and so forth
#' @param controls	 list of control parameters, see \code{\link[nloptr]{nloptr}}
#' @param cl		 cluster object, \code{NULL} (default), of class "\code{MPIcluster}", "\code{SOCKcluster}", "\code{cluster}"
#' @param verbose 	 logical, \code{FALSE} (default) for intermediate output
#' 
#' @return An object of class \code{reml} which consists of a list of named lists
#'  (of elements `\code{model}` and `\code{convergence}`) each storing a fitted covariance model
#'  together with optimization results from a call to \code{\link[nloptr]{nloptr}} as an attribute
#'  named `\code{optres}` if \code{verbose=TRUE}. The default method for estimation is \code{\link[nloptr]{mlsl}} which
#'  uses random starting points and thus produces different results if it is run more than onces. If the results strongly vary,
#'  then the corresponding REML function might have many local minima which precludes the use of this default algorithm and another
#'  one, e.g. `\code{NLOPT_GN_DIRECT}` (see \code{\link[nloptr]{nloptr.print.options}}), might lead to better results. 
#' 
#' @details The function fits a list of covariance models using the REML method. In order to avoid singularities
#'  of the so-called trend matrices make sure to use at least the minimum required number of sample points given by
#'  `\code{Xs}` which depends on trend order, see \code{\link{setCovModel}}. THe use is given an advice if the trend order does
#'  not match the required number of (initial) design points.   
#' 
#' @examples 
#' data(normal)  
#' 
#' # fit 1st statistic and get REML results
#' fitCov(qsd$covT[1],
#'        Xs=as.matrix(qsd$qldata[1:2]),
#'        data=qsd$qldata["mean.T1"],verbose=TRUE)
#'   
#' @seealso \code{\link{setCovModel}} 
#' 
#' @author M. Baaske
#' @rdname fitCov 
#' @export 
fitCov <- function(models, Xs, data, controls = list(),
			     	  cl = NULL, verbose = FALSE) {
		
	if(!is.data.frame(data))
		stop("Expected argument `data` of class `data.frame`.")	
	if(!is.matrix(Xs))
		stop("Expected argument `Xs` to be  a matrix of sample locations.")
			
	if(length(controls)>0L) {		
		opts <- nloptr::nl.opts()
		opts[names(controls)] <- controls
	} else {
		opts <- list("algorithm" = "NLOPT_GN_MLSL",
				"local_opts" = list("algorithm" = "NLOPT_LN_COBYLA","ftol_rel" = 1.0e-6,
						"xtol_rel" = 1.0e-6,"maxeval" = 1000),
				"maxeval" = 200, "xtol_rel" = 1.0e-6, "ftol_rel" = 1.0e-6, "population"=0)	
	}	
	for(i in 1:length(models))
	 models[[i]]$dataT <- as.numeric(data[[i]])
		
 	mods <- doInParallel(models, doREMLfit, Xs=Xs, opts = opts,
			  cl=cl, verbose=verbose)
	  
	if(inherits(mods,"error")) {
		msg <- paste0("REML estimation failed: ",conditionMessage(mods),"\n")
		message(msg)
		return(.qleError(message=msg,
				call=match.call(),error=mods))
	}	
	errId <- which(sapply(mods,function(x) .isError(x)))
	if(verbose) {	  
	  if(any(errId))
		message(paste(c("Failed fitting covariance models with index: ",as.character(errId)), collapse=" ")) 
	  else {
		id <- which(sapply(mods,function(x) x$convergence))
		if(!all(id)) {
		 message(paste(c("REML failed to converge: ",as.character(id)), collapse=" ")) 
		} else
		 message("Successfully fitted covariance parameters.","\n")		
	  }
	}	
	structure(mods,
		opts = opts,
		error = if(length(errId) > 0L) errId else NULL,
		class = "QLFit")
}

#' @name QLmodel
#' 
#' @title Construct the quasi-likelihood approximation model
#' 
#' @description Aggregate and construct the data for quasi-likelihood estimation
#' 
#' @param qldata		data frame of (initial) simulation results (see \code{\link{setQLdata}})
#' @param lb		    numeric vector of lower bounds defining the (hyper)box
#' @param ub 			numeric vector of upper bounds defining the (hyper)box
#' @param obs	    	numeric vector of observed statistics
#' @param mods			list of (fitted) covariance models (see \code{\link{fitSIRFk}}) 
#' @param nfit			number of cycles, \code{nfit=1} (default), after which covariance
#' 						parameters are re-estimated and otherwise only re-used 
#' @param cv.fit 		logical, \code{TRUE} (default), whether to re-fit CV models (re-estimate covariance parameters)	
#' @param var.type  	name of the variance approximation method (see \code{\link{covarTx}})
#' @param useVar    	logical, \code{TRUE} (default), whether to use prediction variances (see details)
#' @param criterion 	global criterion function for sampling and minimization, either "\code{qle}" or "\code{mahal}"				    	
#' @param verbose       logical, \code{FALSE} (default), whether to give intermediate output 
#' 
#' @return An object of class \code{\link{QLmodel}} which stores the data frame of simulation results, bounds on
#'  the parameter space, covariance models for kriging, vector of observed statistics as well as options for
#'  kriging and covariance parameter estimation. 
#' 
#' @details The function aggregates all required information for quasi-likelihood estimation and defines the input object to the 
#'   function \code{\link{qle}}, stores the fitted covariance models of the sample means of the statistics and the type of variance
#'   matrix approximation. For an advanced setup of the estimation procedure and more involved statistical models this function
#'   explicitly offers the data structure to construct individual covariance models for each statistic, see \code{\link{setCovModel}}.
#'   The user has the choice whether or not to make use of kriging prediction variances by `\code{useVar}` to account for the simulation
#' 	 error when constructing the approximation of the variance matrix and the quasi-score function. If \code{useVar=TRUE}, then a kriging
#'   procedure including the computation of prediction variances based on kriging is automatically used. Otherwise the so-called
#'   \emph{dual} approach is employed which has some computational advantage if prediction variances are not required. 
#' 
#' @examples 
#' 
#' data(normal)
#' 
#' # As an example we re-use the stored normal data and fit 
#' # a generalized covariance model to the data using simulation
#' # variances as local variances for REML estimation.
#' mods <- fitSIRFk(qsd$qldata, verbose=TRUE)
#' 
#' # construct QL approximation model
#' qsd <- QLmodel(qsd$qldata,qsd$lower,qsd$upper,
#' 			    c("T1"=2,"T2"=1),mods)
#' 
#' @seealso \code{\link{getQLmodel}}, \code{\link{updateCovModels}} 
#' 
#' @author M. Baaske
#' @rdname QLmodel
#' @export 
QLmodel <- function(qldata, lb, ub, obs, mods, nfit = 1, cv.fit = TRUE,
		    var.type = c("wcholMean","cholMean","wlogMean","logMean","kriging","const"),
				useVar = TRUE, criterion = c("qle","mahal"), verbose = FALSE)
{	
	if(!inherits(qldata,"QLdata"))
	 stop("expected argument `qldata` of class `QLdata`.")
	if(missing(lb) || missing(ub))
	 stop("Arguments `lb` and `ub` are missing.")
 	
 	dx <- attr(qldata,"xdim")
    if(!is.numeric(lb) || !is.numeric(ub) ||
	   length(lb)!=length(ub) || dx != length(ub))
  	  stop("Dimensions of `lb` or `ub` do not match.")
 	
 	obs <- unlist(obs)
	if(anyNA(obs) | any(!is.finite(obs)))
	 stop("`NA`,`NaN` or `Inf`values detected in argument `obs.")
	if(!is.numeric(obs))
	  stop("Argument `obs` must be a (named) numeric vector or list.")
  	stopifnot(!is.null(mods$covT))
	stopifnot(class(mods)=="QLFit")
	
	if(length(mods$covT) != length(obs))
	  stop("Number of covariance models `covT` and length of observations vector `obs` must equal.")
  	var.type <- match.arg(var.type)
	criterion <- match.arg(criterion)
		
	if(is.null(mods$covL) && var.type == "kriging")
	  stop("Covariance models for variance matrix interpolation must be set for argument \'var.type\'.")
	if(!is.numeric(nfit) || length(nfit)>1L )
	 stop("Argument 'nfit must be numeric of length one.")	
	
 	covT <- .extractCovModels(mods$covT,verbose)
	stopifnot(class(covT)=="krige")
	
	covL <- NULL
	if(!is.null(mods$covL)){		
		covL <- .extractCovModels(mods$covL,verbose)
		stopifnot(class(covL)=="krige")
	}
	# reml optimization options 
	opts <- attr(mods,"opts")
	if(is.null(opts) || length(opts) == 0L){
		opts <- list("algorithm" = "NLOPT_GN_MLSL",
				"local_opts" = list("algorithm" = "NLOPT_LN_COBYLA","ftol_rel" = 1.0e-6,
						"xtol_rel" = 1.0e-6,"maxeval" = 1000),
				"maxeval" = 200, "xtol_rel" = 1.0e-6, "ftol_rel" = 1.0e-6, "population"=0)			  
	}
	# minimum required sample size
	minN <- ifelse(min(sapply(covT,	function(x) x$trend)) < 2, dx+2, (dx+1)*(dx+2)/2+1)
	if(nrow(qldata)<minN) {
	 stop(paste0("Choose the size of the initial sample for parameter dimension ",dx,
			" at least of size: ",minN))
 	}	
	
	structure(
	    list("qldata" = qldata,
			 "lower" = lb,
			 "upper" = ub,
			 "covT" = covT,
			 "covL" = covL,			
			 "obs" = obs,		 
			 "var.type" = var.type,			 
			 "krig.type" = if(useVar) "var" else "dual",
			 "criterion" = criterion,			 
			 "minN" = minN,
			 "nfit" = nfit,
			 "cv.fit"=cv.fit),
	  opts = opts,
	  class="QLmodel"
	)
}

# intern
.extractCovModels <- function(covs, verbose = FALSE) {
	if(is.null(covs))
	  return (NULL)
    # which one procudces errors in fitting?
	errId <- which(sapply(covs,function(x) .isError(x)))
	
	if(length(errId)>0L)
	  message(.makeMessage("A total of ",length(errId)," errors detected in fitted covariance models.")) 
	
    structure(
		lapply(covs,
		 function(x) {
			if(verbose)
			 structure(x$model,"optres"=attr(x,"optres"))
			else x$model
		 }
	    ), error = if(length(errId)>0L) errId else NULL,
	  class = "krige"
 	)	
}

#' @name fitSIRFk
#' 
#' @title Estimation of covariance parameters
#' 
#' @description Fit a generalized covariance model to simulation data
#' 
#' @param qldata		object of class \code{QLdata}, a data frame from function \code{\link{setQLdata}}
#' @param set.var 		logical vector of length one or equal to the number of covariance models;
#' 						for values \code{TRUE} (default), set simulation variances as local nugget variances
#' 						for the corresponding covariance model according to its index 
#' @param var.type      name of variance matrix approximation type (see \code{\link{covarTx}})  
#' @param var.opts	    list of arguments passed to \code{\link{setCovModel}}
#' 						(only if `\code{var.type}`="\code{kriging}" and ignored otherwise)
#' @param intrinsic 	logical vector, \code{FALSE} (default), of length one or equal to the number of Cholesky
#' 					    decompositions of variance matrices; as default use an internal nugget variance estimate (see details)
#' 						for kriging the variance matrix of the statistics
#' @param ...			arguments passed to \code{\link{setCovModel}}
#' @param cl			cluster object, \code{NULL} (default), of class "\code{MPIcluster}", "\code{SOCKcluster}", "\code{cluster}"
#' @param controls		list of control parameters passed to \code{\link[nloptr]{nloptr}} for local minimization
#' @param verbose		if \code{TRUE}, show intermediate output
#' 
#' @return A list of fitted covariance models for kriging the sample means of statistics named `\code{covT}` and optionally
#'  the variance matrix of statistics, `\code{covL}`. The object also stores the reml optimization parameters `\code{controls}`. 
#' 
#' @details The function contructs and estimates the parameters of the covariance models by the REML estimatino method for both kriging
#'   the sample means of the statistics and kriging the variance matrix of statistics unless `\code{var.type}`
#'   equals "\code{const}" for the latter. The default covariance model is derived from a (self-similar) intrinsic random function, that is,
#'   the `\code{sirfk}` function of order \eqn{k} (see, e.g. [1]) with \eqn{k=1,2}, for all statistics (including a default quadratic drift term
#'   \eqn{k=2}). The user can also define different covariance models for each statistic separately (see below). Other covariance models can be set
#'   by their name in the argument `\code{model}` which is passed to the function \code{\link{setCovModel}}. Currently, kriging the variance matrix
#'   is done by the `\code{sirfk}` model. 
#'    		
#'   The argument `\code{var.opts}` only sets the options for the covariance models for kriging the variance matrix if this is the users prefered
#'   type of approximation. Further optional arguments, e.g., `\code{var.sim}` used only for the approximatino of the statistics,
#'   `\code{var.opts$var.sim}` for kriging the variance matrix, specify the local vector of \dfn{nugget} values for each sample point depending on
#'   whether or not `\code{set.var}` (which is only used for kriging the statistics) equals \code{TRUE}. Both arguments are passed to
#'   \code{\link{setCovModel}} and must be data frames of lengths (number of columns) corresponding to the number of covariance
#'   models of the statistics and, respectively, to the number of \emph{Cholesky} decomposed terms in case of kriging the variance matrix.
#'   If `\code{set.var=TRUE}` (default), then local nugget variances are estimated by the variance of the sample average of the simulated values of the statistics.
#'   Otherwise the values given in `\code{var.sim}` are used as fixed `nugget` variances and replicated to match the number of sample points.
#'  
#'   The same applies in case of kriging the variance matrix. If `\code{intrinsic=TRUE}`, then local nugget variances
#'   for each of the variance-covariances of the statistics are estimated by a bootstrapping procedure. Otherwise the values given by
#'   `\code{var.opts$var.sim}` (of length one or equal to the number of corresponding sample points) are used directly as local estimates
#'   (which then must exactly match the order of the Cholesky decomposed terms). A global nugget value can be estimated during the REML
#'   estimation which is the default option for both cases unless this parameter is excluded from the covariance parameter estimation
#'   (see \code{\link{setCovModel}}).
#'  
#'   The default optimization algorithm for estimating the covariance parameters is \code{\link[nloptr]{mlsl}} followed by a final local search using
#'   \code{NLOPT_LN_COBYLA}. Note that in this case the estimated parameters may vary when starting the REML procedure several times since starting
#'   points are randomly chosen for \code{\link[nloptr]{mlsl}}. All options for optimization can be modified by the argument `\code{controls}`.
#' 
#'   Note that the returned object can also be constructed manually and passed as an input argument to
#'   \code{\link{QLmodel}} in case the user prefers to set up each covariance model separately. In this case, first use
#'   \code{\link{setCovModel}} to construct the covariance model, then estimate the parameters by \code{\link{fitCov}} and pass a list of
#'   fitted covariance models to function \code{\link{QLmodel}}. The resulting object is the same as obtained by this function. Please see
#'   the function \code{\link{QLmodel}} for an example. 
#' 
#' @seealso \code{\link{setCovModel}}, \code{\link{fitCov}},  \code{\link{QLmodel}}
#' 
#' @author M. Baaske
#' @rdname fitSIRFk
#' @export
fitSIRFk <- function(qldata, set.var = TRUE, var.type = "wcholMean",
						var.opts = list("var.sim"=1e-6), intrinsic = FALSE, ...,
						 controls = list(), cl = NULL, verbose = FALSE)
{	
	args <- list(...)
	stopifnot(is.data.frame(qldata))
	stopifnot(is.logical(set.var) && is.logical(intrinsic))
	
	if(length(args)>0L) {		 
		nms <- formals(setCovModel)
		.checkfun(setCovModel, c(nms,args))			
 	}
	
	xdim <- attr(qldata,"xdim") 									# dimension of parameter to estimate!
	nsim <- attr(qldata,"nsim")										# number of simulation replications
	Xs <- data.matrix(qldata[seq(xdim)])	
	dataT <- qldata[grep("^mean[.]",names(qldata))]					# simulated statistic data
		
	np <- nrow(Xs)	
	nstat <- ncol(dataT)	
	
	# all cov models are equal	
	useVarSim <- !is.null(args$var.sim)	
	# if not used as a fixed nugget: set.var == FALSE	
	set.var <- rep(set.var,length.out=ncol(dataT))
	dfvar <-
	 if(useVarSim) {
	  if(anyNA(args$var.sim))
		 stop("local nugget variance vector has `Na`s for kriging statistics.")
	  rep(as.data.frame(args$var.sim),length.out=ncol(dataT))		
	} else NULL
	
	covT <- 
		lapply(1:ncol(dataT),
		     function(i){
			   args$var.sim <-
				 if(set.var[i]) {					  
					  qldata[[xdim+nstat+i]]/nsim
				 } else if(useVarSim && any(dfvar[[i]]>0)) {				 
					  dfvar[[i]]										# numeric vector of length equal number of locations
				  } else NULL
				 fnargs <- c(list("dataT"=dataT[[i]],	      		    # temporarly add the data				  				  			  
								  "npoints"=np,
								  "type"="statcov"), args)
				 do.call(setCovModel,fnargs)					   
			}
		)

	 covL <- NULL	 
	 if(var.type == "kriging"){
		 # check input
		 args <- var.opts
		 if(length(args)>0L) {
			 nms <- formals(setCovModel)			
			.checkfun(setCovModel, c(nms,args))				 
	 	 }		 	   		   		 
  		 useVarSim <- !is.null(args$var.sim)
		 Lvec <- qldata[grep("^L+",names(qldata))]		 
		 
		 # individually set intrinsic noise terms as local nugget variances
		 # for each covariance model of Cholesky decomposed terms
		 intrinsic <- rep(intrinsic,length.out = ncol(Lvec))		 
		 dfvar <- 
		   if(useVarSim) {
			   if(anyNA(args$var.sim))
				   stop("local nugget variance vector has `Na`s for kriging variance matrix.")
			   rep(as.data.frame(args$var.sim),length.out=ncol(Lvec))			   
		   } else { NULL }
		 
   		 # find number of additional columns in `qldata`
   		 M <- if(attr(qldata,"Nb")>0) (nstat*(nstat+1))/2 else 0		 
		 # first M columns are Cholesky terms,
	     # 2nd are bootstrap variances
		 covL <- lapply(1:(ncol(Lvec)-M),
				  function(i)  {
				     args$var.sim <-
					 if(intrinsic[i] && M>0) {						 							 
						 if(any(Lvec[[i+M]] < 0) || anyNA(Lvec[[i+M]])){	# stored bootstrap variances
							if(verbose)
							  message("Bootstrap variance has negative values or `Na`s. Try a default nugget variance value.")
						    # set small value anyway
							as.numeric(dfvar[[i]]) 
						 } else {
							 # square root of nugget variance as this corresponds
							 # to a single value of the Cholseky decomposition 
							 sqrt(Lvec[[i+M]])
						 }
					 } else if(useVarSim && any(dfvar[[i]]>0)) {						 
						 as.numeric(dfvar[[i]]) 
					 } else NULL
					 	 	    
					 fnargs <- c(list("dataT"=Lvec[[i]],									  		  
									  "npoints"=np,
									  "type"="kriging"), args)							   	  
					 do.call(setCovModel,fnargs)					 
				  }
		 )	 
	 }		 	 
	 # (default) reml optimization options 
	 if(length(controls) > 0L) {		
		 opts <- nloptr::nl.opts()
		 opts[names(controls)] <- controls
	 } else {
		 opts <- list("algorithm" = "NLOPT_GN_MLSL",
				  "local_opts" = list("algorithm" = "NLOPT_LN_COBYLA","ftol_rel" = 1.0e-6,
						 "xtol_rel" = 1.0e-6,"maxeval" = 1000),
				  "maxeval" = 200, "xtol_rel" = 1.0e-6, "ftol_rel" = 1.0e-6, "population"=0)		  
	 }
	 	 
	 # REML fit covariance models (statistics and variance matrices)
	 mods <- doInParallel(c(covT,covL), doREMLfit, Xs=Xs, opts=opts,
			 	cl=cl, verbose=verbose)
		
	 if(inherits(mods,"error")) {
		msg <- paste0("REML estimation failed: ",conditionMessage(mods),"\n")		
		message(msg)
		return(.qleError(message=msg,
				 call=match.call(),error=mods))
	 }
	 errId <- which(sapply(mods,function(x) .isError(x)))
	 if(any(errId)) {
		 msg <- paste(c("Failed fitting covariance parameters: ",
				   as.character(errId)), collapse=" ")
   		 message(msg)
		 return(.qleError(message=msg,error=mods[errId]))
	 } else {
		 if(verbose)
		   cat("Successfully fitted covariance parameters.\n")		 
	 }
	 ret <- structure(
			  list("covT" = mods[1:nstat],
			  	   "var.type" = var.type),
		     opts = opts,
			 error = if(length(errId)>0L) errId else NULL,
			 class = "QLFit")	
	 
	 if(!is.null(covL))
	  ret$covL <- mods[(nstat+1):length(mods)]
	  	
	 return ( ret )	
}

#' @name getQLmodel
#' 
#' @title Setup the quasi-likelihood approximation model all at once  
#' 
#' @description Initial setup of the quasi-likelihood approximation model
#'
#' @param runs   	object of class \code{simQL}, simulation results from \code{\link{simQLdata}}
#' @param lb		lower bounds defining the (hyper)box of the parameter domain for QL estimation
#' @param ub 		upper bounds defining the (hyper)box of the parameter domain for QL estimation
#' @param obs		numeric vector of observed statistics
#' @param X   		matrix of sample locations (model parameters)
#' @param useVar   	logical, \code{TRUE} (default), whether to use prediction variances of any kind
#' @param criterion name of criterion function to be minimized for QL estimation (see \code{\link{qle}})
#' @param ...		arguments passed to \code{\link{fitSIRFk}}, \code{\link{setQLdata}}, \code{\link{setCovModel}} 
#'  				and \code{\link{QLmodel}} for REML estimation of all covariance models
#'  
#' @return Object of class \code{\link{QLmodel}}
#' 
#' @details The function is a wrapper of \code{\link{simQLdata}}, \code{\link{QLmodel}}, \code{\link{fitSIRFk}}
#'  and thus sets up the quasi-likelihood approximation model all at once.
#' 
#' @examples
#' 
#' data(normal)
#' 
#' # simulate model at a minimum of required design points
#' sim <- simQLdata(sim=qsd$simfn,nsim=5,N=8,
#' 			 method="maximinLHS",lb=qsd$lower,ub=qsd$upper)
#' 	 
#' # true and error-free observation
#' obs <- structure(c("T1"=2,"T2"=1), class="simQL")
#' 
#' # construct QL approximation model
#' qsd <- getQLmodel(sim,qsd$lower,qsd$upper,obs,var.type="wcholMean")
#' 
#' 
#' @seealso \code{\link{simQLdata}}, \code{\link{QLmodel}}, \code{\link{fitSIRFk}} 
#' 
#' @author M. Baaske
#' @rdname getQLmodel
#' @export
getQLmodel <- function(runs, lb, ub, obs, X = NULL, useVar = TRUE, criterion = "qle", ...)
{	
	args <- list(...)
	verbose <- isTRUE(args$verbose)
		
	tryCatch({
        if(.isError(runs))
		  stop("Simulations have errors. Please check the input argument `runs`.")
		if(verbose)
		  cat("Collect data for fitting covariance models of statistics.\n")
	  		  
	  	id <- which(is.na(pmatch(names(args),names(formals(setQLdata)))))
		args.tmp <-
		 if(length(id)>0L){
		   if( (c("Nb") %in% names(args)) && !isTRUE(args$intrinsic) )	 
		     args$Nb <- 0  # no bootstrap anyway if not "intrinsic" equals TRUE
		   args[-id]
	     } else NULL 
 		# construct all data
	    qldata <- do.call(setQLdata,c(list(runs,X),args.tmp))
		if(.isError(qldata))
			return(qldata)
		
		# fitting statistics
		if(verbose)
		 cat("Fitting covariance models...\n")	 	
	    id <- which(is.na(pmatch(names(args),names(c(formals(fitSIRFk),formals(setCovModel))))))
	 	args.tmp <-
		 if(length(id)>0L) {
			args[-id]
		 } else args	    
		
		# fitting			 	
	    mods <- do.call(fitSIRFk,c(list(qldata),args.tmp))	
		if(.isError(mods)) {
			message(.makeMessage("Failed fitting covariance models: ","\n"))			
			return(mods)
		}		  
		if(verbose)
		  cat("Setup QL approximation model...\n")
	  	id <- which(is.na(pmatch(names(args),names(formals(QLmodel)))))
	  	if(length(id) > 0L) {
		  args <- args[-id]		  
	 	}
	    do.call(QLmodel,
			 c(list(qldata,
					lb,ub,
			  	    obs,
				    mods,				    
					useVar=useVar,
				    criterion=criterion),
		     args))
		
	}, error = function(e) {
		msg <- .makeMessage("Failed to setup QL model: ",conditionMessage(e))
		message(msg)
		return(.qleError(message=msg,
				 call=match.call(),error=e))	
	   }
	)
}

#' @name updateCovModels
#' 
#' @title Update covariance models
#' 
#' @description The function updates the current covariance models stored in `\code{qsd}`.
#' 
#' @param qsd			object of class \code{\link{QLmodel}} which is to be updated
#' @param nextData		object of class \code{QLdata} which includes new simulation results
#' @param fit 			logical, if \code{TRUE} (default), re-estimate covariance parameters
#' @param cl			cluster object, \code{NULL} (default), of class "\code{MPIcluster}", "\code{SOCKcluster}", "\code{cluster}"
#' @param controls	    list of control parameters passed to \code{\link[nloptr]{nloptr}}
#' @param verbose 		logical, \code{FALSE} (default), whether to show intermediate output
#' 
#' @return Object of class \code{\link{QLmodel}} as a list of updated covariance models
#' 
#' @details The function updates both, the covariance models for kriging the statistics, and, if applicable,
#'  the ones for kriging the variance matrix of statistics based on the new data given by `\code{nextData}`. In practice, the user hardly
#'  needs to call this function except for empirical studies of how additional sample points might influence the overall predictive
#'  quality of the quasi-score and/or criterion function approximations.
#' 
#'  If `\code{fit=TRUE}`, then the function re-estimates the covariance parameters for each statistic separately
#'  each time a total of `\code{qsd$nfit}` new sample points have been added. Thus, we can choose whether to fit the updated
#'  covariance models (by the REML estimation method) each time, e.g. during the estimation by \code{\link{qle}} if `\code{qsd$nfit}`=1, or after
#'  each 2nd, 3rd, and so on newly added point in order to limit the computational overhead. If bootstrapping was used to estimate the nugget variance
#'  of kriging models of the variance matrix, then these are taken from `\code{nextData}`. 
#' 
#' @examples
#' 
#' data(normal)
#' 
#' # old design
#' X <- as.matrix(qsd$qldata[c(1,2)])
#' 
#' # augment design with two additional points
#' Xnew <- multiDimLHS(N=2,qsd$lower,qsd$upper,X=X,
#'            method="augmentLHS",type="matrix")
#' 
#' # new simulations
#' Xsim <- simQLdata(sim=qsd$simfn,nsim=10,X=Xnew)
#' 
#' # prepare data
#' Xdata <- setQLdata(Xsim,Xnew)
#' 
#' # do not re-estimate covariance parameters
#' qsd2 <- updateCovModels(qsd,Xdata,fit=FALSE) 
#'  
#' @seealso \code{\link{setQLdata}}, \code{\link{simQLdata}}, \code{\link{QLmodel}}
#' 
#' @author M. Baaske
#' @rdname updateCovModels
#' @export
updateCovModels <- function(qsd, nextData, fit = TRUE,
						cl = NULL, controls = list(), verbose=FALSE)
{		
	stopifnot(class(qsd) == "QLmodel")
	stopifnot(inherits(nextData,"QLdata"))	
	
	nnew <- NROW(nextData)	
	xdim <- attr(qsd$qldata,"xdim")
	nsim <- attr(qsd$qldata,"nsim")
	nsim.new <- attr(nextData,"nsim")
	
	nstat <- length(qsd$covT)	
	stid <- (xdim+1):(xdim+nstat)
	vid <- c((xdim+nstat+1):(xdim+2*nstat))
	
	vars <- qsd$qldata[vid]
	vars.new <- nextData[vid]
	
	# combine old and new data and
	# check columns also because L+ might not be given
	if(ncol(nextData) != ncol(qsd$qldata))
	  stop("The number of columns of argument `nextData` does not match with `qldata`.")
	
	# merge to new data, one (sample point) added
	qsd$qldata <- rbind(qsd$qldata,nextData)	
	Xs <- as.matrix(qsd$qldata[seq(xdim)])
	np <- nrow(Xs)
		
	if(length(controls)>0L) {		
		# set default optimization controls
		opts <- nloptr::nl.opts()
		opts[names(controls)] <- controls
	} else {
		# use stored optimization controls
		opts <- attr(qsd,"opts")	
	}		 
	# update function 
	fitit <- (fit && !(nrow(Xs) %% qsd$nfit))
	
	#c(xm$fix.nugget, xm$ptol*data[[i]][-(1:(np-nnew))] )
	.update <- function(covT, data, vars.new=NULL){
		mod <- lapply(1:length(covT),		
				function(i) {		   
					xm <- covT[[i]]
					# set starting point
					xm$start <- xm$param[xm$free]	
					if(!is.null(xm$fix.nugget)) {					  
					  xm$fix.nugget <-
						if(!is.null(vars.new)){							
						  c(xm$fix.nugget,vars.new[[i]])
						} else c(xm$fix.nugget,rep(xm$fix.nugget[1],nnew)) # re-use first 
		   			} # else (not using simulation variance for REML)
					if(fitit) {
						# store data for REML (remove afterwards in 'doREMLfit')
						xm$dataT <- data[[i]]										      			
					}
					xm
				})	
		if(fitit) {  				
			res <- doInParallel(mod, doREMLfit, Xs=Xs, opts=opts,
					 cl=cl, verbose=verbose)		
			if(!inherits(res,"error")) {
				.extractCovModels(res,verbose)
			} else {
				msg <- paste0("Error fitting covariance parameters.")
				message(msg)
				structure("message"=msg,"error"=res)				
			}	
		} else structure(mod, class = "krige")
	}
	
	tryCatch({
	  qsd$covT <- .update(qsd$covT,
			              qsd$qldata[stid],
						  nextData[vid]/nsim.new)
				  
	  # update kriging VARIANCE models
 	  # Cholesky terms are the data
	  if(qsd$var.type == "kriging"){
		if(is.null(qsd$covL))
		   stop("A covariance model for kriging the variance matrix must be defined but is `Null`.")
	    qsd$covL <-
	     if(attr(qsd$qldata,"Nb") > 0){ 
			 # if bootstrrapping nugget variances
			 # and only use simulations variances at new points
			 .update(qsd$covL,
					 qsd$qldata[grep("^L[^b]",names(qsd$qldata))],
					 nextData[grep("^Lb",names(qsd$qldata))])
		 } else {
			 .update(qsd$covL,qsd$qldata[grep("^L[^b]",names(qsd$qldata))],NULL) 
		 }	 
  	  }
 	}, error = function(e) {
	     msg <- .makeMessage("Failed to update covariance models: ",
				 conditionMessage(e))		
		return(.qleError(message=msg,call=match.call(),error=e))
	})

	return( qsd )
}
