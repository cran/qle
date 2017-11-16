# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the GPL (>=3).
#
# File: 	krige.R
# Date:  	20/10/2017
# Author: 	Markus Baaske
#
# Monte Carlo based testing procedure using an approximate
# efficient score statistic (here: quasi-deviance) or a generalized
# least squares criterion as the test statistic for a goodness-of-fit test

# collect test results
.qleTest <- function(obj,alpha = 0.05) {	
	degf <- nrow(obj)
	sb <- attr(obj,"sb")
	Sb <- attr(obj,"Sb")	
	pb <- (1+sum(Sb>=sb))/(length(Sb)+1)  # similar to : 1-ecdf(Sb)(sb)
	qt <- try(quantile(Sb,1-alpha,na.rm=TRUE),silent=TRUE)
	if(inherits(qt,"try-error"))
	 stop("Could not compute quantiles of test statistic.")	
	
	ans <- list()
	ans$param <- cbind(obj[c("par","se","rmse","bias","mean")])
	dimnames(ans$param) <- list(row.names(obj),
			c("Estimate", "Std. Error", "RMSE", "Bias", "Mean"))
	
	tname <- attr(obj,"test")
	nm <- ifelse(tname=="mahal", "LS criterion", "Score-test")
	Stest <- noquote(paste("Monte Carlo Test (",nm,"):\n\n  Df: ",degf,sep=""))	
	
	ans$test <- noquote(cbind(sb, pb))# , ifelse(degf>0,pchisq(sb, degf, lower.tail = FALSE),"****")))	
	dimnames(ans$test) <- list("Test E(Q)=0:  ", c("S > 0", "P-value")) #, "Pr(>Chisq)"))
	ans$Stest <- Stest	
	ans$tname <- tname
	attr(ans$test,"Sb") <- Sb
	attr(ans$test,"qt") <- qt
	attr(ans$test,"passed") <- (sb < qt)
	
	return(
	  structure(ans,
	    class=c("qleTest"))
	)	
}

# choose the best acc. to consistency 
# criteria and the smallest value of quasi-score
.evalBestRoot <- function(dm) {
	stopifnot(is.data.frame(dm))	
	nm <- apply(dm,2,which.min)	
	id <- sapply(1:length(nm),function(x) sum(nm[1:x] == nm[x]))
	maxid <- cbind(nm, id)[,2]
	id <- as.numeric(nm[which(maxid == max(maxid),arr.ind=TRUE)])
	if(length(id) > 1L){
		for(i in id)		  	
			dimnames(dm)[[1]][i] <- paste0(c(row.names(dm)[i],"*"),collapse=" ")	   
		mi <- which.min(dm[id,length(dm)])
		dimnames(dm)[[1]][id[mi]] <- paste0(c(row.names(dm)[id[mi]],"*"),collapse="")
	} else {
		dimnames(dm)[[1]][id] <- paste0(c(row.names(dm)[id],"*"),collapse=" ")
	}	
	return(row.names(dm))
}


.evalRoots <- function(QD, par = NULL) {	
	if(.isError(QD)){	  	
		return(.qleError(message=.makeMessage("Evaluation of roots failed."),
					call=sys.call(),error=QD))
 	}
	if(is.null(par)){
	   par <- try(do.call(rbind,lapply(QD,"[[","par")),silent=TRUE)
	   stopifnot(is.matrix(par))
    } 
	
	X <- try(
	      lapply(QD,
			function(qd) {
				qIinv <- 
					if(!.isError(qd)) {
						try(solve(qd$I),silent=TRUE)
					} else NA	  	
				xdim <- ncol(qd$I)
				if(.isError(qIinv) || !is.numeric(qIinv) || anyNA(qIinv)) {		
					msg <- .makeMessage("Failed to invert quasi-information.")
					message(msg)
				   .qleError(message=msg,call=sys.call(),error=qIinv)
				}
				M <- qIinv%*%qd$Iobs				
				structure(c("Value"=qd$value,
							"Minor"=.isPosDef(qd$Iobs),
						   	"|det|"=abs(1-det(M)),
							"|max|"=max(abs(diag(M)-rep(1,xdim))),
							"|trace|"=abs(1-sum(diag(M))/xdim), 
							"|score_max|"=max(abs(qd$score))))
			}),silent=TRUE)
	
	if(inherits(X,"try-error"))
	  return(.qleError(message=.makeMessage("Could not compute consistency checks."),
				call=sys.call(),error=X))
	
	ok <- which(sapply(X,function(x) !.isError(x) && is.numeric(x)))	
	if(length(ok) == 0L){
		msg <- .makeMessage("Could not compute consistency checks.")
		message(msg)
		return(.qleError(message=msg,call=sys.call(),error=X))
  	}		
	
	try({				
	   dm <- NULL
	   M <- as.data.frame(do.call(rbind,X[ok]))
	   row.names(M) <- if(!is.null(par)) row.names(par)
	   if(nrow(M)>1) {		   
	   	   dm <- cbind(par[ok,,drop=FALSE],M)
	   	   row.names(dm) <- .evalBestRoot(M)		   
	   } else {
		   dm <- cbind(par[ok,,drop=FALSE],M)		   
	   }
 	},silent=TRUE)

	isErr <- which(!(1:length(X) %in% ok))
    structure(dm,		
		X = if(length(isErr)>0L) X[isErr] else NULL,
		error = if(length(isErr)>0L) isErr else NULL)	
}


#' @name		checkMultRoot
#' 	
#' @title       Inspect estimated parameters
#' 
#' @description		Check out and compare estimated roots of the quasi-score vector
#' 
#' @param est 		object of class \code{qle}, the estimation results from function \code{\link{qle}}
#' @param par	    list or matrix of estimated parameters as roots of the quasi-score vector
#' @param verbose   logical, \code{TRUE} for intermediate output
#' 
#' @return A data frame with columns named corresponding to each component of the investigated parameter,
#'  `\code{quasi-deviance}`, `\code{Minor}`, `\code{det}`, `\code{max}`,
#'  `\code{trace}` and `\code{score_max}` (see vignette). The second column shows the leading minor of
#'   the observed QI matrix which is not positive definite. If so, then the corresponding parameter estimate
#'   cannot be consistent at least in theory. The rows show the corresponding values for each parameter passed by
#'   `\code{par}`. 
#' 
#' @details Only for the quasi-likelihood approach the function inspects the (numerical) consistency of the found
#'  parameters in `\code{par}` by comparing each observed quasi-information matrix with the expected one.
#'  The degree of dissimilarity of both matrices is measured by certain scalar equivalent criteria (see vignette)
#'  and the parameter for which these are smallest is chosen.
#'  
#' @examples 
#'  data(normal)
#'  # a dummy estimation result
#'  OPT <- qle(qsd,qsd$sim,global.opts=list("maxeval"=0))
#'  
#'  # and just a single root 
#'  checkMultRoot(OPT,verbose=TRUE)
#' 
#' @author M. Baaske
#' @rdname checkMultRoot
#' @export
checkMultRoot <- function(est, par = NULL, verbose = FALSE){			
   if(est$qsd$criterion != "qle")
	  stop("Consistency check of multiple roots only for criterion `qle`.")
   if(.isError(est))
	  stop("The estimation result from function `qle` has errors. Please check the argument `est`.")
  
   # always use estimate from est first
   if(!is.null(par)){
	   par <- .LIST2ROW(par)
	   stopifnot(NCOL(par)==length(est$par))	   
   }
   par <- rbind("par"=est$par,par) 
   info <- attr(est,"optInfo")  			   
   QD <- try(quasiDeviance(par,est$qsd,W=info$W,
			  theta=info$theta,cvm=est$cvm,verbose=verbose),
		  silent=TRUE)
   
   if(.isError(QD)) {
	   msg <- paste0("Failed to get quasi-deviance: ","\n")
	   message(msg)
	   return(.qleError(message=msg,call=match.call(),error=QD))
   }
   dm <- try(.evalRoots(QD,par),silent=TRUE)	
   if(.isError(dm)) {
	   msg <- .makeMessage("Could not check consistency of roots.")
	   message(msg)
	   return(.qleError(message=msg,call=match.call(),error=dm))	   
   }
   return (dm)  
}

# intern use only!
.rootTest <- function(par, value, I, obs, alpha, test, ..., cl = NULL, na.rm = TRUE){	
	aiqm <- NULL
	mScore <- NULL	
	xdim <- length(par)
	hasError <- integer(0)
	  	  
	# re-estimate roots
	opt.args <- list(...)
	RES <- do.call(doInParallel,
			c(list(X=obs,
				FUN=function(obs,...)	{
					searchMinimizer(x0=par,...,obs=obs,check=FALSE)
				},
				cl=cl), opt.args)) 

	if(.isError(RES))
	  return(RES)
	
    # check results again
	ok <- which(sapply(RES,function(x) !.isError(x) && x$convergence>0))
	if(length(ok) == 0L){
		stop(.makeMessage("All re-estimations have errors or did not converge, first errors: ","\n"))						
	} else if(length(ok) < length(RES)){
		message(paste0("A total of ",length(RES)-length(ok)," re-estimations have errors or did not converge."))							
	}	
	# average inverse QI
	invI <- 
		lapply(RES[ok],
			function(x) {						
				try(solve(x$I),silent=TRUE) 
			})
	badInv <- sapply(invI,function(x) inherits(x,"try-error") || anyNA(x))
	if(any(badInv))
	 message(paste0("A total of ",sum(badInv)," inversions of quasi-information matrices failed. Check attribute `info`."))    
	else if(all(badInv)){
	 msg <- paste0("All iversions of quasi-information matrices failed.")
	 message(msg)
	 return(.qleError(message=msg,call=match.call(),error=badInv))	 
	}
	# average matrix of inverse qi matrices
	aiqm <- matrix(
			colMeans(
				do.call(rbind,
					lapply(invI[!badInv],as.numeric)
				)
			),ncol=xdim)				
	
	# estimates	
	mpars <- do.call(rbind,lapply(RES[ok],function(x) x$par))
	mScore <- do.call(rbind,lapply(RES[ok],function(x) x$score))
	has.na <- (rowSums(is.na(cbind(mScore,mpars))) > 0L)	
	
	if(na.rm && any(has.na)) {
		ok <- ok[-which(has.na)]
		mpars <- mpars[ok,,drop=FALSE]
		mScore <- mScore[ok,,drop=FALSE]
		warning("Removed `NA` values from quasi-scores.")
	}
	# average quasi-score
	mScore <- try(colMeans(mScore),silent=TRUE)
	# some (empirical) measures	
	msem <- .MSE(mpars,par)	
	# value of test statistic at re-estimated parameters			
	tvals <- sapply(RES[ok],"[[","value") 
	stopifnot(is.numeric(tvals))
	
	# invert QI for predicted std. error (asymptotic) at estimated theta 
	qi <- try(solve(I),silent=TRUE)
	if(inherits(qi,"try-error") || anyNA(qi))
		message("Inversion of quasi-information matrix failed")
	
	# get efficient score test (with MC parameters)
	B <- structure(
			data.frame(
					cbind("par"=par,
						  "se"=apply(mpars,2,sd),					
						  "rmse"=sqrt(diag(msem)),
						  "bias"=colMeans(t(t(mpars)-par)),
						  "mean"=colMeans(mpars))),
			"sb"=value, "Sb"=tvals,
			"test"=test)
	
	relED <-
	 if(!anyNA(c(msem,qi)) && is.matrix(qi) && is.matrix(msem)) {
		 abs(1-sqrt(diag(msem))/sqrt(diag(qi)))
	} else {
		message("Failed to compute relative difference of empirical and predicted error.")
		NULL
	}
	
	# had errors
	hasError <- which(!(1:length(RES) %in% ok))
	if(length(hasError) > 0L)	
		message(paste0("A total of ",length(hasError)," re-estimations failed."))
		
	# results
	structure(.qleTest(B,alpha),				# test results
			msem=msem,							# mean square error matrix
			aiqm=aiqm,							# average inverse QI (re-estimated parameters)
			qi=qi,								# inverse QI at estimated theta
			relED=relED,
			obs=NULL,							# (MC) observations
			mean.score=mScore,					# average score/gradient
			criterion=test,
		info=list(badInv=which(badInv), 			# inversion errors
					hasNa=which(has.na), 			# indices of NA parameters 
					hasError=hasError,
					iseed=NULL),																
	class=c("qleTest"))	
}

#' @name 	qleTest
#'
#' @title	Monte Carlo testing
#'
#' @description Monte Carlo hypothesis testing 
#'
#' @param est			object of class \code{qle}, the estimation results from function \code{\link{qle}}
#' @param local	    	optional, object of class \code{QSResult}, \code{NULL} (default), local estimation results 
#' @param sim			user supplied simulation function (see \code{\link{qle}})
#' @param ...			arguments passed to the simulation function `\code{sim}`
#' @param nsim			number of model replications to generate the simulated statistics
#' @param obs			optional, \code{NULL} (default), simulated statistics at the hypothesized parameter
#' @param check.root    logical, \code{FALSE} (default), whether to check consistency of estimated parameter (see \code{\link{checkMultRoot}})  						
#' @param alpha			significance level for testing the hypothesis
#' @param na.rm 		logical, \code{TRUE}  (default), whether to remove `NA` values from the matrix of
#' 						re-estimated parameters
#' @param cl			cluster object, \code{NULL} (default), see \code{\link[parallel]{makeCluster}}
#' @param iseed			integer, the seed for initializing the cluster workers for parallel computations 
#' @param verbose   	logical, \code{TRUE} for intermediate output
#'
#' @return An object of class \code{qleTest} as a list of:
#'  \item{param}{ data frame of estimated parameters and error measures}
#' 	\item{test}{ the test result}
#'  \item{Stest}{ name of the test} 
#' and attributes:
#' 	 \item{msem}{ mean square error matrix of re-estimated parameters}
#'   \item{aiqm}{ average inverse quasi-information matrix over all re-estimated parameters}
#' 	 \item{qi}{ inverse quasi-information matrix at the parameter to be tested `\code{est$par}`}
#'   \item{relED}{ relative difference of the empirial and predicted standard error of the parameter to be tested} 
#'   \item{obs}{ list of simulated observed statistics}
#'   \item{optRes}{ results from re-estimating the model parameters for each simulated observation from `\code{obs}`}
#'	 \item{mean.score}{ average quasi-score, respectively, average gradient of the MD at the re-estimated parameters}
#' 	 \item{criterion}{always equal to "\code{qle}"}  
#'   \item{solInfo}{ results of the numerical consistency checks for each re-estimated parameter} 
#' 	 \item{info}{ list of indices of re-estimation results where the inversion of the quasi-information matrix failed,
#'       the re-estimated parameters have \code{NA} values, criterion function minimizations have errors or did not
#'       converge numerically, the integer seed value `\code{iseed}`}
#' 
#'  @details The function tests the null hypothesis \eqn{H_0:\,\hat{\theta}=\theta_0}, that is, whether the statistical
#'  model w.r.t. to the estimated parameter is true, against the alternative \eqn{H_1:\,\hat{\theta}\neq\theta_0} by testing based
#'  on a Monte Carlo approach (see vignette). Due to the approximate nature of the assumed statistical model for the observed data the
#'  exact distribution of the test statistics, that is, the Mahalanobis distance or quasi-deviance, is generally unknown and therefore
#'  its asymptotic distribution might be an unrealistic assumption for the null hypothesis. For this reason, and in order to retrieve an empirical
#'  P-value for testing, we generate (pseudo-)observations from the outcome of the model replications and re-estimate the model
#'  parameter for each realization in the same way as done before when estimating the model parameter. This includes all possible types
#'  of variance approximations (by kriging or average approximation) and types of prediction variance (kriging or the CV-based variance).
#' 
#'  The function expects an estimation result as returned from \code{\link{qle}}. If any generated observations are readily available
#'  at the final parameter estimate, then these can be passed by `\code{obs}`. Otherwise the function first generates those
#'  using `\code{nsim}` model replications at the estimated parameter as part of `\code{est}` or `\code{local}`. The criterion
#'  function approximations are used as it (specified in the object `\code{qsd}`) and will not be further improved by
#'  additional samples during the test.
#'  The value of the test statistic is either chosen as the current criterion function value in `\code{OPT}`
#'  (see  argument `\code{criterion}` in \code{\link{getQLmodel}}) or is taken from the optional argument `\code{local}`. Given the local results
#'  `\code{local}` of class \code{QSResult}, the user can also select a different criterion function as a test statistic than before when
#'  estimating the parameter itself. Apart from the quasi-deviance as a test statistic, in principle, any supported type of a least squares criterion,
#'  more generally, the Mahalanobis distance, can be used depending on the prefered type of variance matrix approximation, see \code{\link{covarTx}}.
#'  Practically, the re-estimations might fail to converge, however, then the user can control the convergence conditions of the local solvers
#'  (including quasi-scoring) by the corresponding control parameters (see \code{\link{searchMinimizer}}). Any failed re-estimations are
#'  excluded from the test results and stored in the attribute `\code{info}`. In addition, as part of the returned data frame `\code{param}`
#'  the empirical standard error, predicted standard error (based on the average inverse quasi-information matrix), the root mean square error,
#'  the bias and sample mean value of the re-estimated parameters are also available.
#' 
#' 	For an example please see the package vignette.
#' 
#' @author M. Baaske
#' @rdname qleTest
#' @export
qleTest <- function(est, local = NULL, sim, ...,
			 		 nsim = 100, obs = NULL, check.root = FALSE, alpha = 0.05,
					  na.rm = TRUE, cl = NULL, iseed = NULL, verbose = FALSE)
{				  
	if(.isError(est))
	  stop("Estimation result has errors. Please see attribute `error`.")    
    # last evaluation of criterion function  	
	if(.isError(est$final))
	  stop("Final criterion function evaluation has errors. Please check attribute `error`.")
		
    args <- list(...)
	# basic checks
	stopifnot(class(est) == "qle")
  	stopifnot(class(est$qsd)=="QLmodel")
	xdim <- attr(est$qsd$qldata,"xdim")
	
	# estimated parameter	 
	if(is.null(local)){		
		local <- est$final
		if(.isError(local) || !attr(est,"optInfo")$minimized)
		 stop("Final optimization result has errors. Please check attribute `final` and `optInfo`.")
	 	else if( local$convergence < 0)
		  warning(paste0("Local search did not converge by method: ",local$method))
    } else {
       stopifnot(class(local)=="QSResult")	 
	   # set current test statistic
	   est$qsd$criterion <- local$criterion		
	}	
	# argument ids
	id <- pmatch(names(args),names(formals(searchMinimizer)))
	# check `searchMinimizer' args 
	opt.args <- args[which(!is.na(id))]	
	
	# check with qsd hidden	
	.checkfun(searchMinimizer,opt.args,
		hide.args=c("x0","qsd"),check.default=FALSE)
		
	# use last Sigma (unless kriging variance matrix)  
	if(est$qsd$var.type != "kriging"){
		info <- attr(est,"optInfo")
		opt.args <- c(opt.args,list(W=info$W,theta=info$theta)) 
  	}
	# use final (local) method by default
	# if not given as an input argument
	if(!("method" %in% names(opt.args)))
	 opt.args$method <- local$method
    
    # MC simulation of observed `data`
	# if no observations supplied	
	if(is.null(obs)){		
		nid <- which(!is.na(id))
		if(length(nid) > 0L)
		 args <- args[-nid]		
		sim <- match.fun(sim)
		# check `sim` input values
		.checkfun(sim,args,remove=TRUE)
		stopifnot(is.numeric(nsim) || nsim > 0)
				
		simFun <- function(x) try(do.call(sim,c(list(x),args)))		
		sim.args <-
			list(sim=simFun,X=local$par,nsim=nsim,
				 mode="list",cl=cl,iseed=iseed,verbose=verbose)
		if(verbose)
		 cat("Simulating observed statistics...","\n")
 		
		obs <-
			tryCatch(
				do.call(simQLdata,sim.args),
				error = function(e) {
					msg <- paste0("Monte Carlo simulating observed statistics failed: ",
							 conditionMessage(e))
					message(msg) 
					.qleError(message=msg,call=match.call(),error=e)		   
				}
			)
		if(.isError(obs)) {			
			message("Simulating the observed statistics failed.")
			return(obs)
	 	}
		
	} else {
		if(.isError(obs))
		  stop("Argument `obs` has errors.")
		else if(class(obs) != "simQL" || !is.list(obs))
		  stop("Argument `obs` must be of class `simQL` and `list` type.")	   
	}
	
	# search minimizers	
	if(verbose)
	 cat("Estimating parameters...","\n")
	
 	RES <- do.call(doInParallel,
			  c(list(X=obs[[1]],
					FUN=function(obs,...)	{
						searchMinimizer(x0=local$par,qsd=est$qsd,...,    
							cvm=est$cvm,obs=obs,info=TRUE,inverted=TRUE,
								check=FALSE,verbose=verbose)
					},
				cl=cl), opt.args))    		
	
	if(inherits(RES,"error")) {
	  msg <- paste0("Could not find MC replicated parameters: ",
			  conditionMessage(RES))
	  message(msg)
	  return(.qleError(message=msg,call=match.call(),error=RES))	
	}
	# check results again
	ok <- which(sapply(RES,function(x) !.isError(x) && x$convergence>0 ))
	if(length(ok) == 0L){
		msg <- paste("All re-estimations have errors or did not converge: ")
		warning(msg)
		return(.qleError(message=msg,call=match.call(),error=RES))	   
	} else if(length(ok) < length(RES))
		message("Errors in re-estimating the parameters. Check attribute `optInfo` and `info`.")
  
	# also check H_0?
	# because sampling MC theta
	# must be done under H0	
	aiqm <- NULL
	mScore <- NULL	    	  
	# average quasi-information matrix
	invI <- 
		lapply(RES[ok],
			function(x) {						
				try(solve(x$I),silent=TRUE) 
			})
	badInv <- sapply(invI,function(x) inherits(x,"try-error") || anyNA(x))
	if(any(badInv))
		message(paste0("A total of ",sum(badInv)," inversions of quasi-information matrices failed. Check attribute `info`."))    
	else if(all(badInv)){
		msg <- paste0("All iversions of quasi-information matrices failed.")
		message(msg)
		return(.qleError(message=msg,call=match.call(),error=badInv))	 
	}
	# average matrix of inverse qi matrices
	aiqm <- matrix(
			 colMeans(
				do.call(rbind,
					lapply(invI[!badInv],as.numeric)
				)
			),ncol=xdim)				
 	
	# estimates	
	mpars <- do.call(rbind,lapply(RES[ok],function(x) x$par))
	mScore <- do.call(rbind,lapply(RES[ok],function(x) x$score))	
	has.na <- (rowSums(is.na(cbind(mScore,mpars))) > 0L)	
	if(na.rm && any(has.na)) {
		ok <- ok[-which(has.na)]
		mpars <- mpars[ok,,drop=FALSE]
		mScore <- mScore[ok,,drop=FALSE]
		warning("Removed `NA` values from quasi-scores.")
  	}
	mScore <- try(colMeans(mScore),silent=TRUE)		
	if(nrow(mpars) < 11L)
	 warning("Only a small number of 10 or less parameters could be re-estimated.")	
    
 	# some (empirical) measures	
	msem <- .MSE(mpars,local$par)
		
	# value of test statistic at re-estimated parameters			
	tvals <- sapply(RES[ok],"[[","value") 
	stopifnot(is.numeric(tvals))
	
	# invert QI for predicted std. error (asymptotic) at estimated theta 
	qi <- try(solve(local$I),silent=TRUE)
	if(inherits(qi,"try-error") || anyNA(qi))
	  message("Inversion of quasi-information matrix failed")
	
	# get efficient score test (with MC parameters)
	B <- structure(
			 data.frame(
			  cbind("par"=local$par,
				    "se"=apply(mpars,2,sd),					
					"rmse"=sqrt(diag(msem)),
					"bias"=colMeans(t(t(mpars)-local$par)),
					"mean"=colMeans(mpars))),
			"sb"=local$value, "Sb"=tvals,
			"test"=est$qsd$criterion)
	
	# had errors
	hasError <- which(!(1:length(RES) %in% ok))
	if(length(hasError) > 0L)	
	  message(paste0("A total of ",length(hasError)," re-estimations failed."))
  
    chk <- NULL
    if(check.root && est$qsd$criterion=="qle") {
		chk <- checkMultRoot(est,verbose=verbose)
		if(.isError(chk))
		 message(.makeMessage("Consistency check for the estimated model parameter failed."))		
	}
  
	relED <-
	 if(!anyNA(c(msem,qi)) && is.matrix(qi) && is.matrix(msem)) {
		 abs(1-sqrt(diag(msem))/sqrt(diag(qi)))
	 } else {
		message("Failed to compute relative difference of empirical and predicted error.")
		NULL
	 }
	
	# results
	structure(.qleTest(B,alpha),				# test results
		    msem=msem,							# mean square error matrix
	    	aiqm=aiqm,							# average inverse QI (re-estimated parameters)
			qi=qi,								# inverse QI at estimated theta
	    	relED=relED,
			obs=if(verbose) obs else NULL,		# (MC) observations
	    	optRes=if(verbose) RES else NULL,	# optimization results
			mean.score=mScore,					# average score/gradient
			criterion=est$qsd$criterion,
			solInfo=chk,						# values of consistency criteria
	  info=list(badInv=which(badInv), 			# inversion errors
			  	hasNa=which(has.na), 			# indices of NA parameters 
		 		hasError=hasError,
				iseed=iseed),																
	 class=c("qleTest"), call=sys.call())	
}

# printing function
#' @name print.qleTest
#' 
#' @title print \code{qleTest} results 
#' 
#' @description print the results from \code{\link{qleTest}}
#' 
#' @param x      object of class \code{qleTest} from \code{\link{qleTest}}
#' @param pl	 ignored
#' @param digits number of (significant) digits
#' @param ... 	 ignored
#' 
#' @rdname print.qleTest
#' @method print qleTest
#' @export 
print.qleTest <- function(x, pl = 1, digits = 5,...) {
	if(.isError(x)){
	 	print(x)
		invisible(return(NULL))
 	}
	
	if(!is.null(attr(x,"call"))){
	  cat("\nCall:\n\n")
	  cat(paste(deparse(attr(x,"call")), sep="\n", collapse = "\n"), "\n\n", sep="")
	}
	# consistency check of solution
	chk <- attr(x,"solInfo")
	if(!is.null(chk)){
		cat("Consistency check - the smaller the better: \n\n")
		print(format(signif(chk,digits=digits)),print.gap=2,right=FALSE,quote=FALSE)
		cat("\n\n")		
	}	
	cat("Coefficients:\n")	
	print(format(x$param, digits=digits),
			print.gap = 2, right=FALSE, quote = FALSE)	
	cat("\n\n")
	cat(x$Stest,"\n\n")
	vals <- c(format(x$test[1], digits=digits),
			   formatC(signif(x$test[2], digits=digits), digits=digits,format="fg", flag="#"))
	names(vals) <- colnames(x$test)
	print(vals, print.gap = 2, right=FALSE, quote = FALSE)
	
	cat("\n\n")
	if(!is.null(attr(x,"mean.score"))) {
	  if(attr(x,"criterion") == "mahal")
		cat("Average gradient: \n\n")
	  else cat("Average quasi-score: \n\n")
	  print(format(attr(x,"mean.score"), digits=digits),
			  print.gap = 2, right=FALSE,quote = FALSE)
	  cat("\n\n")
	  qi <- attr(x,"qi")
	  aiqm <- attr(x,"aiqm")
	  if(!is.null(aiqm) && !.isError(aiqm) && !is.null(qi) && !.isError(qi) && !is.null(attr(x,"relED")))
		  pse <- as.data.frame( cbind(sqrt(diag(aiqm)), sqrt(diag(qi)), attr(x,"relED") ) )
	  	  dimnames(pse) <- list(row.names(x$param),c("Average","Estimate", "|1-RMSE/Estimate|"))
		  cat("Predicted std. errors: \n\n")
		  print(format(pse, digits=digits),
				  print.gap = 2, right=FALSE, quote = FALSE)			  
    }	
	invisible(x)	
}
