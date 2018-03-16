# Copyright (C) 2017 Markus Baaske. All Rights Reserved.
# This code is published under the GPL (>=3).
#
# File: 	helpers.R
# Date:  	12/08/2017
# Author: 	Markus Baaske
#
# Define a set of auxiliary functions 

.isPosDef <- function(X) {
	.C(C_isPositiveDefinite,
			as.matrix(X),
			as.integer(NROW(X)),
			pos=integer(1))$pos
}

.splitList <- function(x, ncl) { 
  lapply(parallel::splitIndices(length(x), ncl), function(i) x[i])
}

.MSE <- function(estimate,param) {	
	stopifnot(is.vector(param))
	stopifnot(is.matrix(estimate))	
	param <- matrix(rep(param,nrow(estimate)),
				ncol=length(param),byrow=TRUE)
	x <- as.matrix(estimate-param)
	(t(x)%*%x)/nrow(x)
}

.COL2LIST  <- function(x) {
 	lapply(seq_len(NCOL(x)), function(i) x[,i])
}

.ROW2LIST  <- function(x) {
	if(is.vector(x))
	 return (list(x))
	else if(is.matrix(x)){
 		nms <- dimnames(x)[[2]]
		lapply(seq_len(NROW(x)), function(i) structure(x[i,], names=nms))
	} else {
		stop("Non supported type of conversion.") 
	}
}

.LIST2ROW <- function(x) {
  if(is.numeric(x))
  	return( rbind(x) )
  else if(is.list(x)) {
  	do.call(rbind,
		lapply(x,function(x)
			 unlist(x,recursive=FALSE))
    )
  } else {
	  stop("Non supported type of conversion.")
  }
}


.NORM <- function(x) {
	crossprod(x)^0.5	
}

.MED <- function(x,y,z) {
	if ( ((x - y) * (z - x)) >= 0 ) return (x)
	else if ( ((y - x) * (z - y)) >= 0 ) return (y)
	else return (z);
}

.PROJMED <-function(x,lb,ub) {
	sapply(seq_len(length(x)),function(i) .MED(x[i],lb[i],ub[i]))
}

.ISO.NORM  <- function(h) {
	h <- data.matrix(h)
	if( ncol(h)>1 )
		h <- c(sqrt(h^2 %*% rep(1,ncol(h))))
	c(abs(h))
}

.WEIGHT.NORM  <- function(h,W=diag(1,ncol(h))) {
	h <- data.matrix(h)	
	sqrt(as.vector(apply(h,1,function(x) crossprod(x,W%*%x))))
}

.distX <- function(X,Y=X,W=diag(1,ncol(X))) {
	idxx <- rep(1:NROW(X),NROW(Y))
	idxy <- rep(1:NROW(Y),rep(NROW(X),NROW(Y)))	
	return (.WEIGHT.NORM(X[idxx,]-Y[idxy,],W))
}


.check.distAll <- function(X,xTol=1e-12) {	
	idxx <- rep(1:NROW(X),NROW(X))
	idxy <- rep(1:NROW(X),rep(NROW(X),NROW(X)))
	H <- matrix(.ISO.NORM( X[idxx,]-X[idxy,]),nrow=NROW(X))
	idx <-which( (H<xTol & col(H)>row(H)),arr.ind=TRUE)
	if(NROW(idx)>0)
		structure(rbind(idx),
				min=min(apply(idx,1,function(a) H[a[1],a[2]])))
	else
		structure(matrix(,0,ncol(X)),
			min=min(H[row(H)<col(H)]))
}

.check.distX <- function(X,x,xTol=1e-8) {	
	idxx <- rep(1:NROW(X),1)
	Y <- matrix(c(x), nrow=NROW(X),ncol=length(c(x)),byrow=TRUE)
	H <- .ISO.NORM(X-Y)
	return (which(H<xTol,arr.ind=TRUE))
}

.min.distX <- function(X,Y=X) {
	idxx <- rep(1:NROW(X),NROW(Y))
	idxy <- rep(1:NROW(Y),rep(NROW(X),NROW(Y)))
	h <- X[idxx,]-Y[idxy,]
	return (min(.ISO.NORM(h)))
}

.min.distIdx <- function(X,Y=X) {
	idxx <- rep(1:NROW(X),NROW(Y))
	idxy <- rep(1:NROW(Y),rep(NROW(X),NROW(Y)))
	h <- X[idxx,]-Y[idxy,]
	return (which.min(.ISO.NORM(h)))
}

#.min.distXY <- function(X,Y,W=NULL) {
#	idxx <- rep(1:NROW(X),NROW(Y))
#	idxy <- rep(1:NROW(Y),rep(NROW(X),NROW(Y)))
#	h <- X[idxx,]-Y[idxy,]	
#	x <- if(!is.null(W)) .WEIGHT.NORM(h,W) else .ISO.NORM(h)	
#	splt <- split(x, ceiling(seq_along(x)/NROW(X)))
#	return (as.vector(unlist(lapply(splt,min))))
#}

.min.distXY <- function(X,Y) {
	idxx <- rep(1:NROW(X),NROW(Y))
	idxy <- rep(1:NROW(Y),rep(NROW(X),NROW(Y)))	
	x <- .ISO.NORM(X[idxx,]-Y[idxy,])	
	splt <- split(x, ceiling(seq_along(x)/NROW(X)))
	return (as.vector(unlist(lapply(splt,min))))
}


.min.distXYIdx <- function(X,Y) {
	idxx <- rep(1:NROW(X),NROW(Y))
	idxy <- rep(1:NROW(Y),rep(NROW(X),NROW(Y)))
	h <- X[idxx,]-Y[idxy,]	
	x <- .ISO.NORM(h)	
	splt <- split(x, ceiling(seq_along(x)/NROW(X)))
	return (as.vector(unlist(lapply(splt,which.min))))
}

.min.distXYTol <- function(X,Y,xTol=1e-10) {
	idxx <- rep(1:NROW(X),NROW(Y))
	idxy <- rep(1:NROW(Y),rep(NROW(X),NROW(Y)))
	h <- X[idxx,]-Y[idxy,]	
	x <- .ISO.NORM(h)	
	splt <- split(x, ceiling(seq_along(x)/NROW(X)))
	do.call(rbind,lapply(splt,function(x) {
		id <- which(x<xTol,arr.ind=TRUE)
		c(id,x[id])
	}))		  
}

# generalized Eigendecomposition
geneigen <- function(A = NULL, B, vl=TRUE, vr=TRUE, only.values = TRUE,
			 check.input = FALSE, verbose = FALSE) {
	# if only one matrix is given
	# just return the eigenvalues
	if(is.null(A)) {
		maxE <- try(eigen(B,only.values=TRUE),silent=TRUE)
		if(inherits(maxE,"try-error")) {
			message("Failed to get Eigen values.")
			return (1)
		}
		return (maxE$values)		
	}	
	# now generalized decomposition
	A <- as.matrix(A)
	B <- as.matrix(B)
	N <- NROW(A)
	if(check.input) {
		if (!N) stop("0 x 0 matrix")
		if (NROW(A)!=NCOL(A) || NROW(B)!=NCOL(B))
			stop("non-square matrix in 'geneigen'")
		if (NROW(A)!=NROW(B) || NCOL(B)!=NCOL(A))
			stop("dimension of matrices do not match in 'geneigen'")
		N <- as.integer(N)
		if(is.na(N))
			stop("invalid nrow(x) in 'geneigen'")
	}

	if(only.values) {
		JOBVL <- "N"
		JOBVR <- "N"
	} else  {
		JOBVL <- if(vl) "V" else "N"
		JOBVR <- if(vr) "V" else "N"
	}

	N <- as.integer(N)
	LWORK <- as.integer(max(1,8*N))

	res <- tryCatch( {
			.Fortran(C_dggev, JOBVL, JOBVR, N, A, N, B, N,
						ALPHAR=numeric(N), ALPHAI=numeric(N),BETA=numeric(N),
						vl=if(JOBVL=="V") matrix(0,nrow=N,ncol=N) else numeric(1), N,
						vr=if(JOBVR=="V") matrix(0,nrow=N,ncol=N) else numeric(1), N,
						double(max(1,LWORK)), LWORK, info=integer(1L))} ,
			error = function(e) {				
				msg <- paste(.makeMessage("'geneigen' failed: "),conditionMessage(e))
				message(msg)
				structure(list(message=msg, call=sys.call()), 
						class=c("error","condition"), error = e)
			}			
	)
	if(verbose && res$info >0)
	  print(paste("'geneigen': lapack info:  ", res$info))
	if( all(res$ALPHAI==0) ) {
        alpha  <- res$ALPHAR
        beta   <- res$BETA
		values <- alpha/beta
    } else {
        alpha <- complex(real=res$ALPHAR,imaginary=res$ALPHAI)
        beta  <- res$BETA
		values <- alpha/beta
    }
	if(only.values) {
		return( values )
	} else {
		vl = if(JOBVL=="V") res$vl else NULL
		vr = if(JOBVR=="V") res$vr else NULL
		
		return(list("values"=values,
					"alpha"=alpha,
					"beta"=beta,"vl"=vl,"vr"=vr))	
	}
    
}

gsiSolve <- function(X,B=diag(1,NROW(X)),cond=1e-12,...,use.solve=TRUE) {
	stopifnot(is.matrix(X))
	if(use.solve){
		x <- try(solve(X,B,...),silent=TRUE)
		if(!inherits(x,"try-error"))
			return (x)
	}
	
	if( NROW(X) > 1 ) {
		SVD <- svd(X)
		if(inherits(SVD,"error"))
			stop("SVD in `gsiSolve`failed.")
		d <- SVD$d
		SVD$v%*% (ifelse(d>cond,1/d,0) * (t(SVD$u)%*%as.matrix(B)))
	} else if(all(is.finite(X))) as.matrix(B)/X
	else stop("Non finite values in matrix `X`.")
}

gsiInv <- function(X, cond=1e-12, ..., use.solve=TRUE) {
	stopifnot(is.matrix(X))
	if(use.solve){
		Xinv <- try(solve(X,...),silent=TRUE)
		if(!inherits(Xinv,"try-error"))
		  return (Xinv)
	}
	
	if(NROW(X) > 1) {
		SVD <- svd(X)
		if(inherits(SVD,"error"))
		  stop("SVD in `gsiInv`failed.")
	    d <- SVD$d
		SVD$v %*% diag(ifelse(d>cond,1/d,0)) %*% t(SVD$u)		
	} else if(all(is.finite(X))) 1/X
	 else stop("Non finite values in matrix `X`.")
}
