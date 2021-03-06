% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/covariance.R
\name{QLmodel}
\alias{QLmodel}
\title{Construct the quasi-likelihood approximation model}
\usage{
QLmodel(qldata, lb, ub, obs, mods, nfit = 1, cv.fit = TRUE,
  var.type = c("wcholMean", "cholMean", "wlogMean", "logMean", "kriging",
  "const"), useVar = TRUE, criterion = c("qle", "mahal"),
  verbose = FALSE)
}
\arguments{
\item{qldata}{data frame of (initial) simulation results (see \code{\link{setQLdata}})}

\item{lb}{numeric vector of lower bounds defining the (hyper)box}

\item{ub}{numeric vector of upper bounds defining the (hyper)box}

\item{obs}{numeric vector of observed statistics}

\item{mods}{list of (fitted) covariance models (see \code{\link{fitSIRFk}})}

\item{nfit}{number of cycles, \code{nfit=1} (default), after which covariance
parameters are re-estimated and otherwise only re-used}

\item{cv.fit}{logical, \code{TRUE} (default), whether to re-fit CV models (re-estimate covariance parameters)}

\item{var.type}{name of the variance approximation method (see \code{\link{covarTx}})}

\item{useVar}{logical, \code{TRUE} (default), whether to use prediction variances (see details)}

\item{criterion}{global criterion function for sampling and minimization, either "\code{qle}" or "\code{mahal}"}

\item{verbose}{logical, \code{FALSE} (default), whether to give intermediate output}
}
\value{
An object of class \code{\link{QLmodel}} which stores the data frame of simulation results, bounds on
 the parameter space, covariance models for kriging, vector of observed statistics as well as options for
 kriging and covariance parameter estimation.
}
\description{
Aggregate and construct the data for quasi-likelihood estimation
}
\details{
The function aggregates all required information for quasi-likelihood estimation and defines the input object to the 
  function \code{\link{qle}}, stores the fitted covariance models of the sample means of the statistics and the type of variance
  matrix approximation. For an advanced setup of the estimation procedure and more involved statistical models this function
  explicitly offers the data structure to construct individual covariance models for each statistic, see \code{\link{setCovModel}}.
  The user has the choice whether or not to make use of kriging prediction variances by `\code{useVar}` to account for the simulation
	 error when constructing the approximation of the variance matrix and the quasi-score function. If \code{useVar=TRUE}, then a kriging
  procedure including the computation of prediction variances based on kriging is automatically used. Otherwise the so-called
  \emph{dual} approach is employed which has some computational advantage if prediction variances are not required.
}
\examples{

data(normal)

# As an example we re-use the stored normal data and fit 
# a generalized covariance model to the data using simulation
# variances as local variances for REML estimation.
mods <- fitSIRFk(qsd$qldata, verbose=TRUE)

# construct QL approximation model
qsd <- QLmodel(qsd$qldata,qsd$lower,qsd$upper,
			    c("T1"=2,"T2"=1),mods)

}
\seealso{
\code{\link{getQLmodel}}, \code{\link{updateCovModels}}
}
\author{
M. Baaske
}
