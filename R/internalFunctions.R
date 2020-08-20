#' @title Compute the sum of squares between the theoretical and observed quantiles of a beta distribution.
#'
#' @description
#' This is a helper function that compute the sum of squares between two theoretical and observed quantiles of a beta distribution (typically the lower and upper bounds of a confidence interval).
#' This function is for internal use to find the best-fit beta distribution for a given confidence interval.
#'
#' @param abPars The shape1 and shape2 parameters of the theoretical beta distribution.
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyBetaPars}}, \code{\link{optim}}, \code{\link{qbeta}}
#'


ssBetaPars<-function(abPars,pLow,pUpp,alpha=0.05){
  res<-(qbeta(alpha/2,abPars[1],abPars[2])-pLow)^2+(qbeta(1-alpha/2,abPars[1],abPars[2])-pUpp)^2
  return(res)
}

#' @title Determine the parameters of the best-fit beta distribution for a given confidence interval for a probability parameter.
#'
#' @description
#' Finds the best-fit beta distribution parameters for a given confidence interval for a probability parameter and returns the shape1, shape2 parameters.
#'
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(50,50).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A vector of length 2 giving the 2 parameter shape1 and shape1 for use with rbeta/dbeta/pbeta/qbeta
#'
#' @seealso
#' \code{\link{ssBetaPars}}, \code{\link{optim}}, \code{\link{dbeta}}
#'


identifyBetaPars<-function(pLow,pUpp,alpha=0.05,initPars=c(50,50),maxiter=1e3){
  if(pLow<0 | pUpp>1 | pLow>pUpp){stop("pLow and pUpp need to be both contained within [0,1] and pLow needs to be lower than pUpp.")}

  res<-suppressWarnings(optim(fn=ssBetaPars,pLow=pLow,pUpp=pUpp,par=initPars,control=list(maxit=maxiter)))
  if(res$convergence!=0){stop("optim() as called by identifyBetaPars() failed to converge.")}

  return(res$par)
}

#' @title Compute the sum of squares between the theoretical and observed quantiles of a normal / Gaussian distribution.
#'
#' @description
#' This is a helper function that compute the sum of squares between two theoretical and observed quantiles of a normal distribution (typically the lower and upper bounds of a confidence interval).
#' This function is for internal use to find the best-fit normal distribution for a given confidence interval.
#'
#' @param muSigPars The mean and standard deviation parameters of the theoretical normal distribution.
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyNormPars}}, \code{\link{optim}}, \code{\link{qnorm}}
#'

ssNormPars<-function(muSigPars,pLow,pUpp,alpha=0.05){
  res<-(qnorm(alpha/2,mean=muSigPars[1],sd=muSigPars[2])-pLow)^2+(qnorm(1-alpha/2,mean=muSigPars[1],sd=muSigPars[2])-pUpp)^2
  return(res)
}


#' @title Determine the parameters of the best-fit normal / Gaussian distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit normal distribution parameters for a given confidence interval and returns the mean and sd parameters.
#'
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(50,50).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A vector of length 2 giving the 2 parameter mean and sd for use with rnorm/dnorm/pnorm/qnorm
#'
#' @seealso
#' \code{\link{ssNormPars}}, \code{\link{optim}}, \code{\link{dnorm}}
#'

identifyNormPars<-function(pLow,pUpp,alpha=0.05,initPars=c(0,1),maxiter=1e3){
  if(pLow>pUpp){stop("pLow needs to be lower than pUpp.")}

  res<-suppressWarnings(optim(fn=ssNormPars,pLow=pLow,pUpp=pUpp,par=initPars,control=list(maxit=maxiter)))
  if(res$convergence!=0){stop("optim() as called by identifyNormPars() failed to converge.")}

  return(res$par)
}


#' @title Compute the sum of squares between the theoretical and observed quantiles of a Poisson distribution.
#'
#' @description
#' This is a helper function that compute the sum of squares between two theoretical and observed quantiles of a normal distribution (typically the lower and upper bounds of a confidence interval).
#' This function is for internal use to find the best-fit normal distribution for a given confidence interval.
#'
#' @param poisPar The rate parameter of the theoretical Poisson distribution.
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyPoisPars}}, \code{\link{optim}}, \code{\link{qpois}}
#'

ssPoisPar<-function(poisPar,pLow,pUpp,alpha=0.05){
  res<-(qpois(alpha/2,lambda=poisPar)-pLow)^2+(qpois(1-alpha/2,lambda=poisPar)-pUpp)^2
  return(res)
}


#' @title Determine the parameters of the best-fit Poisson distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit Poisson distribution parameters for a given confidence interval and returns the rate parameter.
#'
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A single number > 0, giving the initial parameter value to start the optimisation; defaults to 5.
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A single number giving the rate parameter for use with rpois/dpois/ppois/qpois
#'
#' @seealso
#' \code{\link{ssPoisPar}}, \code{\link{optim}}, \code{\link{dpois}}
#'

identifyPoisPars<-function(pLow,pUpp,alpha=0.05,initPars=5,maxiter=1e3){
  if(pLow>pUpp){stop("pLow needs to be lower than pUpp.")}

  res<-suppressWarnings(optim(fn=ssPoisPar,pLow=pLow,pUpp=pUpp,par=initPars,control=list(maxit=maxiter)))
  if(res$convergence!=0){stop("optim() as called by identifyPoisPars() failed to converge.")}

  return(res$par)
}

