#' @title Compute the sum of squares between the theoretical and observed quantiles of a beta distribution.
#'
#' @description
#' This is a helper function that compute the sum of squares between two theoretical and observed quantiles of a beta distribution (typically the lower and upper bounds of a confidence interval).
#' This function is for internal use to find the best-fit beta distribution for a given confidence interval.
#'
#' @param abPars The shape1 and shape2 parameters of the theoretical beta distribution.
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyBetaPars}}, \code{\link{optim}}, \code{\link{qbeta}}
#'

ssBetaPars<-function(abPars,qLow,qUpp,alpha=0.05){
  res<-(qbeta(alpha/2,abPars[1],abPars[2])-qLow)^2+(qbeta(1-alpha/2,abPars[1],abPars[2])-qUpp)^2
  return(res)
}

#' @title Determine the parameters of the best-fit beta distribution for a given confidence interval for a probability parameter.
#'
#' @description
#' Finds the best-fit beta distribution parameters for a given confidence interval for a probability parameter and returns the shape1, shape2 parameters.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(50,50).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A vector of length 2 giving the 2 parameters shape1 and shape1 for use with rbeta/dbeta/pbeta/qbeta.
#'
#' @seealso
#' \code{\link{ssBetaPars}}, \code{\link{optim}}, \code{\link{dbeta}}
#'

identifyBetaPars<-function(qLow,qUpp,alpha=0.05,initPars=c(50,50),maxiter=1e3){
  if(qLow<0 | qUpp>1 | qLow>qUpp){stop("qLow and qUpp need to be both contained within [0,1] and qLow needs to be lower than qUpp.")}

  res<-suppressWarnings(optim(fn=ssBetaPars,qLow=qLow,qUpp=qUpp,par=initPars,control=list(maxit=maxiter)))
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
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyNormPars}}, \code{\link{optim}}, \code{\link{qnorm}}
#'

ssNormPars<-function(muSigPars,qLow,qUpp,alpha=0.05){
  res<-(qnorm(alpha/2,mean=muSigPars[1],sd=muSigPars[2])-qLow)^2+(qnorm(1-alpha/2,mean=muSigPars[1],sd=muSigPars[2])-qUpp)^2
  return(res)
}


#' @title Determine the parameters of the best-fit normal / Gaussian distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit normal distribution parameters for a given confidence interval and returns the mean and sd parameters.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(50,50).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A vector of length 2 giving the 2 parameters mean and sd for use with rnorm/dnorm/pnorm/qnorm.
#'
#' @seealso
#' \code{\link{ssNormPars}}, \code{\link{optim}}, \code{\link{dnorm}}
#'

identifyNormPars<-function(qLow,qUpp,alpha=0.05,initPars=c(0,1),maxiter=1e3){
  if(qLow>qUpp){stop("qLow needs to be lower than qUpp.")}

  res<-suppressWarnings(optim(fn=ssNormPars,qLow=qLow,qUpp=qUpp,par=initPars,control=list(maxit=maxiter)))
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
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyPoisPars}}, \code{\link{optim}}, \code{\link{qpois}}
#'

ssPoisPars<-function(poisPar,qLow,qUpp,alpha=0.05){
  res<-(alpha/2-ppois(qLow,lambda=poisPar))^2+(1-alpha/2-ppois(qUpp,lambda=poisPar))^2
  return(res)
}


#' @title Determine the parameters of the best-fit Poisson distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit Poisson distribution parameters for a given confidence interval and returns the rate parameter.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A single number > 0, giving the initial parameter value to start the optimisation; defaults to 5.
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A single number giving the rate parameter for use with rpois/dpois/ppois/qpois.
#'
#' @seealso
#' \code{\link{ssPoisPars}}, \code{\link{optim}}, \code{\link{dpois}}
#'

identifyPoisPars<-function(qLow,qUpp,alpha=0.05,initPars=5,maxiter=1e3){
  if(qLow>qUpp){stop("qLow needs to be lower than qUpp.")}

  res<-suppressWarnings(optim(fn=ssPoisPars,qLow=qLow,qUpp=qUpp,par=initPars,control=list(maxit=maxiter)))
  if(res$convergence!=0){stop("optim() as called by identifyPoisPar() failed to converge.")}

  return(res$par)
}


#' @title Compute the sum of squares between the theoretical and observed quantiles of a negative binomial distribution.
#'
#' @description
#' This is a helper function that compute the sum of squares between two theoretical and observed quantiles of a negative binomial distribution (typically the lower and upper bounds of a confidence interval).
#' This function is for internal use to find the best-fit negative binomial distribution for a given confidence interval.
#'
#' @param sizeProbPars The size and prob parameters of the theoretical negative binomial distribution.
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyNegBinPars}}, \code{\link{optim}}, \code{\link{qnbinom}}
#'

ssNegBinPars<-function(sizeProbPars,qLow,qUpp,alpha=0.05){
  res<-(alpha/2-pnbinom(qLow,size=sizeProbPars[1],prob=sizeProbPars[2]))^2+(1-alpha/2-pnbinom(qUpp,size=sizeProbPars[1],prob=sizeProbPars[2]))^2
  return(res)
}


#' @title Determine the parameters of the best-fit negative binomial distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit negative binomial distribution parameters for a given confidence interval and returns the size, prob parameters.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(10,0.5).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A vector of length 2 giving the 2 parameters size and prob for use with rnbinom/dnbinom/pnbinom/qnbinom.
#'
#' @seealso
#' \code{\link{ssNegBinPars}}, \code{\link{optim}}, \code{\link{dnbinom}}
#'

identifyNegBinPars<-function(qLow,qUpp,alpha=0.05,initPars=c(10,0.5),maxiter=1e3){
  if(qLow>qUpp){stop("qLow needs to be lower than qUpp.")}

  res<-suppressWarnings(optim(fn=ssNegBinPars,qLow=qLow,qUpp=qUpp,par=initPars,control=list(maxit=maxiter)))
  if(res$convergence!=0){stop("optim() as called by identifyNegBinPars() failed to converge.")}

  return(res$par)
}


#' @title Compute the sum of squares between the theoretical and observed quantiles of a gamma distribution.
#'
#' @description
#' This is a helper function that compute the sum of squares between two theoretical and observed quantiles of a gamma distribution (typically the lower and upper bounds of a confidence interval).
#' This function is for internal use to find the best-fit gamma distribution for a given confidence interval.
#'
#' @param shapeRatePars The shape and rate parameters of the theoretical gamma distribution.
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyGammaPars}}, \code{\link{optim}}, \code{\link{qgamma}}
#'

ssGammaPars<-function(shapeRatePars,qLow,qUpp,alpha=0.05){
  res<-(qgamma(alpha/2,shape=shapeRatePars[1],rate=shapeRatePars[2])-qLow)^2+(qgamma(1-alpha/2,shape=shapeRatePars[1],rate=shapeRatePars[2])-qUpp)^2
  return(res)
}


#' @title Determine the parameters of the best-fit gamma distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit gamma distribution parameters for a given confidence interval and returns the shape, rate parameters.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(1,1).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A vector of length 2 giving the 2 parameters shape and rate for use with rgamma/dgamma/pgamma/qgamma.
#'
#' @seealso
#' \code{\link{ssGammaPars}}, \code{\link{optim}}, \code{\link{dgamma}}
#'

identifyGammaPars<-function(qLow,qUpp,alpha=0.05,initPars=c(1,1),maxiter=1e3){
  if(qLow>qUpp){stop("qLow needs to be lower than qUpp.")}

  res<-suppressWarnings(optim(fn=ssGammaPars,qLow=qLow,qUpp=qUpp,par=initPars,control=list(maxit=maxiter)))
  if(res$convergence!=0){stop("optim() as called by identifyGammaPars() failed to converge.")}

  return(res$par)
}



#' @title Compute the sum of squares between the theoretical and observed quantiles of an exponential distribution.
#'
#' @description
#' This is a helper function that compute the sum of squares between two theoretical and observed quantiles of an exponential distribution (typically the lower and upper bounds of a confidence interval).
#' This function is for internal use to find the best-fit exponential distribution for a given confidence interval.
#'
#' @param ratePar The rate parameter of the theoretical exponential distribution.
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A single number, the sum of squares.
#'
#' @seealso
#' \code{\link{identifyExpPars}}, \code{\link{optim}}, \code{\link{qexp}}
#'

ssExpPars<-function(ratePar,qLow,qUpp,alpha=0.05){
  res<-(qexp(alpha/2,rate=ratePar)-qLow)^2+(qexp(1-alpha/2,rate=ratePar)-qUpp)^2
  return(res)
}


#' @title Determine the parameters of the best-fit exponential distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit exponential distribution parameter for a given confidence interval and returns the rate parameter.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A single number giving the initial parameter value to start the optimisation; defaults to 1.
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A single number giving the rate parameter for use with rexp/dexp/pexp/qexp.
#'
#' @seealso
#' \code{\link{ssExpPars}}, \code{\link{optim}}, \code{\link{dexp}}
#'

identifyExpPars<-function(qLow,qUpp,alpha=0.05,initPars=1,maxiter=1e3){
  if(qLow>qUpp){stop("qLow needs to be lower than qUpp.")}

  res<-suppressWarnings(optim(fn=ssExpPars,qLow=qLow,qUpp=qUpp,par=initPars,control=list(maxit=maxiter)))
  if(res$convergence!=0){stop("optim() as called by identifyExpPar() failed to converge.")}

  return(res$par)
}

