#' @title Find the best-fit beta distribution for a given confidence interval for a probability parameter.
#'
#' @description
#' Finds the best-fit beta distribution for a given confidence interval for a probability parameter; returns the corresponding density, distribution, quantile and sampling functions.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(50,50).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A list with 5 elements:
#' \item{r}{The sampling function.}
#' \item{d}{The density function.}
#' \item{p}{The distribution function.}
#' \item{q}{The quantile function.}
#' \item{pars}{A vector of length 2 giving the two shape parameters for the best-fit beta distribution (\code{shape1} and \code{shape2} as in \code{\link{rbeta}}, \code{\link{dbeta}}, \code{\link{pbeta}}, \code{\link{qbeta}}).}
#'
#' @seealso
#' \code{\link{identifyBetaPars}}, \code{\link{optim}}, \code{\link{dbeta}}
#'
#' @examples
#' b<-getBetaFromCI(qLow=0.1167,qUpp=0.1636,initPars=c(200,800))
#' print(b$pars) # the fitted parameter values
#' b$r(10) # 10 random values from the fitted beta distribution
#' b$d(0.15) # the probability density at x=0.15 for the fitted beta distribution
#' b$p(0.15) # the cumulative density at x=0.15 for the fitted beta distribution
#' b$q(c(0.25,0.5,0.75)) # the 25th, 50th (median) and 75th percentiles of the fitted distribution
#' x<-seq(0,1,length=1e3)
#' y<-b$d(x)
#' plot(x,y,type="l",xlab="",ylab="density") # density plot for the fitted beta distribution
#'
#' @export getBetaFromCI

getBetaFromCI<-function(qLow,qUpp,alpha=0.05,initPars=c(50,50),maxiter=1e3){
  pars<-identifyBetaPars(qLow,qUpp,alpha=alpha,initPars=initPars,maxiter=maxiter)

  rFun<-function(n){rbeta(n,shape1=pars[1],shape2=pars[2])}
  dFun<-function(x){dbeta(x,shape1=pars[1],shape2=pars[2])}
  pFun<-function(q){pbeta(q,shape1=pars[1],shape2=pars[2])}
  qFun<-function(p){qbeta(p,shape1=pars[1],shape2=pars[2])}

  res<-list(r=rFun,d=dFun,q=qFun,p=pFun,pars=pars)
  return(res)
}


#' @title Find the best-fit normal / Gaussian distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit normal distribution for a given confidence interval; returns the corresponding density, distribution, quantile and sampling functions.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values (mean & sd) to start the optimisation; defaults to c(0,1).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A list with 5 elements:
#' \item{r}{The sampling function.}
#' \item{d}{The density function.}
#' \item{p}{The distribution function.}
#' \item{q}{The quantile function.}
#' \item{pars}{A vector of length 2 giving the mean and standard deviation for the best-fit normal distribution (\code{mean} and \code{sd} as in \code{\link{rnorm}}, \code{\link{dnorm}}, \code{\link{pnorm}}, \code{\link{qnorm}}).}
#'
#' @seealso
#' \code{\link{identifyNormPars}}, \code{\link{optim}}, \code{\link{dnorm}}
#'
#' @examples
#' n<-getNormFromCI(qLow=1.08,qUpp=8.92)
#' print(n$pars) # the fitted parameter values (mean & sd)
#' n$r(10) # 10 random values from the fitted normal distribution
#' n$d(6) # the probability density at x=6 for the normal distribution
#' n$p(4.25) # the cumulative density at x=4.25 for the fitted normal distribution
#' n$q(c(0.25,0.5,0.75)) # the 25th, 50th (median) and 75th percentiles of the fitted distribution
#' x<-seq(0,10,length=1e3)
#' y<-n$d(x)
#' plot(x,y,type="l",xlab="",ylab="density") # density plot for the fitted normal distribution
#'
#' @export getNormFromCI

getNormFromCI<-function(qLow,qUpp,alpha=0.05,initPars=c(0,1),maxiter=1e3){
  pars<-identifyNormPars(qLow,qUpp,alpha=alpha,initPars=initPars,maxiter=maxiter)

  rFun<-function(n){rnorm(n,mean=pars[1],sd=pars[2])}
  dFun<-function(x){dnorm(x,mean=pars[1],sd=pars[2])}
  pFun<-function(q){pnorm(q,mean=pars[1],sd=pars[2])}
  qFun<-function(p){qnorm(p,mean=pars[1],sd=pars[2])}

  res<-list(r=rFun,d=dFun,q=qFun,p=pFun,pars=pars)
  return(res)
}


#' @title Find the best-fit Poisson distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit Poisson distribution for a given confidence interval; returns the corresponding probability mass, distribution, quantile and sampling functions.
#' The use of this function within the bootComb package is limited: this is a discrete distribution but since users provide confidence intervals, the corresponding parameters will be best approximated by continuous distributions.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 1 giving the initial parameter value (rate parameter) to start the optimisation; defaults to 5.
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A list with 5 elements:
#' \item{r}{The sampling function.}
#' \item{d}{The probability mass function.}
#' \item{p}{The distribution function.}
#' \item{q}{The quantile function.}
#' \item{pars}{A single number giving the rate parameter for the best-fit Poisson distribution (\code{lambda} as in \code{\link{rpois}}, \code{\link{dpois}}, \code{\link{ppois}}, \code{\link{qpois}}).}
#'
#' @seealso
#' \code{\link{identifyPoisPars}}, \code{\link{optim}}, \code{\link{dpois}}
#'
#' @examples
#' n<-getPoisFromCI(qLow=9,qUpp=22)
#' print(n$par) # the fitted parameter value (lambda)
#' n$r(10) # 10 random values from the fitted Poisson distribution
#' n$d(6) # the probability mass at x=6 for the Poisson distribution
#' n$p(7) # the cumulative probability at x=7 for the fitted Poisson distribution
#' n$q(c(0.25,0.5,0.75)) # the 25th, 50th (median) and 75th percentiles of the fitted distribution
#' x<-0:40
#' y<-n$d(x)
#' barplot(height=y,names.arg=x,xlab="",ylab="probability mass") # bar plot of the fitted Poisson pmf
#'
#' @export getPoisFromCI

getPoisFromCI<-function(qLow,qUpp,alpha=0.05,initPars=5,maxiter=1e3){
  pars<-identifyPoisPars(qLow,qUpp,alpha=alpha,initPars=initPars,maxiter=maxiter)

  rFun<-function(n){rpois(n,lambda=pars)}
  dFun<-function(x){dpois(x,lambda=pars)}
  pFun<-function(q){ppois(q,lambda=pars)}
  qFun<-function(p){qpois(p,lambda=pars)}

  res<-list(r=rFun,d=dFun,q=qFun,p=pFun,par=pars)
  return(res)
}


#' @title Find the best-fit negative binomial distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit negative binomial distribution for a given confidence interval; returns the corresponding probability mass, distribution, quantile and sampling functions.
#' The use of this function within the bootComb package is limited: this is a discrete distribution but since users provide confidence intervals, the corresponding parameters will be best approximated by continuous distributions.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values (size & prob) to start the optimisation; defaults to c(10,0.5).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A list with 5 elements:
#' \item{r}{The sampling function.}
#' \item{d}{The probability mass function.}
#' \item{p}{The distribution function.}
#' \item{q}{The quantile function.}
#' \item{pars}{A vector of length 2 giving the mean and standard deviation for the best-fit negatove binomial distribution (\code{size} and \code{prob} as in \code{\link{rnbinom}}, \code{\link{dnbinom}}, \code{\link{pnbinom}}, \code{\link{qnbinom}}).}
#'
#' @seealso
#' \code{\link{identifyNegBinPars}}, \code{\link{optim}}, \code{\link{dnbinom}}
#'
#' @examples
#' n<-getNegBinFromCI(qLow=1.96,qUpp=19.12)
#' print(n$pars) # the fitted parameter values (size & prob)
#' n$r(10) # 10 random values from the fitted negative binomial distribution
#' n$d(8) # the probability mass at x=8 for the negative binomial distribution
#' n$p(12) # the cumulative probability at x=12 for the fitted negative binomial distribution
#' n$q(c(0.25,0.5,0.75)) # the 25th, 50th (median) and 75th percentiles of the fitted distribution
#' x<-0:30
#' y<-n$d(x)
#' barplot(height=y,names.arg=x,xlab="",ylab="probability mass") # bar plot of the fitted neg. bin. pmf
#'
#' @export getNegBinFromCI

getNegBinFromCI<-function(qLow,qUpp,alpha=0.05,initPars=c(10,0.5),maxiter=1e3){
  pars<-identifyNegBinPars(qLow,qUpp,alpha=alpha,initPars=initPars,maxiter=maxiter)

  rFun<-function(n){rnbinom(n,size=pars[1],prob=pars[2])}
  dFun<-function(x){dnbinom(x,size=pars[1],prob=pars[2])}
  pFun<-function(q){pnbinom(q,size=pars[1],prob=pars[2])}
  qFun<-function(p){qnbinom(p,size=pars[1],prob=pars[2])}

  res<-list(r=rFun,d=dFun,q=qFun,p=pFun,pars=pars)
  return(res)
}


#' @title Find the best-fit gamma distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit gamma distribution for a given confidence interval; returns the corresponding density, distribution, quantile and sampling functions.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values (shape & rate) to start the optimisation; defaults to c(1,1).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A list with 5 elements:
#' \item{r}{The sampling function.}
#' \item{d}{The density function.}
#' \item{p}{The distribution function.}
#' \item{q}{The quantile function.}
#' \item{pars}{A vector of length 2 giving the shape and rate for the best-fit gamma distribution (\code{shape} and \code{rate} as in \code{\link{rgamma}}, \code{\link{dgamma}}, \code{\link{pgamma}}, \code{\link{qgamma}}).}
#'
#' @seealso
#' \code{\link{identifyGammaPars}}, \code{\link{optim}}, \code{\link{dgamma}}
#'
#' @examples
#' n<-getGammaFromCI(qLow=0.82,qUpp=5.14)
#' print(n$pars) # the fitted parameter values (shape & rate)
#' n$r(10) # 10 random values from the fitted gamma distribution
#' n$d(6) # the probability density at x=6 for the gamma distribution
#' n$p(2) # the cumulative density at x=2 for the fitted gamma distribution
#' n$q(c(0.25,0.5,0.75)) # the 25th, 50th (median) and 75th percentiles of the fitted distribution
#' x<-seq(0,8,length=1e3)
#' y<-n$d(x)
#' plot(x,y,type="l",xlab="",ylab="density") # density plot for the fitted gamma distribution
#'
#' @export getGammaFromCI

getGammaFromCI<-function(qLow,qUpp,alpha=0.05,initPars=c(1,1),maxiter=1e3){
  pars<-identifyGammaPars(qLow,qUpp,alpha=alpha,initPars=initPars,maxiter=maxiter)

  rFun<-function(n){rgamma(n,shape=pars[1],rate=pars[2])}
  dFun<-function(x){dgamma(x,shape=pars[1],rate=pars[2])}
  pFun<-function(q){pgamma(q,shape=pars[1],rate=pars[2])}
  qFun<-function(p){qgamma(p,shape=pars[1],rate=pars[2])}

  res<-list(r=rFun,d=dFun,q=qFun,p=pFun,pars=pars)
  return(res)
}


#' @title Find the best-fit exponential distribution for a given confidence interval.
#'
#' @description
#' Finds the best-fit exponential distribution for a given confidence interval; returns the corresponding density, distribution, quantile and sampling functions.
#'
#' @param qLow The observed lower quantile.
#' @param qUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A single number giving the initial rate parameter value to start the optimisation; defaults to 1.
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A list with 5 elements:
#' \item{r}{The sampling function.}
#' \item{d}{The density function.}
#' \item{p}{The distribution function.}
#' \item{q}{The quantile function.}
#' \item{pars}{A single number giving the rate parameter for the best-fit exponential distribution (\code{rate} as in \code{\link{rexp}}, \code{\link{dexp}}, \code{\link{pexp}}, \code{\link{qexp}}).}
#'
#' @seealso
#' \code{\link{identifyExpPars}}, \code{\link{optim}}, \code{\link{dexp}}
#'
#' @examples
#' n<-getExpFromCI(qLow=0.01,qUpp=1.75)
#' print(n$pars) # the fitted rate parameter value
#' n$r(10) # 10 random values from the fitted exponential distribution
#' n$d(2) # the probability density at x=2 for the exponential distribution
#' n$p(1.5) # the cumulative density at x=1.5 for the fitted exponential distribution
#' n$q(c(0.25,0.5,0.75)) # the 25th, 50th (median) and 75th percentiles of the fitted distribution
#' x<-seq(0,5,length=1e3)
#' y<-n$d(x)
#' plot(x,y,type="l",xlab="",ylab="density") # density plot for the fitted exponential distribution
#'
#' @export getExpFromCI

getExpFromCI<-function(qLow,qUpp,alpha=0.05,initPars=1,maxiter=1e3){
  pars<-identifyExpPars(qLow,qUpp,alpha=alpha,initPars=initPars,maxiter=maxiter)

  rFun<-function(n){rexp(n,rate=pars)}
  dFun<-function(x){dexp(x,rate=pars)}
  pFun<-function(q){pexp(q,rate=pars)}
  qFun<-function(p){qexp(p,rate=pars)}

  res<-list(r=rFun,d=dFun,q=qFun,p=pFun,pars=pars)
  return(res)
}
