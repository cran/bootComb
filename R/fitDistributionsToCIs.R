#' @title Find the best-fit beta distribution for a given confidence interval for a probability parameter.
#'
#' @description
#' Finds the best-fit beta distribution for a given confidence interval for a probability parameter; returns the corresponding density, distribution, quantile and sampling functions.
#'
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values to start the optimisation; defaults to c(50,50).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A list with 5 elements:
#' \item{r}{The sampling function.}
#' \item{d}{The density function.}
#' \item{p}{The distribution function.}
#' \item{q}{The quantile function.}
#' \item{pars}{A vector of length 2 giving the two shape parameters for the best-fit normal distribution (\code{shape1} and \code{shape2} as in \code{\link{rbeta}}, \code{\link{dbeta}}, \code{\link{pbeta}}, \code{\link{qbeta}}.}
#'
#' @seealso
#' \code{\link{identifyBetaPars}}, \code{\link{optim}}, \code{\link{dbeta}}
#'
#' @examples
#' b<-getBetaFromCI(pLow=0.1167,pUpp=0.1636,initPars=c(200,800))
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

getBetaFromCI<-function(pLow,pUpp,alpha=0.05,initPars=c(50,50),maxiter=1e3){
  pars<-identifyBetaPars(pLow,pUpp,alpha=alpha,initPars=initPars,maxiter=maxiter)

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
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 2 giving the initial parameter values (mean & sd) to start the optimisation; defaults to c(0,1).
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A list with 5 elements:
#' \item{r}{The sampling function.}
#' \item{d}{The density function.}
#' \item{p}{The distribution function.}
#' \item{q}{The quantile function.}
#' \item{pars}{A vector of length 2 giving the mean and standard deviation for the best-fit normal distribution (\code{mean} and \code{sd} as in \code{\link{rnorm}}, \code{\link{dnorm}}, \code{\link{pnorm}}, \code{\link{qnorm}}.}
#'
#' @seealso
#' \code{\link{identifyNormPars}}, \code{\link{optim}}, \code{\link{dnorm}}
#'
#' @examples
#' n<-getNormFromCI(pLow=1.08,pUpp=8.92)
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

getNormFromCI<-function(pLow,pUpp,alpha=0.05,initPars=c(0,1),maxiter=1e3){
  pars<-identifyNormPars(pLow,pUpp,alpha=alpha,initPars=initPars,maxiter=maxiter)

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
#' Finds the best-fit Poisson distribution for a given confidence interval; returns the corresponding density, distribution, quantile and sampling functions.
#'
#' @param pLow The observed lower quantile .
#' @param pUpp The observed upper quantile.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param initPars A vector of length 1 giving the initial parameter value (rate parameter) to start the optimisation; defaults to 5.
#' @param maxiter Maximum number of iterations for \code{optim}. Defaults to 1e3. Set to higher values if convergence problems are reported.
#'
#' @return A list with 5 elements:
#' \item{r}{The sampling function.}
#' \item{d}{The probability mass function.}
#' \item{p}{The distribution function.}
#' \item{q}{The quantile function.}
#' \item{pars}{A vector of length 1 giving the rate parameter for the best-fit Poisson distribution (\code{lambda} as in \code{\link{rpois}}, \code{\link{dpois}}, \code{\link{ppois}}, \code{\link{qpois}}.}
#'
#' @seealso
#' \code{\link{identifyPoisPars}}, \code{\link{optim}}, \code{\link{dpois}}
#'
#' @examples
#' n<-getPoisFromCI(pLow=9,pUpp=22)
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

getPoisFromCI<-function(pLow,pUpp,alpha=0.05,initPars=5,maxiter=1e3){
  pars<-identifyPoisPars(pLow,pUpp,alpha=alpha,initPars=initPars,maxiter=maxiter)

  rFun<-function(n){rpois(n,lambda=pars)}
  dFun<-function(x){dpois(x,lambda=pars)}
  pFun<-function(q){ppois(q,lambda=pars)}
  qFun<-function(p){qpois(p,lambda=pars)}

  res<-list(r=rFun,d=dFun,q=qFun,p=pFun,par=pars)
  return(res)
}
