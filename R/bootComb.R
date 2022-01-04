#' @title Combine parameter estimates via bootstrap
#'
#' @description
#' This package propagates uncertainty from several estimates when combining these estimates via a function.
#' It does this by using the parametric bootstrap to simulate values from the distribution of each estimate to build up an empirical distribution of the combined parameter.
#' Finally either the percentile method is used or the highest density interval is chosen to derive a confidence interval for the combined parameter with the desired coverage.
#'
#' @param distList If \code{Sigma} is set to NULL, this is a list object where each element of the list is a sampling function for a probability distribution function (i.e. like rnorm, rbeta, ...). If \code{Sigma} is specified, then this needs to be a list of quantile functions for the distributions for each parameter.
#' @param combFun The function to combine the different estimates to a new parameter. Needs to take a single list as input argument, one element of the list for each estimate. This list input argument needs to be a list of same length as distList.
#' @param N The number of bootstrap samples to take. Defaults to 1e6.
#' @param distributions Alternatively to specifying \code{distlist}, the parameters \code{distributions}, \code{qLowVect}, \code{qUppVect} and (optionally) \code{alphaVect} can be specified. The first 3 of these need t be either all specified and be vectors of the same length or all set to NULL. The \code{distributions} parameter needs to be a vector specifying the names of the distributions for each parameter (one of "beta", "exponential", "gamma", "normal", "Poisson" or "NegativeBinomial").
#' @param qLowVect Alternatively to specifying \code{distlist}, the parameters \code{distributions}, \code{qLowVect}, \code{qUppVect} and (optionally) \code{alphaVect} can be specified. The first 3 of these need t be either all specified and be vectors of the same length or all set to NULL. The \code{qLowVect} parameter needs to be a vector specifying the lower confidence interval limits for each parameter.
#' @param qUppVect Alternatively to specifying \code{distlist}, the parameters \code{distributions}, \code{qLowVect}, \code{qUppVect} and (optionally) \code{alphaVect} can be specified. The first 3 of these need t be either all specified and be vectors of the same length or all set to NULL. The \code{qUppVect} parameter needs to be a vector specifying the upper confidence interval limits for each parameter.
#' @param alphaVect Alternatively to specifying \code{distlist}, the parameters \code{distributions}, \code{qLowVect}, \code{qUppVect} and (optionally) \code{alphaVect} can be specified. The first 3 of these need t be either all specified and be vectors of the same length or all set to NULL. The \code{alphaVect} parameter needs to be a vector specifying the alpha level (i.e. 1 minus the coverage) of each confidence interval. Can be specified as a single number if the same for all parameters. Defaults to 0.05.
#' @param Sigma Set to NULL if parameters are assumed to be independent (the default). If specified, this needs to be a valid covariance matrix for a multivariate normal distribution with variances equal to 1 for all variables (in other words, this really is a correlation matrix).
#' @param method The method uses to derive a confidence interval from the empirical distribution of the combined parameter.Needs to be one of 'quantile' (default; uses the percentile method to derive the confidence interval) or hdi' (computes the highest density interval).
#' @param coverage The desired coverage of the resulting confidence interval.Defaults to 0.95.
#' @param doPlot Logical; indicates whether a graph should be produced showing the input distributions and the resulting empirical distribution of the combined estimate together with the reported confidence interval. Defaults to FALSE.
#' @param legPos Legend position (only used if doPlot==TRUE); either NULL (no legend) or one of "top", "topleft", "topright", "bottom", "bottomleft", "bottomright" "left", "right", "center".
#' @param returnBootVals Logical; if TRUE then the parameter values computed from the bootstrapped input parameter values will be returned; values for the individual parameters will be reported as a second list element; defaults to FALSE.
#' @param validRange Optional; if not NULL, a vector of length 2 giving the range within which the values obtained from the bootstrapped input parameters must lie; values outside this range will be discarded. Behaviour that results in the need for this option arises when parameters are not independent. Use with caution.
#' @param seed If desired a random seed can be specified so that the same results can be reproduced.
#'
#' @return A list with 3 elements:
#' \item{conf.int}{A vector of length 2 giving the lower and upper limits of the computed confidence interval.}
#' \item{bootstrapValues}{A vector containing the computed / combined parameter values from the bootstrap samples of the input parameters. (Only non-NULL if \code{returnBootVals} is set to TRUE.)}
#' \item{bootstrapValuesInput}{A list where each element is the vector of the bootstrapped values for the corresponding input parameter. This can be useful to check the dependence structure that was specified. (Only non-NULL if \code{returnBootVals} is set to TRUE.)}
#'
#' @seealso
#' \code{\link[HDInterval]{hdi}}
#'
#' @examples
#' ## Example 1 - product of 2 probability parameters for which only the 95% CIs are reported
#' dist1<-getBetaFromCI(qLow=0.4,qUpp=0.6,alpha=0.05)
#' dist2<-getBetaFromCI(qLow=0.7,qUpp=0.9,alpha=0.05)
#' distListEx<-list(dist1$r,dist2$r)
#' combFunEx<-function(pars){pars[[1]]*pars[[2]]}
#' bootComb(distList=distListEx,
#'          combFun=combFunEx,
#'          doPlot=TRUE,
#'          method="hdi",
#'          N=1e5, # reduced from N=1e6 so that it runs quicker; larger values => more accurate
#'          seed=352)
#'
#' # Alternatively, the same example can be run in just 2 lines of code:
#' combFunEx<-function(pars){pars[[1]]*pars[[2]]}
#' bootComb(distributions=c("beta","beta"),
#'          qLowVect=c(0.4,0.7),
#'          qUppVect=c(0.6,0.9),
#'          combFun=combFunEx,
#'          doPlot=TRUE,
#'          method="hdi",
#'          N=1e5, # reduced from N=1e6 so that it runs quicker; larger values => more accurate
#'          seed=352)
#'
#' ## Example 2 - sum of 3 Gaussian distributions
#' dist1<-function(n){rnorm(n,mean=5,sd=3)}
#' dist2<-function(n){rnorm(n,mean=2,sd=2)}
#' dist3<-function(n){rnorm(n,mean=1,sd=0.5)}
#' distListEx<-list(dist1,dist2,dist3)
#' combFunEx<-function(pars){pars[[1]]+pars[[2]]+pars[[3]]}
#' bootComb(distList=distListEx,combFun=combFunEx,doPlot=TRUE,method="quantile")
#'
#' # Compare with theoretical result:
#' exactCI<-qnorm(c(0.025,0.975),mean=5+2+1,sd=sqrt(3^2+2^2+0.5^2))
#' print(exactCI)
#' x<-seq(-10,30,length=1e3)
#' y<-dnorm(x,mean=5+2+1,sd=sqrt(3^2+2^2+0.5^2))
#' lines(x,y,col="red")
#' abline(v=exactCI[1],col="red",lty=3)
#' abline(v=exactCI[2],col="red",lty=3)
#'
#' ## Example 3 - same as Example 1 but assuming the 2 parameters to be dependent / correlated
#' combFunEx<-function(pars){pars[[1]]*pars[[2]]}
#' bootComb(distributions=c("beta","beta"),
#'          qLowVect=c(0.4,0.7),
#'          qUppVect=c(0.6,0.9),
#'          Sigma=matrix(byrow=2,ncol=2,c(1,0.5,0.5,1)),
#'          combFun=combFunEx,
#'          doPlot=TRUE,
#'          method="hdi",
#'          N=1e5, # reduced from N=1e6 so that it runs quicker; larger values => more accurate
#'          seed=352)
#'
#' @export bootComb

bootComb<-function(distList,combFun,N=1e6,distributions=NULL,qLowVect=NULL,qUppVect=NULL,alphaVect=0.05,Sigma=NULL,method="quantile",coverage=0.95,doPlot=FALSE,legPos="topright",returnBootVals=FALSE,validRange=NULL,seed=NULL){
  if(!is.null(distributions) & !is.null(qLowVect) & !is.null(qUppVect)){
    K<-length(distributions)
    if(K!=length(qLowVect) | K!=length(qUppVect)){
      stop("Arguments 'distributions', 'qLowVect' and 'qUppVect' need to be vectors of the same length.")
    }
    if(length(alphaVect)!=1 & length(alphaVect)!=K){
      stop("Argument 'alphaVect' needs to be either a vector of length 1 (a single number) or a vector of same length as arguments 'distributions', 'qLowVect' and 'qUppVect'.")
    }
    if(sum(alphaVect<0 | alphaVect>1)>0){
      stop("All elements of argument 'alphaVect' need to be values between 0 and 1.")
    }
    if(length(alphaVect)==1){alphaVect<-rep(alphaVect,K)}

    distList<-vector("list",K)

    for(k in 1:K){
      if(tolower(distributions[k])=="beta"){distTmp<-getBetaFromCI(qLow=qLowVect[k],qUpp=qUppVect[k],alpha=alphaVect[k])}
      if(tolower(distributions[k])=="exponential"){distTmp<-getExpFromCI(qLow=qLowVect[k],qUpp=qUppVect[k],alpha=alphaVect[k])}
      if(tolower(distributions[k])=="gamma"){distTmp<-getGammaFromCI(qLow=qLowVect[k],qUpp=qUppVect[k],alpha=alphaVect[k])}
      if(tolower(distributions[k])=="normal"){distTmp<-getNormFromCI(qLow=qLowVect[k],qUpp=qUppVect[k],alpha=alphaVect[k])}
      if(tolower(distributions[k])=="poisson"){distTmp<-getPoisFromCI(qLow=qLowVect[k],qUpp=qUppVect[k],alpha=alphaVect[k])}
      if(tolower(distributions[k])=="negativebinomial"){distTmp<-getNegBinFromCI(qLow=qLowVect[k],qUpp=qUppVect[k],alpha=alphaVect[k])}

      if(is.null(Sigma)){
        distList[[k]]<-distTmp$r
      }else{
        distList[[k]]<-distTmp$q
      }
    }
  }else if(!is.null(distributions) | !is.null(qLowVect) | !is.null(qUppVect)){
    stop("Either all of the arguments 'distributions', 'qLowVect' and 'qUppVect' have to be specified or all have to be set to NULL.")
  }

  pars<-vector("list",length=length(distList))

  if(!is.null(seed)){set.seed(seed)}

  if(is.null(Sigma)){
    for(k in 1:length(distList)){
      pars[[k]]<-distList[[k]](N)
    }
  }else{
    if(sum(diag(Sigma)!=rep(1,nrow(Sigma)))>0){stop("All diagonal elements of the Sigma matrix need to be equal to 1.")}
    parsNorm<-MASS::mvrnorm(n=N,mu=rep(0,k),Sigma=Sigma)
    parsNorm<-pnorm(parsNorm)
    for(k in 1:length(distList)){
      pars[[k]]<-distList[[k]](parsNorm[,k])
    }
    rm(parsNorm)
  }

  combParVals<-combFun(pars)
  if(!is.null(validRange)){
    combParVals<-combParVals[combParVals>=validRange[1] & combParVals<=validRange[2]]
    if(length(combParVals)==0){stop("All sample values are outside the valid range.")}
  }

  if(method=="hdi"){
    if(requireNamespace("HDInterval",quietly=TRUE)){
      ci<-HDInterval::hdi(combParVals,credMass=coverage)
    }else{
      warning("The bootComb argument 'method' was set to 'hdi', but the required R package 'HDInterval is not installed. Falling back on method='quantile'.")
      ci<-quantile(combParVals,probs=c((1-coverage)/2,(1+coverage)/2))
    }
  }else if(method=="quantile"){
    ci<-quantile(combParVals,probs=c((1-coverage)/2,(1+coverage)/2))
  }else{
    stop(paste(sep=" ","The argument 'method' to function 'bootComb' needs to be one of 'hdi', 'quantile' (received <",method,">)."))
  }

  if(doPlot){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    layout(matrix(c(1:length(distList),rep(length(distList)+1,length(distList))),nrow=2,byrow=TRUE))

    for(k in 1:length(distList)){
      hist(pars[[k]],breaks=100,xlab=names(distList)[k],ylab="density",freq=FALSE,main=paste(sep="","Histogram - sampled values for parameter ",k,"."))
    }

    hist(combParVals,breaks=100,freq=FALSE,xlab="combined parameter",ylab="density",main="Empirical distribution of the combined parameter.")
    abline(v=ci[1],lty=2,lwd=2,col="steelblue")
    abline(v=ci[2],lty=2,lwd=2,col="steelblue")
    if(!is.null(legPos)){
      legend(x=legPos,bty="n",horiz=TRUE,lwd=2,lty=2,col=c("steelblue"),legend=c("95% confidence limits"))
    }
  }

  if(!returnBootVals){
    combParVals<-NULL
    pars<-NULL
  }

  return(list(conf.int=ci,bootstrapValues=combParVals,bootstrapValuesInput=pars))
}
