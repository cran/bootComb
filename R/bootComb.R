#' @title Combine parameter estimates via bootstrap
#'
#' @description
#' This package propagates uncertainty from several estimates when combining these estimates via a function.
#' It does this by using the parametric bootstrap to simulate values from the distribution of each estimate to build up an empirical distribution of the combined parameter.
#' Finally either the percentile method is used or the highest density interval is chosen to derive a confidence interval for the combined parameter with the desired coverage.
#'
#' @param distList A list object where each element of the list is a sampling function for a probability distribution function (i.e. like rnorm, rbeta, ...).
#' @param combFun The function to combine the different estimates to a new parameter. Needs to take a single list as input argument, one element of the list for each estimate. This list input argument needs to be a list of same length as distList.
#' @param N The number of bootstrap samples to take. Defaults to 1e6.
#' @param method The method uses to derive a confidence interval from the empirical distribution of the combined parameter.Needs to be one of 'hdi' (default; computes the highest density interval) or 'quantile (uses quantiles to derive the confidence interval).
#' @param coverage The desired coverage of the resulting confidence interval.Defaults to 0.95.
#' @param doPlot Logical; indicates whether a graph should be produced showing the input distributions and the resulting empirical distribution of the combined estimate together with the reported confidence interval. Defaults to FALSE.
#' @param legPos Legend position (only used if doPlot==TRUE); either NULL (no legend) or one of "top", "topleft", "topright", "bottom", "bottomleft", "bottomright" "left", "right", "center".
#' @param returnBootVals Logical; if TRUE then the parameter values computed from the boostrapped input parameter values will be returned; defaults to FALSE.
#' @param validRange Optional; if not NULL, a vector of length 2 giving the range within which the values obtained from the bootstrapped input parameters must lie; values outside this range will be discarded.
#'
#' @return A list with 2 elements:
#' \item{conf.int}{A vector of length 2 giving the lower and upper limits of the computed confidence interval.}
#' \item{bootstrapValues}{A vector containing the computed parameter values from the bootstrap samples of the input parameters. (Only non-NULL if \code{returnBootVals} is set to TRUE)}
#'
#' @seealso
#' \code{\link[HDInterval]{hdi}}
#'
#' @examples
#' ## Example 1 - product of 2 probability parameters for which only the 95% CIs are reported
#' dist1<-getBetaFromCI(pLow=0.4,pUpp=0.6,alpha=0.05)
#' dist2<-getBetaFromCI(pLow=0.7,pUpp=0.9,alpha=0.05)
#' distListEx<-list(dist1$r,dist2$r)
#' combFunEx<-function(pars){pars[[1]]*pars[[2]]}
#' bootComb(distList=distListEx,combFun=combFunEx,doPlot=TRUE)
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
#' @export bootComb

bootComb<-function(distList,combFun,N=1e6,method="hdi",coverage=0.95,doPlot=FALSE,legPos="topright",returnBootVals=FALSE,validRange=NULL){
  pars<-vector("list",length=length(distList))

  for(k in 1:length(distList)){
    pars[[k]]<-distList[[k]](N)
  }

  combParVals<-combFun(pars)
  if(!is.null(validRange)){
    combParVals<-combParVals[combParVals>=validRange[1] & combParVals<=validRange[2]]
    if(length(combParVals)==0){stop("All sample values are outside the valid range.")}
  }

  if(method=="hdi"){
    ci<-HDInterval::hdi(combParVals,credMass=coverage)
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

  if(!returnBootVals){combParVals<-NULL}

  return(list(conf.int=ci,bootstrapValues=combParVals))
}
