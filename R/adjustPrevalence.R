#' @title Adjust a prevalence point estimate for a given assay sensitivity and specificity.
#'
#' @description
#' Given a reported prevalence estimate from an imperfect assay with known sensitivity and specificity, this function will adjust the prevalence point estimate for the assay sensitivity and specificity.
#'
#' @param prevEst The reported prevalence point estimate.
#' @param sens The known assay sensitivity.
#' @param spec The known assay specificity.
#' @param replaceImpossibleValues Logical; not all combinations of prevalence, sensitivity and specificity are possible and it can be that the adjusted prevalence is <0 or >1, so if this parameter is set to TRUE, values below 0 are set to 0, values above 1 to 1. Default to FALSE.
#'
#' @return A vector of the same length as prevEst, returning the adjusted prevalence estimates.
#'
#' @seealso
#' \code{\link{adjPrevSensSpecCI}}, \code{\link{ssBetaPars}}, \code{\link{optim}}, \code{\link{dbeta}}
#'
#' @examples
#' adjPrevSensSpec(prevEst=0.16,sens=0.90,spec=0.95)
#'
#' @export adjPrevSensSpec

adjPrevSensSpec<-function(prevEst,sens,spec,replaceImpossibleValues=FALSE){
  if((length(sens)!=length(prevEst) & length(sens)!=1) | (length(spec)!=length(prevEst) & length(spec)!=1)){stop("sens and spec need to be either vectors of the same length as prevEst or scalars / vectors or length 1.")}

  prevAdj<-(prevEst-(1-spec))/(sens-(1-spec))

  if(replaceImpossibleValues){
    prevAdj[prevAdj<0]<-0
    prevAdj[prevAdj>1]<-1
  }

  return(prevAdj)
}

#' @title Adjust a prevalence point estimate and confidence interval for a given assay sensitivity and specificity (also known only imprecisely).
#'
#' @description
#' This function takes as input a prevalence confidence interval, a sensitivity confidence interval and a specificity confidence interval and returns a confidence interval with the desired coverage of the adjusted prevalence.
#' Optionally the point estimates of prevalence, sensitivity and specificity can also be specified and, if so, these will be returned together with the confidence interval.
#' This function will automatically replace impossible point estimate values with 0 (if estimate <0) or 1 (if estimate >1) and also update the lower, repsectively upper confidence interval limit in this case.
#'
#' @param prevCI A vector of length 2 giving the lower and upper bounds of the confidence interval for the prevalence estimate.
#' @param sensCI A vector of length 2 giving the lower and upper bounds of the confidence interval for the assay sensitivity estimate.
#' @param specCI A vector of length 2 giving the lower and upper bounds of the confidence interval for the assay specificity estimate.
#' @param N A (large) integer giving the number of parametric bootstrap samples to take. Defaults to 1e6.
#' @param method The method uses to derive a confidence interval from the empirical distribution of the combined parameter. Needs to be one of 'hdi' (default; computes the highest density interval) or 'quantile (uses quantiles to derive the confidence interval).
#' @param alpha The desired confidence level; i.e. the returned confidence interval will have coverage 1-alpha.
#' @param doPlot Logical; indicates whether a graph should be produced showing the input estimated distributions for the prevalence, sensitivity and specificity estimates and the resulting empirical distribution of the adjusted prevalence together with the reported confidence interval. Defaults to FALSE.
#' @param prev Optional; if not NULL, and parameters \code{sens} and \code{spec} are also not NULL, then an adjusted point estimate will also be calculated.
#' @param sens Optional; if not NULL, and parameters \code{prev} and \code{spec} are also not NULL, then an adjusted point estimate will also be calculated.
#' @param spec Optional; if not NULL, and parameters \code{prev} and \code{sens} are also not NULL, then an adjusted point estimate will also be calculated.
#' @param ylim Optional; a vector of length 2, giving the vertical limits for the top panel of the produced plot. Only used if \code{doPlot} is set to \code{TRUE}.
#'
#' @return A list object with 2 elements:
#' \item{estimate}{The adjusted prevalence point estimate (only non-NULL if \code{prev}, \code{sens} and \code{spec} are specified).}
#' \item{conf.int}{The confidence interval for the adjusted prevalence.}
#'
#' @seealso
#' \code{\link{bootComb}}, \code{\link{adjPrevSensSpec}}, \code{\link{identifyBetaPars}}, \code{\link{dbeta}}, \code{\link[HDInterval]{hdi}}
#'
#' @examples
#' adjPrevSensSpecCI(
#'   prevCI=binom.test(x=84,n=500)$conf.int,
#'   sensCI=binom.test(x=238,n=270)$conf.int,
#'   specCI=binom.test(x=82,n=88)$conf.int,
#'   doPlot=TRUE,
#'   prev=84/500,
#'   sens=238/270,
#'   spec=82/88)
#'
#' @export adjPrevSensSpecCI

adjPrevSensSpecCI<-function(prevCI,sensCI,specCI,N=1e6,method="hdi",alpha=0.05,doPlot=FALSE,prev=NULL,sens=NULL,spec=NULL,ylim=NULL){

  prevDist<-getBetaFromCI(qLow=prevCI[1],qUpp=prevCI[2],alpha=alpha)
  sensDist<-getBetaFromCI(qLow=sensCI[1],qUpp=sensCI[2],alpha=alpha)
  specDist<-getBetaFromCI(qLow=specCI[1],qUpp=specCI[2],alpha=alpha)

  distList<-list(prevDist$r,sensDist$r,specDist$r)
  combFun<-function(pars){adjPrevSensSpec(prevEst=pars[[1]],sens=pars[[2]],spec=pars[[3]])}

  if(!is.null(prev) & !is.null(sens) & !is.null(spec)){
    adjPrev<-adjPrevSensSpec(prev,sens,spec,replaceImpossibleValues=TRUE)
  }

  adjPrevCI<-bootComb(distList=distList,combFun=combFun,N=N,method=method,coverage=1-alpha,doPlot=FALSE,legPos=NULL,returnBootVals=TRUE,validRange=c(0,1))

  if(doPlot){
    oldpar <- par(no.readonly = TRUE)
    on.exit(par(oldpar))

    par(mfrow=c(2,1))

    x<-seq(0,1,length=1000)
    yPrev<-dbeta(x,prevDist$pars[1],prevDist$pars[2])
    ySens<-dbeta(x,sensDist$pars[1],sensDist$pars[2])
    ySpec<-dbeta(x,specDist$pars[1],specDist$pars[2])
    if(is.null(ylim)){ylim<-c(0,max(c(yPrev,ySens,ySpec)))}

    plot(x,yPrev,xlim=c(0,1),ylim=ylim,xlab="probability parameter",ylab="density",type="l",col="steelblue",lwd=2,main="Estimated densities of prevalence, sensitivity and specificity from their 95% confidence intervals.")
    lines(x,ySens,lwd=2,col="orange")
    lines(x,ySpec,lwd=2,col="salmon")
    legend(x="top",col=c("steelblue","orange","salmon"),legend=c("prevalence","sensitivity","specificity"),lwd=2,horiz=TRUE,bty="n")

    hist(adjPrevCI$bootstrapValues,breaks=100,freq=FALSE,xlab="adjusted prevalence",ylab="density",xlim=c(0,1),main="Histogram of the sensitivity- & specificity-adjusted prevalence")
    box()
    abline(v=adjPrevCI$conf.int[1],lty=2,lwd=2,col="steelblue")
    abline(v=adjPrevCI$conf.int[2],lty=2,lwd=2,col="steelblue")
    if(!is.null(prev) & !is.null(sens) & !is.null(spec)){
      abline(v=adjPrev,lty=2,lwd=2,col="orange")
    }
    legend(x="top",bty="n",horiz=TRUE,lwd=2,lty=2,col=c("steelblue","orange"),legend=c("95% confidence limits","point estimate"))
  }

  if(!is.null(prev) & !is.null(sens) & !is.null(spec)){
    if(adjPrev<adjPrevCI$conf.int[1]){adjPrevCI$conf.int[1]<-adjPrev}
    if(adjPrev>adjPrevCI$conf.int[2]){adjPrevCI$conf.int[2]<-adjPrev}
    res<-list(estimate=adjPrev,conf.int=adjPrevCI$conf.int)
  }else{
    res<-list(estimate=NULL,conf.int=adjPrevCI$conf.int)
  }
    return(res)
}
