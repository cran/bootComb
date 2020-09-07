
#' @title Simulation scenario for the product of two prevlaence estimates.
#'
#' @description
#' This is a simulation to compute the coverage of the confidence interval returned by bootComb() in the case of the product of 2 probability parameter estimates.
#'
#' @param B The number of simulations to run. Defaults to 1e3.
#' @param p1 The true value of the first probability parameter.
#' @param p2 The true value of the second probability parameter.
#' @param nExp1 The size of each simulated experiment to estimate \code{p1}.
#' @param nExp2 TThe size of each simulated experiment to estimate \code{p2}.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#'
#' @return A list with 2 elements:
#'    \item{estimate}{A single number, the proportion of simulations for which the confidence interval contained the true parameter value.}
#'    \item{conf.int}{A 95\% confidence interval for the coverage estimate.}
#'
#' @examples
#' \donttest{
#' simScenProductTwoPrevs(p1=0.35,p2=0.2,nExp1=100,nExp2=1000,B=100)
#'   # B value only for convenience here
#'   # Increase B to 1e3 or 1e4 (be aware this may run for some time).
#'  }
#'
#' @export simScenProductTwoPrevs

simScenProductTwoPrevs<-function(B=1e3,p1,p2,nExp1,nExp2,alpha=0.05){
  trueP<-p1*p2
  res<-rep(0,B)

  for(j in 1:B){
    exp1<-rbinom(1,size=nExp1,prob=p1)
    exp2<-rbinom(1,size=nExp2,prob=p2)

    p1Est<-binom.test(x=exp1,n=nExp1)
    p2Est<-binom.test(x=exp2,n=nExp2)

    betaDist1<-getBetaFromCI(qLow=p1Est$conf.int[1],qUpp=p1Est$conf.int[2],initPars=c(round(nExp1/2),round(nExp1/2)))
    betaDist2<-getBetaFromCI(qLow=p2Est$conf.int[1],qUpp=p2Est$conf.int[2],initPars=c(round(nExp2/2),round(nExp2/2)))

    distListEx<-list(betaDist1$r,betaDist2$r)
    combFunEx<-function(pars){pars[[1]]*pars[[2]]}

    bsOut<-bootComb(distList=distListEx,combFun=combFunEx)

    if(trueP>=bsOut$conf.int[1] & trueP<=bsOut$conf.int[2]){res[j]<-1}
  }

  tmp<-binom.test(x=sum(res),n=length(res))
  coverage<-list(estimate=as.numeric(tmp$estimate),conf.int=tmp$conf.int) # ideally should be 95%

  return(coverage)
}

#' @title Simulation scenario for adjusting a prevalence for sensitivity and specificity.
#'
#' @description
#' This is a simulation to compute the coverage of the confidence interval returned by bootComb() in the case of adjusting a prevalence estimate for estimates of sensitivity and specificity.
#'
#' @param B The number of simulations to run. Defaults to 1e3.
#' @param p The true value of the prevalence parameter.
#' @param sens The true value of the assay sensitivity parameter.
#' @param spec The true value of the assay specificity parameter
#' @param nExp The size of each simulated experiment to estimate \code{p}.
#' @param nExpSens The size of each simulated experiment to estimate \code{sens}.
#' @param nExpSpec The size of each simulated experiment to estimate \code{spec}.
#' @param alpha The confidence level; i.e. the desired coverage is 1-alpha. Defaults to 0.05.
#' @param assumeSensSpecExact Logical; indicates whether coverage should also be computed for the situation where sensitivity and specificity are assumed to be known exactly. Defaults to FALSE.
#'
#' @return A list with 2 or 4 elements, depending whether \code{assumeSensSpecExact} is set to FALSE or TRUE:
#'    \item{estimate}{A single number, the proportion of simulations for which the confidence interval contained the true prevalence parameter value.}
#'    \item{conf.int}{A confidence interval of coverage 1-alpha for the coverage estimate.}
#'    \item{estimate.sensSpecExact}{Returned only if \code{assumeSensSpecExact} is set to TRUE. A single number, the proportion of simulations for which the confidence interval, derived assuming sensitivity and specificity are known exactly, contained the true prevalence parameter value.}
#'    \item{conf.int.sensSpecExact}{Returned only if \code{assumeSensSpecExact} is set to TRUE. A confidence interval of coverage 1-alpha for the coverage estimate in the scenario where sensitivity and specificity are assumed to be known exactly.}
#'
#' @examples
#' \donttest{
#' simScenPrevSensSpec(p=0.15,sens=0.85,spec=0.90,nExp=300,nExpSens=600,nExpSpec=400,B=1000)
#'   # B value only for convenience here
#'   # Increase B to 1e3 or 1e4 (be aware this may run for some time).
#'  }
#'
#' @export simScenPrevSensSpec

simScenPrevSensSpec<-function(B=1e3,p,sens,spec,nExp,nExpSens,nExpSpec,alpha=0.05,assumeSensSpecExact=FALSE){
  trueObsPrev<-p*sens+(1-p)*(1-spec)
  res<-rep(0,B)
  if(assumeSensSpecExact){res2<-rep(0,B)}

  for(j in 1:B){
    expP<-rbinom(1,size=nExp,prob=trueObsPrev)
    expSens<-rbinom(1,size=nExpSens,prob=sens)
    expSpec<-rbinom(1,size=nExpSpec,prob=spec)

    pEst<-binom.test(x=expP,n=nExp)
    sensEst<-binom.test(x=expSens,n=nExpSens)
    specEst<-binom.test(x=expSpec,n=nExpSpec)

    bsOut<-adjPrevSensSpecCI(prevCI=pEst$conf.int,sensCI=sensEst$conf.int,specCI=specEst$conf.int,alpha=alpha,method="hdi")

    if(p>=bsOut$conf.int[1] & p<=bsOut$conf.int[2]){res[j]<-1}

    if(assumeSensSpecExact){
      ciTmp<-c((pEst$conf.int[1]-(1-specEst$estimate))/(sensEst$estimate-(1-specEst$estimate)),(pEst$conf.int[2]-(1-specEst$estimate))/(sensEst$estimate-(1-specEst$estimate)))
      if(p>=ciTmp[1] & p<=ciTmp[2]){res2[j]<-1}
    }
  }

  tmp<-binom.test(x=sum(res),n=length(res))
  coverage<-list(estimate=as.numeric(tmp$estimate),conf.int=tmp$conf.int) # ideally should be 95%

  if(assumeSensSpecExact){
    tmp<-binom.test(x=sum(res2),n=length(res2))
    coverage$estimate.sensSpecExact<-as.numeric(tmp$estimate)
    coverage$conf.int.sensSpecExact<-tmp$conf.int
  }

  return(coverage)
}
