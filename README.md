# bootComb

[![](https://cranlogs.r-pkg.org/badges/bootComb)](https://CRAN.R-project.org/package=bootComb)
[![](https://cranlogs.r-pkg.org/badges/grand-total/bootComb)](https://CRAN.R-project.org/package=bootComb)

Deriving confidence intervals for combinations of independent parameter estimates.

## Description

This is an R package to combine several estimates via a arbitrary function to compute a new parameter.

Uncertainty is propagated through by sampling from probability distributions for each estimate and then deriving an empirical interval containing 95% of the resulting distribution (either via highest density interval or quantiles).

Several specific applications, though the methodology is quite general:
1. Combining several conditional prevalences (see e.g. [Stockdale et al., J. Hepatol. (2020)](https://doi.org/10.1016/j.jhep.2020.04.008))
2. Adjusting a prevalence for assay sensitivity and specificity when the latter are not known exactly and only estimates are available.

The method assumes that the combined estimates are independent. This is generally the case where parameters estimated on different, reasonably large datasets are combined, but in some situations this may not hold (e.g. when parameters themselves are not independent). Future versions of this package plan to include support for some joint distributions.

## Installation

### From CRAN (stable version)

``` r
install.packages('bootComb')
```

### From GitHub (development version)

``` r
# install.packages("devtools")
devtools::install_github("gitMarcH/bootComb")
```

## Example

95% confidence interval for the product of 2 prevalence parameters for which only the 95% confidence intervals are known.

``` r
library(bootComb)

dist1<-getBetaFromCI(qLow=0.2,qUpp=0.8,alpha=0.05)
dist2<-getBetaFromCI(qLow=0.4,qUpp=0.9,alpha=0.05)

distListEx<-list(dist1$r,dist2$r)
combFunEx<-function(pars){pars[[1]]*pars[[2]]}

bootComb(distList=distListEx,combFun=combFunEx,method="hdi")
#> 
#> $conf.int
#>     lower     upper 
#> 0.1049555 0.5872140 
#> attr(,"credMass")
#> [1] 0.95
#> 
#> $bootstrapValues
#> NULL
```
