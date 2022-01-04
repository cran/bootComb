# bootComb

[![](https://cranlogs.r-pkg.org/badges/bootComb)](https://CRAN.R-project.org/package=bootComb)
[![](https://cranlogs.r-pkg.org/badges/grand-total/bootComb)](https://CRAN.R-project.org/package=bootComb)

Deriving confidence intervals for combinations of independent parameter estimates.

## Citation

Please use the following to cite the use of this R package in publications:

*  Henrion, M.Y.R. (2021). bootComb - An R Package to Derive Confidence Intervals for Combinations of Independent Parameter Estimates. International
  Journal of Epidemiology, dyab049. https://doi.org/10.1093/ije/dyab049

## Description

This is an R package to combine several estimates via a arbitrary function to compute a new parameter.

Uncertainty is propagated through by sampling from probability distributions for each estimate and then deriving an empirical interval containing 95% of the resulting distribution (either via highest density interval or quantiles).

Several specific applications, though the methodology is quite general:
1. Combining several conditional prevalences (see e.g. [Stockdale et al., J. Hepatol. (2020)](https://doi.org/10.1016/j.jhep.2020.04.008))
2. Adjusting a prevalence for assay sensitivity and specificity when the latter are not known exactly and only estimates are available.

By default, the method assumes that the combined estimates are independent. This is generally the case where parameters estimated on different, reasonably large datasets are combined, but in some situations this may not hold (e.g. when parameters themselves are not independent).
Where this assumption is likely not met, users can specify a correlation matrix for the input parameters. The package uses a Gaussian copula function to implement the dependence structure.

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

95% confidence interval for the product of 2 prevalence parameters for which only the 95% confidence intervals are known. The `seed` argument is not needed - it just allows exact reproduction of the example.

``` r
library(bootComb)

dist1<-getBetaFromCI(qLow=0.2,qUpp=0.8,alpha=0.05)
dist2<-getBetaFromCI(qLow=0.4,qUpp=0.9,alpha=0.05)

distListEx<-list(dist1$r,dist2$r)
combFunEx<-function(pars){pars[[1]]*pars[[2]]}

bootComb(distList=distListEx,combFun=combFunEx,method="hdi",seed=123)
#> 
#> $conf.int
#>     lower     upper 
#> 0.1034161 0.5854474  
#> attr(,"credMass")
#> [1] 0.95
#> 
#> $bootstrapValues
#> NULL
```

Alternatively, the same result can also be arrived at in just 2 lines of code:

```r
combFunEx<-function(pars){pars[[1]]*pars[[2]]}
bootComb(distributions=c("beta","beta"),qLowVect=c(0.2,0.4),qUppVect=c(0.8,0.9),combFun=combFunEx,method="hdi",seed=123)
#> 
#> $conf.int
#>     lower     upper 
#> 0.1034161 0.5854474  
#> attr(,"credMass")
#> [1] 0.95
#> 
#> $bootstrapValues
#> NULL
```
