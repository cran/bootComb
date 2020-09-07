# bootComb 1.0.0

With the inclusion of additional distributions functions that are supported, this is considered to be the first mature version of the package. Some changes for consistency of coding and better function naming are also minimally not backward compatible, hence the upgrade to version number 1.0.0.

* Added cranlogs badges to README.md
* Changed the way parameters are found for the poisson distribution. Rather than quantile matching, the CDF is matched (as this will be a continuous function rather than a step function in the way the quantile function is).
* Added support for negative binomial, gamma and exonential distributions in addition to the beta, normal and Poisson distributions. Finding the best-fit negative binomial distribution proceeds similar to the procedure for the Poisson distribution as this is another discrete distribution and hence the quantile function is a step function (albeit a two-dimensional step function in this case).
* Renamed the (internal only) function 'ssPoisPar' to 'ssPoisPars' for consistency with other distributions.
* Ranmed arguments 'pLow' and 'pUpp' from getBetaFromCI() to 'qLow' and 'qUpp' for consistency. Updated code from other functions where these parameters were specified.
* Default method for interval computation changed from `hdi` to `quantile` as the latter will always be correct and the former may be wrong if the bootstrap sample of combined parameter values is severely multimodal.
* Dependency on package `HDInterval` turned into a "Suggests" rather than an "Imports" dependency with `bootComb` falling back on `method='quantile'` if `HDInterval` is not available.
* Small update to the example from the README.md file.
* Clarified the names for the functions to run the 2 example simulations.

# bootComb 0.2.0

First version accepted & published on CRAN.

* Various bug fixes.
* Second simulation scenario added (for adjusted prevalence estimation situation).
* Ability to specify number of optimisation iterations for `identifyBetaPars`
* `identifyBetaPars` is no longer an exported function; strictly internal
* Added `getBetaFromCI`, a function that returns beta distribution functions from a given confidence interval (and in particular a random sampling function)
* Added support for normal and Poisson distributions in addition to the beta distribution
* Changed the organisation of code within the source R code files
* Implemented changes requested from the submission of v0.1.0 to CRAN (references in DESCRIPTION file, TRUE/FALSE instead of T/F, restoring user parameters after doing plots)

# bootComb 0.1.0

First version to be submitted to CRAN (but rejected).
