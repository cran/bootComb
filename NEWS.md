# bootComb 1.1.2

* Added the option to specify a dependence structure in `adjPrevSensSpecCI` (i.e. argument `Sigma` allowed and passed through to `bootComb`).
* Added the option to specify a seed value for `adjPrevSensSpecCI`.
* Added the option for bootstrapped values to be returned for `adjPrevSensSpecCI`.
* Added a note to the documentation files for getPoisFromCI and getNegBinFromCI - these two functions are unlikely to be of much use within bootComb (where parameters are necessarily assumed to be continuous).
* Minor fix to one example.

# bootComb 1.1.1

* Reduced the number of bootstrap samples for the example used for the `bootComb` example as CRAN gave an error that examples take over 10s to run.
* Added a seed value for the examples for the `bootComb` function.

# bootComb 1.1.0

* Support added for dependent / correlated parameters (using Gaussian copulas to specify the dependence). This is done by specifying a correlation matrix (new input argument `Sigma`).
* Added an alternative, and simpler, way of calling the main package function `bootComb` (new input arguments `distributions`, `qLowVect`, `qUppVect` and `alphaVect`).
* `returnBootVals` returns now also the values for the individual parameters that get combined. This change allows checking the dependence structure that was specified.
* A random seed can now be specified as an input argument to the `bootComb` function (argument `seed`).
* Update to README.md (example of alternative specification and citation info to avoid hard-coded reference to version number).


# bootComb 1.0.2

* Citation information updated now the IJE paper is out. Please cite https://doi.org/10.1093/ije/dyab049 when using the package.
* Updated the function `adjPrevSensSpecCI` for better behaviour in the case where the provided point estimates for the prevalence, sensitivity and specificity result in an impossible value (adjusted prevalence point estimate <0 or >1). The adjusted point estimates and boundaries of CI have now been updated to avoid impossible values and be consistent with each other (e.g. if the adjusted point estimate would be -0.000123, the lower boundary from the bootstrapping 0.000451 and the upper boundary 0.112345, then the reported point estimate and CI will be 0 [0,0.112345]).

# bootComb 1.0.1

* Fixed a bug with the internal functions `identifyBetaPars`, `identifyNormPars`, `identifyPoisPars`, `identifyNegBinPars`, `identifyGammaPars`, `identifyExpPars` whereby the `alpha` parameter did not get passed through.

# bootComb 1.0.0

With the inclusion of additional distributions functions that are supported, this is considered to be the first mature version of the package. Some changes for consistency of coding and better function naming are also minimally not backward compatible, hence the upgrade to version number 1.0.0.

* Added cranlogs badges to README.md
* Changed the way parameters are found for the Poisson distribution. Rather than quantile matching, the CDF is matched (as this will be a continuous function rather than a step function in the way the quantile function is).
* Added support for negative binomial, gamma and exonential distributions in addition to the beta, normal and Poisson distributions. Finding the best-fit negative binomial distribution proceeds similar to the procedure for the Poisson distribution as this is another discrete distribution and hence the quantile function is a step function (albeit a two-dimensional step function in this case).
* Renamed the (internal only) function 'ssPoisPar' to 'ssPoisPars' for consistency with other distributions.
* Renamed arguments 'pLow' and 'pUpp' from getBetaFromCI() to 'qLow' and 'qUpp' for consistency. Updated code from other functions where these parameters were specified.
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
