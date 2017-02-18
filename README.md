# MSGARCH
Markov-switching GARCH models in R

[![CRAN](http://www.r-pkg.org/badges/version/MSGARCH)](https://cran.r-project.org/package=MSGARCH) [![Downloads](http://cranlogs.r-pkg.org/badges/MSGARCH?color=brightgreen)](http://www.r-pkg.org/pkg/MSGARCH)[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/MSGARCH?color=brightgreen)](http://www.r-pkg.org/pkg/MSGARCH)
[![Pending Pull-Requests](http://githubbadges.herokuapp.com/keblu/MSGARCH/pulls.svg?style=flat)](https://github.com/keblu/MSGARCH/pulls)
[![Github Issues](http://githubbadges.herokuapp.com/keblu/MSGARCH/issues.svg)](https://github.com/keblu/MSGARCH/issues)

## Introduction
Markov-switching GARCH models have become popular to model the structural break in the conditional variance dynamics of financial time series. The R package `MSGARCH` ([Ardia et al., 2016](https://ssrn.com/abstract=2845809), Ardia et al., 2017) implements Markov-switching GARCH-type models very effficiently by using C object-oriented programming techniques. It allows the user to perform simulations as well as Maximum Likelihood and Bayesian estimation of a very large class of Markov-switching GARCH-type models. Risk management tools such as Value-at-Risk and Expected-Shortfall calculations are available. See [Ardia et al. (2016)](https://ssrn.com/abstract=2845809) for further details.

## Contents
* MSGARCH-manual.pdf: This document is the documentation for the MSGARCH package.
* Package: This folder contains the latest developpement version of the package.
* bin: This folder contains the previous stable version of the package.
* test: Test folder with an R file containing all the examples from the documentation for test purpose.

## Installation

The latest stable version of `MSGARCH` is available on CRAN (https://CRAN.R-project.org/package=MSGARCH) and can be installed via:

      R > install.packages("MSGARCH")
  
To install the latest development  version of `MSGARCH` (which may contain bugs!) use these lines:

      R > install.packages("devtools")
      R > require("devtools")
      R > devtools::install_github("keblu/MSGARCH", subdir="Package")

## References

Please cite `MSGARCH` in publications:

Ardia, D., Bluteau, K., Boudt, K., Trottier, D.-A. (2016).  
_Markov-switching GARCH models in R: The MSGARCH package_.  
Working paper.  
https://ssrn.com/abstract=2845809

Ardia, D., Bluteau, K., Boudt, K., Peterson, B., Trottier, D.-A. (2017).    
_MSGARCH package_  

## Acknowledgements

The authors acknowledge Google for financial support via the Google Summer of Code 2016
project "MSGARCH", the International Institute of Forecasters (IIF, https://forecasters.org) and Industrielle-Alliance.