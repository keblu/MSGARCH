[![CRAN](http://www.r-pkg.org/badges/version/MSGARCH)](https://cran.r-project.org/package=MSGARCH) [![Downloads](http://cranlogs.r-pkg.org/badges/MSGARCH?color=brightgreen)](http://www.r-pkg.org/pkg/MSGARCH)[![Downloads](http://cranlogs.r-pkg.org/badges/grand-total/MSGARCH?color=brightgreen)](http://www.r-pkg.org/pkg/MSGARCH)

# MSGARCH
The goal of the MSGARCH package is to implement a package that will give the financial community tools to estimate,
simulate, and test several MSGARCH models used in volatility (i.e., square root of conditional variance) forecasting.
By relying on a hidden/latent variable, these models are able to switch among several processes for the conditional
volatility and therefore, account for structural break in the volatility dynamics. MSGARCH have gained a huge interest 
in the financial risk management community over the recent years as they are better at forecasting volatility and provide
more accurate risk measures. The package follows the structure of rugarch since this is one of the most used packages
for volatility modeling. The core is implemented in C++ while simple R functions will facilitate usage of the package.

## Contents
* MSGARCH-manual.pdf: This document is the documentation for the MSGARCH package.
* Package: This folder contains the latest developpement version of the package.
* bin: This folder contains the previous stable version of the package.
* test: Test folder with an R file containing all the examples from the documentation for test purpose.

## Installation

To install the latest stable version is available on CRAN https://cran.r-project.org/web/packages/MSGARCH/index.html and can be installed via:

      R > install.packages("MSGARCH")
  
To install the latest developpement version of the package (which may contain bugs) use these lines:

      R > install.packages("devtools")
      R > require("devtools")
      R > devtools::install_github("keblu/MSGARCH", subdir="Package")
 
For a full explanation of the package functionallities please read the vignette located at https://ssrn.com/abstract=2845809. 
