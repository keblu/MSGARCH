# MSGARCH
MSGARCH
The goal of the MSGARCH package is to implement a package that will give the financial community tools to estimate,
simulate, and test several MSGARCH models used in volatility (i.e., square root of conditional variance) forecasting.
By relying on a hidden/latent variable, these models are able to switch among several processes for the conditional
volatility and therefore, account for structural break in the volatility dynamics. MSGARCH have gained a huge interest 
in the financial risk management community over the recent years as they are better at forecasting volatility and provide
more accurate risk measures. The package follows the structure of rugarch since this is one of the most used packages
for volatility modeling. The core is implemented in C++ while simple R functions will facilitate usage of the package.

## Contents
* MSGARCH.pdf: This document is the documentation for the MSGARCH package.
* Package: This folder is The latest developpement version of the package.
* bin: This folder continas the previous stable version of the package.
* test: Test folder with an R file containing all the examples from the documentation for test purpose.
* MSGARCH_0.16.tar.gz: Latest stable version of the package.

## Installation

To install the latest stable version of the package download MSGARCH_0.16.tar.gz and use this line:

      R > install.packages("./MSGARCH_0.16.tar.gz", repos = NULL, type = "source")
  
We will eventually upload the latest stable version of the package on CRAN in the coming weeks.
To install the latest developpement version of the package use these lines:

      R > install.packages("devtools")
      R > require("devtools")
      R > devtools::install_github("keblu/MSGARCH", subdir="Package")
 
For a full explanation of the package functionallities please read the vignette located in the vignettes folder of the package.
