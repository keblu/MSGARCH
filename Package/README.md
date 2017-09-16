# MSGARCH
Markov-switching GARCH models in R

## Introduction
Markov-switching GARCH models have become popular to model the structural break in the conditional variance dynamics of financial time series. The R package `MSGARCH` ([Ardia et al., 2016](https://ssrn.com/abstract=2845809), Ardia et al., 20xx) implements Markov-switching GARCH-type models very effficiently by using C object-oriented programming techniques. It allows the user to perform simulations as well as Maximum Likelihood and Bayesian estimation of a very large class of Markov-switching GARCH-type models. Risk management tools such as Value-at-Risk and Expected-Shortfall calculations are available. See [Ardia et al. (2016)](https://ssrn.com/abstract=2845809) for further details. A large-scale empirical study is presented in [Ardia et al. (2017)](https://ssrn.com/abstract=2918413).

## Installation

The latest stable version of `MSGARCH` is available on CRAN (https://CRAN.R-project.org/package=MSGARCH) and can be installed via:

      R > install.packages("MSGARCH")
  
To install the latest development  version of `MSGARCH` (which may contain bugs!) use these lines:

      R > install.packages("devtools")
      R > require("devtools")
      R > devtools::install_github("keblu/MSGARCH", subdir="Package")

## References

Please cite `MSGARCH` in publications:

Ardia, D., Bluteau, K., Boudt, K., Catania, L., Trottier, D.-A. (2016).  
_Markov-switching GARCH models in R: The MSGARCH package_.  
Working paper.  
https://ssrn.com/abstract=2845809

Ardia, D., Bluteau, K., Boudt, K., Catania, L. (2017).    
_Forecasting performance of Markov-switching GARCH models: A large-scale empirical study_.    
Working paper.    
https://ssrn.com/abstract=2918413  

Ardia, D., Bluteau, K., Boudt, K., Catania, L., Peterson, B., Trottier, D.-A. (20xx).    
_MSGARCH package_  

## Acknowledgements

The authors acknowledge Google for financial support via the Google Summer of Code 2016-2017
project "MSGARCH", the International Institute of Forecasters (IIF, https://forecasters.org), Industrielle-Alliance, 
FQRSC (Grant # 2015-NP-179931) and Fonds de Donations at the University of Neuchâtel, Switzerland.