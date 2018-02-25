#' @title The R package MSGARCH
#' @description The \R package \pkg{MSGARCH} implements a comprehensive
#' set of functionalities for Markov-switching GARCH (Haas et al. 2004a) and Mixture of GARCH (Haas et al. 2004b) models,
#' This includes fitting, filtering, forecasting, and simulating.
#' Other functions related to Value-at-Risk and Expected-Shortfall are also available.\cr
#' The main functions of the package are coded
#' in \code{C++} using \pkg{Rcpp} (Eddelbuettel and Francois, 2011)
#' and \pkg{RcppArmadillo} (Eddelbuettel and Sanderson, 2014).\cr
#' \pkg{MSGARCH} focuses on the conditional variance (and higher moments) process.
#' Hence, there is no equation for the mean.
#' Therefore, you must pre-filter via AR(1) before applying the model.\cr
#' The \pkg{MSGARCH} package implements a variety of GARCH specifications together with several conditional distributions.
#' This allows for a rich modeling
#' environment for Markov-switching GARCH models. Each single-regime process
#' is a one-lag process (e.g., GARCH(1,1)).
#' When optimization is performed, we ensure that the variance in each regime is covariance-stationary
#' and strictly positive (refer to the vignette for more information).\cr
#' We refer to Ardia et al. (2017) \url{https://ssrn.com/abstract=2845809} for a detailed
#' introduction to the package and its usage.\cr
#' The authors acknowledge Google for financial support via the Google Summer of Code 2016 & 2017,
#' the International Institute of Forecasters and Industrielle-Alliance.
#' @references Ardia, D. Bluteau, K. Boudt, K. Catania, L. & Trottier, D.-A. (2017).
#' Markov-switching GARCH models in \R: The \pkg{MSGARCH} package.
#' \url{https://ssrn.com/abstract=2845809}
#' @references Eddelbuettel, D. & Francois, R. (2011).
#' \pkg{Rcpp}: Seamless \R and \code{C++} integration.
#' \emph{Journal of Statistical Software}, 40, 1-18.
#' \url{http://www.jstatsoft.org/v40/i08/}
#' @references Eddelbuettel, D. & Sanderson, C. (2014).
#' \pkg{RcppArmadillo}: Accelerating \R with high-performance \code{C++} linear algebra.
#' \emph{Computational Statistics & Data Analysis}, 71, 1054-1063.
#' \url{http://dx.doi.org/10.1016/j.csda.2013.02.005}
#' @references Haas, M. Mittnik, S. & Paolella, MS. (2004).
#' A new approach to Markov-switching GARCH models.
#' \emph{Journal of Financial Econometrics}, 2, 493-530.
#' \url{http://doi.org/10.1093/jjfinec/nbh020} 
#' @references Haas, M. Mittnik, S. & Paolella, M. S. (2004b).
#' Mixed normal conditional heteroskedasticity.
#' \emph{Journal of Financial Econometrics}, 2, 211-250.
#' \url{http://doi.org/10.1093/jjfinec/nbh009}
#' @useDynLib MSGARCH, .registration = TRUE
"_PACKAGE"

#'@name SMI
#'@docType data
#'@aliases SMI
#'@title Swiss market index dataset
#'@description
#'See Mullen et al. (2011) for a description of this dataset.
#'@usage data("SMI")
#'@format \code{zoo} object containing 2,500 observations ranging from 1990-11-12 to 2000-10-20.
#'@source \code{DEoptim} package
#'@references
#'Mullen, K.M, Ardia, D., Gil, D., Windover, D., Cline, J. (2011).
#'\code{DEoptim}: An R Package for Global Optimization by
#'Differential Evolution. \emph{Journal of Statistical Software}, 40(6), 1-26.
#'\url{http://www.jstatsoft.org/v40/i06/}
#'@keywords datasets
NULL

#'@name dem2gbp
#'@docType data
#'@aliases dem2gbp
#'@title DEM/GBP exchange rate log-returns
#'@description
#'The \code{vector} \code{dem2gbp} contains daily observations of the Deutschmark vs British Pound foreign exchange
#'rate log-returns. This dataset has been promoted as an informal benchmark for GARCH
#'time-series software validation. See McCullough and Renfro (1999), and Brooks, Burke, and Persand
#'(2001) for details. The nominal returns are expressed in percent as in Bollerslev and Ghysels
#'(1996). The sample period is from January 3, 1984, to December 31, 1991, for a total of 1974
#'observations.
#'@usage data("dem2gbp")
#'@format \code{vector} of size 1,974.
#'@references Bollerslev T., Ghysels, E. (1996)
#'Periodic autoregressive conditional heteroscedasticity.
#'\emph{Journal of Business and Economic Statistics}, 14, 139-151.
#'@references Brooks C., Burke S. P., Persand G. (2001)
# Benchmarks and the accuracy of GARCH model estimation.
#'\emph{International Journal of Forecasting}, 17, 45-57.
#'@references McCullough B. D., Renfro C. G. (1999)
#'Benchmarks and software standards: A case study of GARCH procedures. 
#'\emph{Journal of Economic and Social Measurement}, 25, 59-71.
#'@keywords datasets
NULL