#' The R package MSGARCH
#'@description The \R package \code{MSGARCH} aims to provide a comprehensive set of functionalities for
#' Markov-switching GARCH processes, including fitting, filtering, forecasting, and simulating.
#' Other functions related to  Value-at-Risk, Expected-Shortfall, and conditional distributions
#' are also available. The main functions of the package are coded in \code{C++} with \code{Rcpp} (Eddelbuettel and Francois, 2011) 
#' and \code{RcppArmadillo} (Eddelbuettel and Sanderson, 2014).
#' In the \R package \code{MSGARCH} there is no equation for the mean as in the \R package \code{rugarch} (Ghalanos, 2015). 
#' This means that we assume that before modeling, the user has filter the mean from their time series.
#' \cr We provide a variety of 
#' single-regime GARCH process and regime-switching process as well as many conditional distributions. 
#' This allows for a rich modeling environment for Markov-switching GARCH models. Each single-regime process in a one-lag process (e.g, GARCH(1,1)).
#' Allowing for only one-lag is sufficient in Finance and reduce models complexity which can become a problem during the optimization process. 
#'  
#'The authors acknowledge Google for financial support via the Google Summer of Code 2016
#' project "MSGARCH"; see \url{https://summerofcode.withgoogle.com/projects/#6497774455488512}
#', the International Institute of Forecasting and Industrielle-Alliance.
#'@references Eddelbuettel, D. & Francois, R. (2011). \code{Rcpp}: Seamless \R and \code{C++} Integration. \emph{Journal of Statistical Software}, 40, pp. 1-18, \url{http://www.jstatsoft.org/v40/i08/}.
#'@references Eddelbuettel, D. & Sanderson, C. (2014). \code{RcppArmadillo}: Accelerating \R with High--Performance \code{C++} Linear Algebra. \emph{Computational Statistics & Data Analysis}, 71, pp. 1054-1063, \url{http://dx.doi.org/10.1016/j.csda.2013.02.005}.
#'@references Ghalanos, A. (2015). \code{rugarch}: Univariate GARCH Models. \url{https://cran.r-project.org/web/packages/rugarch/}.
"_PACKAGE"