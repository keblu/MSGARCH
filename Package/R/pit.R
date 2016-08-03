#'Probability Integral Transform.
#' @description Method returning the predictive Probability integral transform (PIT) in-sample or of a vector of points at \code{t = T + 1}.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param x Vector (of size N) of point at \code{t = T + 1} to be evaluated (used when \code{is.its = FALSE}).
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param do.norm  Boolean indicating if the PIT value are transforms into standard Normal variate. (Default: \code{do.norm = FALSE}).
#' @param is.its  Boolean indicating if the in-sample pit is returned. (Default: \code{is.its = FALSE})
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian Probability integral transform is calculated.
#' If \code{is.its = FALSE}, the points \code{x} are evaluated as \code{t = T + 1} realization and the method uses the variance estimate at \code{t = T + 1}.
#' If \code{is.its = TRUE}, \code{y} is evaluated using their respective variance estimate at each time \code{t}.
#' The \code{do.norm} argument transforms the PIT value into Normal variate so that normality test can be done.
#' @usage pit(spec, theta, y, do.norm = FALSE, is.its = TRUE)
#' pit(fit, do.norm = FALSE, is.its = TRUE) 
#' @examples
#'\dontrun{
#' # load data
#'data("sp500ret")
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)
#'
#'# run pit method in-sample              
#'pit.its = MSGARCH::pit(fit, do.norm = FALSE, is.its = TRUE)                              
#' 
#'plot(pit.its)  
#'                                                                          
#'# generate random draws at T + 1 from model
#'set.seed(123)
#'rnd = MSGARCH::rnd(fit, n = 100000)
#'
#'x = rnd$draws
#'
#'# run pit method on random draws at T + 1 from model
#'pit = MSGARCH::pit(fit, x = x, do.norm = FALSE)
#'
#'plot(pit)
#'}
#' @return A list of class \code{MSGARCH_PIT} containing two components:
#' \itemize{
#' \item \code{pit}:\cr If \code{is.its = FALSE}: Probability integral transform of the points \code{x} at \code{t = T + 1} or Normal variate derived from the Probability integral transform of \code{x} (vector of size N).\cr
#'                   If \code{is.its = TRUE}: In-sample  Probability integral transform or Normal variate derived from the Probability integral transform of \code{y} (vector of size T or matrix of size M x T). 
#' \item \code{x}:\cr If \code{is.its = FALSE}: Vector (of size N) of at point \code{t = T + 1} evaluated.\cr
#'                 If \code{is.its = TRUE}: Vector (of size T) of observations.  
#' }
#'The class \code{MSGARCH_PIT} contains the \code{plot} method.
#' @usage pit(spec, theta, y, do.norm = TRUE, is.its = TRUE)
#' pit(fit, do.norm = TRUE, is.its = TRUE) 
#' @importFrom stats qnorm
#' @export
pit <- function(spec, x, theta, y,  do.norm = FALSE, is.its = FALSE)
{
  UseMethod("pit", spec)
}

#' @export
pit.MSGARCH_SPEC = function(spec, x = NULL, theta, y, do.norm = FALSE, is.its = FALSE) {
  
  y = f.check.y(y)
  
  theta = f.check.theta(spec, theta)
  
  N = nrow(theta)
  
  if(isTRUE(is.its)){
    nx = length(y)
    x = y
  } else {
    nx = length(x)
  }
  tmp = matrix(data = NA, nrow = N, ncol = nx)
  for (i in 1:N) {
    tmp = MSGARCH::cdf(spec = spec, x, theta = theta[i, ], y = y, log = FALSE, is.its = is.its)$cdf
  }
  tmp = colMeans(tmp)
  if (do.norm) {
    tmp = qnorm(tmp, mean = 0, sd = 1)
  }
  if (any(is.nan(tmp))) {
    stop("NaN value in PIT calculation")
  }
  out = list()
  out$pit = tmp
  out$x = x
  out$is.its = is.its
  class(out) = "MSGARCH_PIT"
  return(out)
}

#' @export
pit.MSGARCH_MLE_FIT = function(fit, x = NULL, do.norm = TRUE, is.its = FALSE) {
  
  return(MSGARCH::pit(spec = fit$spec, x =  x, theta = fit$theta, y = fit$y, do.norm = do.norm, is.its = is.its))
  
}

#' @export
pit.MSGARCH_BAY_FIT = function(fit, x = NULL, do.norm = TRUE, is.its = FALSE) {
  
  return(MSGARCH::pit(spec = fit$spec, x =  x, theta = fit$theta, y = fit$y, do.norm = do.norm, is.its = is.its))
  
}