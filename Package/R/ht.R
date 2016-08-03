#' Conditional volatility in each regime.
#' @description Method returning the conditional volatility in each regime.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#'\dontrun{
#'# load data
#'data("sp500ret")
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#' fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)
#'
#'# Compute the conditional volatility
#'ht = MSGARCH::ht(fit)
#'
#'plot(ht)
#'}
#' @return Condititional volatility time serie (array of size (T + 1) x M x K) in each regime.
#' @usage ht(spec, theta, y)
#' ht(fit)
#' @export
ht <- function(spec, theta, y)
{
  UseMethod("ht", spec)
}


#' @export
ht.MSGARCH_SPEC = function(spec, theta, y) {
  y = f.check.y(y)
  
  theta = f.check.theta(spec, theta)
  
  out = spec$rcpp.func$calc_ht(theta, y)
  out = sqrt(out)
  return(out)
}

#' @export
ht.MSGARCH_MLE_FIT = function(fit) {
  
  return(MSGARCH::ht(spec = fit$spec, theta = fit$theta, y = fit$y))
  
}

#' @export
ht.MSGARCH_BAY_FIT = function(fit) {
  
  return(MSGARCH::ht(spec = fit$spec, theta = fit$theta, y = fit$y))
  
}