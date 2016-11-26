#' Conditional variance in each regime.
#' @description Method returning the conditional variance of each regime.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not require when using a fit object).
#' @param y  Vector (of size T) of observations (not require when using a fit object).
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#'# load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#' fit = MSGARCH::fit.mle(spec = spec, y = sp500)
#'
#'# Compute the conditional variance
#'ht = MSGARCH::ht(object = fit)
#'
#'plot(ht)
#' @return Condititional variance (array of size (T + 1) x M x K) in each regime.
#' @export
ht <- function(object, theta, y) {
  UseMethod("ht", object)
}


#' @export
ht.MSGARCH_SPEC <- function(object, theta, y) {
  y     <- f.check.y(y)
  theta <- f.check.theta(object, theta)
  out   <- object$rcpp.func$calc_ht(theta, y)
  class(out) <- "MSGARCH_HT"
  return(out)
}

#' @export
ht.MSGARCH_MLE_FIT <- function(object, theta = NULL, y = NULL) {
  return(MSGARCH::ht(object = object$spec, theta = object$theta, y = object$y))
}

#' @export
ht.MSGARCH_BAY_FIT <- function(object, theta = NULL, y = NULL) {
  return(MSGARCH::ht(object = object$spec, theta = object$theta, y = object$y))
}