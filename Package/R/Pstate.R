#' Filtered state probabilities.
#' @description Method returning the filtered probabilities of the states.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not required when using a fit object) where d must have
#'  the same length as the default parameters of the specification.
#' @param y  Vector (of size T) of observations (not required when using a fit object).
#' @details If a matrix of parameter estimates is given, each parameter estimate (each row) is evaluated individually.
#' @examples 
#' require("MSGARCH")
#'# load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = list(do.init = FALSE))
#'  
#'# compute the filtered state probabilities
#'Pstate  = MSGARCH::Pstate(object = fit)
#'
#'plot(Pstate)
#'@return Filtered state probabilities of class \code{MSGARCH_PSTATE} (array of size (T + 1) x M x K).
#'The class \code{MSGARCH_PSTATE} contains the \code{plot} method.
#' @export
Pstate <- function(object, theta, y) {
  UseMethod("Pstate", object)
}

#' @export
Pstate.MSGARCH_SPEC <- function(object, theta, y) {
  y     <- as.matrix(y)
  theta <- f.check.theta(object, theta)
  out   <- array(dim = c(nrow(y) + 1, nrow(theta), object$K))
  for (i in 1:nrow(theta)) {
    tmp <- object$rcpp.func$get_Pstate_Rcpp(theta[i, ], y, FALSE)
    for (j in 1:object$K) {
      out[, i, j] <- tmp[, j]
    }
  }
  class(out) <- "MSGARCH_PSTATE"
  return(out)
}

#' @export
Pstate.MSGARCH_MLE_FIT <- function(object, theta = NULL, y = NULL) {
  return(MSGARCH::Pstate(object = object$spec, theta = object$theta, y = object$y))
}

#' @export
Pstate.MSGARCH_BAY_FIT <- function(object, theta = NULL, y = NULL) {
  return(MSGARCH::Pstate(object = object$spec, theta = object$theta, y = object$y))
}