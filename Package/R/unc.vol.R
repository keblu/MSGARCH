#' Unconditional volatility of each regime.
#' @description Method returning the unconditional volatility of the process in each state.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not require when using a fit object).
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#' # create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# compute the unconditional volatility in each regime
#'unc.vol = MSGARCH::unc.vol(object = spec, theta = spec$theta0)
#' @return Unconditional volatility (vector of size K or matrix of size M x K) of each regime.
#' @export
unc.vol <- function(object, theta) {
  UseMethod(generic = "unc.vol", object =  object)
}

#' @export
unc.vol.MSGARCH_SPEC <- function(object, theta = NULL) {
  theta <- f.check.theta(object, theta)
  for (i in 1:nrow(theta)) {
    out <- object$rcpp.func$unc_vol_Rcpp(theta, 0)
  }
  
  if(any(is.null(out))){
    out = matrix(data = NA, ncol = object$K, nrow = nrow(theta))
    return(out)
  }
  
  if(any(is.na(out))){
    out = matrix(data = NA, ncol = object$K,nrow = nrow(theta))
    return(out)
  } else {
    if(any(out < 0)) {
      out = matrix(data = NA, ncol = object$K, nrow = nrow(theta))
      return(out)
    }
  }
  
  out <- sqrt(out)
  out = matrix(out, ncol = object$K,nrow = nrow(theta))
  colnames(out) = paste0("State ", 1:object$K)
  return(out)
}

#' @export
unc.vol.MSGARCH_MLE_FIT <- function(object, theta = NULL) {
  return(MSGARCH::unc.vol(object = object$spec, theta = object$theta))
}

#' @export
unc.vol.MSGARCH_BAY_FIT <- function(object, theta = NULL) {
  return(MSGARCH::unc.vol(object = object$spec, theta = object$theta))
}