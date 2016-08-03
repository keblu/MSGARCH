#' Unconditional volatility of each regime.
#' @description Method returning the unconditional volatility of the process in each state.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#' # create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# compute the unconditional volatility in each regime
#'unc.vol = MSGARCH::unc.vol(spec = spec, theta = spec$theta0)
#' @return Unconditional volatility (vector of size K or matrix of size M x K) of each regime. 
#' @usage unc.vol(spec, theta)
#' unc.vol(fit)
#' @export
unc.vol <- function(spec, theta)
{
  UseMethod("unc.vol", spec)
}

#' @export
unc.vol.MSGARCH_SPEC = function(spec, theta) {
  
  if (isTRUE(spec$is.shape.ind)) {
    theta = spec$func$f.do.shape.ind(theta)
  }
  
  if (isTRUE(spec$is.mix)) {
    theta = spec$func$f.do.mix(theta)
  }
  
  theta = f.check.theta(spec, theta)
  
  for(i in 1:nrow(theta)){
    out =  spec$rcpp.func$unc_vol_Rcpp(theta, 0)
  }
  
  out = sqrt(out)
  return(out)
}

#' @export
unc.vol.MSGARCH_MLE_FIT = function(fit) {
  
  return(MSGARCH::unc.vol(spec = fit$spec, theta = fit$theta))
  
}

#' @export
unc.vol.MSGARCH_BAY_FIT = function(fit) {
  
  return(MSGARCH::unc.vol(spec = fit$spec, theta = fit$theta))
  
}