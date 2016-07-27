#' Unconditional variance of each regime.
#' @description Method returning the unconditional variance of the process in each state.
#' @param spec Model specification of class \code{\link{MSGARCH_SPEC}} created with \code{\link{create.spec}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#'spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'unc.vol = MSGARCH::unc.vol(spec = spec, theta = theta)
#' @return Unconditional variance (vector of size K or matrix of size M x K) of each regime. 
#' @export
unc.vol <- function(spec, theta, y = 0)
{
  UseMethod("unc.vol", spec)
}

#' @export
unc.vol.MSGARCH_SPEC = function(spec, theta, y = 0) {
  
  y = as.matrix(y)
  if (isTRUE(spec$is.shape.ind)) {
    theta = spec$func$f.do.shape.ind(theta)
  }
  
  if (isTRUE(spec$is.mix)) {
    theta = spec$func$f.do.mix(theta)
  }
  
  if (is.vector(theta)) {
    theta = matrix(theta, nrow = 1)
  }
  
  for(i in 1:nrow(theta)){
    out =  spec$rcpp.func$unc_vol_Rcpp(theta, y)
  }
  
  out = sqrt(out)
  return(out)
}