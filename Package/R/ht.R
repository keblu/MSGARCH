#' Conditional variance in each regime.
#' @description Method returning the conditional variance of each regime.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'ht = MSGARCH::ht(spec = spec,theta = spec$theta0, y = sp500ret)
#' @return Condititional variance time serie (array of size T + 1 x M x K) for each regime.
#' @export
ht <- function(spec, theta, y )
{
  UseMethod("ht", spec)
}


#' @export
ht.MSGARCH_SPEC = function(spec, theta, y) {
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
  
  out = spec$rcpp.func$calc_ht(theta, y)
  return(out)
}