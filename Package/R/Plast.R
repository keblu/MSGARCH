#' State probabilities at T + 1.
#' @description Method returning the state probabilities at  T + 1.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}. 
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#' @return State probabilities at T + 1 (matrix of size M x K).
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'Plast = MSGARCH::Plast(spec = spec, theta = spec$theta0, y = sp500ret) 
#' @export
Plast <- function(spec, theta, y)
{
  UseMethod("Plast", spec)
}

#' @export
Plast.MSGARCH_SPEC = function(spec, theta, y) {
  
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
  
  out = matrix(data = NA,nrow = nrow(theta), ncol = spec$K)
  for(i in 1:nrow(theta)){
    out[i,] =  spec$rcpp.func$get_Pstate_Rcpp(theta[i,], y, TRUE)
  }
  
  return(out)
}