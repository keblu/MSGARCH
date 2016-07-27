#' State probabilities filtering function.
#' @description Method returning the filtered state probabilities.
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
#'Pstate  = MSGARCH::Pstate(spec = spec, theta = spec$theta0, y = sp500ret)
#'@return Filtered state probabilities (array of size T x M x K).
#' @export
Pstate <- function(spec, theta, y )
{
  UseMethod("Pstate", spec)
}

#' @export
Pstate.MSGARCH_SPEC = function(spec, theta, y) {
  
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
  
  out = array(dim = c(nrow(y), nrow(theta), spec$K))
  for(i in 1:nrow(theta)){
    tmp = spec$rcpp.func$get_Pstate_Rcpp(theta[i,], y, FALSE)
    for(j in 1:spec$K){
      out[,i,j] = tmp[,j]
    }
  }
  
  return(out)
}