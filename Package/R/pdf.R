#' Probability density function at T + 1.
#' @description Method returning the probability density at T + 1 of a vector of points.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param x Vector (of size N) of point to be evaluated
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param log  Boolean indicating if the log-density is returned. (Default: \code{log = TRUE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#'  The \code{\link{pdf}} method uses the last variance estimate by filtering.
#' @examples 
#'# load data
#'data("sp500ret")
#'
#'# create model specification
#'spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'                              
#'set.seed(123)
#'
#'# generate random draws
#'x = rnorm(100)
#'
#'# run pdf method on random draws
#'pdf = MSGARCH::pdf(spec = spec, x = x, theta = spec$theta, y = sp500ret, log = FALSE)
#' @return (Log-)Probability density at T + 1 of the points \code{x} (vector of size N or matrix of size M x N).
#' @export
pdf <- function(spec, x, theta, y, log = TRUE)
{
  UseMethod("pdf", spec)
}

#' @export
pdf.MSGARCH_SPEC = function(spec, x, theta, y, log = TRUE) {
  
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
  
  out = matrix(data = NA,nrow = nrow(theta),ncol = length(x))
  for(i in 1:nrow(theta)){
    out[i,] = spec$rcpp.func$pdf_Rcpp(x, theta[i,], y, log)
  }
  return(out)
}