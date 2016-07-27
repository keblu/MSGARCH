#'Probability Integral Transform at T + 1.
#' @description Method returning the predictive Probability integral transform (PIT).
#' @param spec Model specification of class \code{\link{MSGARCH_SPEC}} created with \code{\link{create.spec}}.
#' @param x Vector (of size N) of point to be evaluated
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param do.norm  Boolean indicating if the PIT value are transforms into standard Normal variate. (\code{do.norm = FALSE})
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian Probability integral transform is calculated.
#'          The \code{do.norm} argument transforms the PIT value into Normal variate so that normality test can be done.
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                            do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'                              
#'set.seed(123)
#'
#'x = rnorm(100)
#'
#'pit = MSGARCH::pit(spec = spec, x = x theta = spec$theta0, y = sp500ret, do.norm = FALSE)
#' @return Probability integral transform of the points \code{x} or Normal variate derived from the Probability integral transform of \code{x} (vector of size N).
#' @export
pit <- function(spec, x, theta, y,  do.norm = FALSE)
{
  UseMethod("pit", spec)
}

#' @export
pit.MSGARCH_SPEC = function(spec, x, theta, y, do.norm = FALSE) {
  
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
  
  N = nrow(theta)
  nx = length(x)
  
  out = matrix(data = NA, nrow = N, ncol = nx)
  for (i in 1:N) {
    out[i, ] = MSGARCH::cdf(spec = spec, x, theta = theta[i, ], y = y, log = FALSE)
  }
  out = colMeans(out)
  if (do.norm) {
    out = qnorm(out, mean = 0, sd = 1)
  }
  if (any(is.nan(out))) {
    stop("NaN value in PIT calculation")
  }
  return(out)
}