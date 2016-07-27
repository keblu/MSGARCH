#' Predictive density function.
#' @description Method returning the predictive probability density of a vector of points.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param x  Vector (of size N) of point to be evaluated
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param log  Boolean indicating if the log-density is returned. (default: \code{log = TRUE})
#' @details  If a matrix of MCMC posterior draws estimates is given, the Bayesian predictive density is calculated.
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
#'pred = MSGARCH::pred(spec = spec, x = x, theta = spec$theta0, y = sp500ret, log = TRUE)
#' @return Predictive density or log-density of \code{x} (vector of size N).
#' @export
pred <- function(spec, x, theta, y, log = TRUE)
{
  UseMethod("pred", spec)
}

#' @export
pred.MSGARCH_SPEC = function(spec, x, theta, y, log = TRUE) {
  
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
    out[i, ] = MSGARCH::pdf(spec, x, theta = theta[i, ], y = y, log = FALSE)
  }
  out = colMeans(out)
  if (log) {
    out = log(out)
  }
  return(out)
}
