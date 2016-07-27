#' Kernel function.
#' @description Method returning the kernel value of a vector of observations.
#' @param spec Model specification of class \code{\link{MSGARCH_SPEC}} created with \code{\link{create.spec}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param log  Boolean indicating if the log kernel is returned. (default: \code{log = TRUE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually.
#'  The kernel is a combination of the prior and the likelihood function. 
#'  The kernel is equal to prior(\eqn{\theta}) + L(y|\eqn{\theta}) where L is the likelihood
#'  of y given the parameter \eqn{\theta}. When doing optimization, the goal is to minimize the negative log-kernel.
#'  \itemize{
#'  \item Details on the prior \cr
#'        The prior is different for each specification. It ensures that the \eqn{\theta} makes the conditional variance process stationary, positive,
#'        and that it respect that the sums of the probabilities in the case of a multiple-regime models are all equal to 1. If any of these three conditions is not respected the prior return \code{-1e10}, meaning that the optimizer or sampler
#'        will know that \eqn{\theta} is not a good candidate.
#'   }
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'kernel = MSGARCH::kernel(spec = spec, theta = theta, y = sp500ret, log = TRUE)
#' @references Hamilton, J. D. (1989) A New Approach to the Economic Analysis of Nonstationary Time Series and the Business Cycle. \emph{Econometrica}, 57, pp.357-38
#' @return Kernel or log-kernel value (scalar or vector of size M) of the vector of observations. 
#' @export
kernel <- function(spec, theta, y, log = TRUE)
{
  UseMethod("kernel", spec)
}

#' @export
kernel.MSGARCH_SPEC = function(spec, theta, y, log = TRUE) {
  
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
  
  lnd = spec$rcpp.func$eval_model(theta, y)
  lnd[is.na(lnd) | is.nan(lnd) | is.infinite(lnd)] = -1e+10
  if (!log) 
    lnd = exp(lnd)
  return(lnd)
}