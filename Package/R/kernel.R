#' Kernel function.
#' @description Method returning the kernel value of a vector of observations given a model specification.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not require when using a fit object).
#' @param y  Vector (of size T) of observations (not require when using a fit object).
#' @param log  Boolean indicating if the log kernel is returned. (Default: \code{log = TRUE})
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
#'# load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#' fit = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = list(do.init = FALSE))
#'
#'# compute the kernel
#'kernel = MSGARCH::kernel(fit, log = TRUE)
#' @references Hamilton, J. D. (1989) A New Approach to the Economic Analysis of Nonstationary Time Series and the Business Cycle. \emph{Econometrica}, 57, pp.357-38
#' @return (Log-)Kernel value (scalar or vector of size M) of the vector of observations.
#' @export
kernel <- function(object, theta, y, log = TRUE) {
  UseMethod("kernel", object)
}

#' @export
kernel.MSGARCH_SPEC <- function(object, theta, y = NULL, log = TRUE) {
  y     <- f.check.y(y)
  theta <- f.check.theta(object, theta)
  lnd   <- object$rcpp.func$eval_model(theta, y)
  lnd[is.na(lnd) | is.nan(lnd) | is.infinite(lnd)] <- -1e+10
  if (!log)
    lnd <- exp(lnd)
  return(lnd)
}

#' @export
kernel.MSGARCH_MLE_FIT <- function(object, theta = NULL, y = NULL, log = TRUE) {
  return(MSGARCH::kernel(object = object$spec, theta = object$theta,
                         y = object$y, log = log))
}

#' @export
kernel.MSGARCH_BAY_FIT <- function(object, theta = NULL, y = NULL, log = TRUE) {
  return(MSGARCH::kernel(object = object$spec, theta = object$theta,
                         y = object$y, log = log))
}