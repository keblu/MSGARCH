#' Predictive function.
#' @description Method returning the predictive probability density in-sample or of a vector of points at \code{t = T + 1}.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param x Vector (of size N) of point at \code{t = T + 1} to be evaluated (used when \code{is.its = FALSE}).
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param log  Boolean indicating if the log-density is returned. (Default: \code{log = TRUE})
#' @param is.its  Boolean indicating if the in-sample predictive is returned. (Default: \code{is.its = FALSE})
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian Probability integral transform is calculated.
#' If \code{is.its = FALSE}, the points \code{x} are evaluated as \code{t = T + 1} realization and the method uses the variance estimate at \code{t = T + 1}.
#' If \code{is.its = TRUE}, \code{y} is evaluated using their respective variance estimate at each time \code{t}.
#' @usage pred(spec, theta, y, log = TRUE, is.its = TRUE)
#' pred(fit, log = TRUE, is.its = TRUE) 
#' @examples 
#'\dontrun{
#'# load data
#'data("sp500ret")
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)
#'                            
#'# run pred method in-sample     
#'pred.its = MSGARCH::pred(fit, log = TRUE, is.its = TRUE)  
#' 
#'plot(pred.its)  
#'                                              
#'# create mesh
#'x = seq(-3,3,0.01)
#'
#'# run pred method on mesh at T + 1
#'pred = MSGARCH::pred(fit, x = x, log = TRUE, is.its = FALSE)
#'
#'plot(pred)
#'}
#' @return A list of class \code{MSGARCH_PRED} containing two components:
#' \itemize{
#' \item \code{pred}:\cr If \code{is.its = FALSE}: (Log-)Predictive of of the points \code{x} at \code{t = T + 1} (vector of size N). \cr
#'                   If \code{is.its = TRUE}: In-sample Predictive of \code{y} (vector of size T or matrix of size M x T). 
#' \item \code{x}:\cr If \code{is.its = FALSE}: Vector (of size N) of point at \code{t = T + 1} evaluated.\cr
#'                 If \code{is.its = TRUE}: Vector (of size T) of observations.
#' }
#'The class \code{MSGARCH_PRED} contains the \code{plot} method.
#' @export
pred <- function(spec, x, theta, y, log = TRUE, is.its = FALSE)
{
  UseMethod("pred", spec)
}

#' @export
pred.MSGARCH_SPEC = function(spec, x = NULL, theta, y, log = TRUE, is.its = FALSE) {
  
  
  y = f.check.y(y)

  theta = f.check.theta(spec, theta)
  
  N = nrow(theta)
  
  if(isTRUE(is.its)){
    nx = length(y)
    x = y
  } else {
    nx = length(x)
  }
  
  tmp = matrix(data = NA, nrow = N, ncol = nx)
  for (i in 1:N) {
    tmp[i, ] = MSGARCH::pdf(spec, x, theta = theta[i, ], y = y, log = FALSE, is.its = is.its)$pdf
  }
  tmp = colMeans(tmp)
  if (log) {
    tmp = log(tmp)
  }
  out = list()
  out$pred = tmp
  out$x = x
  out$is.its = is.its
  class(out) = "MSGARCH_PRED"
  return(out)
}

#' @export
pred.MSGARCH_MLE_FIT = function(fit, x = NULL, log = TRUE, is.its = FALSE) {
  
  return(MSGARCH::pred(spec = fit$spec, x =  x, theta = fit$theta, y = fit$y, log = log, is.its = is.its))
  
}

#' @export
pred.MSGARCH_BAY_FIT = function(fit, x = NULL, log = TRUE, is.its = FALSE) {
  
  return(MSGARCH::pred(spec = fit$spec, x =  x, theta = fit$theta, y = fit$y, log = log, is.its = is.its))
  
}