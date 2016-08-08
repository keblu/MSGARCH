#' Compute Akaike information criterion (AIC).
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @return AIC value.
#' @examples 
#'\dontrun{
#' # load data
#'data("sp500ret")
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#' 
#'# fit the model by MLE                                                         
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)
#'
#'# compute AIC
#'AIC = MSGARCH::AIC(fit)
#'}
#' @details If a matrix of MCMC posterior draws estimates is given, the AIC on the posterior mean is calculated.
#' @references Akaike, H. (1974). A New Look at the Statistical Model Identification. \emph{IEEE Transactions on Automatic Control}, 19, pp. 716-723.
#' @export
AIC <- function(fit)
{
  UseMethod("AIC", fit)
}

#' @export
AIC.MSGARCH_MLE_FIT <- function(fit){
  
  aic = f.AIC(fit$spec, fit$theta, fit$y)
  return(aic)
}

#' @export
AIC.MSGARCH_BAY_FIT <- function(fit){
  
  aic = f.AIC(fit$spec, fit$theta, fit$y)
  return(aic)
}

f.AIC = function(spec, theta, y) {
  
  if (is.vector(theta)) {
    theta = matrix(theta, nrow = 1)
  }
  theta = matrix(colMeans(theta), nrow = 1)
 
  LL =  sum(pred(object = spec, theta = theta, y = y, log = TRUE, is.its = TRUE)$pred, na.rm = TRUE)
  k = dim(theta)[2]
  aic = 2 * k - 2 * LL 
  return(aic)
}