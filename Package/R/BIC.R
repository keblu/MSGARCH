

#' Compute Bayesian information criterion (BIC).
#' @param fit Fit object of type \code{\link{MSGARCH_MLE_FIT}} created with \code{\link{fit.mle}} or \code{\link{MSGARCH_BAY_FIT}} created with \code{\link{fit.bayes}}
#' @references Schwarz, G. (1978). Estimating the dimension of a model. \emph{Annals of Statistics}, 6, pp. 461-464. 
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'                              
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)
#'
#'BIC = MSGARCH::BIC(fit = fit)
#' @details If a matrix of MCMC posterior draws estimates is given, the BIC on the posterior mean is calculated.
#' @return BIC value.
#' @export
BIC <- function(fit)
{
  UseMethod("BIC", fit)
}

#' @export
BIC.MSGARCH_MLE_FIT <- function(fit){
  
  bic = f.BIC(fit$spec, fit$theta, fit$y)
  return(bic)
}

#' @export
BIC.MSGARCH_BAY_FIT <- function(fit){
  
  bic = f.BIC(fit$spec, fit$theta, fit$y)
  return(bic)
}

f.BIC = function(spec, theta, y) {
  
  if (is.vector(theta)) {
    theta = matrix(theta, nr = 1)
  }
  theta = matrix(colMeans(theta), nrow = 1)
  LL = MSGARCH::kernel(spec, theta, y = y, log = TRUE)
  k = dim(theta)[2]
  n = length(y)
  bic = 2 * LL - k * log(n)
  return(bic)
  
}