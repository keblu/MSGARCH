#' Compute Akaike information criterion (AIC).
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @return AIC value.
#' @examples 
#' # load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#' 
#'# fit the model by MLE                                                         
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = list(do.init = FALSE))
#'
#'# compute AIC
#'AIC = MSGARCH::AIC(fit)
#' @details If a matrix of MCMC posterior draws estimates is given, the AIC on the posterior mean is calculated.
#' @references Akaike, H. (1974). A New Look at the Statistical Model Identification. \emph{IEEE Transactions on Automatic Control}, 19, pp. 716-723.
#' @export
AIC <- function(fit) {
  UseMethod(generic = "AIC", object = fit)
}

#' @export
AIC.MSGARCH_MLE_FIT <- function(fit) {
  aic <- f.AIC(spec = fit$spec, theta = fit$theta, y = fit$y)
  return(aic)
}

#' @export
AIC.MSGARCH_BAY_FIT <- function(fit) {
  aic <- f.AIC(spec = fit$spec, theta = fit$theta, y = fit$y)
  return(aic)
}

f.AIC <- function(spec, theta, y) {
  if (is.vector(x = theta)) {
    theta <- matrix(data = theta, nrow = 1)
  }
  theta <- matrix(data = colMeans(x = theta), nrow = 1)
  LL    <- sum(MSGARCH::pdf(object = spec, theta = theta, y = y, log = TRUE,
                            do.its = TRUE)$pdf, na.rm = TRUE)
  k <- dim(x = theta)[2]
  aic <- 2 * k - 2 * LL
  return(aic)
}