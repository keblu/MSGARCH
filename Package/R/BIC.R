#' Compute Bayesian information criterion (BIC).
#' @param fit Fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}
#' @references Schwarz, G. (1978). Estimating the dimension of a model. \emph{Annals of Statistics}, 6, pp. 461-464. 
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
#'# compute BIC
#'BIC = MSGARCH::BIC(fit)
#' @details If a matrix of MCMC posterior draws estimates is given, the BIC on the posterior mean is calculated.
#' @return BIC value.
#' @export
BIC <- function(fit) {
  UseMethod(generic = "BIC", object = fit)
}

#' @export
BIC.MSGARCH_MLE_FIT <- function(fit) {
  bic <- f.BIC(spec = fit$spec, theta = fit$theta, y = fit$y)
  return(bic)
}

#' @export
BIC.MSGARCH_BAY_FIT <- function(fit) {
  bic <- f.BIC(spec = fit$spec, theta = fit$theta, y = fit$y)
  return(bic)
}

f.BIC <- function(spec, theta, y) {
  if (is.vector(x = theta)) {
    theta <- matrix(data = theta, nrow = 1)
  }
  theta <- matrix(colMeans(x = theta), nrow = 1)
  LL <- sum(MSGARCH::pdf(object = spec, theta = theta, y = y, log = TRUE,
                          do.its = TRUE)$pdf, na.rm = TRUE)
  k <- dim(x = theta)[2]
  n <- length(x = y)
  bic <- k * log(x = n) - 2 * LL
  return(bic)
}