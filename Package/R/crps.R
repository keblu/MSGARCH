#'CRPS (continuous ranked probability score) measure.
#' @description Method returning the CRPS at \code{t = T + 1}.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param yn Scalar value to be evaluated at \code{t = T + 1}.
#' @param ctr Control list parameters.
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian CRPS is calculated.
#' @examples
#' require("MSGARCH")
#' # load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = list(do.init = FALSE))
#'
#'# run at T + 1 from model       
#'crps = MSGARCH::crps(object = fit, yn = 0.6)
#'
#' @return A vector with five crps measures
#' @importFrom stats pnorm dnorm
#' @export
crps <- function(object, theta, y, yn, ctr = list(lower = -20, upper = 20, n.mesh = 500, a = 0, b = 1)) {
  UseMethod("crps", object)
}

#' @export
crps.MSGARCH_SPEC = function(object, theta, y, yn, ctr = list(lower = -20, upper = 20, n.mesh = 500, a = 0, b = 1)) {
  
  x   = seq(from = ctr$lower, to = ctr$upper, length.out = ctr$n.mesh)
  n.x = length(x)
  
  w1  = rep(1, n.x) #unif
  w2  = stats::dnorm(x, mean = ctr$a, sd = ctr$b) #center
  w3  = 1 - w2 / dnorm(0, mean = ctr$a, sd = ctr$b) #tails
  w4  = stats::pnorm(x, mean = ctr$a, sd = ctr$b) #tail_r
  w5  = 1 - w4 #tail_l
  W   = cbind(w1, w2, w3, w4, w5)
  
  cdf = MSGARCH::pit(object, x = x, y = y, theta = theta, do.norm = FALSE, do.its = FALSE)$pit
  id  = yn < x
  tmp = (cdf - id)^2
  X   = matrix(data = tmp, nrow = n.x, ncol = 5)
  
  out = apply(W * X, 2, sum) / (ctr$n.mesh - 1)
  names(out) = c("unif", "center", "tails", "tail_r", "tail_l")
  return(out)
}

#' @export
crps.MSGARCH_MLE_FIT <- function(object, theta = NULL, y = NULL, yn = NULL,  ctr = list(lower = -20, upper = 20, n.mesh = 500, a = 0, b = 1)) {
  return(MSGARCH::crps(object = object$spec, theta = object$theta, y = object$y, yn = yn, ctr = ctr))
}

#' @export
crps.MSGARCH_BAY_FIT <- function(object, theta = NULL, y = NULL, yn = NULL,  ctr = list(lower = -20, upper = 20, n.mesh = 500, a = 0, b = 1)) {
  return(MSGARCH::crps(object = object$spec, theta = object$theta, y = object$y, yn = yn, ctr = ctr))
}