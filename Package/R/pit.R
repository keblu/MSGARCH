#'Probability Integral Transform.
#' @description Method returning the predictive probability integral transform (PIT) in-sample or of a vector of points at \code{t = T + 1}.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param x Vector (of size N) of point at \code{t = T + 1} to be evaluated (used when \code{do.its = FALSE}).
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not require when using a fit object).
#' @param y  Vector (of size T) of observations (not require when using a fit object).
#' @param do.norm  Boolean indicating if the PIT value are transforms into standard Normal variate. (Default: \code{do.norm = FALSE}).
#' @param do.its  Boolean indicating if the in-sample pit is returned. (Default: \code{do.its = FALSE})
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian probability integral transform is calculated.
#' If \code{do.its = FALSE}, the points \code{x} are evaluated as \code{t = T + 1} realization and the method uses the variance estimate at \code{t = T + 1}.
#' If \code{do.its = TRUE}, \code{y} is evaluated using their respective variance estimate at each time \code{t}.
#' The \code{do.norm} argument transforms the PIT value into Normal variate so that normality test can be done.
#' @examples
#' # load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = list(do.init = FALSE))
#'
#'# run pit method in-sample              
#'pit.its = MSGARCH::pit(object = fit, do.norm = FALSE, do.its = TRUE)                              
#' 
#'plot(pit.its)  
#'                                                                          
#'# generate random draws at T + 1 from model
#'set.seed(123)
#'sim.ahead = MSGARCH::simahead(object = fit, n = 1, m = 100)
#'
#'x = sim.ahead$draws
#'
#'# run pit method on random draws at T + 1 from model
#'pit = MSGARCH::pit(object = fit, x = x, do.norm = FALSE)
#'
#'plot(pit)
#' @return A list of class \code{MSGARCH_PIT} containing two components:
#' \itemize{
#' \item \code{pit}:\cr If \code{do.its = FALSE}: probability integral transform of the points \code{x} at \code{t = T + 1} or Normal variate derived from the probability integral transform of \code{x} (vector of size N).\cr
#'                   If \code{do.its = TRUE}: In-sample  probability integral transform or Normal variate derived from the probability integral transform of \code{y} (vector of size T or matrix of size M x T). 
#' \item \code{x}:\cr If \code{do.its = FALSE}: Vector (of size N) of at point \code{t = T + 1} evaluated.\cr
#'                 If \code{do.its = TRUE}: Vector (of size T) of observations.  
#' }
#'The class \code{MSGARCH_PIT} contains the \code{plot} method only if \code{do.its = FALSE}.
#' @importFrom stats qnorm
#' @export
pit <- function(object, x, theta, y, do.norm = FALSE, do.its = FALSE) {
  UseMethod("pit", object)
}

#' @export
pit.MSGARCH_SPEC <- function(object, x = NULL, theta = NULL, y = NULL, do.norm = FALSE, do.its = FALSE) {
  y <- f.check.y(y)
  N <- nrow(theta)
  if (isTRUE(do.its)) {
    nx <- length(y)
    x <- y
  } else {
    nx <- length(x)
  }
  tmp <- matrix(data = NA, nrow = N, ncol = nx)
  for (i in 1:N) {
    tmp <- MSGARCH::cdf(object = object, x, theta = theta[i, ], y = y, log = FALSE, do.its = do.its)$cdf
  }
  tmp <- colMeans(tmp, na.rm = TRUE)
  if (do.norm) {
    tmp <- qnorm(tmp, mean = 0, sd = 1)
  }
  
  out <- list()
  out$pit <- tmp
  out$x <- x
  out$do.its <- do.its
  class(out) <- "MSGARCH_PIT"
  return(out)
}

#' @export
pit.MSGARCH_MLE_FIT <- function(object, x = NULL, theta = NULL, y = NULL, do.norm = TRUE,
                               do.its = FALSE) {
  return(MSGARCH::pit(object = object$spec, x = x, theta = object$theta, y = object$y,
                      do.norm = do.norm, do.its = do.its))
}

#' @export
pit.MSGARCH_BAY_FIT <- function(object, x = NULL, theta = NULL, y = NULL, do.norm = TRUE,
                               do.its = FALSE) {
  return(MSGARCH::pit(object = object$spec, x = x, theta = object$theta, y = object$y, 
                      do.norm = do.norm, do.its = do.its))
}