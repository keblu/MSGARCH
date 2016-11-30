#' Step ahead simulation method.
#' @description Method returning step ahead simulation up to time \code{n}.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param n  Mumber of step ahead time step. (Default: \code{n = 1})
#' @param m Number of simulations. (Default: m = 1)
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not require when using a fit object).
#' @param y  Vector (of size T) of observations (not require when using a fit object).

#' @return A list of class \code{MSGARCH_SIM} containing  two components:
#' \itemize{
#' \item \code{draws}:  Matrix (of size m x n) of step ahead simulated draws.
#' \item \code{state}:  Matrix (of size m x n) of step ahead simulated states.
#' }
#' The \code{MSGARCH_SIM} class contains the \code{plot} method.
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually and \code{m = M}.
#' The \code{MSGARCH_SIM} class contains the \code{plot} method. The difference between \code{\link{sim}} and \code{\link{simahead}} is that
#' \code{\link{sim}} starts the simulation a t = 0 creating an entire new process while  \code{\link{simahead}} starts the simulation at t = T + 1
#'  taking in consideration all the information available in the original time serie \code{y}.
#' @examples 
#' # load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#' spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500, ctr = list(do.init = FALSE))
#'  
#'# generate random draws
#'set.seed(123)
#' simahead = MSGARCH::simahead(object = fit, n = 30, m = 100)
#' 
#' plot(simahead)
#' @export
simahead <- function(object, n, m, theta, y) {
  UseMethod("simahead", object)
}

#' @export
simahead.MSGARCH_SPEC <- function(object, n, m, theta = NULL, y = NULL) {
  y <- f.check.y(y)
  theta <- f.check.theta(object, theta)
  if (nrow(theta) == 1) {
    theta <- matrix(theta[rep(1, m), ], ncol = ncol(theta))
  }
  draws <- matrix(data = NA, nrow = nrow(theta), ncol = n)
  state <- matrix(data = NA, nrow = nrow(theta), ncol = n)
  for (j in 1:n) {
    for (i in 1:nrow(theta)) {
      tmp <- object$rcpp.func$rnd_Rcpp(1, theta[i, ], c(y, draws[i, 0:(j - 1)]))
      if (object$K == 1) {
        draws[i, j] <- tmp
        state[i, j] <- 1
      } else {
        draws[i, j] <- tmp$draws
        state[i, j] <- tmp$state
      }
    }
  }
  out <- list()
  out$draws  <- draws
  out$state  <- state
  class(out) <- "MSGARCH_SIM"
  return(out)
}

#' @export
simahead.MSGARCH_MLE_FIT <- function(object, n = 1, m = 1, theta = NULL, y = NULL) {
  return(MSGARCH::simahead(object = object$spec, n = n, m = m, theta = object$theta,
                            y = object$y))
}

#' @export
simahead.MSGARCH_BAY_FIT <- function(object, n = 1, m = 1, theta = NULL, y = NULL) {
  return(MSGARCH::simahead(object = object$spec, n = n, m = m, theta = object$theta,
                            y = object$y))
}
