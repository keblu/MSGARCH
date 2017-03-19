#' Step ahead simulation method.
#' @description Method returning step ahead simulation up to time \code{n}.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param n  Number of steps ahead time steps. (Default: \code{n = 1})
#' @param m Number of simulations. (Default: m = 1)
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not required when using a fit object) where d must have
#'  the same length as the default parameters of the specification.
#' @param y  Vector (of size T) of observations (not required when using a fit object).

#' @return A list of class \code{MSGARCH_SIM} containing  two components:
#' \itemize{
#' \item \code{draws}:  Matrix (of size m x n) of step ahead simulated draws.
#' \item \code{state}:  Matrix (of size m x n) of step ahead simulated states.
#' }
#' The \code{MSGARCH_SIM} class contains the \code{plot} method.
#' @details If a matrix of parameters estimates is given, each parameter estimates is evaluated individually and \code{m = M}.
#' The \code{MSGARCH_SIM} class contains the \code{plot} method. The difference between \code{\link{sim}} and \code{\link{simahead}} is that
#' \code{\link{sim}} starts the simulation a t = 0 creating an entire new process while  \code{\link{simahead}} starts the simulation at t = T + 1
#'  taking in consideration all the information available in the original time series \code{y}.
#' @examples 
#' require("MSGARCH")
#' # load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#' spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation
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
  y <- MSGARCH:::f.check.y(y)
  theta <-  MSGARCH:::f.check.theta(object, theta)
  if (nrow(theta) == 1) {
    P_0 = matrix(MSGARCH::Pstate(object, theta = theta, y =y)[(length(y)+1),,], ncol = object$K)
    P_0 = matrix(rep(P_0,m),ncol =object$K, byrow = TRUE)
    theta <- matrix(theta[rep(1, m), ], ncol = ncol(theta))
  } else {
    P_0 =MSGARCH::Pstate(object, theta = theta, y =y)[(length(y)+1),,1:2]
    m = nrow(theta)
  }
  draws <- matrix(data = NA, nrow = nrow(theta), ncol = n)
  state <- matrix(data = NA, nrow = nrow(theta), ncol = n)
  for (i in 1:m) {
     tmp <- object$rcpp.func$simahead(y=y, n=n, theta = theta[i,],P_0[i,])
      if (object$K == 1) {
        draws[i,] <- tmp
        state[i,] <- rep(1,length(tmp))
      } else {
        draws[i,] <- tmp$draws
        state[i,] <- tmp$state
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
