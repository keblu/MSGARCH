#' Cumulative  function.
#' @description Method returning the cumulative function in-sample or of a vector of points at \code{t = T + 1}.
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param x Vector (of size N) of point at \code{t = T + 1} to be evaluated (used when \code{do.its = FALSE}).
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not require when using a fit object).
#' @param y  Vector (of size T) of observations (not require when using a fit object).
#' @param log  Boolean indicating if the log cumulative is returned. (Default: \code{log = FALSE})
#' @param do.its  Boolean indicating if the in-sample cdf is returned. (Default: \code{do.its = FALSE})
#' @details If a matrix of parameter estimates is given, each parameter estimates is evaluated individually. 
#' If \code{do.its = FALSE}, the points \code{x} are evaluated as \code{t = T + 1} realization and the method uses the variance estimate at \code{t = T + 1}.
#' If \code{do.its = TRUE}, \code{y} is evaluated using their respective variance estimate at each time \code{t}.
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
#'# run pdf method in-sample
#'cdf.its = MSGARCH::cdf(object = fit, log = FALSE, do.its = TRUE)
#'
#'# create mesh
#'x = seq(-3,3,0.01)
#'
#'# run cdf method on mesh at T + 1
#'cdf = MSGARCH::cdf(object = fit, x = x, log = FALSE, do.its = FALSE)
#'
#'plot(cdf)
#' @return A list of class \code{MSGARCH_CDF} containing two components:
#' \itemize{
#' \item \code{cdf}:\cr  If \code{do.its = FALSE}: (Log-)Cumulative of the points \code{x} at \code{t = T + 1} (vector of size N or matrix of size M x N).\cr
#'                   If \code{do.its = TRUE}: In-sample (Log-)Cumulative of \code{y} (vector of size T or matrix of size M x T). 
#' \item \code{x}: \cr If \code{do.its = FALSE}: Vector (of size N) of point at \code{t = T + 1} evaluated.\cr
#'                 If \code{do.its = TRUE}: Vector (of size T) of observations.
#' }
#'The class \code{MSGARCH_CDF} contains the \code{plot} method only if \code{do.its = FALSE}.
#' @export
cdf <- function(object, x, theta, y, log = FALSE, do.its = FALSE) {
  UseMethod("cdf", object)
}

#' @export
cdf.MSGARCH_SPEC <- function(object, x = NULL, theta, y, log = FALSE, do.its = FALSE) {
  y <- f.check.y(y)
  if (is.vector(theta)) {
    theta <- matrix(theta, nrow = 1)
  }
 
  theta_check <- f.check.theta(object, theta)
  if (isTRUE(do.its)) {
    x <- y
    tmp <- matrix(data = NA, nrow = nrow(theta), ncol = length(y) - 1)
    for (i in 1:nrow(theta)) {
      tmp2 <- matrix(data = NA, nrow = length(y) - 1, ncol = object$K)
      if (object$K == 1) {
        tmp2 <- object$rcpp.func$cdf_Rcpp_its(theta_check[i, ], y, FALSE)
        tmp[i, ] <- tmp2
      } else {
        Pstate <- MSGARCH::Pstate(object = object, theta = theta[i, ], y = y)
        Pstate.tmp <- matrix(data = NA, nrow = dim(Pstate)[1],
                             ncol = dim(Pstate)[3])
        for (j in 1:dim(Pstate)[3]) {
          Pstate.tmp[, j] <- Pstate[, , j]
        }
        tmp2 <- object$rcpp.func$cdf_Rcpp_its(theta_check[i, ], y, FALSE)
        tmp[i, ] <- rowSums(tmp2 * Pstate.tmp[2:(nrow(Pstate.tmp) - 1), ])
      }
    }
    tmp <- cbind(NA, tmp)
  } else {
    tmp <- matrix(data = NA, nrow = nrow(theta), ncol = length(x))
    for (i in 1:nrow(theta)) {
      tmp[i, ] <- object$rcpp.func$cdf_Rcpp(x, theta_check[i, ], y, FALSE)
    }
  }
  out <- list()
  if (log) {
    tmp <- log(tmp)
  }
  out$cdf    <- tmp
  out$x      <- x
  out$do.its <- do.its
  class(out) <- "MSGARCH_CDF"
  return(out)
}

#' @export
cdf.MSGARCH_MLE_FIT <- function(object, x = NULL, theta = NULL, y = NULL,
                                log = FALSE, do.its = FALSE) {
  return(MSGARCH::cdf(object = object$spec, x = x, theta = object$theta,
                      y = object$y, log = log, do.its = do.its))
}

#' @export
cdf.MSGARCH_BAY_FIT <- function(object, x = NULL, theta = NULL, y = NULL,
                                log = FALSE, do.its = FALSE) {
  return(MSGARCH::cdf(object = object$spec, x = x, theta = object$theta,
                      y = object$y, log = log, do.its = do.its))
}
