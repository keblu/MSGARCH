#'Probability Integral Transform.
#' @description Method returning the predictive probability integral transform (PIT) in-sample or of a vector of points  consider as one step ahead draws (\code{t = T + 1}).
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param x Vector (of size N) of pointsevaluated at \code{t = T + 1} (used when \code{do.its = FALSE}).
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not required when using a fit object) where d must have
#'  the same length as the default parameters of the specification.
#' @param y  Vector (of size T) of observations (not required when using a fit object).
#' @param do.norm  Boolean indicating if the PIT values are transformed into standard Normal variate. (Default: \code{do.norm = FALSE}).
#' @param do.its  Boolean indicating if the in-sample pit is returned. (Default: \code{do.its = FALSE})
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian probability integral transform is calculated.
#' @details If a matrix of parameter estimates is given, each parameter estimate (each row) is evaluated individually. 
#' If \code{do.its = FALSE}, the vector \code{x} are evaluated as \code{t = T + 1} realization and the method uses the variance estimate at \code{t = T + 1}.
#' If \code{do.its = TRUE} and  each column of  \code{x} is evaluated a their respective time \code{t} indicated by their column index.
#' Finally if \code{x = NULL} the vector \code{y} is evaluated using their respective variance estimate at each time \code{t}.
#' The \code{do.norm} argument transforms the PIT value into Normal variates so that normality test can be done.
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
#' \item \code{do.its}: Orinigal user inputed \code{do.its} for reference.
#' }
#' @importFrom stats qnorm
#' @export
pit <- function(object, x, theta, y, do.norm = FALSE, do.its = FALSE) {
  UseMethod("pit", object)
}

#' @export
pit.MSGARCH_SPEC <- function(object, x = NULL, theta = NULL, y = NULL, do.norm = FALSE, do.its = FALSE) {
  y <- MSGARCH:::f.check.y(y)
  if (is.vector(theta)) {
    theta <- matrix(theta, nrow = 1)
  }
  theta_check <- MSGARCH:::f.check.theta(object, theta)
  if (isTRUE(do.its)) {
    if(is.null(x)){
      x = matrix(data = y, ncol = length(y))
    } else {
      x = matrix(x)
      if(ncol(x) == 1){
        x = matrix(x,ncol = length(y), nrow = nrow(x))
      } else {
          stop("x have more than 1 column:x must be a vector, NULL, or a matrix of size N x 1")
      }
    }
    tmp <- matrix(data = 0, nrow = nrow(x), ncol = length(y))
    for (i in 1:nrow(theta)) {
      if (object$K == 1) {
        tmp2 <- object$rcpp.func$cdf_Rcpp_its(theta_check[i,], y, x, FALSE)
        tmp  <-  tmp + tmp2[,,1]
      } else {
        Pstate <- MSGARCH::Pstate(object = object, theta = theta[i, ], y = y)
        Pstate.tmp <- matrix(data = NA, nrow = dim(Pstate)[1], ncol = dim(Pstate)[3])
        for (j in 1:dim(Pstate)[3]) {
          Pstate.tmp[, j] <- Pstate[, , j]
        }
        tmp2 <- object$rcpp.func$cdf_Rcpp_its(theta_check[i,], y, x, FALSE)
        for(k in 1:object$K){
          tmp = tmp + tmp2[,,k] * matrix(Pstate.tmp[1:(nrow(Pstate.tmp) - 1), k],ncol  = length(y), nrow = nrow(x),byrow = TRUE)
        }
      }
      
    }
    tmp = tmp/nrow(theta)
  } else {
    x = matrix(x)
    if(ncol(x) !=1){
      stop("x must be a vector or a matrix of size N x 1")
    }
    tmp <- matrix(data = 0, nrow = nrow(x), ncol = 1)
    for (i in 1:nrow(theta)) {
      # DA we need to check the inputs
      str = "CDF FAIL IN CPP"
      is.ok = tryCatch({
        tmp <- tmp + object$rcpp.func$cdf_Rcpp(x, theta_check[i,], y, FALSE)
        TRUE
      }, warning = function(warn) {
        MSGARCH:::f.error(str)
        browser()
      }, error = function(err) {
        MSGARCH:::f.error(str)
        browser()
      })
    }
    tmp = tmp/nrow(theta)
  }
  if (do.norm) {
    tmp <- qnorm(tmp, mean = 0, sd = 1)
  }
  out <- list()
  out$pit <- tmp
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