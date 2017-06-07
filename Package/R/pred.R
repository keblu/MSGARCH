#' Predictive function.
#' @description Method returning the predictive probability density in-sample or of a vector of points  consider as one step ahead draws (\code{t = T + 1}).
#' @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param x Vector (of size N).
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not required when using a fit object) where d must have
#'  the same length as the default parameters of the specification.
#' @param y  Vector (of size T) of observations (not required when using a fit object).
#' @param log  Boolean indicating if the log-density is returned. (Default: \code{log = FALSE})
#' @param do.its  Boolean indicating if the in-sample predictive is returned. (Default: \code{do.its = FALSE})
#' @details If a matrix of parameter estimates is given, each parameter estimate (each row) is evaluated individually. 
#' If \code{do.its = FALSE}, the vector \code{x} are evaluated as \code{t = T + 1} realization and the method uses the variance estimate at \code{t = T + 1}.
#' If \code{do.its = TRUE} and  each column of  \code{x} is evaluated a their respective time \code{t} indicated by their column index.
#' Finally if \code{x = NULL} the vector \code{y} is evaluated using their respective variance estimate at each time \code{t}.
#' @examples 
#'require("MSGARCH")
#'# load data
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
#'# run pred method in-sample     
#'pred.its = MSGARCH::pred(object = fit, log = TRUE, do.its = TRUE)  
#'
#'sum(pred.its$pred, na.rm = TRUE)
#'                                              
#'# create mesh
#'x = seq(-3,3,0.01)
#'
#'# run pred method on mesh at T + 1
#'pred = MSGARCH::pred(object = fit, x = x, log = TRUE, do.its = FALSE)
#'
#'plot(pred)
#' @return A list of class \code{MSGARCH_PRED} containing two components:
#' \itemize{
#' \item \code{pred}:\cr If \code{do.its = FALSE}: (Log-)Predictive of of the points \code{x} at \code{t = T + 1} (vector of size N). \cr
#'                   If \code{do.its = TRUE}: In-sample Predictive of \code{y} (vector of size T or matrix of size M x T). 
#' \item \code{do.its}: Original user inputed \code{do.its} for reference.
#' }
#' @export
pred <- function(object, x, theta, y, log = FALSE, do.its = FALSE) {
  UseMethod("pred", object)
}

#' @export
pred.MSGARCH_SPEC <- function(object, x = NULL, theta = NULL, y = NULL, log = FALSE, do.its = FALSE) {
  y <- MSGARCH:::f.check.y(y)
  if (is.vector(theta)) {
    theta <- matrix(theta, nrow = 1)
  }
  theta_check <- MSGARCH:::f.check.theta(object, theta)
  if (isTRUE(do.its)) {
    if(is.null(x)){
      x = matrix(data = y,ncol = length(y))
    } else {
      x = matrix(x)
      if(ncol(x) == 1){
        x = matrix(x,ncol = length(y), nrow = nrow(x))
      } else {
          stop("x have more than 1 column: x must be a vector, NULL, or a matrix of size N x 1")
      }
    }
    tmp <- matrix(data = 0, nrow = nrow(x), ncol = length(y))
    for (i in 1:nrow(theta)) {
      if (object$K == 1) {
        tmp2 <- object$rcpp.func$pdf_Rcpp_its(theta_check[i,], y, x, FALSE)
        tmp  <- tmp + tmp2[,,1]
      } else {
        Pstate <- MSGARCH::Pstate(object = object, theta = theta[i, ], y = y)
        Pstate.tmp <- matrix(data = NA, nrow = dim(Pstate)[1], ncol = dim(Pstate)[3])
        for (j in 1:dim(Pstate)[3]) {
          Pstate.tmp[, j] <- Pstate[, , j]
        }
        tmp2 <- object$rcpp.func$pdf_Rcpp_its(theta_check[i,], y, x, FALSE)
        for(k in 1:object$K){
          tmp = tmp + tmp2[,,k] * matrix(Pstate.tmp[1:(nrow(Pstate.tmp) - 1), k],ncol  = length(y), nrow = nrow(x),byrow = TRUE)
        }
      }
    }
    tmp = tmp/nrow(theta)
  } else {
    if(is.null(x)){
      stop("x is NULL: x must be a vector or a matrix of size N x 1")
    }
    x = matrix(x)
    if(ncol(x) !=1){
      stop("x have more than 1 column: x must be a vector or a matrix of size N x 1")
    }
    tmp <- matrix(data = 0, nrow = nrow(x), ncol = 1)
    for (i in 1:nrow(theta)) {
      # DA we need to check the inputs
      str = "PDF FAIL IN CPP"
      is.ok = tryCatch({
        tmp <- tmp + object$rcpp.func$pdf_Rcpp(x, theta_check[i,], y, FALSE)
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
  if (log) {
    tmp <- log(tmp)
  }
  out <- list()
  out$pred    <- tmp
  out$do.its <- do.its
  class(out) <- "MSGARCH_PRED"
  return(out)
}

#' @export
pred.MSGARCH_MLE_FIT <- function(object, x = NULL, theta = NULL, y = NULL, log = FALSE,
                                do.its = FALSE) {
  return(MSGARCH::pred(object = object$spec, x = x, theta = object$theta, y = object$y,
                       log = log, do.its = do.its))
}

#' @export
pred.MSGARCH_BAY_FIT <- function(object, x = NULL, theta = NULL, y = NULL, log = FALSE,
                                do.its = FALSE) {
  return(MSGARCH::pred(object = object$spec, x = x, theta = object$theta, y = object$y,
                       log = log, do.its = do.its))
}