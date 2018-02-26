#' @title Volatility filtering.
#' @description Method returning the in-sample conditional volatility.
#' @param object Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}} or fit object of type \code{MSGARCH_ML_FIT}
#' created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}.
#' @param par Vector (of size d) or matrix (of size \code{nmcmc} x d) of parameter
#' estimates where d must have the same length as the default parameters of the specification.
#' @param data  Vector (of size T) of observations.
#' @param newdata Vector (of size T*) of new observations. (Default \code{newdata = NULL})
#' @param ... Not used. Other arguments to \code{Volatility}.
#' @return In-sample condititional volatility (vector of size T + T*) of class \code{MSGARCH_CONDVOL}.\cr
#' The \code{MSGARCH_CONDVOL} class contains the \code{plot} method.
#' @details If a matrix of MCMC posterior draws is given, the
#' Bayesian predictive conditional volatility is calculated.
#' @examples
#' # create specification
#' spec <- CreateSpec()
#' 
#' # load data
#' data("SMI", package = "MSGARCH")
#' 
#' # in-sample volatility from specification
#' par <- c(0.1, 0.1, 0.8, 0.2, 0.1, 0.8, 0.99, 0.01)
#' vol <- Volatility(object = spec, par = par, data = SMI)
#' head(vol)
#' plot(vol)
#' 
#' # in-sample volatility from ML fit
#' fit <- FitML(spec = spec, data = SMI)
#' vol <- Volatility(object = fit)
#' head(vol)
#' plot(vol)
#' 
#' \dontrun{
#' # in-sample volatility from MCMC fit
#' set.seed(1234)
#' fit <- FitMCMC(spec = spec, data = SMI)
#' vol <- Volatility(object = fit)
#' head(vol)
#' plot(vol)
#' }
#' @export
Volatility <- function(object, ...) {
  UseMethod(generic = "Volatility", object)
}

#' @rdname Volatility
#' @export
Volatility.MSGARCH_SPEC <- function(object, par, data, ...) {
  out  <- f_CondVol(object = object, par = par, data = data,
                    do.its = TRUE, ctr = list())
  return(out$vol)
}

#' @rdname Volatility
#' @export
Volatility.MSGARCH_ML_FIT <- function(object, newdata = NULL, ...) {
  data <- c(object$data, newdata)
  if(is.ts(object$data)){
    if(is.null(newdata)){
      data = zoo::zooreg(data, order.by =  c(zoo::index(data)))
    } else {
      data = zoo::zooreg(data, order.by =  c(zoo::index(data),zoo::index(data)[length(data)]+(1:length(newdata))))
    }
    data = as.ts(data)
  }
  out  <- f_CondVol(object = object$spec, par = object$par, data = data,
                  do.its = TRUE, ctr = list())
  return(out$vol)
}

#' @rdname Volatility
#' @export
Volatility.MSGARCH_MCMC_FIT <- function(object, newdata = NULL, ...) {
  data <- c(object$data, newdata)
  if(is.ts(object$data)){
    if(is.null(newdata)){
      data = zoo::zooreg(data, order.by =  c(zoo::index(data)))
    } else {
      data = zoo::zooreg(data, order.by =  c(zoo::index(data),zoo::index(data)[length(data)]+(1:length(newdata))))
    }
    data = as.ts(data)
  }
  out  <- f_CondVol(object = object$spec, par = object$par, data = data,
                  do.its = TRUE, ctr = list())
  return(out$vol)
}