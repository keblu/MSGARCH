#' @rdname predict
#' @title predict method.
#' @description Method returning conditional volatility forecasts and density forecasts  of the process.
#' @param object Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}} or fit object of type \code{MSGARCH_ML_FIT}
#' created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}.
#' @param newdata Vector (of size T*) of new observations. (Default \code{newdata = NULL})
#' @param nahead  Scalar indicating the number of step-ahead evaluation.
#' @param do.return.draw  Are the sampled simulation draws returned? (Default \code{do.return.draw = FALSE})
#' @param par Vector (of size d) or matrix (of size \code{nmcmc} x d) of
#' parameter estimates where d must have
#' the same length as the default parameters of the specification.
#' @param ctr A list of control parameters:
#'        \itemize{
#'        \item \code{nsim} (integer >= 0):
#'        Number indicating the number of simulation done for the
#'        conditional volatlity forecast at \code{nahead > 1}. (Default: \code{nsim = 10000L})
#'        }
#' @param ... Not used. Other arguments to \code{Forecast}.
#' @return A list of class \code{MSGARCH_CONDVOL} with the following elements:
#' \itemize{
#'  \item \code{vol}: Condititional volatility Forecast (vector of size \code{nahead}).
#'  \item \code{draw}: If \code{do.return.draw = TRUE}:\cr
#'  Draws sample from the predictive distributions  (matrix of size \code{nahead} x \code{nsim}).\cr
#'  If \code{do.return.draw = FALSE}: \code{NULL}
#'  }
#' The \code{MSGARCH_FORECAST} class contains the \code{plot} method.
#' @details If a matrix of MCMC posterior draws is given, the
#' Bayesian predictive conditional volatility forecasts are calculated.
#' @examples
#' # load data
#' data("SMI", package = "MSGARCH")
#'
#' # create model specification
#' # MS(2)-GARCH(1,1)-Normal (default)
#' spec <- CreateSpec()
#'
#' # fit the model on the data by ML
#' fit <- FitML(spec = spec, data = SMI)
#'
#' # compute the In-sample conditional volatility from the fitted model
#' forecast <- predict(object = fit, nahead = 5L)
#' plot(forecast)

#' @rdname predict
#' @export
predict.MSGARCH_SPEC <- function(object, newdata = NULL, nahead = 1L, do.return.draw = FALSE, par = NULL, ctr = list(), ...) {
  data <- c(object$data, newdata)
  out  <- f_CondVol(object = object, par = par, data = newdata, nahead = nahead,
                    do.its = FALSE, ctr = ctr)
  if(!isTRUE(do.return.draw)){
    out$draw = NULL
  }
  class(out) <- "MSGARCH_FORECAST"
  return(out)
}

#' @rdname predict
#' @export
predict.MSGARCH_ML_FIT <- function(object, newdata = NULL, nahead = 1L, do.return.draw = FALSE, ctr = list(), ...) {
  data <- c(object$data, newdata)
  out  <- f_CondVol(object = object$spec, par = object$par, data = data, nahead = nahead,
                    do.its = FALSE, ctr = ctr)
  if(!isTRUE(do.return.draw)){
    out$draw = NULL
  }
  class(out) <- "MSGARCH_FORECAST"
  return(out)
}

#' @rdname predict
#' @export
predict.MSGARCH_MCMC_FIT <- function(object, newdata = NULL, nahead = 1L, do.return.draw = FALSE, ctr = list(), ...) {
  data <- c(object$data, newdata)
  out  <- f_CondVol(object = object$spec, par = object$par, data = data, nahead = nahead,
                    do.its = FALSE, ctr = ctr)
  if(!isTRUE(do.return.draw)){
    out$draw = NULL
  }
  class(out) <- "MSGARCH_FORECAST"
  return(out)
}