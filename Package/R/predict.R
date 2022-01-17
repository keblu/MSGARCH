#' @rdname predict
#' @title predict method.
#' @description Conditional volatility (and predictive distribution) forecasts.
#' @param object Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}} or fit object of type \code{MSGARCH_ML_FIT}
#' created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}.
#' @param newdata Vector (of size T*) of new observations. (Default \code{newdata = NULL})
#' @param nahead  Scalar indicating the number of step-ahead evaluation.
#' @param do.return.draw  Are simulation draws from the predictive distribution 
#' returned? (Default \code{do.return.draw = FALSE})
#' @param par Vector (of size d) or matrix (of size \code{nmcmc} x d) of
#' parameter estimates where d must have
#' the same length as the default parameters of the specification.
#' @param do.cumulative Logical indicating if the conditional volatility 
#' prediction is computed on the cumulative simulations (typically log-returns, as they can be aggregated).
#' (Default: \code{do.cumulative = FALSE})
#' @param ctr A list of control parameters:
#'        \itemize{
#'        \item \code{nsim} (integer >= 0):
#'        Number indicating the number of simulation done for the
#'        conditional volatlity forecast at \code{nahead > 1}. (Default: \code{nsim = 10000L})
#'        }
#' @param ... Not used. Other arguments to \code{predict}.
#' @return A list of class \code{MSGARCH_FORECAST} with the following elements:
#' \itemize{
#'  \item \code{vol}: Condititional volatility forecast (vector of size \code{nahead}).
#'  \item \code{draw}: If \code{do.return.draw = TRUE}:\cr
#'  Draws sampled from the predictive distributions (matrix of size \code{nahead} x \code{nsim}).\cr
#'  If \code{do.return.draw = FALSE}:\cr
#'  \code{NULL}
#'  }
#' The \code{MSGARCH_FORECAST} class contains the \code{plot} method.
#' @details If a matrix of MCMC posterior draws is given, the
#' Bayesian predictive conditional volatility (and predictive distribution) 
#' forecasts are returned.
#' @examples
#' # create specification
#' spec <- CreateSpec()
#' 
#' # load data
#' data("SMI", package = "MSGARCH")
#' 
#' # predict from specification
#' par <- c(0.1, 0.1, 0.8, 0.2, 0.1, 0.8, 0.99, 0.01)
#' set.seed(1234)
#' pred <- predict(object = spec, par = par, newdata = SMI, nahead = 5L)
#' head(pred)
#' plot(pred)
#' 
#' # predict from ML fit
#' fit <- FitML(spec = spec, data = SMI)
#' set.seed(1234)
#' pred <- predict(object = fit, nahead = 5L, do.return.draw = TRUE)
#' head(pred)
#' plot(pred)
#' 
#' \dontrun{
#  # predict from MCMC fit
#' set.seed(1234)
#' fit <- FitMCMC(spec = spec, data = SMI)
#' pred <- predict(object = fit, nahead = 5L, do.return.draw = TRUE)
#' plot(pred)
#' }
#' @importFrom stats predict
#' @export
predict.MSGARCH_SPEC <- function(object, newdata = NULL, nahead = 1L, 
                                 do.return.draw = FALSE, par = NULL, 
                                 do.cumulative = FALSE, ctr = list(), ...) {
  out  <- f_CondVol(object = object, par = par, data = newdata, nahead = nahead,
                    do.its = FALSE, do.cumulative = do.cumulative, ctr = ctr)
  if(!isTRUE(do.return.draw)){
    out$draw = NULL
  }
  
  class(out) <- "MSGARCH_FORECAST"
  return(out)
}

#' @rdname predict
#' @export
predict.MSGARCH_ML_FIT <- function(object, newdata = NULL, 
                                   nahead = 1L, do.return.draw = FALSE, 
                                   do.cumulative = FALSE, ctr = list(), ...) {
  data <- c(object$data, newdata)
  if(is.ts(object$data)){
    if(is.null(newdata)){
      data = zoo::zooreg(data, order.by =  c(zoo::index(data)))
    } else {
      data = zoo::zooreg(data, order.by =  c(zoo::index(data),zoo::index(data)[length(data)]+(1:length(newdata))))
    }
    data = as.ts(data)
  }
  out  <- f_CondVol(object = object$spec, par = object$par, data = data, nahead = nahead,
                    do.its = FALSE, do.cumulative = do.cumulative, ctr = ctr)
  if(!isTRUE(do.return.draw)){
    out$draw = NULL
  }
  class(out) <- "MSGARCH_FORECAST"
  return(out)
}

#' @rdname predict
#' @export
predict.MSGARCH_MCMC_FIT <- function(object, newdata = NULL, nahead = 1L, 
                                     do.return.draw = FALSE, do.cumulative = FALSE, ctr = list(), ...) {
  data <- c(object$data, newdata)
  if(is.ts(object$data)){
    if(is.null(newdata)){
      data = zoo::zooreg(data, order.by =  c(zoo::index(data)))
    } else {
      data = zoo::zooreg(data, order.by =  c(zoo::index(data),zoo::index(data)[length(data)]+(1:length(newdata))))
    }
    data = as.ts(data)
  }
  out  <- f_CondVol(object = object$spec, par = object$par, data = data, nahead = nahead,
                    do.its = FALSE, do.cumulative = do.cumulative, ctr = ctr)
  if(!isTRUE(do.return.draw)){
    out$draw = NULL
  }
  class(out) <- "MSGARCH_FORECAST"
  return(out)
}