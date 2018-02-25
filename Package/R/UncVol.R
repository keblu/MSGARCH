#' @title Unconditional volatility.
#' @description Method returning the unconditional volatility of the process.
#' @param object Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}} or fit object of type \code{MSGARCH_ML_FIT}
#' created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT}
#' created with \code{\link{FitMCMC}}.
#' @param par Vector (of size d) or matrix (of size \code{nmcmc} x d) of parameter
#' estimates where d must have
#' the same length as the default parameters of the specification.
#' @param ... Not used. Other arguments to \code{UncVol}.
#' @param ctr A list of control parameters:
#'        \itemize{
#'        \item \code{nsim} (integer >= 0) :
#'        Number of simulations used for the estimation of the
#'        unconditional volatility. (Default: \code{nsim = 250L})
#'        \item \code{nahead} (integer >= 0) :
#'        Number of step ahead performed to estimate the
#'        unconditional volatility .(Default: \code{nahead = 5000L})
#'        \item \code{nburn} (integer >= 0) :
#'        Number of discarded step to estimate the
#'        unconditional volatility. (Default: \code{nburn = 1000L})
#'        }
#' @return A \code{scalar} of unconditional volatility.
#' @details If a matrix of MCMC posterior draws is given, the
#' Bayesian unconditional volatility is calculated.
#'  The unconditional volatility is estimated by first simulating \code{nsim}
#'  paths up to \code{nburn + nahead},
#'  calculating a forecast of the conditional volatility at each step ahead,
#'  discarding the first \code{nburn} step ahead conditional volatilities forecasts,
#'  and computing the mean of the remaining \code{nahead - nburn} conditional
#'  volatilites forecasts. This method is based on the fact that
#'  the conditional volatility forecast will converge to the unconditional volatilty
#'  the further the forecast his from the starting point.
#'  We take the average as a way to remove the noise that comes with the simulation process.
#'  Overall, this method allows to compute the unconditional volatilty complex models.
#' @examples
#' # create model specification
#' # MS(2)-GARCH(1,1)-Normal (default)
#' spec <- CreateSpec()
#'
#' \dontrun{
#'   # compute the unconditional volatility of the process
#'   par <- c(0.1, 0.1, 0.8, 0.2, 0.1, 0.8, 0.99, 0.01)
#'   UncVol(object = spec, par = par)
#'   
#'   # load data
#'   data("SMI", package = "MSGARCH")
#'   
#'   # fit the model on the data by ML
#'   fit <- FitML(spec = spec, data = SMI)
#'   UncVol(object = fit)
#' }
#' @export
UncVol <- function(object, ...) {
  UseMethod(generic = "UncVol", object = object)
}

#' @rdname UncVol
#' @export
UncVol.MSGARCH_SPEC <- function(object, par = NULL, ctr = list(), ...) {
  object <- f_check_spec(object)

  if (is.vector(par)) {
    par <- t(as.matrix(par))
  }
  if (nrow(par) == 1) {
    ctr   <- f_process_ctr(ctr, type = 2)
    nsim <- ctr$nsim
  } else {
    if(is.null(ctr$nsim)){
      nsim = 1
    } else {
      nsim = ctr$nsim
    }
  }
  ctr    <- f_process_ctr(ctr, type = 2)
  tmp <- f_CondVol(object = object,
                 par = par,
                 data = c(1,1),
                 do.its = FALSE,
                 nahead = ctr$nburn + ctr$nahead,
                 ctr = list(nsim = nsim))$vol
  out <- mean(tmp[ctr$nburn:ctr$nahead])
  return(out)
}

#' @rdname UncVol
#' @export
UncVol.MSGARCH_ML_FIT <- function(object, ctr = list(), ...) {
  out <- UncVol(object = object$spec, par = object$par, ctr = ctr)
  return(out)
}

#' @rdname UncVol
#' @export
UncVol.MSGARCH_MCMC_FIT <- function(object, ctr = list(), ...) {
  out <- UncVol(object = object$spec, par = object$par, ctr = ctr)
  return(out)
}
