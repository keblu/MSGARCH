#' @title  Bayesian information criterion (BIC).
#' @description Method which computes the Bayesian information criterion (BIC) from a fit object of type
#' \code{MSGARCH_ML_FIT} created with \code{\link{FitML}}
#' or \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}.
#' @param fit Fit object of type \code{MSGARCH_ML_FIT} created with \code{\link{FitML}}
#' or \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}
#' @return BIC value.
#' @details Computes the Bayesian information criterion (BIC) based on the work of Schwarz (Schwarz, 1978).
#' If a matrix of MCMC posterior draws is given, the BIC on the posterior mean is calculated.
#' @references Schwarz, G. (1978).
#' Estimating the dimension of a model.
#' \emph{Annals of Statistics}, 6, 461-464.
#' @examples
#' # load data
#' data("SMI", package = "MSGARCH")
#'
#' # create model specification
#' # MS(2)-GARCH(1,1)-Normal (default)
#' spec <- CreateSpec()
#'
#' # fit the model on data by ML
#' fit <- FitML(spec = spec, data = SMI)
#'
#' # compute BIC
#' BIC(fit)
#' @export
BIC <- function(fit) {
  UseMethod(generic = "BIC", object = fit)
}

#' @rdname BIC
#' @export
BIC.MSGARCH_ML_FIT <- function(fit) {
  out <- f_BIC(spec = fit$spec, par = fit$par, data = fit$data)
  return(out)
}

#' @rdname BIC
#' @export
BIC.MSGARCH_MCMC_FIT <- function(fit) {
  out <- f_BIC(spec = fit$spec, par = fit$par, data = fit$data)
  return(out)
}

f_BIC <- function(spec, par, data) {
  spec <- f_check_spec(spec)
  data <- f_check_y(data)
  if (is.vector(x = par)) {
    par <- matrix(data = par, nrow = 1L)
  }
  par <- matrix(colMeans(x = par), nrow = 1L)
  LL  <- Kernel(object = spec, par = par, data = data, log = TRUE, do.prior = FALSE)
  k   <- dim(x = par)[2L]
  n   <- length(data)
  out <- k * log(x = n) - 2 * LL
  return(out)
}
