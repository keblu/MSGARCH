#' @title  Akaike information criterion (AIC).
#' @description Method which computes the Akaike information criterion (AIC) from a fit object of type
#' \code{MSGARCH_ML_FIT} created with \code{\link{FitML}}
#' or \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}.
#' @param fit Fit object of type \code{MSGARCH_ML_FIT} created with \code{\link{FitML}}
#' or \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}.
#' @return AIC value.
#' @details Computes the Akaike information criterion (AIC) based on the work of Akaike (Akaike, 1974).
#' If a matrix of MCMC posterior draws is given, the AIC on the posterior mean is calculated.
#' @references Akaike, H. (1974).
#' A new look at the statistical model identification.
#' \emph{IEEE Transactions on Automatic Control}, 19, 716-723.
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
#' # compute AIC
#' AIC(fit)
#' @export
AIC <- function(fit) {
  UseMethod(generic = "AIC", object = fit)
}

#' @rdname AIC
#' @export
AIC.MSGARCH_ML_FIT <- function(fit) {
  out <- f_AIC(spec = fit$spec, par = fit$par, data = fit$data)
  return(out)
}

#' @rdname AIC
#' @export
AIC.MSGARCH_MCMC_FIT <- function(fit) {
  out <- f_AIC(spec = fit$spec, par = fit$par, data = fit$data)
  return(out)
}

f_AIC <- function(spec, par, data) {
  spec <- f_check_spec(spec)
  data <- f_check_y(data)
  if (is.vector(x = par)) {
    par <- matrix(data = par, nrow = 1L)
  }
  par <- matrix(data = colMeans(x = par), nrow = 1L)
  LL  <- Kernel(object = spec, par = par, data = data, log = TRUE, do.prior = FALSE)
  k   <- dim(x = par)[2L]
  out <- 2 * k - 2 * LL
  return(out)
}
