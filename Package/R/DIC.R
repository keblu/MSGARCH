#' @title  Deviance Information Criterion (DIC).
#' @description Method which computes the Deviance Information Criterion (DIC) from a fit object of type
#' \code{MSGARCH_MCMC_FIT} created with \code{\link{FitMCMC}}.
#' @param fit Fit object of type \code{MSGARCH_MCMC_FIT} created
#' with \code{\link{FitMCMC}}.
#' @return A list with the following elements:
#'        \itemize{
#'        \item \code{DIC}: Deviance Information Criterion.
#'        \item \code{IC}: Bayesian Predictive Information Criterion (IC = 2 * pV + D.bar).
#'        \item \code{pV}: Effective number of parameters (pV = var(D)/2).
#'        \item \code{D.bar}: Expected value of the deviance over the posterior.
#'        }
#' @details Computes the Deviance information criterion of Spiegelhalter et al. (2002).
#' @references Spiegelhalter, David J., et al. (2002).
#' Bayesian measures of model complexity and fit.
#' \emph{Journal of the Royal Statistical Society: Series B}, 64, 583-639
#' @examples
#' # load data
#' data("SMI", package = "MSGARCH")
#'
#' # create model specification
#' # MS(2)-GARCH(1,1)-Normal (default)
#' spec <- CreateSpec()
#'
#' # fit the model on data by MCMC
#' set.seed(123)
#' fit <- FitMCMC(spec = spec, data = SMI, ctr = list(n.burn = 500L, n.mcmc = 500L))
#'
#' # compute DIC
#' DIC(fit)
#' @importFrom stats var
#' @export
DIC <- function(fit) {
  UseMethod(generic = "DIC", object = fit)
}

#' @rdname DIC
#' @export
DIC.MSGARCH_MCMC_FIT <- function(fit) {
  out <- f_DIC(spec = fit$spec, par = fit$par, data = fit$data)
  return(out)
}

f_DIC <- function(spec, par, data) {
  spec <- f_check_spec(spec)
  data <- f_check_y(data)
  if (is.vector(x = par)) {
    par <- matrix(data = par, nrow = 1L)
  }
  LL <- vector(mode = "numeric", length = nrow(par))
  for (i in 1:nrow(par)) {
    LL[i] <- Kernel(object = spec, par = par[i, ], data = data, log = TRUE, do.prior = FALSE)
  }
  D.bar <- -2 * mean(x = LL)
  pV    <- stats::var(x = -2 * LL)/2
  out   <- list(DIC = pV + D.bar, IC = 2 * pV + D.bar, pV = pV, D.bar = D.bar)
  return(out)
}
