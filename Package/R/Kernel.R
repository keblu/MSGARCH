# @title Kernel function.
# @description Method returning the kernel value of a vector of observations given a model specification.
# @param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{CreateSpec}}
# or fit object of type \code{MSGARCH_ML_FIT} created with \code{\link{FitML}} or \code{MSGARCH_MCMC_FIT}
# created with \code{\link{FitMCMC}}.
# @param par Vector (of size d) or matrix (of size \code{nmcmc} x d) of parameter
# estimates (not required when using a fit object) where d must have the same length
# as the default parameters of the specification.
# @param data  Vector (of size T) of observations (not required when using a fit object).
# @param log  Logical indicating if the log kernel is returned. (Default: \code{log = TRUE})
# @param do.prior  Logical indicating if the prior is evaluated. (Default: \code{do.prior = TRUE})
# @return (Log-)kernel value (scalar or vector of size \code{nmcmc}) of the vector of observations.
# @details If a matrix of parameter estimates is given, each parameter estimate (each row) is evaluated individually.
#  The kernel is a combination of the prior and the likelihood function.
#  The kernel is equal to LP(\eqn{\psi}) + LL(data|\eqn{\psi}) where LP is the log-prior of \eqn{\psi}
#  and LL is the log-likelihood of \code{data} given the parameter \eqn{\psi}.\cr
#  The prior is different for each specification. It ensures that the \eqn{\psi} makes
#  the conditional variance process covariance-stationary, positive,
#  and that the sums of the probabilities in the case of a
#  multiple-regime model are all equal to one. If any of these three conditions is
#  not respected the prior returns \code{-1e10}, meaning that the optimizer or the sampler
#  will know that \eqn{\psi} is not a good candidate.
# @examples
# # load data
# data("SMI", package = "MSGARCH")
#
# # create model specification
# # MS(2)-GARCH(1,1)-Normal (default)
# spec <- CreateSpec()
#
# # fit the model on the data by ML
# fit <- FitML(spec = spec, data = SMI)
#
# # compute the kernel
# Kernel(fit, log = TRUE)
# @export
#Kernel <- function(object, par, data, log = TRUE, do.prior = FALSE) {
#  UseMethod(generic = "Kernel", object)
#}

# @export
Kernel <- function(object, par, data = NULL, log = TRUE, do.prior = FALSE) {
  object <- f_check_spec(object)
  data   <- f_check_y(data)
  par    <- f_check_par(object, par)
  lnd    <- object$rcpp.func$eval_model(par, data, do.prior)
  
  if (isTRUE(do.prior)) {
    prior.mean = object$rcpp.func$get_mean()
    prior.sd = object$rcpp.func$get_sd()
    names(prior.mean) = names(prior.sd) = object$label[1:length(prior.mean)]
    if (isTRUE(object$fixed.pars.bool)) {
      lnd = lnd - sum(dnorm(par[, names(object$fixed.pars)],
                            prior.mean[names(object$fixed.pars)],
                            prior.sd[names(object$fixed.pars)],
                            log = TRUE))
    }
    
    if (isTRUE(object$regime.const.pars.bool) & object$K > 1) {
      lnd = lnd - sum(dnorm(par[, paste0(object$regime.const.pars, "_", (2:object$K))],
                            prior.mean[paste0(object$regime.const.pars, "_", (2:object$K))],
                            prior.sd[paste0(object$regime.const.pars, "_", (2:object$K))], log = TRUE))
    }
  }
  
  lnd[is.na(lnd) | is.nan(lnd) | is.infinite(lnd)] <- -1e+10
  if (!log)
    lnd <- exp(lnd)
  return(lnd)
}
