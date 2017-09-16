#' @title Maximum Likelihood estimation.
#' @description Method that performs Maximum Likelihood estimation
#' of a \code{MSGARCH_SPEC} object on a set of observations.
#' @param spec Model specification created with \code{\link{CreateSpec}}.
#' @param data Vector (of size T) of observations.
#' @param ctr A list of control parameters:
#'        \itemize{
#'        \item \code{par0}: Vector (of size d) where d must have
#'         the same length as the default parameters of the specification.
#'         It is the starting value for the optimizer (if empty the
#'         the method automatically set starting parameters; see *Details*).
#'        \item \code{do.se} Logical. Should standard errors be computed?
#'        (Default: \code{do.se = TRUE}).
#'        \item \code{do.plm} Logical. If \code{do.plm = FALSE}, parameter transformation
#'        during the optimization step is performed without ensuring stationarity
#'        for the volatility processes. For combinations of parameters that do not
#'        imply stationarity the likelihood value is fixed at -1e10. If
#'        \code{fixed} is defined in the  \code{list} \code{contraint.spec}
#'        of \link{CreateSpec}, \code{do.plm = TRUE}
#'        is used. (Default: \code{do.plm = FALSE})
#'        \item \code{OptimFUN}: Custom optimization function (see *Details*).
#'        }
#' @return A list of class \code{MSGARCH_ML_FIT} with the following elements:
#'        \itemize{
#'        \item \code{par}: Vector (of size d) of optimal parameters.
#'        \item \code{loglik}: Log-likelihood of \code{y} given the optimal parameters.
#'        \item \code{Inference}: \code{list} with elements \code{MatCoef} and \code{Hessian}.
#'         \code{MatCoef} is a matrix (of size d x 4) with optimal parameter estimates, standard errors, t-stats, and p-values.
#'         \code{Hessian} is the Hessian (matrix of size d x d) of the negative log-likelihood function
#'         evaluated at the optimal parameter estimates \code{par}.
#'        \item \code{spec}: Model specification of class \code{MSGARCH_SPEC}
#'        created with \code{\link{CreateSpec}}.
#'        \item \code{data}: Vector (of size T) of observations.
#'        \item \code{ctr}: \code{list} of the control used for the fit.
#'        }
#' The \code{MSGARCH_ML_FIT} with the following methods:
#' \itemize{
#' \item \code{\link{AIC}}: Compute Akaike information criterion (AIC).
#' \item \code{\link{BIC}}: Compute Bayesian information criterion (BIC).
#' \item \code{\link{Volatility}}: In-sample conditional volatility filterting of the overall process.
#' \item \code{\link{Forecast}}: Forecast of the conditional volatility of the overall process.
#' \item \code{\link{UncVol}}: Unconditional volatility in each regime and the overall process.
#' \item \code{\link{Pred}}: Predictive method.
#' \item \code{\link{PIT}}: Probability Integral Transform.
#' \item \code{\link{Risk}}: Value-at-Risk and Expected-Shortfall methods.
#' \item \code{\link{Sim}}: Simulation method.
#' \item \code{\link{State}}: State probabilities methods.
#' \item \code{\link{ExtractStateFit}}: Single-regime model extractor.
#' \item \code{summary}: Summary of the fit.
#' }
#' @details By default, \code{OptimFUN} is set such that optimization is done via the well known Broyden-
#' Fletcher-Goldfarb-Shanno (BFGS) algorithm using the \code{optim} function with \code{method =
#' "BFGS"}.
#' Starting values when \code{par0} is not provided are chosen automatically
#' before optimization (see Ardia et al. (2016) for more details)\cr
#' \code{OptimFUN} allows for a custom optimizer to be used. The function must take
#' the form: \cr \code{function(vPw, f_nll, spec, data, do.plm)}, \cr
#' where \code{vPw} are starting parameters (transformed), \code{f_nll} is the function
#' to be minimize, \code{spec} is the specification, \code{data} is the data,
#' and \code{do.plm} the originally inputed or default \code{do.plm}.
#' The inputs \code{spec}, \code{data}, and \code{do.plm}
#' must be passed as inputs in the optimizer (see *Examples*).
#' It must output a list with the following elements:
#' \itemize{
#' \item \code{value}: Optimal negative log-likelihood.
#' \item \code{par}: Optimal parameters.
#' }
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
#' summary(fit)
#'
#' # custom optimizer example
#'f_custom_optim <- function(vPw, f_nll, spec, data, do.plm){
#'  out <- stats::optim(vPw, f_nll, spec = spec, data = data,
#'                      do.plm = do.plm, method = "Nelder-Mead")
#'  return(out)
#'}
#'
#' set.seed(123)
#' fit <- FitML(spec, data = SMI, ctr = list(OptimFUN = f_custom_optim))
#' summary(fit)
#' @importFrom stats runif
#' @export
FitML <- function(spec, data, ctr = list()) {
  UseMethod(generic = "FitML", spec)
}

#' @export
FitML.MSGARCH_SPEC <- function(spec, data, ctr = list()) {
  
  time.start <- Sys.time()
  spec <- f_check_spec(spec)
  data <- f_check_y(data)
  ctr  <- f_process_ctr(ctr)
  
  if ((isTRUE(spec$fixed.pars.bool)) || (isTRUE(spec$regime.const.pars.bool))) {
    f_check_fixedpars(spec$fixed.pars, spec)
    ctr$do.plm <- TRUE
  }
  
  if (is.null(ctr$par0)) {
    vPw  <- f_StargingValues(data, spec, ctr)
    par0 <- matrix(f_mapPar(vPw, spec, ctr$do.plm), nrow = 1L, dimnames = list(NULL, names(vPw)))
  } else {
    par0 = ctr$par0
    vPw <- f_unmapPar(par0, spec, ctr$do.plm)
    if (isTRUE(spec$regime.const.pars.bool)) {
      vPw <- f_remove_regimeconstpar(vPw, spec$regime.const.pars, spec$K)
    }
    if (isTRUE(spec$fixed.pars.bool)) {
      vPw <- f_substitute_fixedpar(vPw, spec$fixed.pars)
    }
  }
  optimizer <- ctr$OptimFUN(vPw, f_nll, spec, data, ctr$do.plm)
  
  llk <- -optimizer$value
  
  if (llk == 1e+10) {
    str <- "FitML -> Error during optimization"
    f_error(str)
    stop()
  }
  
  vPw <- optimizer$par
  vPn <- f_mapPar(vPw, spec, ctr$do.plm)
  np <- length(vPw)
  
  if (isTRUE(spec$fixed.pars.bool)) {
    vPn <- f_add_fixedpar(vPn, spec$fixed.pars)
    vPn <- vPn[colnames(spec$par0)]
  }
  
  if (isTRUE(spec$regime.const.pars.bool)) {
    vPn <- f_add_regimeconstpar(vPn, spec$K, spec$label)
  }
  
  par <- matrix(vPn, nrow = 1L, dimnames = list(NULL, names(vPn)))
  par <- f_sort_par(spec, par)
  par <- as.vector(par)
  names(par) <- spec$label
  vPww <- f_unmapPar(par, spec, ctr$do.plm)
  
  if (isTRUE(spec$regime.const.pars.bool)) {
    vPww <- f_remove_regimeconstpar(vPww, spec$regime.const.pars, spec$K, for.se = TRUE)
  }
  if (isTRUE(spec$fixed.pars.bool)) {
    vPww <- f_substitute_fixedpar(vPww, spec$fixed.pars)
  }
  
  elapsed.time <- Sys.time() - time.start
  
  if (isTRUE(ctr$do.se)) {
    Inference <- f_InferenceFun(vPww, data, spec, do.plm = ctr$do.plm)
  } else {
    Inference <- NULL
  }
  
  out <- list(par = par, loglik = llk, spec = spec, data = data,
              Inference = Inference, ctr = ctr)
  
  class(out) <- "MSGARCH_ML_FIT"
  return(out)
}
