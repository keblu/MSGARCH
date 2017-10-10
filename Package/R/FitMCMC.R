#' @title MCMC/Bayesian estimation.
#' @description Method that performs MCMC/Bayesian estimation of
#' a \code{MSGARCH_SPEC} object on a set of observations.
#' @param spec Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}}.
#' @param data  Vector (of size T) of observations.
#' @param ctr A list of control parameters:
#'        \itemize{
#'        \item \code{par0}: Vector (of size d) where d must have
#'        the same length as the default parameters of the specification.
#'        It is the starting value for the chain (if empty the
#'        the method automatically set starting parameters; see *Details*).
#'        \item \code{n.burn} (integer >= 0): Number of discarded draws.
#'        (Default: \code{n.burn = 5000L})
#'        \item \code{n.mcmc} (integer > 0): Number of draws.
#'        (Default: \code{n.mcmc = 10000L})
#'        \item \code{n.thin} (integer > 0): Thinning factor (every \code{N.thin}
#'        draws are kept). (Default: \code{n.thin = 10L})
#'        \item \code{SamplerFUN}: Custom MCMC sampler (see *Details*).
#'        }
#' @return A list of class \code{MSGARCH_MCMC_FIT} with the following elements:
#' \itemize{
#' \item \code{par}: The MCMC chain (matrix from the \R package
#' \code{coda} (Plummer et al., 2006) of size \code{N.mcmc} / \code{N.thin} x d).
#' \item \code{accept}: Acceptation rate of the sampler.
#' \item \code{spec}:  Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}}.
#' \item \code{data}:  Vector (of size T) of observations.
#' \item \code{ctr}: \code{list} of the control used for the fit.
#' }
#' The \code{MSGARCH_MCMC_FIT} with the following methods:
#' \itemize{
#' \item \code{\link{AIC}}: Compute Akaike information criterion (AIC).
#' \item \code{\link{BIC}}: Compute Bayesian information criterion (BIC).
#' \item \code{\link{DIC}}: Compute Deviance Information Criterion (DIC).
#' \item \code{\link{Volatility}}: In-sample conditional volatility filterting of the overall process.
#' \item \code{\link{Forecast}}: Forecast of the conditional volatility of the overall process.
#' \item \code{\link{UncVol}}: Unconditional volatility in each regime and the overall process.
#' \item \code{\link{Pred}}: Predictive method.
#' \item \code{\link{PIT}}: Probability integral transform.
#' \item \code{\link{Risk}}: Value-at-Risk And Expected-Shortfall methods.
#' \item \code{\link{Sim}}: Simulation method.
#' \item \code{\link{State}}: State probabilities methods.
#' \item \code{\link{ExtractStateFit}}: Single-regime model extractor.
#' \item \code{summary}: Summary of the fit.
#' }
#' @details The total number of draws is equal to \code{n.mcmc / n.thin}.
#' The MCMC/Bayesian estimation relies on an \pkg{Rcpp} implementation of the adaptive sampler of Vihola (2012). 
#' The implementation is based on the R package \pkg{adaptMCMC} (Andreas, 2012).
#' Starting values when \code{par0} is not provided are chosen automatically
#' before sampling (see Ardia et al. (2016) for more details).\cr
#' \code{SamplerFUN} allows for a custom sampler to be used. The function
#' must take the form: \cr \code{function(f_posterior, data, spec, par0, ctr)}, \cr
#' where  \code{f_posterior} is the function to optimize, \code{data} is
#' the data, \code{spec} is the specification,
#' \code{par0} are the starting parameters, and \code{ctr} are the control
#' parameters. The inputs \code{spec} and \code{data},
#' must be passed as inputs in the sampler (see *Examples*).
#' The custom sampler must output a matrix containing the MCMC chain. \cr
#' @references Andreas, S. (2012).
#' \code{adaptMCMC}: Implementation of a generic adaptive Monte Carlo Markov chain sampler.
#' \url{https://cran.r-project.org/package=adaptMCMC}
#' @references Ardia, D. Bluteau, K. Boudt, K. Catania, L. & Trottier, D.-A. (2016).
#' Markov-switching GARCH models in \R: The MSGARCH package.
#' \url{https://ssrn.com/abstract=2845809}
#' @references MacDonald, I.L., Zucchini, W. (1997).
#' Hidden Markov and other models for discrete-valued time series.
#' \emph{CRC press}.
#' @references Plummer, M. Best, N. Cowles, K. & Vines, K. (2006).
#' \code{coda}: Convergence diagnosis and output analysis for MCMC.
#' \emph{R News}, 6, 7-11.
#' \url{https://cran.r-project.org/package=coda}
#' @references Vihola, M. (2012).
#' Robust adaptive Metropolis algorithm with coerced acceptance rate.
#' \emph{Statistics and Computing}, 22, 997-1008.
#' @examples
#' # load data
#' data("SMI", package = "MSGARCH")
#'
#' # create model specification
#' # MS(2)-GARCH(1,1)-Normal (default)
#' spec <- CreateSpec()
#'
#' # fit the model on the data by MCMC
#' set.seed(123)
#' fit <- FitMCMC(spec = spec, data = SMI, ctr = list(n.burn = 500L, n.mcmc = 500L, n.thin = 1L))
#' summary(fit)
#'
#' # custom sampler example
#' \dontrun{
#' library("mcmc")
#' f_MCMC <- function(f_posterior, data, spec, par0, ctr){
#'   par <- mcmc::metrop(f_posterior, initial = par0, nbatch = ctr$n.mcmc + ctr$n.burn,
#'                         data = data, spec = spec)$batch
#'   return(par)
#' }
#'
#' set.seed(123)
#' fit <- FitMCMC(spec, data = SMI, ctr  = list(SamplerFUN = f_MCMC,
#'                                              n.burn = 500L, n.mcmc = 500L, n.thin = 1L))
#' summary(fit)
#' }
#' @import coda
#' @export
FitMCMC <- function(spec, data, ctr = list()) {
  UseMethod(generic = "FitMCMC", spec)
}

#' @method FitMCMC MSGARCH_SPEC
#' @export
FitMCMC.MSGARCH_SPEC <- function(spec, data, ctr = list()) {

  time.start <- Sys.time()
  spec <- f_check_spec(spec)
  data <- f_check_y(data)
  ctr  <- f_process_ctr(ctr)
  ctr$do.plm <- TRUE

  if (isTRUE(spec$fixed.pars.bool)) {
    f_check_fixedpars(spec$fixed.pars, spec)
  }

  if (is.null(ctr$par0)) {
    par0 <- f_StargingValues(data, spec, ctr)
  } else {
    par0 = ctr$par0
    if (isTRUE(spec$regime.const.pars.bool)) {
      par0 <- f_remove_regimeconstpar(par0, spec$regime.const.pars, spec$K)
    }
    if (isTRUE(spec$fixed.pars.bool)) {
      par0 <- f_substitute_fixedpar(par0, spec$fixed.pars)
    }
    par0 <- f_unmapPar(par0, spec, do.plm = TRUE)
  }
  par    <- ctr$SamplerFUN(f_posterior = f_posterior, data = data, spec = spec, par0 = par0, ctr = ctr)
  np     <- length(par0)
  accept <- 1 - mean(duplicated(par))
  by     <- seq(from = (ctr$n.burn + 1L), to = (ctr$n.burn + ctr$n.mcmc), by = ctr$n.thin)

  par    <- par[by, , drop = FALSE]
  vLower <- spec$lower
  vUpper <- spec$upper

  names(vLower) <- spec$label
  names(vUpper) <- spec$label

  vLower <- vLower[names(par0)]
  vUpper <- vUpper[names(par0)]

  par  <- f_map(par, vLower, vUpper)
  par0 <- f_map(par0, vLower, vUpper)

  vpar_full <- par[1L, ]

  if (isTRUE(spec$fixed.pars.bool)) {
    vpar_full <- f_add_fixedpar(vpar_full, spec$fixed.pars)
    vpar_full <- vpar_full[spec$label]
    names(vpar_full) <- spec$label
  }

  if (isTRUE(spec$regime.const.pars.bool)) {
    vpar_full <- f_add_regimeconstpar(vpar_full, spec$K, spec$label)
  }

  if (length(vpar_full) != ncol(par)) {

    vAddedPar <- vpar_full[!names(vpar_full) %in% colnames(par)]

    if (isTRUE(spec$fixed.pars.bool)) {
      par <- cbind(par, matrix(rep(vAddedPar, nrow(par)), nrow = nrow(par),
                               byrow = TRUE, dimnames = list(NULL, names(vAddedPar))))
      par <- par[, spec$label]
    }

    if (isTRUE(spec$regime.const.pars.bool)) {
      par <- f_add_regimeconstpar_matrix(par, spec$K, spec$label)
    }
  }

  par <- f_sort_par(spec, par)
  par <- coda::mcmc(par)
  ctr$par0 <- par0
  elapsed.time <- Sys.time() - time.start

  out <- list(par = par, accept = accept, data = data, spec = spec, ctr = ctr)
  class(out) <- "MSGARCH_MCMC_FIT"
  return(out)
}
