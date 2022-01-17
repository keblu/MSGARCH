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
#'        \item \code{nburn} (integer >= 0): Number of discarded draws.
#'        (Default: \code{nburn = 5000L})
#'        \item \code{nmcmc} (integer > 0): Number of draws.
#'        (Default: \code{nmcmc = 10000L})
#'        \item \code{nthin} (integer > 0): Thinning factor (every \code{nthin}
#'        draws are kept). (Default: \code{nthin = 10L})
#'        \item \code{do.sort} (bool): Logical indicating if the MCMC draws are post-sorted
#'         following Geweke (2007). By default, \code{do.sort = TRUE}, such that the
#'         MCMC draws are ordered to ensure that unconditional variance is an 
#'         increasing function of the regime (identification constraint). If the user sets
#'          \code{do.sort = FALSE}, no sorting is imposed, and label switching can occur (see *Details*).
#'        \item \code{SamplerFUN}: Custom MCMC sampler (see *Details*).
#'        }
#' @return A list of class \code{MSGARCH_MCMC_FIT} with the following elements:
#' \itemize{
#' \item \code{par}: The MCMC chain (matrix from the \R package
#' \code{coda} (Plummer et al., 2006) of size \code{nmcmc} / \code{nthin} x d).
#' \item \code{accept}: Acceptance rate of the sampler.
#' \item \code{spec}:  Model specification of class \code{MSGARCH_SPEC}
#' created with \code{\link{CreateSpec}}.
#' \item \code{data}:  Vector (of size T) of observations.
#' \item \code{ctr}: \code{list} of the control used for the fit.
#' }
#' The \code{MSGARCH_MCMC_FIT} with the following methods:
#' \itemize{
#' \item \code{\link{DIC}}: Deviance Information Criterion (DIC).
#' \item \code{simulate}: Simulation.
#' \item \code{\link{Volatility}}: In-sample conditional volatility.
#' \item \code{predict}: Forecast of the conditional volatility (and predictive distribution).
#' \item \code{\link{UncVol}}: Unconditional volatility.
#' \item \code{\link{PredPdf}}: Predictive density (pdf).
#' \item \code{\link{PIT}}: Probability Integral Transform.
#' \item \code{\link{Risk}}: Value-at-Risk and Expected-Shortfall.
#' \item \code{\link{State}}: State probabilities (smoothed, filtered, predictive, Viterbi).
#' \item \code{\link{ExtractStateFit}}: Single-regime model extractor.
#' \item \code{summary}: Summary of the fit.
#' }
#' @details The total number of draws is equal to \code{nmcmc / nthin}.
#' The MCMC/Bayesian estimation relies on an \pkg{Rcpp} implementation of the adaptive sampler of Vihola (2012). 
#' The implementation is based on the R package \pkg{adaptMCMC} (Andreas, 2012).
#' Starting values when \code{par0} is not provided are chosen automatically
#' before sampling (see Ardia et al. (2019) for more details).\cr
#' \code{SamplerFUN} allows for a custom sampler to be used. The function
#' must take the form: \cr \code{function(f_posterior, data, spec, par0, ctr)}, \cr
#' where  \code{f_posterior} is the function to optimize, \code{data} is
#' the data, \code{spec} is the specification,
#' \code{par0} are the starting parameters, and \code{ctr} are the control
#' parameters. The inputs \code{spec} and \code{data},
#' must be passed as inputs in the sampler (see *Examples*).
#' The custom sampler must output a matrix containing the MCMC chain. \cr
#' When \code{do.sort = TRUE}, sorting of each MCMC draw conditional on the unconditional variance is done across homogeneous regime specification.\cr
#' 
#' @references Andreas, S. (2012).
#' \code{adaptMCMC}: Implementation of a generic adaptive Monte Carlo Markov chain sampler.
#' \url{https://cran.r-project.org/package=adaptMCMC}
#' 
#' @references Ardia, D. Bluteau, K. Boudt, K. Catania, L. Trottier, D.-A. (2019).
#' Markov-switching GARCH models in \R: The \pkg{MSGARCH} package.
#' \emph{Journal of Statistical Software}, 91(4), 1-38.
#' \doi{10.18637/jss.v091.i04}
#' 
#' @references Geweke J (2007).
#' Interpretation and Inference in Mixture Models: Simple MCMC Works.
#' \emph{Computational Statistics & Data Analysis}, 51(7), 3529-3550.
#' \doi{10.1016/j.csda.2006.11.026}
#' 
#' @references MacDonald, I.L., Zucchini, W. (1997).
#' \emph{Hidden Markov and other models for discrete-valued time series}.
#' CRC press.
#' 
#' @references Plummer, M. Best, N. Cowles, K. & Vines, K. (2006).
#' \code{coda}: Convergence diagnosis and output analysis for MCMC.
#' \emph{R News}, 6, 7-11.
#' \url{https://cran.r-project.org/package=coda}
#' 
#' @references Vihola, M. (2012).
#' Robust adaptive Metropolis algorithm with coerced acceptance rate.
#' \emph{Statistics and Computing}, 22, 997-1008.
#' \doi{10.1007/s11222-011-9269-5}
#' 
#' @examples
#' # create model specification
#' spec <- CreateSpec()
#' 
#' # load data
#' data("SMI", package = "MSGARCH")
#'
#' # fit the model on the data by MCMC
#' set.seed(123)
#' fit <- FitMCMC(spec = spec, data = SMI, ctr = list(nburn = 500L, nmcmc = 500L, nthin = 1L))
#' summary(fit)
#'
#' # custom sampler example
#' \dontrun{
#' library("mcmc")
#' f_MCMC <- function(f_posterior, data, spec, par0, ctr){
#'   par <- mcmc::metrop(f_posterior, initial = par0, nbatch = ctr$nmcmc + ctr$nburn,
#'                         data = data, spec = spec)$batch
#'   colnames(par) = names(par0)
#'   return(par)
#' }
#'
#' set.seed(123)
#' fit <- FitMCMC(spec, data = SMI, ctr  = list(SamplerFUN = f_MCMC,
#'                                              nburn = 500L, nmcmc = 500L, nthin = 1L))
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
  data_ <- f_check_y(data)
  ctr  <- f_process_ctr(ctr)
  ctr$do.plm <- TRUE

  if (isTRUE(spec$fixed.pars.bool)) {
    f_check_fixedpars(spec$fixed.pars, spec)
  }

  if (is.null(ctr$par0)) {
    par0 <- f_StargingValues(data_, spec, ctr)
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
  par    <- ctr$SamplerFUN(f_posterior = f_posterior, data = data_, spec = spec, par0 = par0, ctr = ctr)
  np     <- length(par0)
  accept <- 1 - mean(duplicated(par))
  by     <- seq(from = (ctr$nburn + 1L), to = (ctr$nburn + ctr$nmcmc), by = ctr$nthin)

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
  if(isTRUE(ctr$do.sort)){
    par <- f_sort_par(spec, par)
  }
  par <- coda::mcmc(par)
  ctr$par0 <- par0
  elapsed.time <- Sys.time() - time.start

  out <- list(par = par, accept = accept, data = data, spec = spec, ctr = ctr)
  class(out) <- "MSGARCH_MCMC_FIT"
  return(out)
}
