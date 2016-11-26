#' Bayesian estimation.
#' @description Method that performs Bayesian estimation of a \code{MSGARCH_SPEC} object on a set of observations.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param y  Vector (of size T) of observations.
#' @param ctr  A list of control parameters. \cr
#'        The control parameters have three components:
#'        \itemize{
#'        \item \code{N.burn} (integer >= 0): Number of discarded draws. (Default: \code{N.burn = 5000})
#'        \item \code{N.mcmc} (integer > 0) : Number of draws. (Default: \code{N.mcmc = 10000})
#'        \item \code{N.thin} (integer > 0) : Thinning factor (every \code{N.thin} draws are kept). (Default: \code{N.thin = 10})
#'        \item \code{theta0} : Starting value for the chain (if empty the specification default value are used).
#'        \item \code{do.enhance.theta0} : Boolean indicating if the default parameters value are enhance using \code{y} variance. (Default: \code{do.enhance.theta0 = FALSE})
#'        }
#' @return  A list of class \code{MSGARCH_BAY_FIT} containing four components:
#' \itemize{
#' \item \code{theta} : The MCMC chain (matrix from the \R package \code{coda} (Plummer et al., 2006) of size \code{N.mcmc / N.thin}).
#' \item \code{accept} : Acceptation rate of the sampler.
#' \item \code{y} :  Vector (of size T) of observations.
#' \item \code{spec} :  Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' }
#' The \code{MSGARCH_BAY_FIT} contains these methods:
#' \itemize{
#' \item \code{\link{AIC}} : Compute Akaike information criterion (AIC).
#' \item \code{\link{BIC}} : Compute Bayesian information criterion (BIC).
#' \item \code{\link{DIC}} : Compute Deviance Information Criterion (DIC).
#' \item \code{\link{ht}}  : Conditional volatility in each regime.
#' \item \code{\link{kernel}} : Kernel method.
#' \item \code{\link{unc.vol}} : Unconditional volatility in each regime.
#' \item \code{\link{pred}} : Predictive method.
#' \item \code{\link{pit}} : Probability Integral Transform.
#' \item \code{\link{risk}} : Value-at-Risk And Expected-Shortfall methods.
#' \item \code{\link{simahead}} : Step ahead simulation method.
#' \item \code{\link{sim}} : Simulation method.
#' \item \code{\link{pdf}} : Probability density function.
#' \item \code{\link{cdf}} : Cumulative function.
#' \item \code{\link{Pstate}} : State probabilities filtering method.
#' \item \code{summary} : Summary of the fit.
#' }
#' @details The total number of draws is equal to \code{N.mcmc / N.thin}.
#' The Bayesian estimation uses the \R package \code{adaptMCMC} (Andreas, 2012) which  
#' implements the adaptive sampler of Vihola (2012). 
#'  The starting parameters are the specification default parameters.
#' The argument \code{do.enhance.theta0}
#'  uses the volatilities of rolling windows of \code{y} and adjust the default parameter of
#'  the specification so that the unconditional volatility of each regime
#'  is set to different quantiles of the volatilities of the rolling windows of \code{y}.
#' @examples 
#' # load data
#' data("sp500")
#' sp500 = sp500[1:1000]
#' 
#' # create model specification
#' spec = MSGARCH::create.spec() 
#' 
#' # fit the model on the data with Bayesian estimation
#' set.seed(123)
#' fit = MSGARCH::fit.bayes(spec = spec, y = sp500, 
#'                          ctr = list(N.burn = 500, N.mcmc = 1000, N.thin = 1))
#'                          
#'summary(fit)
#' @references Andreas, S. (2012). \code{adaptMCMC}: Implementation of a Generic Adaptive Monte Carlo Markov Chain Sampler. \url{https://cran.r-project.org/package=adaptMCMC}.
#' @references Metropolis, N.; Rosenbluth, A. W.; Rosenbluth, M. N.; Teller, A. H. & Teller, E. (1953). Equation of State Calculations by Fast Computing Machines. \emph{Journal of Chemical Physics}, 21, pp. 1087-1092.
#' @references Plummer, M. Best, N. Cowles, K. & Vines, K. (2006). \code{CODA}: Convergence Diagnosis and Output Analysis for MCMC. \emph{R News}, 6, pp.7-11. \url{https://cran.r-project.org/package=coda}.
#' @references Vihola, M. (2012). Robust Adaptive Metropolis Algorithm with Coerced Acceptance Rate. \emph{Statistics and Computing}, 22, pp. 997-1008.
#' @import adaptMCMC
#' @export
fit.bayes <- function(spec, y, ctr = list()) {
  UseMethod("fit.bayes", spec)
}

#' @method fit.bayes MSGARCH_SPEC
#' @export
fit.bayes.MSGARCH_SPEC <- function(spec, y, ctr = list()) {
  y <- f.check.y(y)
  l.ctr <- f.process.ctr(ctr)
  if (is.null(l.ctr$theta0)) {
    l.ctr$theta0 <- spec$theta0
  }
  if (isTRUE(ctr$do.enhance.theta0)) {
    l.ctr$theta0 <- f.enhance.theta(spec = spec, theta = l.ctr$theta0, y = y)
  }
  
  f.kernel <- function(x, log = TRUE) {
    name <- spec$name
    unique.spec <- unique(name, FALSE)
    for (i in 1:length(unique.spec)) {
      idx <- name == unique.spec[i]
      options(warn = -1)
      unc.vol <- MSGARCH::unc.vol(spec, x)
      options(warn = 0)
      unc.vol <- unc.vol[idx]
      if (any(is.na(unc.vol))) {
        return(-1e+10)
      } else {
        if (length(unc.vol) != 1) {
          for (j in 1:(length(unc.vol) - 1)) {
            if (unc.vol[j] > unc.vol[j + 1]) {
              return(-1e+10)
            }
          }
        }
      }
    }
    out <- MSGARCH::kernel(spec, x, y, log = log)
    return(out)
  }
  outmh <- adaptMCMC::MCMC(p = f.kernel, n = l.ctr$N.burn + l.ctr$N.mcmc,
                           init = l.ctr$theta0, adapt = TRUE, acc.rate = 0.4)
  by <- seq(from = (l.ctr$N.burn + 1), to = (l.ctr$N.burn + l.ctr$N.mcmc),
            by = l.ctr$N.thin)
  outmh$samples <- outmh$samples[by, , drop = FALSE]
  accept <- 1 - mean(duplicated(outmh$samples))
  theta <- adaptMCMC::convert.to.coda(outmh)
  colnames(theta) <- spec$label
  l.out <- list(theta = theta, accept = accept, y = y, spec = spec)
  class(l.out) <- "MSGARCH_BAY_FIT"
  return(l.out)
}