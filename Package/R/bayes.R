#' Bayesian estimation.
#' @description Method that performs Bayesian estimation of a \code{\link{MSGARCH_SPEC}} object on a set of observations.
#' @param spec Model specification of class \code{\link{MSGARCH_SPEC}} created with \code{\link{create.spec}}.
#' @param y  Vector (of size T) of observations.
#' @param ctr  A list of control parameters. \cr
#'        The control parameters have three components:
#'        \itemize{
#'        \item \code{N.burn} (integer >= 0): Number of discarded draws. (default: \code{N.burn = 1000})
#'        \item \code{N.mcmc} (integer > 0) : Number of draws. (default: \code{N.mcmc = 5000})
#'        \item \code{N.thin} (integer > 0) : Thinning factor (every \code{N.thin} draws are kept). (default: \code{N.thin = 10})
#'        }
#' @return  A list of class \code{\link{MSGARCH_BAY_FIT}} containing four variables:
#' \itemize{
#' \item \code{theta} : The Bayesian chain (matrix from the \R package \code{coda} (Plummer et al., 2006) of size \code{N.mcmc / N.thin}).
#' \item \code{accept} : Acceptation rate of the sampler.
#' \item \code{y} : Initial vector (of size T) of observations.
#' \item \code{spec} : Initial specification.
#' }
#' The \code{\link{MSGARCH_BAY_FIT}} possess this method:
#' \itemize{
#' \item \code{\link{AIC}} : Compute Akaike information criterion (AIC).
#' \item \code{\link{BIC}} : Compute Bayesian information criterion (BIC).
#' \item \code{\link{DIC}} : Compute Deviance Information Criterion (DIC).
#' }
#' @details The total number of draws is equal to \code{N.burn + N.mcmc}.
#' The Bayesian estimation uses the \R package \code{adaptMCMC} (Andreas, 2012) which  
#' implements the adaptive sampler of Vihola (2012).
#' @examples 
#' data("sp500ret")
#' 
#' spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                               do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#' 
#' fit = MSGARCH::fit.bayes(spec = spec, y = sp500ret, 
#'                          l.ctr = list(N.burn = 100,N.mcmc = 1000, N.thin = 1))
#'                          
#' @references Andreas, S. (2012). \code{adaptMCMC}: Implementation of a Generic Adaptive Monte Carlo Markov Chain Sampler. \url{https://cran.r-project.org/web/packages/adaptMCMC/}.
#' @references Metropolis, N.; Rosenbluth, A. W.; Rosenbluth, M. N.; Teller, A. H. & Teller, E. (1953). Equation of State Calculations by Fast Computing Machines. \emph{Journal of Chemical Physics}, 21, pp. 1087-1092.
#' @references Plummer, M. Best, N. Cowles, K. & Vines, K. (2006). \code{CODA}: Convergence Diagnosis and Output Analysis for MCMC. \emph{R News}, 6, pp.7-11. \url{https://cran.r-project.org/web/packages/coda/}.
#' @references Vihola, M. (2012). Robust Adaptive Metropolis Algorithm with Coerced Acceptance Rate. \emph{Statistics and Computing}, 22, pp. 997-1008.
#' @import adaptMCMC
#' @export
fit.bayes <- function(spec, y, l.ctr)
{
  UseMethod("fit.bayes", spec)
}

#' @method fit.bayes MSGARCH_SPEC
#' @export
fit.bayes.MSGARCH_SPEC = function(spec, y, ctr) {
  
  l.ctr = f.process.ctr(ctr)
  # For Identification problem we make sure that the first
  # models of each type of spec in the MSGARCH list of model
  # have the lowest uncondtional volatility
  f.kernel = function(x, log = TRUE) {
    name = spec$name
    unique.spec = unique(name, FALSE)
 
    for (i in 1:length(unique.spec)) {
      idx = name == unique.spec[i]
      unc.vol = MSGARCH::unc.vol(spec = spec,theta = x)
      unc.vol = unc.vol[idx]
      
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
    out = MSGARCH::kernel(spec = spec, theta = x, y = y, log = log)
    return(out)
  }
  
  ## ==> MCMC estimation
  outMH = adaptMCMC::MCMC(p = f.kernel, n = l.ctr$N.burn + 
      l.ctr$N.mcmc, init = spec$theta0, adapt = TRUE, acc.rate = 0.4)
  
  by = seq(from = (l.ctr$N.burn + 1), to = (l.ctr$N.burn + l.ctr$N.mcmc), by = l.ctr$N.thin)
  nby = length(by)
  outMH$samples = outMH$samples[by, , drop = FALSE]
  accept = 1 - mean(duplicated(outMH$samples))
  theta = adaptMCMC::convert.to.coda(outMH)
  colnames(theta) = spec$label
  
  l.out = list(theta = theta, accept = accept, y = y, spec = spec)
  class(l.out) = "MSGARCH_BAY_FIT"
  return(l.out)
}


