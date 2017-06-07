#' Compute Deviance Information Criterion (DIC).
#' @param fit Fit object of type \code{MSGARCH_BAY_FIT} created with \code{\link{fit.bayes}}.
#' @references Spiegelhalter, David J., et al. (2002). Bayesian measures of model complexity and fit. \emph{Journal of the Royal Statistical Society: Series B (Statistical Methodology)}
#' @examples 
#'require("MSGARCH")
#'# load data
#'data("sp500")
#'sp500 = sp500[1:1000]
#'
#'# create model specification
#'spec = MSGARCH::create.spec()
#'
#'# fit the model by Bayesian estimation 
#'set.seed(123)                                                           
#'fit = MSGARCH::fit.bayes(spec = spec, y = sp500)
#'
#'# compute DIC
#'DIC = MSGARCH::DIC(fit)
#' @return A list containing four variables:
#'        \itemize{
#'        \item \code{DIC} : Deviance Information Criterion.
#'        \item \code{IC}  : Bayesian Predictive Information Criterion (IC = 2 * pV + D.bar).
#'        \item \code{pV}  : Effective number of parameters (pV = var(D)/2)
#'        \item \code{D.bar}: Expected value of the deviance over the posterior
#'        }
#' @details Compute the Deviance information criterion of Spiegelhalter, David J., et al. (2002). We define the deviance as: \deqn{D(\theta) = -2LLH(\mathbf{y}|\theta),} where \eqn{\mathbf{y}} are the data, \eqn{\theta} 
#'  are the parameters, and LLH() is the log-likelihood function.
#'  The expectation \deqn{\bar{D} = {\mathbf{E}}^{\theta}[D(\theta)],} where \eqn{{\mathbf{E}}^{\theta}} is the expectation over all theta in a MCMC chain,
#'  is a measure of how well the model fits the data. The larger this expectation is, the worse is the fit. The effective number of parameters of the model can be defined as
#' \deqn{p_{V} = {\frac{1}{2}}\widehat{var}\left(D(\theta)\right),} where \eqn{\widehat{var}} is the the population variance estimator.
#'   The larger the effective number of parameters is, the easier it is for the model to fit the data, and so the deviance needs to be penalized.
#'  Finally DIC is defined as: \deqn{\mathit{DIC} = p_{V}+\bar{D}.}
#' @importFrom stats var
#' @export
DIC <- function(fit) {
  UseMethod(generic = "DIC", object = fit)
}

#' @export
DIC.MSGARCH_BAY_FIT <- function(fit) {
  DIC <- f.DIC(spec = fit$spec, theta = fit$theta, y = fit$y)
  return(DIC)
}

f.DIC <- function(spec, theta, y) {
  if (is.vector(x = theta)) {
    theta <- matrix(data = theta, nrow = 1)
  }
  LL = vector(mode = "numeric",length = nrow(theta))
  for(i in 1:nrow(theta)){
      LL[i]  <- sum(MSGARCH:::pred(object = spec, theta = theta[i,], y = y,
                                 log = TRUE, do.its = TRUE)$pred)
  }
  D.bar     <- -2 * mean(x = LL)
  pV        <- var(x = -2 * LL) / 2
  out       <- list(DIC = pV + D.bar, IC = 2 * pV + D.bar, pV = pV, D.bar = D.bar)
  return(out)
}