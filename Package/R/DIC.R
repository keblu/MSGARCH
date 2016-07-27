#' Compute Deviance Information Criterion (DIC).
#' @param fit Fit object of type \code{\link{MSGARCH_BAY_FIT}} created with \code{\link{fit.bayes}}.
#' @references Gelman, A. Carlin, J. B. Stern, H. S. & Rubin, D. B. (2003). Bayesian Data Analysis. \emph{Chapman and Hall/CRC}
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'                              
#'fit = MSGARCH::fit.bay(spec = spec, y = sp500ret)
#'
#'DIC = MSGARCH::DIC(fit = fit)
#' @return A list containing six variables:
#'        \itemize{
#'        \item \code{DIC} : Deviance Information Criterion.
#'        \item \code{IC}  : Bayesian Predictive Information Criterion.
#'        \item \code{pD}  : Effective number of parameters (pD = Dbar - Dhat)
#'        \item \code{pV}  : Effective number of parameters (pV = var(D)/2)
#'        \item \code{D.bar}: Expected value of the deviance over the posterior
#'        \item \code{D.hat}: Deviance at the mean posterior estimate
#'        }
#' @export
DIC <- function(fit)
{
  UseMethod("DIC", fit)
}

#' @export
DIC.MSGARCH_BAY_FIT <- function(fit){
  
  DIC = f.DIC(fit$spec, fit$theta, fit$y)
  return(DIC)
}

f.DIC = function(spec, theta, y) {
  
  if (is.vector(theta)) {
    theta = matrix(theta, nr = 1)
  }
  
  LL = MSGARCH::kernel(spec, theta, y = y, log = TRUE)
  D.bar = -2*mean(LL)
  theta.bar <- colMeans(theta)
  D.hat = -2*MSGARCH::kernel(spec, theta.bar, y = y, log = TRUE)
  pD <- D.bar - D.hat
  pV <- var(-2*LL)/2
  out = list(DIC= pD + D.bar,IC= 2*pD + D.bar, pD = pD, pV = pV, D.bar = D.bar, D.hat = D.hat)
  return(out)
}