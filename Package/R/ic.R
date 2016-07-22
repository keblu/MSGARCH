
#' Compute Bayesian information criterion (BIC).
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param spec Model specification created with \code{\link{f.create.spec}}.
#' @references Schwarz, G. (1978). Estimating the dimension of a model. \emph{Annals of Statistics}, 6, pp. 461-464. 
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'bic = spec$f.bic(theta = spec$theta0, y = sp500ret, spec = spec)
#' @return BIC values (scalar or vector of size M).
#' @export
f.bic = function(theta, y, spec) {
  
  if (is.vector(theta)) {
    theta = matrix(theta, nr = 1)
  }
  LL = spec$f.kernel(theta, y = y, log = TRUE)
  k = dim(theta)[2]
  n = length(y)
  bic = 2 * LL - k * log(n)
  return(bic)
  
}

#' Compute Akaike information criterion (AIC).
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param spec Model specification created with \code{\link{f.create.spec}}.
#' @return AIC values (scalar or vector of size M).
#' @examples 
#'data("sp500ret")
#'
#'spec = MSGARCH::f.create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'
#'aic = spec$f.aic(theta = spec$theta0, y = sp500ret, spec = spec)
#' @references Akaike, H. (1974). A New Look at the Statistical Model Identification. \emph{IEEE Transactions on Automatic Control}, 19, pp. 716-723.
#' @export
f.aic = function(theta, y, spec) {
  
  if (is.vector(theta)) {
    theta = matrix(theta, nr = 1)
  }
  LL = spec$f.kernel(theta, y = y, log = TRUE)
  k = dim(theta)[2]
  aic = 2 * LL - 2 * k
  return(aic)
  
}