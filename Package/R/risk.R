#' Value-at-Risk And Expected-shortfall functions at T + 1.
#' @description Method returning the Value-at-Risk and Expected-shortfall at T + 1.
#' @param spec Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates.
#' @param y  Vector (of size T) of observations.
#' @param level Vector (of size R) of Value-at-risk and Expected-shortfall levels.\cr
#'  (Default: \code{level = c(0.95,0.99)})
#' @param ES  Boolean indicating if Expected-shortfall is also calculated. (Default: \code{ES = TRUE})
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian Value-at-Risk and Expected-shortfall are calculated.
#' @return A list containing one or two components:
#' \itemize{
#' \item \code{VaR} : Value-at-Risk at T + 1 at the choosen levels (vector of size R).
#' \item \code{ES}  : Expected-shortfall at T + 1 at the choosen levels (vector of size R).
#' }
#' @examples 
#'# load data
#'data("sp500ret")
#'
#'# create model specification
#'spec = MSGARCH::create.spec(model = c("sGARCH","sGARCH"), distribution = c("norm","norm"),
#'                              do.skew = c(FALSE,FALSE), do.mix = FALSE, do.shape.ind = FALSE) 
#'set.seed(123)
#'                              
#'# compute the Value-at-Risk and Expected-shortfall 
#'risk = MSGARCH::risk(spec = spec, theta = spec$theta0, y = sp500ret,
#'                     level = c(0.95,0.99), ES = TRUE)
#' @importFrom stats integrate sd uniroot                    
#' @export
risk <- function(spec, theta, y, level = c(0.95,0.99), ES = TRUE)
{
  UseMethod("risk", spec)
}

#' @export
risk.MSGARCH_SPEC = function(spec, theta, y, level = c(0.95,0.99), ES = TRUE) {
  
  y = as.matrix(y)
  if (isTRUE(spec$is.shape.ind)) {
    theta = spec$func$f.do.shape.ind(theta)
  }
  
  if (isTRUE(spec$is.mix)) {
    theta = spec$func$f.do.mix(theta)
  }
  
  if (is.vector(theta)) {
    theta = matrix(theta, nrow = 1)
  }
  
  p = 1 - level
  n = length(y)
  N = nrow(theta)
  np = length(p)
  xmin = min(y) - sd(y)
  xmax = 0
  itermax = 100
  tol = 1e-05
  
  tmp.VaR = NULL
  
  f.pdf = function(x) {
    out = MSGARCH::pred(spec, x, theta, y, log = FALSE)
    return(out)
  }
  
  f.fun = function(x, pi) {
    out = integrate(f.pdf, lower = xmin, upper = x)$value - pi
    return(out)
  }
  
  # gross approximation for VaR
  out = list()
  for (i in 1:np) {
    tmp.VaR[i] = uniroot(f.fun, lower = xmin, upper = xmax, pi = p[i])$root
  }
  
  # add precision by Newton-Raphson
  out$VaR = vector("double", np)
  for (i in 1:np) {
    p_i = p[i]
    level_i = level[i]
    calc.step = function(V) {
      PDF = MSGARCH::pred(spec, x = V, theta = theta, y = y, log = F)
      CDF = MSGARCH::pit(spec, x = V, theta = theta, y = y)
      lPDF = log(PDF)
      err = p_i - CDF
      step = err * exp(-lPDF)
      return(list(step = step, err = abs(err)))
    }
    # Starting values
    VaR = tmp.VaR[i]
    # VaR calculation
    newStep = calc.step(VaR)
    delta = newStep$step
    VaR = VaR + delta
    covERR = newStep$err
    for (j in 1:itermax) {
      if (covERR < tol) {
        break
      }
      newStep = calc.step(VaR)
      delta = newStep$step
      VaR = VaR + delta
      covERR = newStep$err
    }
    out$VaR[i] = VaR
  }
  
  if (isTRUE(ES)) {
    out$ES = vector("double", np)
    for (i in 1:np) {
      
      f.condMean = function(x) {
        out = x * MSGARCH::pred(spec, x, theta, y, log = FALSE)
        return(out)
      }
      out$ES[i] = integrate(f.condMean, lower = -Inf, upper = out$VaR[i], stop.on.error = FALSE)$value/p[i]
    }
  }
  return(out)
}