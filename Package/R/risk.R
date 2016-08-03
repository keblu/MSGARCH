#' Value-at-Risk And Expected-shortfall.
#' @description Method returning the Value-at-Risk and Expected-shortfall in-sample or at \code{t = T + 1} based on the predictive density.
#'@param object Model specification of class \code{MSGARCH_SPEC} created with \code{\link{create.spec}}
#' or fit object of type \code{MSGARCH_MLE_FIT} created with \code{\link{fit.mle}} or \code{MSGARCH_BAY_FIT}
#' created with \code{\link{fit.bayes}}.
#' @param theta Vector (of size d) or matrix (of size M x d) of parameter estimates (not require when using a fit object).
#' @param y  Vector (of size T) of observations (not require when using a fit object).
#' @param level Vector (of size R) of Value-at-risk and Expected-shortfall levels.\cr
#'  (Default: \code{level = c(0.95,0.99)})
#' @param ES  Boolean indicating if Expected-shortfall is also calculated. (Default: \code{ES = TRUE})
#' @param is.its  Boolean indicating if the in-sample risk estimator are returned. (Default: \code{is.its = FALSE})
#' @details If a matrix of MCMC posterior draws estimates is given, the Bayesian Value-at-Risk and Expected-shortfall are calculated.
#' If \code{is.its = FALSE}, \code{x} the risk estimator at \code{t = T + 1}, the method uses the variance estimated at \code{t = T + 1}.
#' If \code{is.its = TRUE}, The in-sample risk estimator are calculated.
#' @return A list containing of class \code{MSGARCH_RISK} containing two or three components:
#' \itemize{
#' \item \code{VaR} : \cr If \code{is.its = FALSE}: Value-at-Risk at \code{t = T + 1} at the choosen levels (vector of size R).\cr
#'                    If \code{is.its = TRUE}: In-sample Value-at-Risk at the choosen levels (Matrix of size T x R).
#' \item \code{ES}  :\cr If \code{is.its = FALSE}: Expected-shortfall at \code{t = T + 1} at the choosen levels (vector of size R).\cr
#'                    If \code{is.its = TRUE}: In-sample Expected-shortfall at the choosen levels (Matrix of size T x R).
#' \item \code{y}  :\cr Vector (of size T) of observations.
#' }
#' @examples 
#'\dontrun{
#'# load data
#'data("sp500ret")
#'
#'# create model specification
#'spec = MSGARCH::create.spec() 
#'
#'# fit the model on the data with ML estimation using DEoptim intialization
#' set.seed(123)
#'fit = MSGARCH::fit.mle(spec = spec, y = sp500ret)
#' 
#'# compute the Value-at-Risk and Expected-shortfall 
#'# Risk estimation in-sample 
#'risk.its = MSGARCH::risk(object = fit, level = c(0.95,0.99), ES = TRUE, is.its = TRUE)
#'                     
#'# Risk estimation at T + 1                     
#'risk = MSGARCH::risk(object = fit, level = c(0.95,0.99), ES = TRUE, is.its = FALSE)
#'}
#' @importFrom stats integrate sd uniroot                    
#' @export
risk <- function(object, theta, y, level = c(0.95,0.99), ES = TRUE, is.its = FALSE)
{
  UseMethod("risk", object)
}

#' @export
risk.MSGARCH_SPEC = function(object, theta, y, level = c(0.95,0.99), ES = TRUE, is.its = FALSE) {
  y = c(0,y)
  y = f.check.y(y)
  out = list()
  
  theta = f.check.theta(object, theta)
  
  ny =  nrow(y)
  if(!isTRUE(is.its)){
    start = ny 
    step = -(ny - 1)
    end  = ny
  } else {
    start = 2
    step = -1
    end  = ny - 1
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
  out$VaR = matrix(data = NA, nrow = end - start + 1, ncol = np)
  if (isTRUE(ES)) {
    out$ES = matrix(data = NA, nrow = end - start + 1, ncol = np)
  }
  for(v in start:end){
      
    f.pdf = function(x) {
      out = MSGARCH::pred(object, x, theta, y[1:v], log = FALSE)$pred
      return(out)
    }
    
    f.fun = function(x, pi) {
      out = integrate(f.pdf, lower = xmin, upper = x)$value - pi
      return(out)
    }
    
    # gross approximation for VaR
  
    for (i in 1:np) {
      tmp.VaR[i] = uniroot(f.fun, lower = xmin, upper = xmax, pi = p[i])$root
    }
    
    # add precision by Newton-Raphson
    
    for (i in 1:np) {
      p_i = p[i]
      level_i = level[i]
      calc.step = function(V) {
        PDF = MSGARCH::pred(object, x = V, theta = theta, y = y[1:v], log = F)$pred
        CDF = MSGARCH::pit(object, x = V, theta = theta, y = y[1:v])$pit
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
      out$VaR[v+step,i] = VaR
    }
    
    if (isTRUE(ES)) {
      
      for (i in 1:np) {
        
        f.condMean = function(x) {
          out = x * MSGARCH::pred(object, x, theta, y[1:v], log = FALSE)$pred
          return(out)
        }
        out$ES[v+step,i] = integrate(f.condMean, lower = -Inf, upper = out$VaR[v+step,i], stop.on.error = FALSE)$value/p[i]
      }
    }
  }
  out$y = y[2:length(y)]
  colnames(out$VaR) = level
  colnames(out$ES) = level
  out$VaR = rbind(NA, out$VaR) 
  out$ES = rbind(NA, out$ES)
  class(out) = "MSGARCH_RISK"
  return(out)
}

#' @export
risk.MSGARCH_MLE_FIT = function(object, theta = NULL, y = NULL, level = c(0.95,0.99), ES = TRUE, is.its = FALSE) {
  
  return(MSGARCH::risk(object = object$spec,  theta = object$theta, y = object$y, level = level, ES = ES, is.its = is.its))
  
}

#' @export
risk.MSGARCH_BAY_FIT = function(object, theta = NULL, y = NULL, level = c(0.95,0.99), ES = TRUE, is.its = FALSE) {
  
  return(MSGARCH::risk(object = object$spec,  theta = object$theta, y = object$y, level = level, ES = ES, is.its = is.its))
  
}